#!/usr/bin/env python
"""
@file depth_honesty.py
@brief Does the atom finder's score measure DEPTH? The same found atoms scored as reported, with their depths randomised,
       and moved half a lattice period; plus a strict score and a finder-free depth contrast from the volume.

Why (2026-10-08): validate.match_found_to_gt matches a found atom to any unclaimed GT atom within 0.6 A in plane and
2.0 A in depth, whatever its species. Along a Pb column (an atom every 3.9 A) every depth is within 1.95 A of a Pb, and
the depth rms is taken over matched atoms only, so it cannot exceed ~2/sqrt(3) = 1.15 A. Recall and "depth error" can
therefore look respectable with no depth information at all. This scores, per run (manifest keys of make_arm_page):

  reported    validate.finder_report as analyse_sweep ran it (tol 0.6 A / 2.0 A, any species)
  random z    the same atoms, depths drawn uniformly over the slab (mean of 20 draws) -- no depth information
  shift 1.95  the same atoms, every depth moved by half the Pb period -- deliberately wrong depths
  strict      tol 0.6 A / 0.5 A AND the found species equals the GT species
  frame0      'reported' in the GT frame registered on the volume only (align.register), not refined on these atoms
  strict null 'strict' for the random-depth atoms: what strict scoring gives with no depth information
  (finder-free depth measures were tried and dropped: a correlation with the atom template failed its own control --
  the half-period-shifted template scored the same, the slab envelope -- and a 3.9 A lock-in on each column's depth
  profile did not separate good runs from bad. The depth blur is wider than the atom spacing, so depth lives in the
  kernel fit, not in visible peaks: see results/2026-W41/atomfind_honesty/ for the kernel-shift test that does work.)
  vacuum      mean |phase| in the vacuum layers over the slab's: how much the reconstruction dumps outside the specimen

The GT frame is the one the reported score used (refined on the run's own found atoms), so the nulls are scored in
exactly the same frame.

    ~/hyperspy-bundle/bin/python analysis/depth_honesty.py [--runs known_a080 fixed_1e7 ...] [--csv out.csv]
"""
from __future__ import annotations

import argparse
import copy
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_arm_page as mp                                    # noqa: E402  (Runs, manifest, region GT cache)
import figdata                                                # noqa: E402

DEFAULT = ["known_a050", "known_a060", "known_a065", "known_a070", "known_a080", "known_a090", "known_a100",
           "win_070", "win_052", "kick_a060", "long_m3", "long_m6", "fixed_1e7", "fixed_1e8", "hail_1e7"]


def frame(info):
    """Volume, GT, found atoms and the alignment exactly as the reported score built them."""
    from atomfind import align
    V, meta = figdata.load(info["fd"])
    dx, dz, centre = meta["dx"], meta["dz"], meta["centre"]
    cfg = figdata.config_for(os.path.join(info["fd"], "phase_vol.npy"), mp.gt_dir(int(round(2 * centre[0]))), dz, dx, centre)
    pos, Z = align.load_gt(cfg)
    found = np.load(os.path.join(info["af"], "found_atoms.npy"))
    al = align.refine_with_atoms(align.register(V, dx, pos, Z, cfg), found, pos, Z, cfg)
    return V, cfg, pos, Z, found, al


def score(found, pos, Z, al, cfg, tol_z=None, species=False):
    """recall per species (bulk), precision, depth rms -- validate.finder_report, optionally strict."""
    from atomfind import validate
    c = copy.copy(cfg)
    if tol_z is not None:
        c.match_tol_z_A = tol_z
    if not species:
        r, _ = validate.finder_report(found, pos, Z, al, c)
        return dict(Pb=r["Pb"]["recall_bulk"], Ti=r["Ti"]["recall_bulk"], O=r["O"]["recall_bulk"], prec=r["precision"], z=r["z_rms_A"])
    # species-required: score each found species against GT of that species only, then pool
    out, nmatch, nfound, dzs = {}, 0, len(found), []
    from atomfind import align as _al
    win = _al.in_window(pos, c)
    for zz, nm in ((82, "Pb"), (22, "Ti"), (8, "O")):
        f = found[found["species"] == zz]
        keep = (Z == zz)
        r, m = validate.finder_report(f, pos[keep], Z[keep], al, c) if len(f) else (None, None)
        # recall_bulk in finder_report is per GT species; with one species in play it is that species' recall
        out[nm] = r[nm]["recall_bulk"] if r else 0.0
        if r:
            nmatch += r["n_matched"]
            dzs += list(m["match_dz"][m["match_gi"] >= 0])
    out["prec"] = nmatch / max(nfound, 1)
    out["z"] = float(np.sqrt(np.mean(np.square(dzs)))) if dzs else np.nan
    return out


def columns(pos, win, tol=0.6):
    """Column id per atom: atoms within tol A of each other in plane share a column (single linkage), so off-centre Ti
    stays in its B-site column with the apical O. -1 outside the window."""
    idx = np.where(win)[0]
    xy = pos[idx, :2]
    lab = -np.ones(len(pos), int)
    parent = np.arange(len(idx))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    order = np.argsort(xy[:, 0])
    for a_ in range(len(order)):
        i = order[a_]
        for b_ in range(a_ + 1, len(order)):
            j = order[b_]
            if xy[j, 0] - xy[i, 0] > tol:
                break
            if np.hypot(*(xy[j] - xy[i])) <= tol:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[rj] = ri
    lab[idx] = [find(i) for i in range(len(idx))]
    return lab


def vacuum_ratio(V, cfg):
    nl = V.shape[0]
    z = (np.arange(nl) + 0.5) * cfg.dz
    vac = (z < cfg.trim_z_A[0]) | (z > cfg.trim_z_A[1])
    a = np.abs(V - np.median(V))
    return float(a[vac].mean() / a[~vac].mean())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", default=DEFAULT)
    ap.add_argument("--csv", default=os.path.join(mp.REPO, "aberration_experiment", "results", "2026-W41", "depth_honesty.csv"))
    a = ap.parse_args()
    import json
    R = mp.Runs(json.load(open(mp.MANIFEST)))
    rng = np.random.default_rng(0)
    rows = []
    hdr = (f"{'run':12s} | {'reported Pb/Ti/O prec z':27s} | {'frame not fitted to them':27s} | {'random depths':27s} | {'shifted 1.95 A':27s} | "
           f"{'strict 0.5 A + species':21s} | {'strict, random depths':21s} | vacuum")
    print(hdr); print("-" * len(hdr))
    for key in a.runs:
        info = R.find(key)
        if not (info["found"] and info["fd"] and info["af"]):
            print(f"{key:12s} | (no volume or found atoms)"); continue
        V, cfg, pos, Z, found, al = frame(info)
        rep = score(found, pos, Z, al, cfg)
        from atomfind import align as _al
        al0 = _al.register(V, cfg.dx if hasattr(cfg, "dx") else al.dx, pos, Z, cfg)   # frame NOT fitted to these atoms
        rep0 = score(found, pos, Z, al0, cfg)
        lo, hi = cfg.trim_z_A[0] / cfg.dz, cfg.trim_z_A[1] / cfg.dz
        rnd, rnds = [], []
        for _ in range(20):
            f = found.copy(); f["layer"] = rng.uniform(lo, hi, len(f))
            rnd.append(score(f, pos, Z, al, cfg)); rnds.append(score(f, pos, Z, al, cfg, tol_z=0.5, species=True))
        mean = lambda L: {k: float(np.nanmean([r[k] for r in L])) for k in L[0]}
        rnd, rnds = mean(rnd), mean(rnds)
        f = found.copy(); f["layer"] = f["layer"] + 1.95 / cfg.dz
        sh = score(f, pos, Z, al, cfg)
        st = score(found, pos, Z, al, cfg, tol_z=0.5, species=True)
        vr = vacuum_ratio(V, cfg)
        fmt = lambda r: f"{100*r['Pb']:3.0f}/{100*r['Ti']:3.0f}/{100*r['O']:3.0f} {r['prec']:4.2f} {r['z']:4.2f}"
        fms = lambda r: f"{100*r['Pb']:3.0f}/{100*r['Ti']:3.0f}/{100*r['O']:3.0f} {r['prec']:4.2f}"
        print(f"{key:12s} | {fmt(rep):27s} | {fmt(rep0):27s} | {fmt(rnd):27s} | {fmt(sh):27s} | {fms(st):21s} | {fms(rnds):21s} | "
              f"{vr:4.2f}", flush=True)
        rows.append(dict(run=key, tag=info["tag"], **{f"rep_{k}": v for k, v in rep.items()}, **{f"rnd_{k}": v for k, v in rnd.items()},
                         **{f"shift_{k}": v for k, v in sh.items()}, **{f"strict_{k}": v for k, v in st.items()},
                         **{f"frame0_{k}": v for k, v in rep0.items()}, **{f"strictnull_{k}": v for k, v in rnds.items()},
                         vacuum_ratio=vr))
    os.makedirs(os.path.dirname(a.csv), exist_ok=True)
    with open(a.csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"wrote {a.csv}")


if __name__ == "__main__":
    main()
