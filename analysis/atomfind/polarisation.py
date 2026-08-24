#!/usr/bin/env python
"""@file polarisation.py
@brief Ti-O6 off-centring (local polarisation proxy) from the blind found atoms, with
uncertainty propagated from the calibrated conformal intervals.

Closes the loop from located atoms to the physical observable: for each located Ti,
delta = r_Ti - centroid(its six nearest located O) -- the same B-site off-centring proxy
used on the ground-truth model in analysis/figures/pol_vortex.py. Uncertainty on delta is
propagated by Monte Carlo from the per-atom 95% conformal half-widths exported by
run_atomfind.py, i.e. from quantities available WITHOUT ground truth. Ground truth enters
only to score the result.

Run:  python atomfind/polarisation.py [--out DIR] [--data-dir DIR]   (cwd analysis/)
Outputs <out_dir>/polarisation.npz + a printed report.
"""
from __future__ import annotations
import os, sys, csv, json
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from atomfind import config, align

CAGE_CUT_A = 2.8       # O-neighbour cutoff for the octahedral cage (Ti-O ~1.95 A)
N_CAGE = 6             # a complete octahedron
N_MC = 400             # Monte-Carlo draws for the propagated interval
MATCH_TOL_A = 0.6      # found-Ti <-> GT-Ti correspondence gate


def load_found(path):
    """@brief Read found_atoms.csv -> (element, xyz, 95% half-widths, guided flag)."""
    rows = list(csv.DictReader(open(path)))
    el = np.array([r["element"] for r in rows])
    xyz = np.array([[float(r["X_A"]), float(r["Y_A"]), float(r["Z_A"])] for r in rows])
    hw = np.array([[float(r["halfwidth95_x_A"]), float(r["halfwidth95_y_A"]),
                    float(r["halfwidth95_z_A"])] for r in rows])
    guided = np.array([int(r["guided"]) for r in rows])
    return el, xyz, hw, guided


def offcentring(ti, ox, cut=CAGE_CUT_A, need=N_CAGE):
    """@brief delta = r_Ti - centroid(`need` nearest O within `cut`).
    @return (delta [n,3] with NaN where the cage is incomplete, ok mask, neighbour idx)."""
    d, idx = cKDTree(ox).query(ti, k=need, distance_upper_bound=cut)
    ok = np.all(np.isfinite(d), axis=1)
    delta = np.full((len(ti), 3), np.nan)
    delta[ok] = ti[ok] - ox[idx[ok]].mean(axis=1)
    return delta, ok, idx


def analyse(cfg=None, csv_path=None, out_dir=None, seed=0):
    """@brief Full pipeline: cages from found atoms, MC uncertainty, scoring vs GT."""
    cfg = (cfg or config.preset("NL70_coherent")).resolve()
    out_dir = out_dir or cfg.out_dir
    csv_path = csv_path or os.path.join(out_dir, "found_atoms.csv")
    zlo, zhi = cfg.bulk_z_A
    cx, cy = cfg.scan_center_xy
    h = cfg.scan_window_A / 2.0
    R = {}

    # -- ground truth, same window and depth band ---------------------------
    pos, Z = align.load_gt(cfg)
    inwin = (np.abs(pos[:, 0] - cx) < h) & (np.abs(pos[:, 1] - cy) < h)
    gti = pos[(Z == 22) & inwin & (pos[:, 2] > zlo) & (pos[:, 2] < zhi)]
    gdel, gok, _ = offcentring(gti, pos[Z == 8])     # cage from the whole cell
    R["gt_n_Ti"], R["gt_n_caged"] = len(gti), int(gok.sum())
    R["gt_mean_abs_delta_A"] = float(np.linalg.norm(gdel[gok], axis=1).mean())

    # -- reconstructed ------------------------------------------------------
    el, xyz, hw, _ = load_found(csv_path)
    fti_all, fti_hw_all = xyz[el == "Ti"], hw[el == "Ti"]
    fox, fox_hw = xyz[el == "O"], hw[el == "O"]
    keep = (fti_all[:, 2] > zlo) & (fti_all[:, 2] < zhi)
    fti, fti_hw = fti_all[keep], fti_hw_all[keep]
    fdel, fok, _ = offcentring(fti, fox)
    R["rec_n_Ti"], R["rec_n_caged"] = len(fti), int(fok.sum())
    R["rec_cage_completeness"] = float(fok.mean())

    # -- correspondence -----------------------------------------------------
    d, j = cKDTree(gti[gok]).query(fti[fok], distance_upper_bound=MATCH_TOL_A)
    m = np.isfinite(d)
    ti_m = fti[fok][m]
    A, B = fdel[fok][m], gdel[gok][j[m]]             # reconstructed / true off-centring
    R["n_matched"] = int(m.sum())

    # -- Monte-Carlo propagated uncertainty (no GT) -------------------------
    rng = np.random.default_rng(seed)
    sTi = fti_hw[fok][m] / 1.96                      # conformal 95% half-width -> sigma
    _, i6 = cKDTree(fox).query(ti_m, k=N_CAGE, distance_upper_bound=CAGE_CUT_A)
    sO = fox_hw[i6] / 1.96
    draws = np.empty((N_MC, len(ti_m), 3))
    for s in range(N_MC):
        draws[s] = ((ti_m + rng.normal(size=sTi.shape) * sTi)
                    - (fox[i6] + rng.normal(size=sO.shape) * sO).mean(axis=1))
    sig = draws.std(axis=0)

    # -- scoring ------------------------------------------------------------
    err = A - B
    for k, ax in enumerate("xyz"):
        e = np.abs(err[:, k])
        R[f"delta_{ax}"] = dict(
            rms=float(np.sqrt((err[:, k] ** 2).mean())), median=float(np.median(e)),
            p90=float(np.percentile(e, 90)), p99=float(np.percentile(e, 99)),
            bias=float(err[:, k].mean()), gt_sd=float(B[:, k].std()),
            corr=float(np.corrcoef(A[:, k], B[:, k])[0, 1]),
            sigma_mc=float(np.median(sig[:, k])),
            coverage95=float((np.abs(err[:, k]) <= 1.96 * sig[:, k]).mean()))
    nA, nB = np.linalg.norm(A[:, :2], axis=1), np.linalg.norm(B[:, :2], axis=1)
    ang = np.degrees(np.arccos(np.clip((A[:, :2] * B[:, :2]).sum(1) / (nA * nB + 1e-12), -1, 1)))
    ip = np.linalg.norm(err[:, :2], axis=1)
    R["in_plane"] = dict(gt_mean=float(nB.mean()), rec_mean=float(nA.mean()),
                         rms=float(np.sqrt((ip ** 2).mean())), median=float(np.median(ip)),
                         p90=float(np.percentile(ip, 90)),
                         ang_median=float(np.median(ang)), ang_p90=float(np.percentile(ang, 90)),
                         frac_within_30deg=float((ang < 30).mean()))
    # the tail: Ti whose cage composition is wrong (a missing/misassigned O moves the centroid)
    bad = ip > 0.05
    R["tail"] = dict(frac=float(bad.mean()), n=int(bad.sum()),
                     median_err_good=float(np.median(ip[~bad])),
                     median_err_bad=float(np.median(ip[bad])),
                     cov95_good=float((np.abs(err[~bad, 0]) <= 1.96 * sig[~bad, 0]).mean()),
                     cov95_bad=float((np.abs(err[bad, 0]) <= 1.96 * sig[bad, 0]).mean()),
                     corr_x_excl_tail=float(np.corrcoef(A[~bad, 0], B[~bad, 0])[0, 1]))

    np.savez(os.path.join(out_dir, "polarisation.npz"), ti=ti_m, delta=A, delta_gt=B,
             sigma=sig, tail=bad)
    with open(os.path.join(out_dir, "polarisation.json"), "w") as f:
        json.dump(R, f, indent=2)
    return R


def report(R):
    """@brief Human-readable summary of analyse()."""
    L = [f"[GT ] Ti in window+bulk {R['gt_n_Ti']}, complete O6 cage {R['gt_n_caged']}, "
         f"mean |delta| {R['gt_mean_abs_delta_A']:.3f} A",
         f"[REC] Ti in bulk {R['rec_n_Ti']}, complete cage from FOUND O {R['rec_n_caged']} "
         f"({100*R['rec_cage_completeness']:.0f}%), matched to GT {R['n_matched']}"]
    for ax in "xyz":
        d = R[f"delta_{ax}"]
        L.append(f"  delta_{ax}: median|e| {d['median']:.3f}  p90 {d['p90']:.3f}  "
                 f"RMS {d['rms']:.3f} A | GT sd {d['gt_sd']:.3f} | r {d['corr']:.3f} | "
                 f"MC sigma {d['sigma_mc']:.3f} | 95% cov {100*d['coverage95']:.0f}%")
    p = R["in_plane"]
    L.append(f"  IN-PLANE: |delta| GT {p['gt_mean']:.3f} vs rec {p['rec_mean']:.3f} A; "
             f"vector error median {p['median']:.3f} p90 {p['p90']:.3f} RMS {p['rms']:.3f} A; "
             f"direction median {p['ang_median']:.1f} deg, {100*p['frac_within_30deg']:.0f}% <30 deg")
    t = R["tail"]
    L.append(f"  TAIL (cage-composition errors): {100*t['frac']:.0f}% of Ti, median error "
             f"{t['median_err_bad']:.3f} A vs {t['median_err_good']:.3f} A for the rest; "
             f"95% interval covers {100*t['cov95_good']:.0f}% of the good but only "
             f"{100*t['cov95_bad']:.0f}% of the tail (intervals are conditional on correct detection)")
    return "\n".join(L)


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Ti-O6 off-centring from the located atoms.")
    ap.add_argument("--preset", default="NL70_coherent")
    ap.add_argument("--out", default=None,
                    help="directory holding found_atoms.csv, and where the outputs go "
                         "(default $ATOMFIND_OUT, else ./atomfind_out)")
    ap.add_argument("--data-dir", default=None, help="see run_atomfind.py --data-dir")
    ap.add_argument("--seed", type=int, default=0, help="Monte-Carlo seed")
    a = ap.parse_args()
    if a.data_dir:
        config.set_data_dir(a.data_dir)
    cfg = config.preset(a.preset)
    if a.out:
        cfg.out_dir = a.out
    print(report(analyse(cfg=cfg, out_dir=cfg.out_dir, seed=a.seed)))


if __name__ == "__main__":
    main()
