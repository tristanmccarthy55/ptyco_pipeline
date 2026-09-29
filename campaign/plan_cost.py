#!/usr/bin/env python
"""@file plan_cost.py
@brief What a large-probe leg would cost to simulate and reconstruct: the geometry options and their price.

For each sweep row (read from campaign/ceos_sweep.tsv, commented-out planned rows included) the probe is built on a
box that holds it, and every geometry the pipeline can make is enumerated:

  region side S   the in-plane cut of the tiled crystal (simulate_4dstem --region-side); 210 A exists (GT cache on
                  Blythe and locally), any other side needs a new simulation region and make_gt_cache --region-side S
  BIN b           the recon window is W = S / b, an integer division (the sim bins its fine detector by b)
  detector        recorded to 1.2 alpha (the presolve guard's minimum) or 1.66 alpha (what 80 mrad ran at)
  positions       1600 (20 A at 0.5 A, the default), 900 (15 A), 625 (12.5 A): the step is what matters, not the field

Probe sizes are the planner's (read from the row's note: d90/d99 about the probe's brightest point, the measure every
plan used; for these speckled probes it moves by up to ~15 % with the grid, so it is taken from one place). The probe is
built here only for the window loss, measured about the beam axis, which does not move with the grid.
A geometry is kept if the window holds the probe -- it loses no more intensity than LOSS_OK, the most a leg has lost
and still reconstructed (round_a070 from its 17.5 A window) -- and the region keeps scan + d99 + 20 A inside its side
(the seam rule of plan_probe.region_geometry). Its price, from measured runs:

  pattern N       2 theta_max W / lambda
  data            positions x N^2 x 4 bytes
  host peak       4 x data (a 1422 px recon peaks at ~4x its data; the node gives one job ~184 GB)
  GROUPING        16 (1426 / N)^2 on a 48 GB L40 (16 fits at 1426); below 1 it does not fit
  recon time      4.5 h (N / 1422)^2 (NL / 18) (positions / 1600) (1 + f^2) / (1 + 0.72^2): the CEOS 80 mrad legs,
                  1419 px, NL 18, 1600 positions, took 3.9-5.2 h; the 110 mrad BIN 1 leg (1426 px, NL 34) took ~7 h
                  against 8.5 predicted. f is the presolve's share of the pattern width: half, but widened to 1.2 alpha
                  where half would cut the aperture (run_synthetic_recon_ML.m) -- 0.72 for the 80 mrad reference at
                  1.66 alpha, the whole pattern (1.0) at a 1.2 alpha detector
  unproven        N above 1426 px, the largest pattern the engine has run
  sim grid        2 x 200 mrad x S / lambda: the simulator samples the potential to 200 mrad whatever it records, so its
                  grid grows with the region (4264 px at 210 A, the only side run so far)

    ~/hyperspy-bundle/bin/python campaign/plan_cost.py --labels arm_a070 arm_a080 ceosopt_a090 ceosopt_a100
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plan_probe import LAM, LOSS_OK, sizes, window_loss     # noqa: E402

TSV = os.path.join(HERE, "ceos_sweep.tsv")
BOX_Z = 27.525                    # the full reconstructed box: 5 cells + 2 x 4 A vacuum
SIDES = (210.0, 280.0, 350.0, 420.0, 560.0)
DET = (1.2, 1.66)
POSITIONS = ((1600, 20.0), (900, 15.0), (625, 12.5), (400, 10.0))
N_CAP = 12288                     # the largest probe grid built here (16 GB laptop); past it the leg is judged analytically
HOST_GB, WALL_H, PROVEN_N = 184.0, 48.0, 1426
T_REF_H, N_REF, NL_REF, POS_REF = 4.5, 1422, 18, 1600


def rows_all(path=TSV):
    out = {}
    with open(path) as f:
        lines = [l[1:] if l.startswith(("#arm_", "#ceosopt_")) else l for l in f]
    for r in csv.DictReader((l for l in lines if not l.startswith("#")), delimiter="\t"):
        out[r["label"]] = r
    return out


def nyquist_nl(alpha):
    return max(1, int(BOX_Z * 2 * (alpha / 1000.0) ** 2 / LAM + 0.5))


def probe(r, ext, n):
    import abtem
    ab = {k: float(v) for k, v in json.loads(r["aber_json"]).items()} if r["aber_json"] not in ("-", "") else \
        {"C30": float(r["c3"]), "C50": float(r["c5"])}
    return np.asarray(abtem.Probe(energy=300e3, semiangle_cutoff=float(r["alpha"]), extent=ext, gpts=n,
                                  defocus=float(r["c1"]), aberrations=ab).build(lazy=False).array)


def planned_sizes(note):
    """(d50, d90, d99) from a planner note: plan_arm writes 'd50/d90/d99 a/b/c A', plan_ceos 'd90 X A, d99 Y A'."""
    import re
    m = re.search(r"d50/d90/d99 ([0-9.]+)/([0-9.]+)/([0-9.]+)", note)
    if m:
        return tuple(float(x) for x in m.groups())
    m = re.search(r"d90 ([0-9.]+) A, d99 ([0-9.]+) A", note)
    return (float("nan"), float(m.group(1)), float(m.group(2))) if m else None


def cost(label, r):
    a = float(r["alpha"])
    nl = nyquist_nl(a)
    got = planned_sizes(r.get("note", ""))
    ext = float(np.ceil(2.6 * got[2] / 50) * 50) if got else 300.0
    # the most favourable geometry: window = d99, detector 1.2 alpha, the fewest positions -- if even that cannot fit
    # the node, no geometry can, and nothing is built
    if got:
        n_min = 2 * 1.2 * a / 1000.0 * got[2] / LAM
        d_min = POSITIONS[-1][0] * n_min ** 2 * 4 / 1e9
        if 4 * d_min > HOST_GB or 16 * (PROVEN_N / n_min) ** 2 < 1:
            return dict(label=label, alpha=a, d50=got[0], d90=got[1], d99=got[2], nl=nl, box=None,
                        best={str(f): None for f in DET}, options=[],
                        why=f"even a {got[2]:.0f} A window at 1.2 alpha is N {n_min:.0f} px: {POSITIONS[-1][0]} positions "
                            f"are {d_min:.0f} GB of data, {4 * d_min:.0f} GB peak")
    for _ in range(4):                                    # a box that holds the probe: d99 < 0.4 of it
        over = min(1.3, N_CAP * LAM / (2 * a / 1000.0 * ext))
        n = int(np.floor(2 * over * a / 1000.0 * ext / LAM / 64) * 64)
        P = probe(r, ext, n)
        _, d50, d90, d99 = sizes(P, ext)
        if d99 < 0.4 * ext:
            break
        ext = float(np.ceil(2.6 * d99 / 50) * 50)
    if got:
        d50, d90, d99 = got
    opts = []
    for S in SIDES:
        for b in range(1, 13):
            W = S / b
            if W > ext:
                continue
            loss = window_loss(P, ext, W)
            if loss > LOSS_OK:
                continue
            for npos, scan in POSITIONS:
                if scan + d99 + 20 > S:
                    continue
                for f in DET:
                    N = int(2 * round(f * a / 1000.0 * W / LAM))
                    data = npos * N * N * 4 / 1e9
                    grp = 16 * (PROVEN_N / N) ** 2
                    fp = min(1.0, max(0.5, 1.2 / f))              # presolve width / N (f = detector / alpha)
                    t = T_REF_H * (N / N_REF) ** 2 * (nl / NL_REF) * (npos / POS_REF) * (1 + fp ** 2) / (1 + 0.72 ** 2)
                    ok = 4 * data <= HOST_GB and grp >= 1 and t <= 0.85 * WALL_H
                    opts.append(dict(side=S, bin=b, window=round(W, 1), loss=round(100 * loss, 3), det=f,
                                     det_mrad=round(f * a), N=N, positions=npos, scan=scan, data_GB=round(data, 1),
                                     peak_GB=round(4 * data), grouping=int(min(32, grp)), hours=round(t, 1),
                                     new_region=S != 210.0, unproven=N > PROVEN_N, fits=ok,
                                     sim_px=int(round(2 * 0.2 * S / LAM))))
    # per detector: the geometry that fits with the most positions, then in the existing 210 A region, then fewest hours
    best = {}
    for f in DET:
        fit = sorted((o for o in opts if o["fits"] and o["det"] == f),
                     key=lambda o: (-o["positions"], o["new_region"], o["hours"]))
        best[str(f)] = fit[0] if fit else None
    return dict(label=label, alpha=a, d50=round(d50, 1), d90=round(d90, 1), d99=round(d99, 1), nl=nl,
                box=ext, best=best, options=opts, why="" if any(best.values()) else "host memory, GPU or walltime")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--labels", nargs="+", required=True)
    ap.add_argument("--out", default=None, help="write every option as JSON here")
    a = ap.parse_args()
    try:
        import abtem
        abtem.config.set({"local_diagnostics.progress_bar": False})
    except Exception:
        pass
    R = rows_all()
    res = []
    for lab in a.labels:
        c = cost(lab, R[lab]); res.append(c)
        print(f"{lab:16s} {c['alpha']:>4.0f} mrad  d90 {c['d90']:6.1f}  d99 {c['d99']:6.1f} A  NL {c['nl']:2d}", flush=True)
        for f, b in c["best"].items():
            if b is None:
                print(f"    detector {f} alpha: no geometry fits ({c['why'] or 'host memory, GPU or walltime'})", flush=True)
                continue
            print(f"    detector {f} alpha ({b['det_mrad']} mrad): side {b['side']:.0f} / BIN {b['bin']} = window "
                  f"{b['window']} A (loses {b['loss']} %), N {b['N']} px, {b['positions']} pos: data {b['data_GB']} GB, "
                  f"peak {b['peak_GB']} GB, GROUPING {b['grouping']}, ~{b['hours']} h"
                  + (f"  [new region + GT cache; sim grid {b['sim_px']} px]" if b["new_region"] else "")
                  + ("  [N unproven]" if b["unproven"] else ""),
                  flush=True)
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(res, fh, indent=1)
        print("wrote", a.out)


if __name__ == "__main__":
    main()
