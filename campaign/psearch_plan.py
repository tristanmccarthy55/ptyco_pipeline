#!/usr/bin/env python
"""
@file psearch_plan.py
@brief Plan a STAGED outer search from a start drawn at the corrector's quoted precision (C1, C3, A1).

The experiment-like route (the user, 2026-10-16): the operator's corrector reports C1, C3 and A1 with 95 % intervals; the
start probe is the truth plus errors DRAWN from those intervals (sigma = ci95 / 1.96, a seeded draw; A1 as a vector
a = C12 (cos 2 phi12, sin 2 phi12), each component drawn); the outer search is bounded by the same intervals:
  stage 1  C1 x C3 on an n1 x n1 grid, +- ci95 about the drawn start, A1 at the drawn start
  stage 2  the A1 vector on an n2 x n2 grid, +- ci95 per component about the drawn start, at stage 1's C1 / C3
  stage 3  (--stages 3, the default since 2026-10-20) C1 x C3 again on stage 1's grid, at stage 2's A1 -- the first
           run showed why: with A1 still ~40 A wrong, stage 1's C1 x C3 minimum slid to the grid corner (C1 32 A off),
           while stage 2 found A1 to ~4 A anyway. A1 is robust to a C1 error; C1 is not to an A1 error.
Higher orders stay at their true values (stated wherever the result is quoted). Writes into --out: start.json (the truth
for the record, the draw, the start, the grids), stage1.tsv and stage2.tsv (one row per trial: idx, the two search
coordinates, c1_A, c3_A, c12_A, phi12, recon_dir; "-" = taken in the job from the previous stage's answer: stage 2's
c1_A / c3_A, stage 3's c12_A / phi12); stage3.tsv with --stages 3.

    python campaign/psearch_plan.py --name armf_a080_lab_fs2p5_src0p1 --nl 18 --c1=-958.6 --c3=-406600 --c12 12.4 \\
        --phi12 1.467822 --ci95 "39.5 8040 59" --seed 1 --out search_af_<leg>
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np


def a1_vector(c12, phi12):
    return np.array([c12 * np.cos(2 * phi12), c12 * np.sin(2 * phi12)])


def a1_polar(v):
    return float(np.hypot(*v)), float(np.arctan2(v[1], v[0]) / 2)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True, help="the lab leg's recon name (trial dirs recon_af_<name>_s<k>t<i>_NL<nl>)")
    ap.add_argument("--nl", type=int, required=True)
    ap.add_argument("--c1", type=float, required=True, help="true C1 (abTEM defocus) [A]")
    ap.add_argument("--c3", type=float, required=True, help="true C3 [A]")
    ap.add_argument("--c12", type=float, default=0.0, help="true A1 magnitude (abTEM C12) [A]")
    ap.add_argument("--phi12", type=float, default=0.0, help="true A1 angle (abTEM phi12) [rad]")
    ap.add_argument("--ci95", default="39.5 8040 59", help="the corrector's 95 %% intervals of C1, C3, A1 [A]")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n1", type=int, default=7)
    ap.add_argument("--n2", type=int, default=5)
    ap.add_argument("--stages", type=int, choices=[2, 3], default=3)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    ci1, ci3, cia = (float(x) for x in a.ci95.split())
    rng = np.random.default_rng(a.seed)
    d1, d3 = rng.standard_normal(2) * np.array([ci1, ci3]) / 1.96
    da = rng.standard_normal(2) * cia / 1.96
    true_a = a1_vector(a.c12, a.phi12)
    s1, s3, sa = a.c1 + d1, a.c3 + d3, true_a + da
    s12, sphi = a1_polar(sa)
    g1 = s1 + np.linspace(-ci1, ci1, a.n1); g3 = s3 + np.linspace(-ci3, ci3, a.n1)
    gx = sa[0] + np.linspace(-cia, cia, a.n2); gy = sa[1] + np.linspace(-cia, cia, a.n2)
    os.makedirs(a.out, exist_ok=True)
    rows1 = [(i, x, y) for i, (x, y) in enumerate((x, y) for x in g1 for y in g3)]
    with open(os.path.join(a.out, "stage1.tsv"), "w") as fh:
        fh.write("idx\tx\ty\tc1_A\tc3_A\tc12_A\tphi12\trecon_dir\n")
        for i, x, y in rows1:
            fh.write(f"{i}\t{x:.4f}\t{y:.2f}\t{x:.4f}\t{y:.2f}\t{s12:.4f}\t{sphi:.6f}\trecon_af_{a.name}_s1t{i}_NL{a.nl}\n")
    with open(os.path.join(a.out, "stage2.tsv"), "w") as fh:
        fh.write("idx\tx\ty\tc1_A\tc3_A\tc12_A\tphi12\trecon_dir\n")
        for i, (x, y) in enumerate((x, y) for x in gx for y in gy):
            m, ph = a1_polar(np.array([x, y]))
            fh.write(f"{i}\t{x:.4f}\t{y:.4f}\t-\t-\t{m:.4f}\t{ph:.6f}\trecon_af_{a.name}_s2t{i}_NL{a.nl}\n")
    if a.stages == 3:
        with open(os.path.join(a.out, "stage3.tsv"), "w") as fh:
            fh.write("idx\tx\ty\tc1_A\tc3_A\tc12_A\tphi12\trecon_dir\n")
            for i, x, y in rows1:
                fh.write(f"{i}\t{x:.4f}\t{y:.2f}\t{x:.4f}\t{y:.2f}\t-\t-\trecon_af_{a.name}_s3t{i}_NL{a.nl}\n")
    inside = dict(c1=bool(abs(d1) <= ci1), c3=bool(abs(d3) <= ci3), a1x=bool(abs(da[0]) <= cia), a1y=bool(abs(da[1]) <= cia))
    rec = dict(seed=a.seed, ci95=dict(c1_A=ci1, c3_A=ci3, a1_A=cia), sigma="ci95 / 1.96 (A1: per vector component)",
               truth_record_only=dict(c1_A=a.c1, c3_A=a.c3, c12_A=a.c12, phi12=a.phi12, a1_vec=true_a.tolist()),
               draw=dict(dc1_A=float(d1), dc3_A=float(d3), da1_vec=da.tolist()),
               start=dict(c1_A=float(s1), c3_A=float(s3), c12_A=s12, phi12=sphi, a1_vec=sa.tolist()),
               truth_inside_grid=inside, grid1=dict(c1_A=g1.tolist(), c3_A=g3.tolist()),
               grid2=dict(a1x_A=gx.tolist(), a1y_A=gy.tolist()), stages=a.stages,
               n_trials=[len(rows1), a.n2 * a.n2] + ([len(rows1)] if a.stages == 3 else []))
    json.dump(rec, open(os.path.join(a.out, "start.json"), "w"), indent=1)
    print(f"[psearch] seed {a.seed}: start C1 {s1:.1f} (true {a.c1:g}, d {d1:+.1f}), C3 {s3:.0f} (d {d3:+.0f}), "
          f"A1 {s12:.1f} A @ {sphi:.3f} rad (true {a.c12:g} @ {a.phi12:.3f}; d ({da[0]:+.1f}, {da[1]:+.1f})); "
          f"{len(rows1)} + {a.n2 * a.n2}{f' + {len(rows1)}' if a.stages == 3 else ''} trials; truth inside the grids: {inside}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
