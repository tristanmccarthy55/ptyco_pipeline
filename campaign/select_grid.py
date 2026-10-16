#!/usr/bin/env python
"""
@file select_grid.py
@brief Pick one stage of a staged outer search (campaign/psearch_plan.py): the minimum of the trials' final error over
       the stage's two coordinates (C1 x C3, or the A1 vector).

Same score and same rule as campaign/select_trial.py (whose final_error and quadratic_minimum it reuses): each trial's
final Fourier error, the minimum of a quadratic fitted in the 3 x 3 block about the best grid point when that quadratic
has a minimum inside the block, else the best grid point; a best point on the grid's edge is flagged. With fewer than
half the trials finished it falls back to the grid centre (the drawn start) and says so. Always writes, into --out,
stage<k>_best.json / stage<k>_best.txt ("<C1 A> <C3 A> <C12 A> <phi12 rad>") and stage<k>_trials.csv; with --final also
search_best.json / search_best.txt, read by run_recon_synthetic_ML.slurm PROBE_FROM. Exits 0.
  stage 1  --axes c1: coordinates are C1 / C3 [A]; A1 from the manifest (the drawn start)
  stage 2  --axes a1: coordinates are the A1 vector (C12 cos 2 phi12, C12 sin 2 phi12) [A]; C1 / C3 from --base
The truth (start.json's record) never enters the choice.

    python campaign/select_grid.py --manifest search_af_<leg>/stage1.tsv --stage 1 --axes c1 --out search_af_<leg>
    python campaign/select_grid.py --manifest search_af_<leg>/stage2.tsv --stage 2 --axes a1 \\
        --base search_af_<leg>/stage1_best.txt --final --out search_af_<leg>
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from select_trial import final_error, quadratic_minimum   # noqa: E402  (same score, same fit as the C1/C3 search)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--stage", type=int, required=True)
    ap.add_argument("--axes", choices=["c1", "a1"], required=True)
    ap.add_argument("--base", default=None, help="stage 2: the stage-1 answer (C1 C3 C12 phi12)")
    ap.add_argument("--final", action="store_true", help="also write search_best.{json,txt}")
    ap.add_argument("--min-done", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    root = os.path.dirname(os.path.dirname(os.path.abspath(a.manifest)))      # the repo: search_af_<leg>/ sits in it
    rows = list(csv.DictReader(open(a.manifest), delimiter="\t"))
    for r in rows:
        r["x"], r["y"] = float(r["x"]), float(r["y"])
        r["error"] = final_error(os.path.join(root, r["recon_dir"]))
    ux, uy = sorted({r["x"] for r in rows}), sorted({r["y"] for r in rows})
    E = np.full((len(ux), len(uy)), np.nan)
    for r in rows:
        E[ux.index(r["x"]), uy.index(r["y"])] = r["error"]
    done = int(np.isfinite(E).sum())
    rec = dict(stage=a.stage, axes=a.axes, n_trials=len(rows), n_finished=done, grid=dict(x=ux, y=uy))
    centre = np.array([np.median(ux), np.median(uy)])
    if done < max(1, a.min_done * len(rows)):
        xy = centre
        rec.update(status=f"FALLBACK: only {done} of {len(rows)} trials finished -- the drawn start", method="grid centre")
    else:
        i, j = np.unravel_index(np.nanargmin(E), E.shape)
        rec["best_grid"] = dict(x=ux[i], y=uy[j], error=float(E[i, j]))
        rec["edge"] = dict(x=bool(i in (0, len(ux) - 1)), y=bool(j in (0, len(uy) - 1)))
        xy = np.array([ux[i], uy[j]], dtype=float); method = "best grid point"
        blk = [(ux[ii], uy[jj], E[ii, jj]) for ii in range(max(0, i - 1), min(len(ux), i + 2))
               for jj in range(max(0, j - 1), min(len(uy), j + 2)) if np.isfinite(E[ii, jj])]
        # fitted in coordinates normalised to the grid (C3 sits near -4e5 A: a raw fit is ill-conditioned)
        sc = np.array([max(np.ptp(ux), 1e-12), max(np.ptp(uy), 1e-12)])
        nb = [((x - centre[0]) / sc[0], (y - centre[1]) / sc[1], e) for x, y, e in blk]
        q = quadratic_minimum(np.array(nb), [0, 1]) if len(nb) >= 6 else None
        if q is not None:
            xy = centre + q * sc; method = "quadratic about the best grid point"
        warn = [k for k, v in rec["edge"].items() if v]
        rec.update(status="ok" + (f" -- best grid point on the {' and '.join(warn)} EDGE: the truth may lie outside" if warn else ""),
                   method=method)
    if a.axes == "c1":
        c1, c3 = float(xy[0]), float(xy[1])
        c12, phi = float(rows[0]["c12_A"]), float(rows[0]["phi12"])
    else:
        b = open(a.base).read().split()
        c1, c3 = float(b[0]), float(b[1])
        c12, phi = float(np.hypot(*xy)), float(np.arctan2(xy[1], xy[0]) / 2)
    rec.update(x=float(xy[0]), y=float(xy[1]), c1_A=c1, c3_A=c3, c12_A=c12, phi12=phi)
    os.makedirs(a.out, exist_ok=True)
    line = f"{c1:.4f} {c3:.2f} {c12:.4f} {phi:.6f}\n"
    names = [f"stage{a.stage}_best"] + (["search_best"] if a.final else [])
    for nm in names:
        json.dump(rec, open(os.path.join(a.out, nm + ".json"), "w"), indent=1)
        open(os.path.join(a.out, nm + ".txt"), "w").write(line)
    with open(os.path.join(a.out, f"stage{a.stage}_trials.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["idx", "x", "y", "c1_A", "c3_A", "c12_A", "phi12", "final_error", "recon_dir"])
        for r in rows:
            w.writerow([r["idx"], r["x"], r["y"], r["c1_A"], r["c3_A"], r["c12_A"], r["phi12"], r["error"], r["recon_dir"]])
    print(json.dumps({k: rec[k] for k in ("stage", "status", "method", "c1_A", "c3_A", "c12_A", "phi12", "n_finished", "n_trials")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
