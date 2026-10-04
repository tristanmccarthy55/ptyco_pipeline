#!/usr/bin/env python
"""
@file select_trial.py
@brief Pick an outer search's start probe: the minimum of the trials' final reconstruction error over (C1, C3).

campaign/run_thin_atomfind.sh SEARCH_C1F / SEARCH_C3F reconstructs the same data once per (C1, C3) pair with that trial
probe FIXED (a job array; search_af_<leg>/trials.tsv lists idx, f_c1, f_c3, c1_A, c3_A, recon_dir, the factors being on
the kicked start) and runs this once every trial has ended. The score is each trial's final Fourier error, the last row
of its *_error_trace.csv -- the solver's own residual, as in the step-1 C1 search, where 50 iterations ranked the trials
as 200 did. The answer is the minimum of a quadratic fitted to the trials in the 3 x 3 block about the best grid point,
when that quadratic has a true minimum inside the block; otherwise the best grid point itself. A best grid point on the
edge of the grid is flagged: the search range was too small on that axis.

Always writes, into --out: search_best.json (the record), search_best.txt ("<C1 A> <C3 A>", read by
run_recon_synthetic_ML.slurm PROBE_FROM), search_trials.csv and search_surface.png; and exits 0. With fewer than half the
trials finished it falls back to the kicked start and says so in "status", so the refined legs and the pack downstream
still run. The simulation's own C1 / C3 (a trial's aberrations.json) are recorded for the analysis only; they never enter
the choice.

    python campaign/select_trial.py --manifest search_af_<leg>/trials.tsv --out search_af_<leg>
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os

import numpy as np


def final_error(recon_dir):
    """The last error in the trial's newest *_error_trace.csv, or nan when the trial left none (crashed, timed out)."""
    traces = sorted(glob.glob(os.path.join(recon_dir, "analysis", "**", "*_error_trace.csv"), recursive=True),
                    key=os.path.getmtime)
    if not traces:
        return np.nan
    last = np.nan
    with open(traces[-1]) as fh:
        for line in fh:
            if line[:1].isdigit():
                try:
                    last = float(line.strip().split(",")[1])
                except (IndexError, ValueError):
                    pass
    return last


def quadratic_minimum(pts, axes):
    """Stationary point of a least-squares quadratic through pts (rows: f_c1, f_c3, error) over the varying axes, if it is
    a minimum inside the points' bounding box; else None. axes: which of (c1, c3) vary (a 1-D search has one)."""
    X = pts[:, :2][:, axes]; E = pts[:, 2]
    if X.shape[1] == 1:
        A = np.c_[np.ones(len(X)), X[:, 0], X[:, 0] ** 2]
        if len(X) < 3:
            return None
        a, b, c = np.linalg.lstsq(A, E, rcond=None)[0]
        if c <= 0:
            return None
        s = np.array([-b / (2 * c)])
    else:
        x, y = X[:, 0], X[:, 1]
        A = np.c_[np.ones(len(x)), x, y, x * x, x * y, y * y]
        if len(x) < 6 or np.linalg.matrix_rank(A) < 6:
            return None
        a, b, c, d, e, g = np.linalg.lstsq(A, E, rcond=None)[0]
        H = np.array([[2 * d, e], [e, 2 * g]])
        if np.any(np.linalg.eigvalsh(H) <= 0):
            return None
        s = np.linalg.solve(H, -np.array([b, c]))
    if np.any(s < X.min(0) - 1e-12) or np.any(s > X.max(0) + 1e-12):
        return None
    out = np.ones(2); out[np.array(axes)] = s
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, help="search_af_<leg>/trials.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-done", type=float, default=0.5, help="fraction of trials that must have finished")
    a = ap.parse_args()
    root = os.path.dirname(os.path.dirname(os.path.abspath(a.manifest)))      # the repo: search_af_<leg>/ sits in it
    rows = list(csv.DictReader(open(a.manifest), delimiter="\t"))
    for r in rows:
        r["f_c1"], r["f_c3"] = float(r["f_c1"]), float(r["f_c3"])
        r["c1_A"], r["c3_A"] = float(r["c1_A"]), float(r["c3_A"])
        r["error"] = final_error(os.path.join(root, r["recon_dir"]))
    kick_c1, kick_c3 = rows[0]["c1_A"] / rows[0]["f_c1"], rows[0]["c3_A"] / rows[0]["f_c3"]
    u1, u3 = sorted({r["f_c1"] for r in rows}), sorted({r["f_c3"] for r in rows})
    E = np.full((len(u1), len(u3)), np.nan)
    for r in rows:
        E[u1.index(r["f_c1"]), u3.index(r["f_c3"])] = r["error"]
    done = int(np.isfinite(E).sum())
    rec = dict(n_trials=len(rows), n_finished=done, kicked=dict(c1_A=kick_c1, c3_A=kick_c3),
               grid=dict(f_c1=u1, f_c3=u3))
    sim = None
    for r in rows:                                       # the truth, for the record only
        p = os.path.join(root, r["recon_dir"], "01", "aberrations.json")
        if os.path.exists(p):
            ab = json.load(open(p))
            c1s, c3s = float(ab["defocus_A"]), float(ab["aberrations_Cnm_A_phi_rad"].get("C30", np.nan))
            sim = dict(c1_A=c1s, c3_A=c3s, f_c1=c1s / kick_c1, f_c3=c3s / kick_c3)
            break
    rec["sim_record_only"] = sim
    if done < max(1, a.min_done * len(rows)):
        f = np.ones(2)
        rec.update(status=f"FALLBACK: only {done} of {len(rows)} trials finished -- starting from the kicked probe",
                   method="kicked start")
    else:
        i, j = np.unravel_index(np.nanargmin(E), E.shape)
        rec["best_grid"] = dict(f_c1=u1[i], f_c3=u3[j], error=float(E[i, j]))
        rec["edge"] = dict(c1=bool(len(u1) > 1 and i in (0, len(u1) - 1)), c3=bool(len(u3) > 1 and j in (0, len(u3) - 1)))
        axes = [k for k, n in enumerate((len(u1), len(u3))) if n >= 3]
        f = np.array([u1[i], u3[j]], dtype=float); method = "best grid point"
        if axes:
            blk = [(u1[ii], u3[jj], E[ii, jj]) for ii in range(max(0, i - 1), min(len(u1), i + 2))
                   for jj in range(max(0, j - 1), min(len(u3), j + 2)) if np.isfinite(E[ii, jj])]
            q = quadratic_minimum(np.array(blk), axes)
            if q is not None:
                f[axes] = q[axes]; method = "quadratic about the best grid point"
        warn = [k for k, v in rec["edge"].items() if v]
        rec.update(status="ok" + (f" -- best grid point on the {' and '.join(warn)} EDGE: widen the search there" if warn else ""),
                   method=method)
    rec.update(f_c1=float(f[0]), f_c3=float(f[1]), c1_A=float(kick_c1 * f[0]), c3_A=float(kick_c3 * f[1]))
    os.makedirs(a.out, exist_ok=True)
    json.dump(rec, open(os.path.join(a.out, "search_best.json"), "w"), indent=1)
    with open(os.path.join(a.out, "search_best.txt"), "w") as fh:
        fh.write(f"{rec['c1_A']:.4f} {rec['c3_A']:.2f}\n")
    with open(os.path.join(a.out, "search_trials.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["idx", "f_c1", "f_c3", "c1_A", "c3_A", "final_error", "recon_dir"])
        for r in rows:
            w.writerow([r["idx"], r["f_c1"], r["f_c3"], r["c1_A"], r["c3_A"], r["error"], r["recon_dir"]])
    try:                                                 # a figure must never cost the answer
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6.4, 5.2))
        if len(u1) > 1 and len(u3) > 1:
            d1, d3 = np.diff(u1).mean() / 2, np.diff(u3).mean() / 2
            im = ax.imshow(E, origin="lower", cmap="viridis", aspect="auto",
                           extent=[u3[0] - d3, u3[-1] + d3, u1[0] - d1, u1[-1] + d1])
            fig.colorbar(im, ax=ax, label="final error (fixed-probe trial)")
            ax.plot(rec["f_c3"], rec["f_c1"], "w+", ms=14, mew=2, label="chosen start")
            if sim:
                ax.plot(sim["f_c3"], sim["f_c1"], "rx", ms=10, mew=2, label="simulation (record only)")
            ax.set_xlabel("C3 / kicked C3"); ax.set_ylabel("C1 / kicked C1")
        else:
            x = u1 if len(u1) > 1 else u3
            ax.plot(x, E.ravel(), "o-"); ax.axvline(rec["f_c1"] if len(u1) > 1 else rec["f_c3"], color="k", label="chosen")
            if sim:
                ax.axvline(sim["f_c1"] if len(u1) > 1 else sim["f_c3"], color="r", ls="--", label="simulation (record only)")
            ax.set_xlabel(("C1" if len(u1) > 1 else "C3") + " / kicked value"); ax.set_ylabel("final error")
        ax.legend(fontsize=8, loc="upper right")
        ax.set_title(f"outer search: {done}/{len(rows)} trials; {rec['method']}", fontsize=9)
        fig.tight_layout(); fig.savefig(os.path.join(a.out, "search_surface.png"), dpi=110); plt.close(fig)
    except Exception as e:
        print(f"surface figure failed ({e})")
    print(json.dumps({k: rec[k] for k in ("status", "method", "f_c1", "f_c3", "c1_A", "c3_A", "n_finished", "n_trials")}))
    if sim:
        print(f"record only -- the simulation's C1/C3 sit at f_c1 {sim['f_c1']:.4f}, f_c3 {sim['f_c3']:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
