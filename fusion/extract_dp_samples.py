#!/usr/bin/env python
"""@file extract_dp_samples.py
@brief Pull a few MB of diffraction patterns out of a multi-GB simulation, for figures.

The full `data_dp.hdf5` is 1.5-6 GB and has no business being copied to a laptop to draw one panel.
This reads it in chunks where it lives and keeps what a detector-plane figure needs: the mean
pattern, the mean inside each domain's interior, a few single patterns at named probe positions,
and the angle map the hollow mask is cut from. CPU-only and I/O-bound -- fine on a login node.

Patterns are returned in the SIMULATION convention: `save_outputs` stores the transpose for MATLAB,
which is undone here.

    <abtem env>/python extract_dp_samples.py --run fusion/runs/fusion --out fig_dp.npz
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True, help="simulation dir containing 01/ and hollow_budget.json")
    ap.add_argument("--truth", default=os.path.join(HERE, "sample", "toy_truth.npz"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--margin", type=float, default=6.0, help="clearance from a domain wall [A]")
    args = ap.parse_args(argv)

    import h5py
    import simulate_fusion as SF
    import sign_test as ST

    b = json.load(open(os.path.join(args.run, "hollow_budget.json")))
    truth = np.load(args.truth, allow_pickle=True)
    with h5py.File(os.path.join(args.run, "01", "data_position.hdf5"), "r") as f:
        pos = np.asarray(f["probe_positions_0"][...], dtype=float).T
    sel = ST.interior(pos, truth, b, args.margin)
    names = [str(n) for n in truth["names"]]

    # named single patterns: nearest scan positions to the three distinct columns of one cell --
    # Pb at (0,0), Ti+O_ap at (a/2,a/2), O_eq at (a/2,0) -- inside domain A
    a = float(truth["a"])
    iA = np.where(sel["A"])[0]
    if len(iA) == 0:                       # a small test scan may not reach any domain interior
        print("  no scan positions inside domain A's interior -- picking single patterns from "
              "the whole scan instead")
        iA = np.arange(len(pos))
    cA = pos[iA].mean(0)
    base = np.floor(cA / a) * a
    targets = {"on Pb column": base, "on Ti-O column": base + a / 2,
               "on O_eq column": base + np.array([a / 2, 0.0])}
    pick = {k: int(np.argmin(np.hypot(*(pos - v).T))) for k, v in targets.items()}

    dp_path = os.path.join(args.run, "01", "data_dp.hdf5")
    with h5py.File(dp_path, "r") as f:
        dp = f["dp"]
        n_pos, nb = dp.shape[0], dp.shape[1]
        tot = np.zeros((nb, nb))
        dom = {n: np.zeros((nb, nb)) for n in names}
        cnt = {n: 0 for n in names}
        for s in range(0, n_pos, 500):
            blk = np.asarray(dp[s:s + 500], dtype=np.float64).transpose(0, 2, 1)
            tot += blk.sum(0)
            for n in names:
                w = sel[n][s:s + 500]
                if w.any():
                    dom[n] += blk[w].sum(0)
                    cnt[n] += int(w.sum())
        singles = {k: np.asarray(dp[i], dtype=np.float64).T for k, i in pick.items()}

    theta = SF.detector_axes(nb, float(b["d_alpha_mrad"]))
    out = dict(mean=tot / n_pos, theta_mrad=theta, d_alpha_mrad=float(b["d_alpha_mrad"]),
               convergence_mrad=float(b["convergence_mrad"]),
               scan_step_A=float(b["scan_step_A"]), n_pos=n_pos,
               single_labels=np.array(list(singles)),
               singles=np.stack(list(singles.values())),
               single_xy=np.stack([pos[i] for i in pick.values()]))
    for n in names:
        out[f"mean_{n}"] = dom[n] / max(cnt[n], 1)
        out[f"count_{n}"] = cnt[n]
    np.savez_compressed(args.out, **out)
    print(f"wrote {args.out}  ({os.path.getsize(args.out)/1e6:.1f} MB) from {n_pos} patterns, "
          f"{nb}x{nb} px; domain interiors " + " ".join(f"{n}={cnt[n]}" for n in names))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
