#!/usr/bin/env python
"""@file sign_test.py
@brief Read the sign of P_z straight from the hollow 4D-STEM data, without a depth reconstruction.

The depth-layer reconstruction turned out to be a lossy route to the sign: it converged to a
smeared object that fits the data ~6.8% worse than the truth (`depth_constraint.py`), so the
information was there and the optimiser did not use it. This takes the direct route instead.

The fusion makes it a TWO-hypothesis problem, which is what makes this cheap. EELS supplies
|delta_z| and projected ptychography supplies delta_xy, so each domain is known up to a single
binary choice -- the stacking order along the beam. Simulating those two candidates and asking
which the measured patterns match is then a likelihood ratio, not a reconstruction.

Two things make it affordable. Each domain is internally uniform and periodic, so its diffraction
pattern depends only on the probe's position WITHIN the unit cell; and the candidates need
simulating only at the intra-cell positions the measurement actually visited. That is a few
hundred positions rather than the 6400 of the full scan.

Only detector pixels OUTSIDE the hole are compared, so this is a genuine hollow-detector
measurement: the electrons inside the hole have gone to the spectrometer.

GPU only, and only through SLURM:

    bash fusion/run_gpu.sh sign_test --data fusion/runs/fusion
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sign_encoding as SE      # noqa: E402
import simulate_fusion as SF    # noqa: E402
import toy_sample as T          # noqa: E402


def load_measured(data_dir: str, scan: str = "01"):
    """@brief The simulated experiment: flux-normalised patterns and their probe positions."""
    import h5py
    d = os.path.join(data_dir, scan)
    with h5py.File(os.path.join(d, "data_position.hdf5"), "r") as f:
        pos = np.asarray(f["probe_positions_0"][...], dtype=float).T      # (Npos, 2) x, y
    budget = json.load(open(os.path.join(data_dir, "hollow_budget.json")))
    return pos, os.path.join(d, "data_dp.hdf5"), budget


def interior(pos, truth, budget, margin_A: float):
    """@brief Positions whose probe stays inside ONE domain.

    A uniform candidate cannot describe a column near a wall, and the probe is not a point: at
    100 mrad with 20 A of overfocus it illuminates a ~4 A disc and spreads further through the
    slab, so the margin has to cover that, not just the cell.
    """
    a, n_lat = float(truth["a"]), int(truth["n_lat"])
    mid = float(truth["box_A"]) / 2.0
    x, y = pos[:, 0], pos[:, 1]
    sel = {}
    for name, (qx, qy) in zip([str(n) for n in truth["names"]], truth["quads"]):
        mx = (x > margin_A) & (x < mid - margin_A) if qx == 0 else (x > mid + margin_A) & (x < 2 * mid - margin_A)
        my = (y > margin_A) & (y < mid - margin_A) if qy == 0 else (y > mid + margin_A) & (y < 2 * mid - margin_A)
        sel[name] = mx & my
    return sel


def read_patterns(dp_path: str, idx) -> np.ndarray:
    """@brief Flux-normalised measured patterns at the given scan indices, in the SIM convention.

    `simulate_4dstem.save_outputs` writes `dp = A.transpose(0, 2, 1)` so that MATLAB, which
    reverses axes on read, sees dp(dy, dx, k). Reading it back in Python therefore gives the
    TRANSPOSE of the simulated pattern, and it has to be undone before comparing against freshly
    simulated candidates -- the detector mask is radially symmetric and would hide the error, so a
    silent transpose here would just look like a 50% success rate.
    """
    import h5py
    idx = np.asarray(sorted(idx))
    with h5py.File(dp_path, "r") as f:
        dp = f["dp"]
        out = np.empty((len(idx), dp.shape[1], dp.shape[2]), dtype=np.float64)
        for k in range(0, len(idx), 256):                    # h5py wants increasing selections
            j = idx[k:k + 256]
            out[k:k + len(j)] = np.asarray(dp[j], dtype=np.float64)
    out = out.transpose(0, 2, 1)                             # back to the simulation's convention
    return out / out.reshape(len(out), -1).sum(1)[:, None, None]


def log_likelihood_ratio(meas, p_a, p_b, keep, eps: float = 1e-14) -> np.ndarray:
    """@brief Per-pattern Poisson log-likelihood ratio for hypothesis A over B.

    Positive favours A. Only pixels in `keep` (outside the hole) enter, so the statistic uses
    exactly the electrons hollow ptychography actually records.
    """
    m = meas[:, keep]
    a = np.maximum(p_a[:, keep], eps)
    b = np.maximum(p_b[:, keep], eps)
    a = a / a.sum(1, keepdims=True)
    b = b / b.sum(1, keepdims=True)
    m = m / m.sum(1, keepdims=True)
    return (m * (np.log(a) - np.log(b))).sum(1)


def run(data_dir: str, truth_npz: str, n_pos: int = 200, margin_A: float = 6.0,
        hsa_frac: float = 0.75, bin_factor: int = 4, slice_thickness_A: float = 1.0,
        device: str = "gpu", seed: int = 0, out: str | None = None) -> dict:
    """@brief Decide sign(delta_z) per domain from the measured patterns and two candidates each."""
    SF.require_gpu_job(device)
    s4 = SF._import_sim4d()
    truth = np.load(truth_npz, allow_pickle=True)
    pos, dp_path, budget = load_measured(data_dir)
    alpha = float(budget["convergence_mrad"])
    hsa = hsa_frac * alpha
    n_lat, n_z = int(truth["n_lat"]), int(truth["n_z"])
    names = [str(n) for n in truth["names"]]
    rng = np.random.default_rng(seed)

    sel = interior(pos, truth, budget, margin_A)
    print(f"sign test: {len(pos)} measured positions, hole {hsa:.0f} mrad "
          f"(only pixels OUTSIDE it are used)")
    print(f"domain interiors at {margin_A:.0f} A margin: "
          + "  ".join(f"{n} {int(sel[n].sum())}" for n in names) + "\n")

    results = {}
    for k, name in enumerate(names):
        idx_all = np.where(sel[name])[0]
        if len(idx_all) < 8:
            print(f"  {name}: too few interior positions, skipped")
            continue
        idx = np.sort(rng.choice(idx_all, size=min(n_pos, len(idx_all)), replace=False))
        p_meas = read_patterns(dp_path, idx)
        n_b = p_meas.shape[-1]
        theta = SF.detector_axes(n_b, float(budget["d_alpha_mrad"]))
        keep = theta >= hsa                                   # what ptychography records

        # the two candidates: EELS fixes |delta_z|, projection fixes delta_xy, sign is the unknown
        d = np.asarray(truth["deltas"][k], float)
        cand = {}
        for tag, dz in (("up", +abs(d[2])), ("down", -abs(d[2]))):
            atoms, _ = T.build(n_lat, n_z,
                               domains=T.uniform_domains([d[0], d[1], dz]))
            cand[tag] = SE.patterns_from_atoms(atoms, pos[idx], alpha, bin_factor,
                                               slice_thickness_A, device)
        llr = log_likelihood_ratio(p_meas, cand["up"], cand["down"], keep)
        want = "up" if d[2] > 0 else "down"
        got = "up" if llr.sum() > 0 else "down"
        frac = float(np.mean((llr > 0) == (d[2] > 0)))
        results[name] = dict(n=int(len(idx)), llr_total=float(llr.sum()),
                             llr_per_pattern=float(llr.mean()),
                             per_pattern_correct=frac, truth=want, decided=got,
                             correct=bool(got == want))
        print(f"  {name}: truth {want:>4}, decided {got:>4}  "
              f"{'OK ' if got == want else 'WRONG'}  |  LLR {llr.sum():+.4g} "
              f"({llr.mean():+.3g}/pattern), {frac*100:.0f}% of individual patterns correct")

    ok = sum(r["correct"] for r in results.values())
    print(f"\n  {ok}/{len(results)} domains correct from the hollow data alone, "
          f"no depth reconstruction")
    res = {"domains": results, "hsa_mrad": hsa, "n_pos": n_pos, "margin_A": margin_A}
    if out:
        json.dump(res, open(out, "w"), indent=2)
        print(f"[sign-test] wrote {out}")
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", default=os.path.join(HERE, "runs", "fusion"),
                    help="the simulated experiment (contains 01/ and hollow_budget.json)")
    ap.add_argument("--truth", default=os.path.join(HERE, "sample", "toy_truth.npz"))
    ap.add_argument("--n-pos", type=int, default=200, help="measured positions per domain")
    ap.add_argument("--margin", type=float, default=6.0, help="clearance from a domain wall [A]")
    ap.add_argument("--hsa", type=float, default=0.75)
    ap.add_argument("--bin-factor", type=int, default=4)
    ap.add_argument("--slice-thickness", type=float, default=1.0)
    ap.add_argument("--device", default="gpu", choices=["gpu", "cpu"])
    ap.add_argument("--out", default=os.path.join(HERE, "runs", "sign_test.json"))
    args = ap.parse_args(argv)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    run(args.data, args.truth, args.n_pos, args.margin, args.hsa, args.bin_factor,
        args.slice_thickness, args.device, out=args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
