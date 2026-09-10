#!/usr/bin/env python
"""@file check_recon.py
@brief Did the reconstruction actually recover DEPTH? One command, before any physics is read off.

The first four-domain run resolved the lattice beautifully in projection and yet carried no depth
structure at all, which invalidated everything downstream. That failure is invisible in the engine
log, so it needs its own check. This reports, for a PtychoShelves multislice reconstruction:

  1. convergence   -- the data-fit trace, sampled at the iterations where it was actually recorded.
                      NOTE the engine's `fourier_error_out` scores EVERY detector pixel, while the
                      hollow mask is applied only inside modulus_constraint.m, so for a hollow run
                      the metric counts pixels the solver is told to ignore and RISES as the model
                      drifts there. It is a valid monitor only at HSA = 0; treat it as bookkeeping
                      otherwise.
  2. lattice       -- is the in-plane structure there at all, and can Pb be told from Ti?
  3. DEPTH         -- power at the lattice period along the beam, after detrending, relative to the
                      mean band power. This is the number that decides whether a depth-resolved
                      readout is possible: ~1 means no atomic planes along z, however good the
                      projection looks.

    ~/hyperspy-bundle/bin/python check_recon.py --recon <Niter*.mat> [--budget hollow_budget.json]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyze_fusion as AF     # noqa: E402


def error_trace(path: str):
    """@brief The data-fit trace at the iterations where it was recorded (the rest are unfilled)."""
    import h5py
    with h5py.File(path, "r") as f:
        if "fourier_error_out" not in f:
            return None, None, {}
        e = np.asarray(f["fourier_error_out"]).ravel().astype(float)
        par = {k: float(np.asarray(f["par"][k]).ravel()[0])
               for k in ("Nlayers", "regularize_layers", "beta_LSQ", "probe_modes",
                         "number_iterations")
               if k in f.get("par", {})}
    it = np.where(np.isfinite(e) & (e > 0))[0]
    return it + 1, e[it], par


def depth_power(V, budget, truth):
    """@brief Power at the lattice period along the beam vs the mean band power, per column type.

    A resolved atomic stacking gives >> 1. About 1 means the depth axis holds no lattice at all,
    which is the failure the projection cannot show you.
    """
    c = float(truth["c"])
    cells, xy, z = AF.column_maps(V, budget, truth)
    dz = z[1] - z[0]
    out = {}
    for kind in AF.KINDS:
        P = np.array([cells[i][kind] for i in range(len(cells))])
        t = np.polynomial.polynomial.polyvander(z, 3)          # remove a smooth envelope
        R = P - (t @ np.linalg.lstsq(t, P.T, rcond=None)[0]).T
        f = np.fft.rfftfreq(len(z), d=dz)
        S = np.abs(np.fft.rfft(R, axis=1)) ** 2
        k = int(np.argmin(np.abs(f - 1.0 / c)))
        band = (f > 0.05) & (f < 0.5 / dz)
        out[kind] = float(S[:, k].mean() / S[:, band].mean())
    return out, len(cells), dz


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--recon", required=True)
    ap.add_argument("--budget", default=None)
    ap.add_argument("--truth", default=os.path.join(HERE, "sample", "toy_truth.npz"))
    args = ap.parse_args(argv)

    truth = np.load(args.truth, allow_pickle=True)
    cand = [args.budget] if args.budget else []
    cand += [os.path.join(HERE, "runs", d, "hollow_budget.json") for d in ("fusion", "probe_test")]
    bpath = next((p for p in cand if p and os.path.exists(p)), None)
    if bpath is None:
        raise SystemExit("no hollow_budget.json found (pass --budget)")
    budget = json.load(open(bpath))

    print(f"reconstruction : {args.recon}")
    print(f"budget         : {os.path.relpath(bpath, HERE)}")

    it, err, par = error_trace(args.recon)
    if it is not None and len(it):
        print(f"\n1. convergence   {', '.join(f'{k} {v:g}' for k, v in par.items())}")
        step = max(1, len(it) // 8)
        print("   " + "  ".join(f"{i}:{e:.1f}" for i, e in zip(it[::step], err[::step])))
        drop = 100 * (err[0] - err[-1]) / err[0]
        rising = err[-1] > err[len(err) // 2]
        print(f"   {err[0]:.1f} -> {err[-1]:.1f}  ({drop:+.1f}%);"
              f" {'RISING in the second half' if rising else 'monotone-ish decrease'}")
        print("   (scores all pixels; for a hollow run it counts the ignored ones -- bookkeeping only)")

    V = AF.load_volume(args.recon)
    print(f"\n2. lattice       object {V.shape}")
    try:
        res, n_cells, dz = depth_power(V, budget, truth)
    except SystemExit as e:
        print(f"   FAILED: {e}")
        return 1

    print(f"\n3. DEPTH         {n_cells} cells, {V.shape[0]} layers x {dz:.2f} A")
    print(f"   power at the {float(truth['c']):.3f} A lattice period / mean band power:")
    for k, v in res.items():
        print(f"     {k:>4}: {v:6.2f}  {'#' * int(min(v, 20))}")
    best = max(res.values())
    print(f"\n   VERDICT: " + ("depth structure IS resolved -- a depth-resolved readout is possible"
                              if best > 3 else
                              "NO depth structure (best %.2f, need >>1). The projection may still be"
                              " excellent; a depth-resolved readout is not possible from this run."
                              % best))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
