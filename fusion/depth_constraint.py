#!/usr/bin/env python
"""@file depth_constraint.py
@brief Do the diffraction data constrain the DEPTH distribution at all, or only the projection?

The reconstruction resolved the lattice in projection and returned a smooth entrance-weighted
envelope along the beam, with no atomic planes (`check_recon.py`). Two explanations are possible
and they demand opposite responses:

  (a) the optimiser failed -- the data DO distinguish the true stacking from a smeared one, and a
      better-conditioned reconstruction would find it; or
  (b) the data are indifferent -- a smeared solution fits as well as the true one, so no
      reconstruction can recover the depth and the specimen or the optics must change.

This decides between them without running a reconstruction at all. Two structures are compared:

  TRUE    the membrane as built;
  SMEARED every column keeps its atoms, their species and their (x, y) -- only their z is
          replaced by a uniform spread through the slab.

The SMEARED structure therefore has an IDENTICAL projected potential and no atomic planes along
the beam: it is the reconstruction's answer in caricature. If its diffraction patterns differ
strongly from the true ones, the data know the difference and (a) holds. If they barely differ,
(b) holds. The difference is reported in the same units as `sign_encoding.py`, so the two can be
read side by side against the 4.9 % sign signal.

GPU only, and only through SLURM (the login node has no CUDA driver):

    bash fusion/run_depth_constraint.sh
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


def smear_z(atoms, seed: int = 0):
    """@brief Same atoms, same species, same (x, y) -- z spread uniformly through the slab.

    Preserving the count per (x, y, species) is what keeps the PROJECTED potential identical, so
    the only thing that changes is the depth distribution. Column membership is taken from the
    in-plane position rounded to 1e-3 A, which is exact here because the toy lattice is built on
    a grid.
    """
    out = atoms.copy()
    z = out.get_positions()[:, 2]
    z0, z1 = z.min(), z.max()
    pos = out.get_positions()
    key = [(round(x, 3), round(y, 3), s)
           for (x, y), s in zip(pos[:, :2], out.get_chemical_symbols())]
    groups: dict = {}
    for i, k in enumerate(key):
        groups.setdefault(k, []).append(i)
    rng = np.random.default_rng(seed)
    for k, idx in groups.items():
        n = len(idx)
        znew = z0 + (np.arange(n) + 0.5) / n * (z1 - z0)      # uniform, no planes
        pos[np.array(idx), 2] = znew
    out.set_positions(pos)
    return out


def check_projection_identical(a, b, tol: float = 1e-9) -> float:
    """@brief The two structures must project identically, or the comparison means nothing."""
    ka = sorted((round(x, 3), round(y, 3), s)
                for (x, y), s in zip(a.get_positions()[:, :2], a.get_chemical_symbols()))
    kb = sorted((round(x, 3), round(y, 3), s)
                for (x, y), s in zip(b.get_positions()[:, :2], b.get_chemical_symbols()))
    if ka != kb:
        raise SystemExit("the smeared structure does not project identically -- comparison invalid")
    return float(np.abs(np.sort(a.get_positions()[:, 2]) - np.sort(b.get_positions()[:, 2])).mean())


def run(thicknesses, n_lat: int = 4, n_scan: int = 8, convergence_mrad: float = 100.0,
        hsa_frac: float = 0.75, bin_factor: int = 4, slice_thickness_A: float = 1.0,
        theta_deg: float = 20.0, device: str = "gpu", out: str | None = None) -> dict:
    """@brief Sweep thickness: how strongly do the data distinguish true stacking from smeared?"""
    SF.require_gpu_job(device)
    s4 = SF._import_sim4d()
    w, d_ti = T.polar_mode()
    t = np.radians(theta_deg)
    delta = abs(d_ti) * np.array([np.sin(t), 0.0, np.cos(t)])
    a = float(T.REF.a_tet)
    hsa = hsa_frac * convergence_mrad
    pos = SE.scan_positions(a, n_scan, 0.0)
    rows = []
    print("depth constraint: TRUE stacking vs SMEARED (identical projection, no planes along z)")
    print(f"alpha {convergence_mrad:.0f} mrad, hole {hsa:.0f} mrad, {n_scan}x{n_scan} positions "
          f"over one {a:.3f} A cell.  NOISELESS.\n")
    for n_z in thicknesses:
        atoms, truth = T.build(n_lat, n_z, domains=T.uniform_domains(delta))
        sm = smear_z(atoms)
        moved = check_projection_identical(atoms, sm)
        p_true = SE.patterns_from_atoms(atoms, pos, convergence_mrad, bin_factor,
                                        slice_thickness_A, device)
        p_smear = SE.patterns_from_atoms(sm, pos, convergence_mrad, bin_factor,
                                         slice_thickness_A, device)
        n_b = p_true.shape[-1]
        box = a * n_lat
        theta = SF.detector_axes(n_b, (bin_factor / box) * s4.wavelength_a() * 1e3)
        D = SE.discriminability(p_true, p_smear, theta, hsa)
        row = dict(n_z=int(n_z), thickness_A=n_z * float(T.REF.c_tet),
                   mean_z_move_A=moved, **D)
        rows.append(row)
        print(f"  {n_z:>3} cells ({row['thickness_A']:5.1f} A): depth signal "
              f"{100*D['frac_total']:6.3f}%  (outside the hole {100*D['frac_outside']:6.3f}%)"
              f"   [atoms moved {moved:.2f} A in z on average]")
    res = {"thicknesses": rows, "convergence_mrad": convergence_mrad, "hsa_mrad": hsa,
           "theta_deg": theta_deg, "n_lat": n_lat, "n_scan": n_scan}
    best = max(r["frac_total"] for r in rows)
    print(f"\n  VERDICT: " + (
        "the data DO constrain depth (%.2f%% at best) -- the smeared solution is a WORSE fit, so "
        "the reconstruction failed to find the true one and better conditioning should help."
        % (100 * best) if best > 0.005 else
        "the data barely distinguish them (%.3f%% at best) -- a smeared solution fits about as "
        "well, so no reconstruction can recover the depth at this thickness/optics."
        % (100 * best)))
    if out:
        json.dump(res, open(out, "w"), indent=2)
        print(f"[depth] wrote {out}")
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--thickness", type=int, nargs="+", default=[3, 5, 10, 20])
    ap.add_argument("--n-lat", type=int, default=4)
    ap.add_argument("--n-scan", type=int, default=8)
    ap.add_argument("--convergence", type=float, default=100.0)
    ap.add_argument("--hsa", type=float, default=0.75)
    ap.add_argument("--theta", type=float, default=20.0)
    ap.add_argument("--bin-factor", type=int, default=4)
    ap.add_argument("--slice-thickness", type=float, default=1.0)
    ap.add_argument("--device", default="gpu", choices=["gpu", "cpu"])
    ap.add_argument("--out", default=os.path.join(HERE, "runs", "depth_constraint.json"))
    args = ap.parse_args(argv)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    run(args.thickness, args.n_lat, args.n_scan, args.convergence, args.hsa,
        args.bin_factor, args.slice_thickness, args.theta, args.device, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
