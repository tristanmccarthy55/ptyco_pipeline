#!/usr/bin/env python
"""@file sign_encoding.py
@brief Is the SIGN of P_z encoded in the diffraction data at all, and how thick must the specimen be?

The first four-domain reconstruction resolved the lattice beautifully in projection but returned no
depth structure whatsoever, so the sign readout had nothing to work with. Before tuning any
reconstruction, this asks the prior question: does the 4D-STEM data contain the sign, and if so at
what thickness and at what dose.

Projected potential is IDENTICAL for an up and a down cell, atom for atom, so any difference is
purely dynamical -- it grows with multiple scattering, hence with thickness. That is the quantity
measured here.

Method. Two UNIFORM membranes differing only in sign(delta_z), scanned at IDENTICAL probe
positions covering one unit cell. Comparing uniform specimens at matched positions is what makes
this valid: a domain-to-domain comparison inside one scan is dominated by intra-cell probe-position
sampling, which is far larger than the structural difference and swamps it.

Metric. For flux-normalised patterns p_up, p_dn the Poisson discriminability per electron is
    D = sum_k (p_up - p_dn)^2 / p_avg
so M positions at N electrons each give SNR = sqrt(M N D), and the dose needed for SNR 3 is
N = 9/(M D). D is also split by scattering angle, which answers the question the hollow geometry
raises: if the sign lives at low angles it goes down the hole to the spectrometer and hollow
ptychography is blind to it however good the reconstruction.

GPU only, and only through SLURM -- the login node has no CUDA driver:

    bash fusion/run_sign_encoding.sh
    THICKNESS="5 10 20 40 60" N_SCAN=12 bash fusion/run_sign_encoding.sh
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import simulate_fusion as SF     # noqa: E402
import toy_sample as T           # noqa: E402


def scan_positions(a: float, n: int, origin: float) -> np.ndarray:
    """@brief An n x n grid over ONE unit cell -- uniform intra-cell sampling, identical for both."""
    u = (np.arange(n) + 0.5) / n * a + origin
    X, Y = np.meshgrid(u, u, indexing="ij")
    return np.stack([X.ravel(), Y.ravel()], axis=1)


def patterns_from_atoms(atoms, pos, convergence_mrad: float, bin_factor: int,
                        slice_thickness_A: float, device: str):
    """@brief Flux-normalised diffraction patterns of ANY structure at the given probe positions.

    The scattering is sim/simulate_4dstem.py exactly as the production pipeline uses it; only the
    structure and the positions change. Flux normalisation makes two structures directly
    comparable pattern by pattern.
    """
    import abtem
    s4 = SF._import_sim4d()
    s4.DEVICE = device
    s4.BIN_FACTOR = bin_factor
    s4.SLICE_THICKNESS_A = slice_thickness_A
    s4.CONVERGENCE_MRAD = convergence_mrad
    pot = s4.build_potential(atoms)
    probe = s4.build_probe(pot)
    arr = np.asarray(s4.run_scan_binned(probe, atoms, abtem.CustomScan(pos)), dtype=np.float64)
    return arr / arr.reshape(len(arr), -1).sum(1)[:, None, None]


def patterns(delta, n_lat: int, n_z: int, pos, convergence_mrad: float, bin_factor: int,
             slice_thickness_A: float, device: str):
    """@brief Flux-normalised diffraction patterns of one uniform membrane at the given positions."""
    atoms, truth = T.build(n_lat, n_z, domains=T.uniform_domains(delta))
    return patterns_from_atoms(atoms, pos, convergence_mrad, bin_factor,
                               slice_thickness_A, device), truth


def assert_distinct(p_a, p_b, what: str) -> None:
    """@brief Two DIFFERENT structures must not produce bit-identical patterns.

    One GPU sweep returned exactly 0.0 for the 3-cell up/down pair while the CPU gives ~3% at the
    identical settings and the two potentials differ by ~96%. Physics cannot produce an exact zero,
    so an identical pair means the computation, not the specimen, and it must stop the run rather
    than enter the results as "no signal".
    """
    if np.array_equal(p_a, p_b):
        raise SystemExit(f"{what}: the two structures gave BIT-IDENTICAL patterns although they "
                         f"differ -- a computation defect, not a physical null. Re-run this "
                         f"thickness, or cross-check with --device cpu.")


def discriminability(p_up, p_dn, theta_mrad, hsa_mrad, eps: float = 1e-14):
    """@brief Poisson discriminability per electron, total and split by the hollow hole.

    @return dict with D_total, D_outside (kept for ptychography) and D_inside (sent to EELS).
    """
    d2 = (p_up - p_dn) ** 2 / np.maximum(0.5 * (p_up + p_dn), eps)
    out = theta_mrad >= hsa_mrad
    ad = np.abs(p_up - p_dn)
    tot = 0.5 * (p_up.sum(axis=(1, 2)) + p_dn.sum(axis=(1, 2)))
    return {"D_total": float(d2.sum(axis=(1, 2)).mean()),
            "D_outside": float(d2[:, out].sum(axis=1).mean()),
            "D_inside": float(d2[:, ~out].sum(axis=1).mean()),
            # NOISELESS separability: the fraction of the pattern that differs at all. This is the
            # quantity that matters before dose enters -- a reconstruction has to be sensitive to
            # THIS, whatever the electron count.
            "frac_total": float((ad.sum(axis=(1, 2)) / tot).mean()),
            "frac_outside": float((ad[:, out].sum(axis=1) / tot).mean())}


def run(thicknesses, n_lat: int = 4, n_scan: int = 8, convergence_mrad: float = 100.0,
        hsa_frac: float = 0.75, bin_factor: int = 4, slice_thickness_A: float = 1.0,
        theta_deg: float = 20.0, device: str = "gpu", out: str | None = None) -> dict:
    """@brief Sweep thickness and report the dose the sign costs, in total and outside the hole."""
    SF.require_gpu_job(device)
    s4 = SF._import_sim4d()
    w, d_ti = T.polar_mode()
    t = np.radians(theta_deg)
    delta = abs(d_ti) * np.array([np.sin(t), 0.0, np.cos(t)])
    a = float(T.REF.a_tet)
    hsa = hsa_frac * convergence_mrad
    rows = []
    print(f"sign encoding: |delta| {abs(d_ti):.3f} A at {theta_deg:.0f} deg from the beam "
          f"(delta_z = {delta[2]:+.3f} A)")
    print(f"alpha {convergence_mrad:.0f} mrad, hollow semi-angle {hsa:.0f} mrad, "
          f"{n_scan}x{n_scan} positions over one {a:.3f} A cell")
    print("NOISELESS throughout: exact patterns are compared, no shot noise anywhere. The dose "
          "line is\na derived aside, not part of the measurement.\n")
    for n_z in thicknesses:
        pos = scan_positions(a, n_scan, 0.0)
        flip = np.array([1.0, 1.0, -1.0])          # flip delta_z ONLY: same in-plane, same |delta|
        p_up, truth = patterns(delta, n_lat, n_z, pos, convergence_mrad, bin_factor,
                               slice_thickness_A, device)
        p_dn, _ = patterns(delta * flip, n_lat, n_z, pos, convergence_mrad, bin_factor,
                           slice_thickness_A, device)
        assert_distinct(p_up, p_dn, f"sign encoding at {n_z} cells")
        n_b = p_up.shape[-1]
        lam = s4.wavelength_a()
        box = float(T.REF.a_tet) * n_lat
        theta = SF.detector_axes(n_b, (bin_factor / box) * lam * 1e3)
        D = discriminability(p_up, p_dn, theta, hsa)
        M = len(pos)
        thick_A = n_z * float(T.REF.c_tet)
        row = dict(n_z=int(n_z), thickness_A=thick_A, M=M, **D)
        for tag, key in (("total", "D_total"), ("outside", "D_outside")):
            N = 9.0 / (M * D[key]) if D[key] > 0 else np.inf     # electrons/pattern for SNR 3
            row[f"dose_{tag}_e_per_A2"] = N / (a / n_scan) ** 2
        rows.append(row)
        print(f"  {n_z:>3} cells ({thick_A:5.1f} A):  NOISELESS pattern difference "
              f"{100*D['frac_total']:6.3f}%  (outside the hole {100*D['frac_outside']:6.3f}%, "
              f"{100*D['D_outside']/max(D['D_total'],1e-30):4.1f}% of the information)")
        print(f"                        [if shot noise were added: SNR 3 at "
              f"{row['dose_outside_e_per_A2']:.1e} e/A^2]")
    # NULL CONTROL: with delta purely in-plane, flipping delta_z changes nothing, so the two
    # structures are literally identical and D must come out at the numerical floor. Anything
    # else would mean the measured D is an artefact of the build/centring, not the polarisation.
    n_z0 = int(thicknesses[0])
    d_in = abs(d_ti) * np.array([1.0, 0.0, 0.0])
    pos = scan_positions(a, n_scan, 0.0)
    q_up, _ = patterns(d_in, n_lat, n_z0, pos, convergence_mrad, bin_factor,
                       slice_thickness_A, device)
    q_dn, _ = patterns(d_in * np.array([1.0, 1.0, -1.0]), n_lat, n_z0, pos, convergence_mrad,
                       bin_factor, slice_thickness_A, device)
    n_b = q_up.shape[-1]
    box = float(T.REF.a_tet) * n_lat
    theta_ax = SF.detector_axes(n_b, (bin_factor / box) * s4.wavelength_a() * 1e3)
    null = discriminability(q_up, q_dn, theta_ax, hsa)
    ratio = rows[0]["D_total"] / max(null["D_total"], 1e-30)
    print(f"\n  NULL (delta in-plane, so the flip is a no-op): D_total {null['D_total']:.3e}")
    print(f"  signal / null at {n_z0} cells = {ratio:.3g}  -> "
          + ("the measured D is the polarisation" if ratio > 100 else
             "*** D is contaminated by a build/centring artefact -- do not trust it ***"))

    res = {"thicknesses": rows, "null": null, "signal_over_null": float(ratio),
           "convergence_mrad": convergence_mrad, "hsa_mrad": hsa,
           "theta_deg": theta_deg, "delta": delta.tolist(), "n_lat": n_lat, "n_scan": n_scan}
    if out:
        json.dump(res, open(out, "w"), indent=2)
        print(f"\n[sign] wrote {out}")
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--thickness", type=int, nargs="+", default=[5, 10, 20, 40],
                    help="slab thicknesses in unit cells")
    ap.add_argument("--n-lat", type=int, default=4, help="in-plane cells (periodic; 4 is plenty)")
    ap.add_argument("--n-scan", type=int, default=8, help="probe positions per cell edge")
    ap.add_argument("--convergence", type=float, default=100.0)
    ap.add_argument("--hsa", type=float, default=0.75, help="hollow semi-angle / alpha")
    ap.add_argument("--theta", type=float, default=20.0, help="polar tilt from the beam [deg]")
    ap.add_argument("--bin-factor", type=int, default=4)
    ap.add_argument("--slice-thickness", type=float, default=1.0)
    ap.add_argument("--device", default="gpu", choices=["gpu", "cpu"])
    ap.add_argument("--out", default=os.path.join(HERE, "runs", "sign_encoding.json"))
    args = ap.parse_args(argv)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    run(args.thickness, args.n_lat, args.n_scan, args.convergence, args.hsa,
        args.bin_factor, args.slice_thickness, args.theta, args.device, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
