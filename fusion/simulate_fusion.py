#!/usr/bin/env python
"""@file simulate_fusion.py
@brief One 4D-STEM simulation of the toy membrane, split by a hollow detector into ptychography
       + EELS + simultaneous HAADF.

The experiment of Song et al. / Lei & Wang (arXiv:2506.22352) as a forward simulation: a single
pixelated detector with a HOLE in it. Electrons inside the hollow semi-angle (HSA) pass through to
the spectrometer; everything outside is recorded for multislice hollow ptychography (MHP); the
outer annulus is a simultaneous virtual HAADF. Nothing is duplicated -- the scattering is
`sim/simulate_4dstem.py` exactly as the validated PTO/STO pipeline uses it (same probe, same
Lobato potential, same binning, same PtychoShelves output contract), driven here with the toy
membrane and its own scan. That module is imported, never modified: the labyrinth's scan region
and rotation stay untouched.

The hole is applied as a MASK, not by deleting data: PtychoShelves' engine takes a global `mask1`
(1 = use this pixel, 0 = ignore) which `+engines/+GPU/private/modulus_constraint.m` turns into the
excluded set. So one simulation serves every hollow semi-angle -- only the mask changes, which is
what makes the HSA sweep cheap.

Outputs (into <out-dir>/):
  01/data_dp.hdf5, data_position.hdf5, probe_initial.mat, sim_meta.mat   PtychoShelves inputs
  01/mask_hsa<f>.mat                                                     hollow masks, one per HSA
  hollow_budget.json                                                     measured dose split + HAADF
  haadf.npy                                                              simultaneous HAADF image

    <blythe abtem env>/python simulate_fusion.py --test          # tiny local shape check
    <blythe abtem env>/python simulate_fusion.py --out-dir ...   # production (GPU)
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_HSA = (0.0, 0.50, 0.75, 0.95)      # hollow semi-angle as a fraction of the convergence


def require_gpu_job(device: str) -> None:
    """@brief Refuse to start a GPU run outside SLURM: the login node has no CUDA driver.

    abtem on the login node fails deep inside cupy with cudaErrorInsufficientDriver, several
    frames from anything recognisable, after the potential has already been built. Fail here
    instead, with the command to use.
    """
    if device == "gpu" and not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit(
            "refusing to run a GPU job outside SLURM -- the Blythe login node has no CUDA "
            "driver.\n  submit it:   bash fusion/run_sign_encoding.sh"
            "   (or fusion/run_fusion_sim.slurm)\n  or force CPU: --device cpu   (small tests only)")


def _import_sim4d():
    """@brief Import sim/simulate_4dstem.py so ALL the scattering is the validated code."""
    sim_dir = os.path.abspath(os.path.join(HERE, "..", "sim"))
    if sim_dir not in sys.path:
        sys.path.insert(0, sim_dir)
    import simulate_4dstem as s4          # noqa: E402
    return s4


# ---------------------------------------------------------------- the hollow detector
def detector_axes(n_b: int, d_alpha_mrad: float) -> np.ndarray:
    """@brief Radial scattering angle [mrad] of every binned detector pixel.

    The binning in `simulate_4dstem._crop_center_to_multiple` keeps DC centred, so the zero-order
    beam sits at index n_b/2 -- the same pixel PtychoShelves is told about via
    `p.ctr = [fix(Ndpx/2)+1, ...]`. Everything downstream (masks, HAADF, dose split) hangs off this.
    """
    i = np.arange(n_b) - n_b / 2.0
    return np.hypot(*np.meshgrid(i, i, indexing="ij")) * d_alpha_mrad


def hollow_mask(theta_mrad: np.ndarray, hsa_mrad: float) -> np.ndarray:
    """@brief PtychoShelves `mask1`: 1 = keep for ptychography, 0 = inside the hole (-> EELS)."""
    return (theta_mrad >= hsa_mrad).astype(np.float32)


def dose_budget(dp, theta_mrad: np.ndarray, hsa_list_mrad, haadf_inner_mrad: float,
                det_max_mrad: float) -> dict:
    """@brief Where the electrons actually go, measured on the simulated patterns.

    Not estimated from (HSA/alpha)^2: the scattered intensity is redistributed by the specimen, so
    the split is read off the diffraction patterns themselves, averaged over the scan.
    """
    mean_dp = np.asarray(dp, dtype=np.float64).reshape(-1, *theta_mrad.shape).mean(0)
    tot = mean_dp.sum()
    out = {"total": 1.0, "hsa": {}}
    for hsa in hsa_list_mrad:
        f = float(mean_dp[theta_mrad < hsa].sum() / tot)
        out["hsa"][f"{hsa:.1f}"] = {"eels_fraction": f, "ptycho_fraction": 1.0 - f}
    out["haadf_fraction"] = float(
        mean_dp[(theta_mrad >= haadf_inner_mrad) & (theta_mrad <= det_max_mrad)].sum() / tot)
    return out


def virtual_image(dp, theta_mrad: np.ndarray, inner: float, outer: float, shape) -> np.ndarray:
    """@brief Virtual annular-detector image from the pixelated data (BF, ABF or HAADF)."""
    m = (theta_mrad >= inner) & (theta_mrad <= outer)
    arr = np.asarray(dp, dtype=np.float64).reshape(-1, *theta_mrad.shape)
    return (arr * m).sum(axis=(1, 2)).reshape(shape)


# ---------------------------------------------------------------- the run
def run(out_dir: str, vasp: str, truth_npz: str, scan_window_A: float, scan_step_A: float,
        convergence_mrad: float, hsa_fracs=DEFAULT_HSA, haadf_inner_mrad: float = 100.0,
        bin_factor: int = 4, slice_thickness_A: float = 1.0, device: str = "gpu",
        phonons: int = 0, test: bool = False, recon_full_box: bool = True) -> str:
    """@brief Simulate the toy membrane and emit everything the two channels need."""
    import ase.io
    require_gpu_job(device)
    s4 = _import_sim4d()

    atoms = ase.io.read(vasp)
    truth = np.load(truth_npz, allow_pickle=True)
    box_a = float(atoms.cell.lengths()[0])
    box_z = float(atoms.cell.lengths()[2])
    ctr = box_a / 2.0                                    # the four domains meet at the box centre

    # drive the VALIDATED module through its own globals -- no edit to sim/simulate_4dstem.py
    s4.DEVICE = device
    s4.BIN_FACTOR = bin_factor
    s4.SLICE_THICKNESS_A = slice_thickness_A
    s4.CONVERGENCE_MRAD = convergence_mrad
    s4.SCAN_CENTER_X_A, s4.SCAN_CENTER_Y_A = ctr, ctr
    s4.SCAN_WINDOW_A = scan_window_A
    s4.SCAN_STEP_A = scan_step_A
    s4.N_PHONONS = phonons
    s4.PER_SPECIES_SIGMA = bool(phonons)                 # per-species RT B factors when TDS is on

    print("=" * 68)
    print(f"[fusion] toy membrane {len(atoms)} atoms | box {box_a:.3f} x {box_a:.3f} x {box_z:.3f} A")
    print(f"[fusion] slab {float(truth['slab_thickness_A']):.2f} A "
          f"({int(truth['n_z'])} cells, column-homogeneous) | alpha = {convergence_mrad:.0f} mrad "
          f"| detector {s4.DETECTOR_MAX_ANGLE_MRAD:.0f} mrad")
    print(f"[fusion] scan {scan_window_A} A window at {scan_step_A} A step, centred on the "
          f"four-domain corner ({ctr:.2f}, {ctr:.2f})")
    print("=" * 68)

    pot = s4.build_potential(atoms, announce=True)
    probe = s4.build_probe(pot)
    scan, pos_xy, ny, _ = s4.make_scan(test)
    arr = s4.run_scan_binned(probe, atoms, scan)         # (M, N_b, N_b), y-fastest
    n_b = arr.shape[-1]

    scan_dir = os.path.join(out_dir, "01")
    dp_path, pos_path, A = s4.save_outputs(arr, pos_xy, __import__("pathlib").Path(scan_dir), box_a)
    s4.selftest_ordering(dp_path, pos_path, A, pos_xy, ny)
    beam_thickness = box_z if recon_full_box else float(truth["slab_thickness_A"])
    s4.write_driver_geometry(n_b, box_a, beam_thickness, __import__("pathlib").Path(scan_dir))

    # ---- the hollow detector: masks, measured dose split, simultaneous HAADF -------------
    from scipy.io import savemat
    lam = s4.wavelength_a()
    d_alpha_mrad = (bin_factor / box_a) * lam * 1e3
    theta = detector_axes(n_b, d_alpha_mrad)
    nx = len(pos_xy) // ny

    hsa_mrad = [f * convergence_mrad for f in hsa_fracs]
    for f, h in zip(hsa_fracs, hsa_mrad):
        m = hollow_mask(theta, h)
        p = os.path.join(scan_dir, f"mask_hsa{f:.2f}.mat")
        savemat(p, {"mask": m})
        print(f"[hollow] HSA {f:4.2f}a = {h:6.1f} mrad -> {p}  "
              f"({int((m == 0).sum())} of {n_b*n_b} px sent to EELS)")

    budget = dose_budget(A, theta, hsa_mrad, haadf_inner_mrad, s4.DETECTOR_MAX_ANGLE_MRAD)
    budget.update(dict(convergence_mrad=convergence_mrad, d_alpha_mrad=d_alpha_mrad, Ndpx=n_b,
                       detector_max_mrad=s4.DETECTOR_MAX_ANGLE_MRAD,
                       haadf_inner_mrad=haadf_inner_mrad, hsa_fracs=list(hsa_fracs),
                       hsa_mrad=hsa_mrad, box_A=box_a, box_z_A=box_z,
                       slab_thickness_A=float(truth["slab_thickness_A"]),
                       scan_shape=[nx, ny], scan_step_A=scan_step_A,
                       scan_window_A=scan_window_A, scan_center_A=[ctr, ctr],
                       beam_thickness_A=beam_thickness, n_phonons=phonons))
    with open(os.path.join(out_dir, "hollow_budget.json"), "w") as fh:
        json.dump(budget, fh, indent=2)

    print("\n== where the electrons go (measured on the simulated patterns) ==")
    for f, h in zip(hsa_fracs, hsa_mrad):
        b = budget["hsa"][f"{h:.1f}"]
        print(f"  HSA {f:4.2f}a ({h:6.1f} mrad):  EELS {b['eels_fraction']*100:5.1f}%   "
              f"ptychography {b['ptycho_fraction']*100:5.1f}%")
    print(f"  simultaneous HAADF [{haadf_inner_mrad:.0f}, {s4.DETECTOR_MAX_ANGLE_MRAD:.0f}] mrad: "
          f"{budget['haadf_fraction']*100:.2f}% of the beam")

    np.save(os.path.join(out_dir, "haadf.npy"),
            virtual_image(A, theta, haadf_inner_mrad, s4.DETECTOR_MAX_ANGLE_MRAD, (nx, ny)))
    np.save(os.path.join(out_dir, "bf.npy"),
            virtual_image(A, theta, 0.0, convergence_mrad, (nx, ny)))
    print(f"[fusion] wrote haadf.npy / bf.npy  ({nx} x {ny}) and hollow_budget.json")
    print(f"\nDone -> {out_dir}")
    return out_dir


def rebudget(out_dir: str, hsa_fracs, haadf_inner_mrad: float = 100.0) -> dict:
    """@brief Re-derive the dose split (and the masks) for new hollow angles, no re-simulation.

    The hole is applied at reconstruction time, so a finished dataset answers "what if the hole
    were this big" directly. Used for the HSA sweep and for the detector panel of the figure.
    """
    import h5py
    from scipy.io import savemat
    scan_dir = os.path.join(out_dir, "01")
    with h5py.File(os.path.join(scan_dir, "data_dp.hdf5"), "r") as f:
        dp = f["dp"][...]                                  # (Npos, N_b_x, N_b_y)
    b = json.load(open(os.path.join(out_dir, "hollow_budget.json")))
    n_b = int(b["Ndpx"])
    theta = detector_axes(n_b, float(b["d_alpha_mrad"]))
    alpha = float(b["convergence_mrad"])
    hsa_mrad = [f * alpha for f in hsa_fracs]
    for f, h in zip(hsa_fracs, hsa_mrad):
        savemat(os.path.join(scan_dir, f"mask_hsa{f:.2f}.mat"), {"mask": hollow_mask(theta, h)})
    nb = dose_budget(dp, theta, hsa_mrad, haadf_inner_mrad, float(b["detector_max_mrad"]))
    b["hsa"].update(nb["hsa"])
    b["hsa_fracs"] = sorted(set(list(b["hsa_fracs"]) + list(hsa_fracs)))
    b["hsa_mrad"] = [f * alpha for f in b["hsa_fracs"]]
    json.dump(b, open(os.path.join(out_dir, "hollow_budget.json"), "w"), indent=2)
    print(f"[hollow] re-derived the split for {len(hsa_fracs)} hollow angles (no re-simulation)")
    for f, h in zip(hsa_fracs, hsa_mrad):
        e = b["hsa"][f"{h:.1f}"]["eels_fraction"]
        print(f"  HSA {f:4.2f}a ({h:6.1f} mrad):  EELS {e*100:5.1f}%   ptychography {(1-e)*100:5.1f}%")
    return b


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default=os.path.join(HERE, "runs", "fusion"))
    ap.add_argument("--vasp", default=os.path.join(HERE, "sample", "toy_membrane.vasp"))
    ap.add_argument("--truth", default=os.path.join(HERE, "sample", "toy_truth.npz"))
    ap.add_argument("--scan-window", type=float, default=24.0, help="scan window [A]")
    ap.add_argument("--scan-step", type=float, default=0.3, help="scan step [A]")
    ap.add_argument("--convergence", type=float, default=100.0, help="probe semi-angle [mrad]")
    ap.add_argument("--bin-factor", type=int, default=4)
    ap.add_argument("--slice-thickness", type=float, default=1.0,
                    help="multislice slice [A]; 1.0 resolves the 2.08 A PbO/TiO2 plane spacing")
    ap.add_argument("--phonons", type=int, default=0, help="frozen-phonon configs (0 = coherent)")
    ap.add_argument("--haadf-inner", type=float, default=100.0)
    ap.add_argument("--hsa", type=float, nargs="+", default=list(DEFAULT_HSA),
                    help="hollow semi-angles as fractions of the convergence angle")
    ap.add_argument("--device", default="gpu", choices=["gpu", "cpu"])
    ap.add_argument("--test", action="store_true", help="3x3 scan shape check")
    ap.add_argument("--rebudget", action="store_true",
                    help="re-derive the dose split + masks for --hsa on an EXISTING dataset "
                         "(no re-simulation: the hole is applied at reconstruction time)")
    args = ap.parse_args(argv)

    if args.rebudget:
        rebudget(args.out_dir, tuple(args.hsa), args.haadf_inner)
        return 0

    os.makedirs(args.out_dir, exist_ok=True)
    run(args.out_dir, args.vasp, args.truth, args.scan_window, args.scan_step,
        args.convergence, tuple(args.hsa), args.haadf_inner, args.bin_factor,
        args.slice_thickness, args.device, args.phonons, args.test)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
