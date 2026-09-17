#!/usr/bin/env python
"""
@file make_probe.py
@brief Write probe_initial.mat for a trial defocus C1, through the simulator's own
       build_initial_probe, so the probe grid is the data's by construction.

The C1 defocus search (aberration_experiment/NEXT_PHASE.md step 1) holds C3/C5 at their corrector
values and reconstructs the same data through a ladder of fixed probes that differ only in C1. Every
trial probe must sit on exactly the grid the sim wrote the data on: the binned detector's N_b pixels,
extent = box / BIN, the data's semi-angle and energy. That geometry is read from the sim's own
sim_meta.mat, the aberration tableau from its aberrations.json, and the probe is built by
simulate_4dstem.build_initial_probe with those globals set -- the very function that wrote the true
probe. Checked 2026-09-17 against the packed true probes (unit-normalised): overlap 1.000000 and
max|dP|/max|P| = 2.0e-7 (a70, BIN 4, 356 px) and 2.8e-6 (a90, BIN 2, 712 px) -- float32 rounding.

Normalisation mirrors the sim: ||P||^2 = mean total counts per pattern, estimated from a strided
subset of data_dp.hdf5 (the engine then refines the scale itself, initial_probe_rescaling).

Self-check: when probe_initial_true.mat is present its overlap with the trial probe is printed, and
when the trial C1 equals the sim's C1 the two must agree (shape-normalised) or the writer exits
non-zero -- a writer that cannot reproduce the truth must not feed a search.

    python sim/make_probe.py --sim-dir sim_out_af_a70_lab/01 --c1=-40 --out recon_x/01/probe_initial.mat
    python sim/make_probe.py --sim-dir <packed sim>/01 --c1=-60 --norm-from true --out /tmp/p.mat   # no data locally

Negative values: use --c1=-40 (with '=') so argparse does not read them as flags.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

HERE = Path(__file__).resolve().parent
SELF_CHECK_MAX_REL = 2e-3     # max|P/||P|| - T/||T|||, relative to max|T|; measured <= 2.8e-6 locally
                              # (loose on purpose: Blythe's numpy/BLAS may round differently)
SELF_CHECK_MIN_OVERLAP = 0.9999


def _load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def mean_pattern_total(dp_path: Path, n_patterns: int) -> tuple[float, int]:
    """Mean total counts per pattern over ~n_patterns evenly strided positions."""
    import h5py
    with h5py.File(dp_path, "r") as f:
        ds = f["dp"]                                  # (Npos, N_b, N_b)
        npos = ds.shape[0]
        idx = np.unique(np.linspace(0, npos - 1, min(n_patterns, npos)).round().astype(int))
        tot = [float(np.asarray(ds[i], dtype=np.float64).sum()) for i in idx]
    return float(np.mean(tot)), len(idx)


def shape_compare(P, T):
    """(overlap, max relative difference) of two probes after normalising each to unit norm."""
    p = P / np.linalg.norm(P); t = T / np.linalg.norm(T)
    return float(abs(np.vdot(t.ravel(), p.ravel()))), float(np.abs(p - t).max() / np.abs(t).max())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sim-dir", type=Path, required=True,
                    help="the sim's 01/ dir (or a recon 01/ that links to it): sim_meta.mat, aberrations.json")
    ap.add_argument("--c1", type=float, required=True, help="trial abTEM defocus [A] (the sim's sign convention)")
    ap.add_argument("--c3", type=float, default=None, help="override C30 [A] (default: aberrations.json)")
    ap.add_argument("--c5", type=float, default=None, help="override C50 [A] (default: aberrations.json)")
    ap.add_argument("--out", type=Path, required=True, help="probe_initial.mat to write")
    ap.add_argument("--norm-from", choices=["data", "true"], default="data",
                    help="||P||^2 from the data's mean pattern total (default) or from probe_initial_true.mat")
    ap.add_argument("--norm-patterns", type=int, default=100, help="patterns sampled for --norm-from data")
    a = ap.parse_args(argv)

    sd = a.sim_dir
    meta = loadmat(str(sd / "sim_meta.mat"), simplify_cells=True)["meta"]
    ab_path = sd / "aberrations.json"
    if not ab_path.exists():
        raise SystemExit(f"[probe] {ab_path} missing -- the C3/C5 tableau must come from the sim, not be typed")
    ab = json.loads(ab_path.read_text())

    sim = _load_module("simulate_4dstem", HERE / "simulate_4dstem.py")
    n_b, box, binf = int(meta["Ndpx"]), float(meta["box_A"]), int(meta["bin_factor"])
    alpha, e_kev = float(meta["convergence_mrad"]), float(meta["energy_kev"])
    if abs(sim.ENERGY_EV / 1e3 - e_kev) > 1e-6 or abs(float(ab["energy_eV"]) - sim.ENERGY_EV) > 1e-3:
        raise SystemExit(f"[probe] energy mismatch: sim module {sim.ENERGY_EV} eV, sim_meta {e_kev} keV, "
                         f"aberrations.json {ab['energy_eV']} eV")
    if abs(float(ab["convergence_mrad"]) - alpha) > 1e-9:
        raise SystemExit(f"[probe] alpha mismatch: sim_meta {alpha} vs aberrations.json {ab['convergence_mrad']}")

    aberr_sim = {k: float(v) for k, v in ab["aberrations_Cnm_A_phi_rad"].items()}
    aberr = dict(aberr_sim)
    for key, val in (("C30", a.c3), ("C50", a.c5)):
        if val is not None:
            if key in aberr and abs(aberr[key] - val) > 1e-6 * max(1.0, abs(val)):
                print(f"[probe] NOTE {key} override {val:g} A differs from the sim's {aberr[key]:g} A")
            aberr[key] = float(val)
    c1_true = float(ab["defocus_A"])

    # the sim's globals, set exactly as its main() sets them for this data
    sim.CONVERGENCE_MRAD = alpha
    sim.BIN_FACTOR = binf
    sim.ABERRATIONS = aberr
    sim.DEFOCUS_A = float(a.c1)
    P = sim.build_initial_probe(n_b, box, aberrated=True).astype(np.complex128)
    assert P.shape == (n_b, n_b), f"probe {P.shape} != DP size {(n_b, n_b)}"

    true_path = sd / "probe_initial_true.mat"
    T = loadmat(str(true_path))["probe"].astype(np.complex128) if true_path.exists() else None
    if a.norm_from == "true":
        if T is None:
            raise SystemExit(f"[probe] --norm-from true but {true_path} missing")
        itot, n_used = float(np.linalg.norm(T) ** 2), 0
    else:
        if not (sd / "data_dp.hdf5").exists():
            raise SystemExit(f"[probe] {sd / 'data_dp.hdf5'} missing (use --norm-from true without the data)")
        itot, n_used = mean_pattern_total(sd / "data_dp.hdf5", a.norm_patterns)
    P *= np.sqrt(itot) / (np.linalg.norm(P) + 1e-30)
    P32 = P.astype(np.complex64)

    if a.out.is_symlink():
        raise SystemExit(f"[probe] {a.out} is a symlink (a known-probe link?) -- refusing to write through it")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    savemat(str(a.out), {"probe": P32, "p": {"binning": False, "detector": {"binning": False}}})

    info = dict(alpha_mrad=alpha, energy_eV=sim.ENERGY_EV, c1_A=float(a.c1), c1_sim_A=c1_true,
                dc1_A=float(a.c1) - c1_true, aberrations_Cnm_A_phi_rad=aberr, bin_factor=binf, n_b=n_b,
                box_A=box, extent_A=box / binf, norm_from=a.norm_from, norm_patterns=n_used,
                mean_pattern_total=itot, sim_dir=str(sd.resolve()))
    try:
        import abtem
        info["abtem_version"] = abtem.__version__
    except Exception:
        pass
    try:
        pp = _load_module("plan_probe", HERE.parent / "campaign" / "plan_probe.py")
        fwhm, d50, d90, d99 = pp.sizes(P, box / binf)
        info.update(fwhm_A=float(fwhm), d50_A=float(d50), d90_A=float(d90), d99_A=float(d99))
    except Exception as e:                         # sizes are provenance, never a reason to fail
        info["sizes_error"] = str(e)

    ok = True
    if T is not None:
        ov, rel = shape_compare(P, T)
        info.update(overlap_with_sim_probe=ov, max_rel_diff_vs_sim_probe=rel,
                    norm_ratio_vs_sim_probe=float(np.linalg.norm(P) / np.linalg.norm(T)))
        print(f"[probe] vs probe_initial_true.mat (C1 {c1_true:g} A): overlap {ov:.6f}, "
              f"max rel diff {rel:.2e}, norm ratio {info['norm_ratio_vs_sim_probe']:.5f}")
        if abs(float(a.c1) - c1_true) < 1e-9 and aberr == aberr_sim:     # the sim's own probe, exactly
            ok = ov >= SELF_CHECK_MIN_OVERLAP and rel <= SELF_CHECK_MAX_REL
            info["self_check"] = "PASS" if ok else "FAIL"
            print(f"[probe] SELF-CHECK {info['self_check']}: trial C1 = sim C1, so the writer must reproduce "
                  f"the sim's probe (overlap >= {SELF_CHECK_MIN_OVERLAP}, max rel diff <= {SELF_CHECK_MAX_REL:g})")
    json_path = a.out.with_suffix(".json")
    if json_path.is_symlink():
        json_path.unlink()
    json_path.write_text(json.dumps(info, indent=2))
    print(f"[probe] wrote {a.out} ({n_b}x{n_b}, extent {box / binf:.3f} A, alpha {alpha:g} mrad, "
          f"C1 {a.c1:g} A = sim {c1_true:+g} {float(a.c1) - c1_true:+g}; C30 {aberr.get('C30', 0):g}, "
          f"C50 {aberr.get('C50', 0):g}; d90 {info.get('d90_A', float('nan')):.2f} A) + {json_path.name}")
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
