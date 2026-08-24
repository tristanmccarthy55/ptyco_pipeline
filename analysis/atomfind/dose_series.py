#!/usr/bin/env python
"""@file dose_series.py
@brief Portability harness: run the finder on the DOSE SERIES (NL105, 0.1 A, coherent).

The portability benchmark -- a different reconstruction geometry (105 layers, dz 0.666,
405x405), ~half NL70's phase amplitude, real recon noise -- run with the SAME config as
NL70 (nothing tuned per-dose). Uses the Pb rev2 kernel for both species (the rev2 Ti/O
kernels are broken; a data-derived kernel is worse). Measured results: RESULTS.md §4.

Usage:
    python atomfind/dose_series.py                       # all doses + the baselines
    python atomfind/dose_series.py 1e10 --no-baselines   # one dose, v3 only (fast)
    python atomfind/dose_series.py --json dose.json      # machine-readable table
"""
from __future__ import annotations
import os, sys
import numpy as np
from dataclasses import replace

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from atomfind import (config, align, psf as psfmod, deconv, find, validate,
                      uncertainty)

DOSE_ROOT = os.path.expanduser("~/Desktop/dose_series")
REV2_KERNEL = "psf_Pb_rev2_d1e10_vol.npy"     # resolved via config.data_path


def load_recon_mat(path):
    """PtychoShelves Niter*.mat -> complex64 (nL, Ny, Nx). object_roi is a cell array of
    per-layer complex arrays stored column-major, so each layer is transposed."""
    import h5py
    with h5py.File(path, "r") as f:
        g = f["outputs"]
        refs = g["object_roi"][:, 0]
        layers = []
        for r in refs:
            a = f[r][:]
            if a.dtype.names:
                a = a["real"] + 1j * a["imag"]
            layers.append(a.T)
        return np.array(layers).astype(np.complex64)


def dose_cfg(dose_tag):
    """Config for one dose volume. NL105 geometry; single shared kernel; no per-dose tuning."""
    return replace(
        config.preset("NL70_coherent"),
        name=f"dose{dose_tag}",
        recon_vol=os.path.join(DOSE_ROOT, f"dose{dose_tag}", "Niter200.mat"),
        dz=0.666,
        dose_e_per_A2=float(dose_tag),
        single_atom_vol=REV2_KERNEL,
        ti_kernel_vol=None,              # rev2 Ti kernel is broken -> Pb shape for both
        trim_z_A=(2.0, 68.0),            # 105 layers x 0.666 = 70 A
        zmax_show_A=68.0,
        exit_band_z_A=58.0,
        bulk_z_A=(10.0, 58.0),
    ).resolve()


def load_dose_volume(cfg):
    """Phase volume + dx for a dose-series .mat, mirroring align.load_phase."""
    vol = load_recon_mat(cfg.recon_vol)
    V = np.angle(vol).astype(float)
    V -= np.median(V, axis=(1, 2), keepdims=True)
    dx = cfg.dx if cfg.dx is not None else cfg.scan_window_A / V.shape[2]
    return V, dx


def _baselines(V, dx, al, pos, Z, olab, cfg, verbose=True):
    """@brief The deconvolve-then-peak-pick family, scored on the SAME volume as v3.

    Mirrors run_atomfind.py exactly (same kernels, same relative floors) so the dose
    comparison is like-for-like with the noiseless one in RESULTS.md sec.2. This is the
    comparison the report calls discriminating: on noiseless data v3's oxygen margin is
    real but modest, and the claim is that it WIDENS once there is noise to survive.

    @return {name: finder_report} for raw / RL / MEM peak-picking.
    """
    psfs = psfmod.build_psfs(V, dx, al, pos, Z, cfg)
    default_psf = psfs[psfs["_default"]]
    dec, dinfo = deconv.richardson_lucy_3d(V, default_psf, cfg)
    mem, minfo = deconv.mem_3d(V, default_psf, cfg)
    if verbose:
        print(f"  deconv RL {dinfo}")
        print(f"  deconv MEM {minfo}")
    out = {}
    for name, vol, kw in (("peaks3d_raw", V,   dict(rel_floor=0.02)),
                          ("peaks3d_rl",  dec, dict(rel_floor=0.02)),
                          ("peaks3d_mem", mem, dict(rel_floor=0.001, max_atoms=4000))):
        f = find.peaks3d(vol, cfg, dx, **kw)
        out[name], _ = validate.finder_report(f, pos, Z, al, cfg, olabel=olab)
    return out


def _o_split(rep):
    """(overlapped-O recall, isolated-O recall) -- the population where the method differs.

    Axially overlapped oxygen sits 1.95 A from a Ti along the beam, inside the axial response,
    and is the whole of the method's margin; in-plane-isolated oxygen is easy for everyone.
    Aggregate O recall averages the two and hides the effect.
    """
    ov = rep.get("O_axial_overlap", {}).get("recall", float("nan"))
    iso = rep.get("O_inplane_isolated", {}).get("recall", float("nan"))
    return ov, iso

def run_one(dose_tag, verbose=True, baselines=True):
    cfg = dose_cfg(dose_tag)
    if not os.path.exists(cfg.recon_vol):
        print(f"[dose{dose_tag}] MISSING {cfg.recon_vol}")
        return None
    V, dx = load_dose_volume(cfg)
    pos, Z = align.load_gt(cfg)
    al = align.register(V, dx, pos, Z, cfg)
    kernels = psfmod.species_kernels(cfg, dx)
    found, seeds = find.find_atoms_v3(V, cfg, dx, kernels)
    al_r = align.refine_with_atoms(al, found, pos, Z, cfg)
    olab = validate.classify_oxygen(pos, Z, cfg)
    rep, m = validate.finder_report(found, pos, Z, al_r, cfg, olabel=olab)
    base = _baselines(V, dx, al_r, pos, Z, olab, cfg, verbose=verbose) if baselines else {}
    if verbose:
        print(f"\n=== dose {dose_tag} ===  phase p99.9={np.percentile(np.clip(V,0,None),99.9):.3f}")
        print(f"  align  before: OFF {al.OFF:+.2f} A   after: OFF {al_r.OFF:+.2f} A  "
              f"mZ {al_r.mZ*100:+.2f}%  (corr_depth {al.corr_depth:.2f})")
        print(f"  found {rep['n_found']} ({(found['guided']==1).sum()} guided)  "
              f"prec {rep['precision']:.2f}  xy-RMS {rep.get('xy_rms_A', float('nan')):.3f} A  "
              f"z-RMS {rep['z_rms_A']:.2f} A")
        for s in ("Pb", "Ti", "O"):
            r = rep[s]
            print(f"    {s}: recall {r['recall']:.0%} (bulk {r['recall_bulk']:.0%})  "
                  f"z-RMS {r['z_rms_A']:.2f} A")
        print(f"  {validate.confusion_line(rep)}")
        if base:
            ov, iso = _o_split(rep)
            print(f"  {'detector':>20}  {'O bulk':>7} {'O overlap':>10} {'O isolated':>11}")
            print(f"  {'v3 (this work)':>20}  {rep['O']['recall_bulk']:>6.0%} "
                  f"{ov:>10.0%} {iso:>11.0%}")
            for k, lbl in (("peaks3d_raw", "peak-pick (raw)"),
                           ("peaks3d_rl", "RL + peak-pick"),
                           ("peaks3d_mem", "MEM + peak-pick")):
                r = base[k]
                bov, biso = _o_split(r)
                print(f"  {lbl:>20}  {r['O']['recall_bulk']:>6.0%} {bov:>10.0%} {biso:>11.0%}")
        # UQ: model sigma (single kernel here -> kernel-mismatch term is 0) + conformal.
        # This exercises the TRANSFER path: calibrate on THIS volume's own matched atoms.
        try:
            qtab = uncertainty.calibrate(found, m, cfg, alphas=cfg.uq_alphas,
                                         min_n=cfg.uq_min_stratum)
            for ln in uncertainty.uncertainty_report(found, m, cfg, qtab, alphas=cfg.uq_alphas):
                print("  " + ln)
        except Exception as e:
            print(f"  UQ: {e}")
    return dict(cfg=cfg, al=al, al_refined=al_r, found=found, rep=rep, V=V, dx=dx,
                pos=pos, Z=Z, baselines=base)


def main():
    import argparse, json
    ap = argparse.ArgumentParser(description="Dose-series portability benchmark + the "
                                             "deconvolve-then-peak-pick head-to-head.")
    ap.add_argument("doses", nargs="*", default=None,
                    help="e.g. 1e10 1e8 (default: all four)")
    ap.add_argument("--no-baselines", action="store_true",
                    help="v3 only; skips the two 3-D deconvolutions (much faster)")
    ap.add_argument("--json", default=None, help="write the full per-dose table here")
    ap.add_argument("--data-dir", default=None, help="see run_atomfind.py --data-dir")
    a = ap.parse_args()
    if a.data_dir:
        config.set_data_dir(a.data_dir)
    tags = a.doses or ["1e10", "1e8", "1e6", "1e4"]

    rows, blob = [], {}
    for t in tags:
        out = run_one(t, baselines=not a.no_baselines)
        if not out:
            continue
        rep = out["rep"]
        ov, iso = _o_split(rep)
        rows.append((t, rep["n_found"], rep["precision"], rep["Pb"]["recall"],
                     rep["Ti"]["recall"], rep["O"]["recall"],
                     validate.confusion_rate(rep), rep["z_rms_A"], ov, iso,
                     {k: _o_split(v)[0] for k, v in out["baselines"].items()}))
        blob[t] = {"v3": rep, "baselines": out["baselines"]}

    if rows:
        print("\n" + "=" * 78)
        print(f"{'dose':>6} {'found':>6} {'prec':>6} {'Pb':>6} {'Ti':>6} {'O':>6} "
              f"{'conf':>7} {'z-RMS':>7}")
        for r in rows:
            print(f"{r[0]:>6} {r[1]:>6} {r[2]:>6.2f} {r[3]:>5.0%} {r[4]:>5.0%} {r[5]:>5.0%} "
                  f"{r[6]:>6.1%} {r[7]:>6.2f}A")

        if not a.no_baselines:
            # The discriminating comparison: recall of AXIALLY OVERLAPPED oxygen (1.95 A from
            # a Ti along the beam) as a function of dose. The aggregate O column above averages
            # this together with the easy in-plane-isolated oxygen and hides the effect.
            print("\n" + "=" * 78)
            print("AXIALLY OVERLAPPED OXYGEN -- recall vs dose (the discriminating population)")
            print(f"{'dose':>6} {'v3':>7} {'raw pk':>8} {'RL pk':>8} {'MEM pk':>8} {'margin':>8}")
            for r in rows:
                b = r[10]
                best = max([v for v in b.values() if v == v], default=float("nan"))
                print(f"{r[0]:>6} {r[8]:>6.0%} {b.get('peaks3d_raw', float('nan')):>7.0%} "
                      f"{b.get('peaks3d_rl', float('nan')):>7.0%} "
                      f"{b.get('peaks3d_mem', float('nan')):>7.0%} "
                      f"{r[8] - best:>+7.0%}")
            print("(margin = v3 minus the BEST baseline at that dose; the report's claim is "
                  "that it widens as dose falls)")

    if a.json and blob:
        with open(a.json, "w") as f:
            json.dump(blob, f, indent=2, default=float)
        print(f"\nwrote {a.json}")


if __name__ == "__main__":
    main()
