#!/usr/bin/env python
"""Portability harness: run the finder on the DOSE SERIES (NL105, 0.1 A, coherent).

This is the dataset the portability fixes are for -- a different reconstruction geometry
(105 layers, dz 0.666, 405x405), roughly half NL70's phase amplitude, and real
reconstruction noise. Nothing here is tuned per-dose: the point is that ONE configuration
works across volumes.

Kernel: the Pb rev2 kernel for BOTH species (the rev2 Ti/O kernels are broken -- argmax in
the corner, SNR ~2.5 -- and a data-derived kernel is worse: axial FWHM 4.0 A vs 2.0 A).

Usage:
    ~/hyperspy-bundle/bin/python atomfind/dose_series.py            # all doses
    ~/hyperspy-bundle/bin/python atomfind/dose_series.py 1e10       # one dose
"""
from __future__ import annotations
import os, sys
import numpy as np
from dataclasses import replace

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from atomfind import config, align, psf as psfmod, find, validate

DOSE_ROOT = os.path.expanduser("~/Desktop/dose_series")
REV2_KERNEL = os.path.expanduser("~/Desktop/psf_Pb_rev2_d1e10_vol.npy")


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
    )


def load_dose_volume(cfg):
    """Phase volume + dx for a dose-series .mat, mirroring align.load_phase."""
    vol = load_recon_mat(cfg.recon_vol)
    V = np.angle(vol).astype(float)
    V -= np.median(V, axis=(1, 2), keepdims=True)
    dx = cfg.dx if cfg.dx is not None else cfg.scan_window_A / V.shape[2]
    return V, dx


def run_one(dose_tag, verbose=True):
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
        if "sigma_coverage_1s" in rep:
            cv = rep["sigma_coverage_1s"]
            print(f"  sigma coverage 1s: x {cv['x']:.0%} y {cv['y']:.0%} z {cv['z']:.0%}")
    return dict(cfg=cfg, al=al, al_refined=al_r, found=found, rep=rep, V=V, dx=dx,
                pos=pos, Z=Z)


def main():
    tags = sys.argv[1:] or ["1e10", "1e8", "1e6", "1e4"]
    rows = []
    for t in tags:
        out = run_one(t)
        if out:
            rep = out["rep"]
            rows.append((t, rep["n_found"], rep["precision"], rep["Pb"]["recall"],
                         rep["Ti"]["recall"], rep["O"]["recall"],
                         validate.confusion_rate(rep), rep["z_rms_A"]))
    if rows:
        print("\n" + "="*78)
        print(f"{'dose':>6} {'found':>6} {'prec':>6} {'Pb':>6} {'Ti':>6} {'O':>6} {'conf':>7} {'z-RMS':>7}")
        for t, n, p, pb, ti, o, cf, zr in rows:
            print(f"{t:>6} {n:>6} {p:>6.2f} {pb:>5.0%} {ti:>5.0%} {o:>5.0%} {cf:>6.1%} {zr:>6.2f}A")


if __name__ == "__main__":
    main()
