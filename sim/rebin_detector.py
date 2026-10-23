#!/usr/bin/env python
"""
@file rebin_detector.py
@brief Record a finished sim on a smaller camera: crop its patterns to a smaller angle and/or sum B x B pixels into one,
       written as a new tagged dataset (aberration_experiment/detector256/README.md).

POST-processing, in the style of sim/add_poisson_noise.py: the multislice sim is not re-run, the original is untouched.
Two independent steps, both symmetric about the zero-angle pixel (the sim writes it at N/2 of an even N):
  crop  keep +-theta (--crop-mrad), i.e. the central n_f pixels, n_f/2 a multiple of B; the angle per pixel is unchanged,
        so the window the solver models (lambda / d_alpha) is unchanged and the object pixel grows (105 A / n_f).
  bin   SUM B x B blocks (an idealised camera: B^2 point samples per pixel). The blocks start ON the zero-angle pixel, which
        is exactly how the GPU engine bins its fine model when p.detector.upsampling = log2(B) (utils.binning_2D on the
        FFT-ordered array, zero angle at index 1): proved against the engine's code by
        aberration_experiment/detector256/test_upsampling_frame.py. NOT the centred blocks a real camera has -- the engine
        models these, and a fit with mismatched blocks gives the residual floor + phase ramp of 2026-09-23.

Probes, all the TRUE probe (probe_initial_true.mat) on the new geometry:
  probe_initial_true.mat        the fine grid (n_f px, the same angle per pixel as the sim): what an upsampled recon
                                (k = log2 B) models -- the true probe cropped in FOURIER space to +-theta. Exact while the
                                aperture fits inside the crop (80 mrad in 96: the power kept is printed); the reference
                                run's own presolve did the same crop (rescale_inputs, interpolateFT).
  probe_initial_true_up<k>.mat  k = 0 .. log2 B: the probe on the grid a recon with p.detector.upsampling = k models,
                                n_c 2^k px at the same object pixel = the fine probe cropped in REAL space (a smaller window;
                                the share of the probe outside it is printed). up<log2 B> is the fine probe itself.
sim_meta.mat: Ndpx = n_c, d_alpha_rad x B, rbf, dx_object_A, detector_record_mrad, bin_factor x B (so box_A / bin_factor is
still the data grid's window, as sim/make_probe.py assumes), plus rebin_* fields. data_position.hdf5 (A) and
aberrations.json are linked. rebin_detector.json records every number.

    python sim/rebin_detector.py --in-dir sim_out_af_armf_a080_lab_fs2p5_src0p1 --crop-mrad 96 --bin 4
      -> sim_out_af_armf_a080_lab_fs2p5_src0p1_detc96b4/01/ (1420 px at +-133.1 mrad -> 1024 at +-96.0 -> 256 px)
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import h5py
import numpy as np
from scipy.io import loadmat, savemat

CHUNK = 100                   # patterns per read (100 x 1420^2 float32 = 0.8 GB)
MAX_PROBE_LOSS = 1e-3         # the Fourier crop must keep the aperture: refuse if it loses more probe power than this


def det_tag(crop_mrad: float | None, b: int) -> str:
    """The data dir tag, as campaign/run_thin_atomfind.sh writes it: _detc<crop>[b<B>] ('.' -> 'p'), _detb<B> uncropped."""
    t = "_det"
    if crop_mrad:
        t += "c" + f"{crop_mrad:g}".replace(".", "p")
    if b > 1:
        t += f"b{b}"
    return t


def plan(n_in: int, d_alpha_mrad: float, crop_mrad: float | None, b: int) -> tuple[int, int, int]:
    """(start, n_f, n_c): crop rows/cols [start, start + n_f) of an n_in-pixel pattern whose zero-angle pixel is n_in // 2,
    then B x B blocks -> n_c = n_f / B. n_f / 2 is a multiple of B (the zero-angle pixel starts a block, and n_c is even, so
    it lands on n_c / 2: the engine's centre convention): the fewest whole blocks that keep +-crop_mrad (to 1e-3 px: the
    sims' angle per pixel is 0.18749989 mrad, so 96 mrad is 512.00003 px), capped at the pattern (no crop: all that fits).
    The same 'keep at least the angle' rule as the driver's 1.2 alpha presolve guard, without its float round-up."""
    if n_in % 2:
        raise ValueError(f"pattern width {n_in} is odd: the sim always writes an even width with zero angle at N/2")
    if b < 1 or b & (b - 1):
        raise ValueError(f"bin {b}: a power of 2 (the engine upsamples by 2^k)")
    cap = (n_in // 2 // b) * b
    half = cap if not crop_mrad else min(cap, b * int(np.ceil(crop_mrad / d_alpha_mrad / b - 1e-3 / b)))
    if half < b:
        raise ValueError(f"crop {crop_mrad} mrad at {d_alpha_mrad:.4f} mrad/px leaves no pixels")
    n_f = 2 * half
    return n_in // 2 - half, n_f, n_f // b


def rebin(dp: np.ndarray, start: int, n_f: int, b: int) -> np.ndarray:
    """(M, n, n) patterns -> (M, n_f / b, n_f / b): crop, then SUM b x b blocks from the crop's first pixel."""
    c = np.asarray(dp, dtype=np.float64)[:, start:start + n_f, start:start + n_f]
    if b > 1:
        n_c = n_f // b
        c = c.reshape(c.shape[0], n_c, b, n_c, b).sum(axis=(2, 4))
    return c.astype(np.float32)


def fourier_crop(P: np.ndarray, n_f: int) -> tuple[np.ndarray, float]:
    """The probe on an n_f grid at the SAME angle per pixel (same window, coarser object pixel): keep the central n_f
    Fourier coefficients (zero frequency at N/2 <-> n_f/2, the data's own crop). Returns (probe, power fraction kept)."""
    n = P.shape[0]
    F = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(P)))
    s = n // 2 - n_f // 2
    Fc = F[s:s + n_f, s:s + n_f]
    kept = float((np.abs(Fc) ** 2).sum() / (np.abs(F) ** 2).sum())
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(Fc))), kept


def real_crop(P: np.ndarray, n: int) -> tuple[np.ndarray, float]:
    """The central n x n of a probe (a smaller window at the same object pixel). Returns (probe, power fraction outside)."""
    s = P.shape[0] // 2 - n // 2
    c = P[s:s + n, s:s + n]
    return c, 1.0 - float((np.abs(c) ** 2).sum() / (np.abs(P) ** 2).sum())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", required=True, type=Path, help="sim output dir (the one CONTAINING 01/)")
    ap.add_argument("--crop-mrad", type=float, default=None, help="keep +- this angle (default: no crop)")
    ap.add_argument("--bin", type=int, default=1, help="sum B x B pixels (1, 2, 4, 8; default 1)")
    ap.add_argument("--out-dir", type=Path, default=None, help="default: <in-dir><tag> alongside the input")
    ap.add_argument("--scan", default="01")
    a = ap.parse_args(argv)
    if not a.crop_mrad and a.bin == 1:
        raise SystemExit("nothing to do: give --crop-mrad and/or --bin > 1")

    ind = a.in_dir / a.scan
    for f in ("data_dp.hdf5", "data_position.hdf5", "sim_meta.mat", "probe_initial_true.mat"):
        if not (ind / f).exists():
            raise SystemExit(f"ERROR: missing {ind / f} -- point --in-dir at a finished sim with its true probe")
    meta = loadmat(ind / "sim_meta.mat", simplify_cells=True)["meta"]
    d_alpha_mrad = float(meta["d_alpha_rad"]) * 1e3
    lam = 12.2643 / np.sqrt(float(meta["energy_kev"]) * 1e3 * (1 + 0.97845e-6 * float(meta["energy_kev"]) * 1e3))
    out = a.out_dir or a.in_dir.parent / f"{a.in_dir.name}{det_tag(a.crop_mrad, a.bin)}"
    outd = out / a.scan
    outd.mkdir(parents=True, exist_ok=True)
    (outd / "rebin_detector.json").unlink(missing_ok=True)      # rewritten last, when the copy is whole again

    with h5py.File(ind / "data_dp.hdf5", "r") as fin:
        dp = fin["dp"]
        npos, n_in, n_in2 = dp.shape
        if n_in != n_in2 or n_in != int(meta["Ndpx"]):
            raise SystemExit(f"ERROR: data {dp.shape} vs sim_meta Ndpx {meta['Ndpx']}")
        start, n_f, n_c = plan(n_in, d_alpha_mrad, a.crop_mrad, a.bin)
        tot_in = tot_out = 0.0
        with h5py.File(outd / "data_dp.hdf5", "w") as fout:
            dset = fout.create_dataset("dp", shape=(npos, n_c, n_c), dtype=np.float32)
            for i in range(0, npos, CHUNK):
                blk = dp[i:i + CHUNK]
                r = rebin(blk, start, n_f, a.bin)
                dset[i:i + CHUNK] = r
                tot_in += float(blk.astype(np.float64).sum()); tot_out += float(r.astype(np.float64).sum())
    itot = tot_out / npos
    theta = n_f // 2 * d_alpha_mrad
    print(f"[rebin] {npos} patterns {n_in} px (+-{n_in // 2 * d_alpha_mrad:.1f} mrad) -> crop {n_f} px "
          f"(rows {start}..{start + n_f - 1}, +-{theta:.2f} mrad) -> {a.bin}x{a.bin} sums -> {n_c} px; "
          f"counts kept {tot_out / tot_in:.6f}")

    # probes on the new geometry (savemat as the sim does: MATLAB load() + complex, same axis order as the input file)
    P = np.asarray(loadmat(ind / "probe_initial_true.mat")["probe"])
    if P.ndim != 2 or P.shape != (n_in, n_in):
        raise SystemExit(f"ERROR: probe_initial_true.mat holds {P.shape}, expected ({n_in}, {n_in})")
    P = P.astype(np.complex128)
    Pf, kept = fourier_crop(P, n_f)
    if 1.0 - kept > MAX_PROBE_LOSS:
        raise SystemExit(f"ERROR: the +-{theta:.1f} mrad crop cuts the probe's aperture ({1 - kept:.2e} of its power)")
    pinfo = dict(fourier_crop_power_kept=kept, grids={})
    k_max = int(np.log2(a.bin))
    pmeta = {"binning": False, "detector": {"binning": False}}
    for k in range(k_max + 1):
        Pk, outside = real_crop(Pf, n_c * 2 ** k)
        Pk = Pk * np.sqrt(itot) / (np.linalg.norm(Pk) + 1e-30)
        savemat(str(outd / f"probe_initial_true_up{k}.mat"), {"probe": Pk.astype(np.complex64), "p": pmeta})
        win = n_c * 2 ** k * lam / (n_f * float(meta["d_alpha_rad"]))
        pinfo["grids"][f"up{k}"] = dict(px=int(n_c * 2 ** k), window_A=win, power_outside=outside)
        print(f"[probe] up{k}: {n_c * 2 ** k} px, window {win:.2f} A, {100 * outside:.3f} % of the probe outside")
    link = outd / "probe_initial_true.mat"
    if link.exists() or link.is_symlink():
        link.unlink()
    os.symlink(f"probe_initial_true_up{k_max}.mat", link)
    print(f"[probe] Fourier crop to +-{theta:.2f} mrad keeps {kept:.8f} of the probe power")

    dx = lam / (n_f * float(meta["d_alpha_rad"]))
    new = dict(meta)
    new.update(Ndpx=int(n_c), d_alpha_rad=float(meta["d_alpha_rad"]) * a.bin, d_alpha_mrad=d_alpha_mrad * a.bin,
               rbf=float(meta["convergence_mrad"]) / (d_alpha_mrad * a.bin), dx_object_A=float(dx),
               detector_record_mrad=float(theta), bin_factor=int(meta["bin_factor"]) * a.bin,
               rebin_crop_px=int(n_f), rebin_bin=int(a.bin), rebin_source=str(ind.resolve()),
               rebin_crop_mrad_requested=float(a.crop_mrad or 0.0))
    savemat(str(outd / "sim_meta.mat"), {"meta": new})
    rec = dict(source=str(ind.resolve()), n_in=n_in, d_alpha_mrad_in=d_alpha_mrad, crop_mrad_requested=a.crop_mrad,
               crop_start=start, n_f=n_f, theta_max_mrad=theta, bin=a.bin, n_c=n_c, d_alpha_mrad_out=d_alpha_mrad * a.bin,
               object_pixel_A=dx, counts_kept=tot_out / tot_in, mean_pattern_total=itot, probes=pinfo,
               block_alignment="blocks start on the zero-angle pixel (engine binning_2D on FFT-ordered data)")
    for name in ("data_position.hdf5", "aberrations.json"):
        src = ind / name
        if src.exists():
            dst = outd / name
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            os.symlink(os.path.relpath(src, outd), dst)
    # LAST: the record marks the copy whole (campaign/run_thin_atomfind.sh reuses a copy only when it is there)
    (outd / "rebin_detector.json").write_text(json.dumps(rec, indent=2))
    print(f"[rebin] wrote {outd}  (object pixel {dx:.4f} A; positions linked)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
