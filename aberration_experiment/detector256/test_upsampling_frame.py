#!/usr/bin/env python
"""
@file test_upsampling_frame.py
@brief Unit test, no MATLAB: sim/rebin_detector.py bins the data exactly the way the GPU engine bins its own model when
       p.detector.upsampling = 2 -- same blocks, same array frame -- traced through the engine's code path, ported line by
       line to numpy (MATLAB's column-major reshapes, its fftshift_2D index rule, crop_pad, interpolateFT).

Why: data summed in blocks that differ from the model's blocks by one fine pixel give a residual floor and a phase ramp
(the campaign paid for that on 2026-09-23). Upsampling models the summing, so data and model must agree to the pixel.

The engine's path (file : line, as of 2026-10-23), each step ported below and its source line asserted to still be there:
  data   load_data.m (matlab_aps) h5read '/dp'  -> MATLAB reverses the axes   (dp[n, a, b] -> D(b, a, n))
         binning.m       unbinning_2D(data, 2^k) / 4^k          (each detector pixel spread over a 2^k x 2^k block)
         prep_data_matlab.m (matlab_ps; both loaders call it)   fmag = fftshift_2D(sqrt(data)) * renorm
         load_from_p.m   diffraction = (fmag / renorm)^2; custom_data_flip [0 0 1]: fftshift_2D -> permute -> ifftshift_2D
         rescale_inputs.m (presolve only)                       fftshift_2D -> crop_pad -> ifftshift_2D
         init_solver.m   binning_2D(Diffraction, 2^k) * 4^k     (back to the detector's pixels)
  model  get_reciprocal_model.m   aPsi2 = binning_2D(|Psi|^2, 2^k)   (block MEAN; the probe's iteration-0 rescale absorbs
                                  the 4^k: LSQML.m probe_amp_corr = sqrt(sum modF^2 / sum aPsi^2))
         modulus_constraint.m     modF, aPsi unbinned back to the fine grid with unbinning_2D

The claim tested (T3): if the engine's fine model equals the R1 data in its own frame (the full-pixel reconstruction's
target, which the reference shows the engine reaches), then its binned model equals the R2 data up to that one constant.
R1 and R2 share the fine grid, the object pixel, the probe file and the flips, so R2's fine forward model IS R1's.

    ~/hyperspy-bundle/bin/python aberration_experiment/detector256/test_upsampling_frame.py
"""
from __future__ import annotations

import importlib.util
import re
import sys
import tempfile
from math import ceil, floor
from pathlib import Path

import h5py
import numpy as np
from scipy.io import loadmat, savemat

REPO = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("rebin_detector", REPO / "sim" / "rebin_detector.py")
rb = importlib.util.module_from_spec(spec); spec.loader.exec_module(rb)

N_IN, D_ALPHA_MRAD, ALPHA, NPOS = 1420, 0.18749989426407307, 80.0, 2   # the kept sim's own sim_meta.mat (fs 2.5 + src 0.1)
LAM = 0.0196877
FAILS = []


def check(name, ok, detail=""):
    print(f"{'ok  ' if ok else 'FAIL'}  {name}{('  -- ' + detail) if detail else ''}")
    if not ok:
        FAILS.append(name)


# ------------------------------------------------------------------ the engine's helpers, ported literally
def fftshift_2D(x):                                  # +math/fftshift_2D.m: idx = [p+1:m 1:p], p = ceil(m/2)
    for k in (0, 1):
        m = x.shape[k]; p = ceil(m / 2)
        x = np.take(x, np.r_[p:m, 0:p], axis=k)
    return x


def ifftshift_2D(x):                                 # +math/ifftshift_2D.m: p = floor(m/2)
    for k in (0, 1):
        m = x.shape[k]; p = floor(m / 2)
        x = np.take(x, np.r_[p:m, 0:p], axis=k)
    return x


def mreshape(x, shape):                              # MATLAB reshape (column-major)
    return np.reshape(x, shape, order="F")


def binning_2D(x, b, centered=False):                # +utils/binning_2D.m (Npix >= binning branch)
    x = x if x.ndim == 3 else x[:, :, None]
    if centered:                                     # x(ceil(b/2):end-ceil(b/2)-1, ...)
        c = ceil(b / 2)
        x = x[c - 1:x.shape[0] - c - 1, c - 1:x.shape[1] - c - 1, :]
    n1, n2, n3 = x.shape
    if n1 % b or n2 % b:
        n1, n2 = (n1 // b) * b, (n2 // b) * b; x = x[:n1, :n2, :]
    x = mreshape(x, (b, n1 // b, b, n2 // b, n3)).sum(axis=(0, 2)) / (b * b)
    if centered:                                     # padarray(x, [1,1], 'replicate', 'post')
        x = np.concatenate([x, x[-1:, :, :]], axis=0); x = np.concatenate([x, x[:, -1:, :]], axis=1)
    return x


def unbinning_2D(x, u):                              # +utils/unbinning_2D.m
    x = x if x.ndim == 3 else x[:, :, None]
    n1, n2, n3 = x.shape
    x = mreshape(x, (1, n1, 1, n2, n3))
    x = np.tile(x, (u, 1, u, 1, 1))
    return mreshape(x, (n1 * u, n2 * u, n3))


def crop_pad(img, outsize, fill=0):                  # +utils/crop_pad.m (1-based centre floor(N/2)+1)
    nin = np.array(img.shape[:2]); nout = np.array(outsize[:2])
    if np.all(nin == nout):
        return img
    out = np.zeros(tuple(nout) + img.shape[2:], dtype=img.dtype) + fill
    co = (nout // 2 + 1) - (nin // 2 + 1)
    o0 = np.maximum(co, 0); o1 = np.minimum(co + nin, nout)
    i0 = np.maximum(-co, 0); i1 = np.minimum(-co + nout, nin)
    out[o0[0]:o1[0], o0[1]:o1[1]] = img[i0[0]:i1[0], i0[1]:i1[1]]
    return out


def interpolateFT(im, outsize):                      # +utils/interpolateFT.m
    imFT = fftshift_2D(np.fft.fft2(im, axes=(0, 1)))
    out = crop_pad(imFT, outsize)
    return np.fft.ifft2(ifftshift_2D(out), axes=(0, 1)) * (outsize[0] * outsize[1] / (im.shape[0] * im.shape[1]))


def engine_data(h5_path, k=0, presolve=None):
    """The data array the engine compares its (binned) model with, from a data_dp.hdf5, following the path above."""
    with h5py.File(h5_path, "r") as f:
        dp = f["dp"][...].astype(np.float64)
    data = np.transpose(dp, (2, 1, 0))                         # h5read reverses the axes: D(i,j,n) = dp[n,j,i]
    if k:
        data = unbinning_2D(data, 2 ** k) / 4 ** k             # binning.m (upsampling branch)
    renorm = np.sqrt(1 / (data.sum(axis=(0, 1)).max() / (data.shape[0] * data.shape[1])))
    fmag = fftshift_2D(np.sqrt(data)) * renorm                 # prep_data_matlab.m (fmask = true)
    diff = (fmag / renorm) ** 2                                # load_from_p.m
    diff = ifftshift_2D(np.transpose(fftshift_2D(diff), (1, 0, 2)))   # custom_data_flip [0 0 1]
    if presolve is not None:                                   # rescale_inputs.m (presolve), size in FINE pixels
        diff = ifftshift_2D(crop_pad(fftshift_2D(diff), (presolve, presolve)))
    if k:
        diff = binning_2D(diff, 2 ** k) * 4 ** k               # init_solver.m
    return diff


def rel(a, b):
    return float(np.abs(a - b).max() / np.abs(b).max())


def overlap(a, b):
    return float(abs(np.vdot(a, b)) / (np.linalg.norm(a) * np.linalg.norm(b)))


# ------------------------------------------------------------------ the engine's source still says what is ported
def require(rel_path, *snippets):
    src = re.sub(r"\s+", " ", (REPO / rel_path).read_text(errors="replace"))
    for s in snippets:
        check(f"source {rel_path}: {s[:60]}", re.sub(r"\s+", " ", s) in src)


def test_sources():
    P = "ptycho/"
    require(P + "+detector/+prep_data/+matlab_aps/load_data.m", "data = h5read(files{p.scanID},'/dp');")
    require(P + "+detector/+prep_data/+matlab_aps/binning.m", "data = utils.unbinning_2D(data, upsample) / upsample^2;")
    require(P + "+detector/+prep_data/+matlab_aps/matlab_aps.m",
            "p = detector.prep_data.matlab_ps.prep_data_matlab(p, data, fmask);")
    require(P + "+detector/+prep_data/+matlab_ps/prep_data_matlab.m", "p.fmag = fftshift_2D(sqrt(data)) * p.renorm;")
    require(P + "+engines/+GPU/+initialize/load_from_p.m",
            "self.diffraction = (single(p.fmag .* p.fmask) / single(p.renorm) ).^2;",
            "self.diffraction = fftshift_2D(self.diffraction);",
            "self.diffraction = permute(self.diffraction, [2,1,3]);",
            "self.diffraction = ifftshift_2D(self.diffraction);")
    require(P + "+engines/+GPU/+initialize/rescale_inputs.m", "self.diffraction = fftshift_2D(self.diffraction);",
            "self.diffraction = crop_pad(self.diffraction, self.Np_p);", "self.diffraction = ifftshift_2D(self.diffraction);")
    require(P + "+engines/+GPU/+initialize/init_solver.m",
            "Diffraction = utils.binning_2D(Diffraction, 2^par.upsampling_data_factor) * (2^(2*par.upsampling_data_factor));")
    require(P + "+engines/+GPU/private/get_reciprocal_model.m", "aPsi2 = utils.binning_2D(aPsi2, 2^par.upsampling_data_factor);")
    require(P + "+engines/+GPU/private/modulus_constraint.m", "modF = utils.unbinning_2D(modF, 2^par.upsampling_data_factor );",
            "aPsi = utils.unbinning_2D(aPsi, 2^par.upsampling_data_factor );")
    require("+utils/binning_2D.m", "centered = false;")
    require(P + "run_synthetic_recon_ML.m", "custom_data_flip  = [0,0,1];", "prepare.data_preparator        = 'matlab_aps';")
    # the 2026-10-23 changes this experiment depends on
    require(P + "+core/initial_checks.m", "if ~strcmpi(p.prepare.data_preparator, 'matlab_aps')")
    require(P + "+core/ptycho_model_probe.m", "size(probe,2)] == asize*2^p.detector.upsampling")
    require(P + "+core/apply_binning.m", "utils.crop_pad(p.object{ii},p.object_size(ii,:), mean(p.object{ii}(:)))")
    require(P + "+engines/+GPU/LSQML.m", "zeros([par.Np_p_presolve / 2^par.upsampling_data_factor, par.p.numpos])")


# ------------------------------------------------------------------ a fake sim, run through the real script
def true_probe(n):
    """An aberrated, NON-symmetric probe on the sim's grid (80 mrad aperture), centred at n/2 as abTEM builds it."""
    k = (np.arange(n) - n // 2) * D_ALPHA_MRAD * 1e-3 / LAM          # 1/A, centred frame
    kx, ky = np.meshgrid(k, k, indexing="ij"); a = np.hypot(kx, ky) * LAM; phi = np.arctan2(ky, kx)
    chi = 2 * np.pi / LAM * (-958.6 * a ** 2 / 2 - 4.066e5 * a ** 4 / 4 + 12.4 * a ** 2 / 2 * np.cos(2 * (phi - 1.47))
                             + 4287.6 * a ** 3 / 3 * np.cos(phi - 5.6) + 164.3 * a ** 3 / 3 * np.cos(3 * (phi - 0.68)))
    F = (a <= ALPHA * 1e-3) * np.exp(-1j * chi)
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F)))


def fake_sim(root: Path):
    d = root / "sim_out_fake" / "01"; d.mkdir(parents=True)
    rng = np.random.default_rng(7)
    dp = rng.random((NPOS, N_IN, N_IN)).astype(np.float32) + 0.1     # no symmetry at all: any frame error shows
    with h5py.File(d / "data_dp.hdf5", "w") as f:
        f.create_dataset("dp", data=dp)
    with h5py.File(d / "data_position.hdf5", "w") as f:
        f.create_dataset("probe_positions_0", data=np.zeros((2, NPOS)))
    meta = dict(Ndpx=N_IN, d_alpha_rad=D_ALPHA_MRAD * 1e-3, d_alpha_mrad=D_ALPHA_MRAD, energy_kev=300.0,
                convergence_mrad=ALPHA, bin_factor=2, box_A=210.0, region_side_A=210.0, beam_thickness_A=27.525,
                detector_sampling="point", dx_object_A=LAM / (N_IN * D_ALPHA_MRAD * 1e-3))
    savemat(str(d / "sim_meta.mat"), {"meta": meta})
    P = true_probe(N_IN).astype(np.complex64)
    savemat(str(d / "probe_initial_true.mat"), {"probe": P, "p": {"binning": False, "detector": {"binning": False}}})
    return d.parent, P


def run_rebin(sim_dir, crop, b):
    args = ["--in-dir", str(sim_dir), "--bin", str(b)] + (["--crop-mrad", str(crop)] if crop else [])
    assert rb.main(args) == 0
    return Path(f"{sim_dir}{rb.det_tag(crop, b)}") / "01"


def main() -> int:
    test_sources()
    with tempfile.TemporaryDirectory() as tmp:
        sim, P = fake_sim(Path(tmp))
        full = sim / "01" / "data_dp.hdf5"
        r1 = run_rebin(sim, 96, 1)          # R1: +-96 mrad, full pixels
        r2 = run_rebin(sim, 96, 4)          # R2/R3: +-96 mrad, 4 x 4 sums
        r4 = run_rebin(sim, 133, 4)         # R4: +-133 mrad (1416), 4 x 4 sums -> 354
        r512 = run_rebin(sim, 96, 2)        # C512: +-96 mrad, 2 x 2 sums -> 512 (upsampling 1)
        r128 = run_rebin(sim, 96, 8)        # C128: +-96 mrad, 8 x 8 sums -> 128 (upsampling 3)

        # T1 the script's geometry
        m1 = loadmat(r1 / "sim_meta.mat", simplify_cells=True)["meta"]
        m2 = loadmat(r2 / "sim_meta.mat", simplify_cells=True)["meta"]
        m4 = loadmat(r4 / "sim_meta.mat", simplify_cells=True)["meta"]
        check("T1 plan: R1 and R2 crop the same 1024 px (96 mrad = 512.00003 px), R4 1416 (+-132.75 mrad)",
              rb.plan(N_IN, D_ALPHA_MRAD, 96, 1)[1:] == (1024, 1024) and rb.plan(N_IN, D_ALPHA_MRAD, 96, 4)[1:] == (1024, 256)
              and rb.plan(N_IN, D_ALPHA_MRAD, 133, 4)[1:] == (1416, 354) and rb.plan(N_IN, D_ALPHA_MRAD, 133, 1)[1] == 1420
              and rb.plan(N_IN, D_ALPHA_MRAD, 95.9, 1)[1] == 1024 and rb.plan(N_IN, D_ALPHA_MRAD, 96.1, 1)[1] == 1026)
        check("T1 R1 1420 -> 1024 px at 0.1875 mrad, object pixel 0.1025 A",
              m1["Ndpx"] == 1024 and abs(m1["d_alpha_rad"] - D_ALPHA_MRAD * 1e-3) < 1e-12 and abs(m1["dx_object_A"] - 0.10253) < 1e-4,
              f"{m1['Ndpx']} px, {m1['d_alpha_rad'] * 1e3:.4f} mrad, {m1['dx_object_A']:.5f} A")
        check("T1 R2 -> 256 px at 0.75 mrad, same object pixel; box/bin_factor = the 26 A data window",
              m2["Ndpx"] == 256 and abs(m2["d_alpha_rad"] - 4 * D_ALPHA_MRAD * 1e-3) < 1e-12
              and abs(m2["dx_object_A"] - m1["dx_object_A"]) < 1e-12 and abs(m2["box_A"] / m2["bin_factor"] - 26.25) < 0.01,
              f"{m2['Ndpx']} px, {m2['d_alpha_rad'] * 1e3:.3f} mrad, window {m2['box_A'] / m2['bin_factor']:.2f} A")
        check("T1 R4 1420 -> 1416 -> 354 px", m4["Ndpx"] == 354 and m4["rebin_crop_px"] == 1416, f"{m4['Ndpx']} px")

        # T2 R1 = the full data cropped by the engine's own presolve crop (rescale_inputs): the zero-angle pixel is centred
        # in the engine's convention. (The reference's presolve itself was 1026 px: the guard's ceil of 512.00003.)
        E_ref_pre = engine_data(full, 0, presolve=1024)
        E1 = engine_data(r1 / "data_dp.hdf5", 0)
        check("T2 R1 data == the engine's presolve crop of the full data (1420 -> 1024 in the engine frame)",
              E1.shape == E_ref_pre.shape and rel(E1, E_ref_pre) < 1e-6, f"max rel diff {rel(E1, E_ref_pre):.1e}")

        # T3 THE TEST: R2 data == the engine's binned model of the R1 truth, after the iteration-0 probe rescale
        E2 = engine_data(r2 / "data_dp.hdf5", 2)
        model = binning_2D(E1, 4)                          # get_reciprocal_model.m on a perfect fine model
        c = E2.sum() / model.sum()                         # LSQML.m probe_amp_corr^2 (= 16: sum vs mean)
        check("T3 R2 data == engine binning of the R1-frame truth (k = 2)", E2.shape == (256, 256, NPOS)
              and abs(c - 16) < 1e-6 and rel(E2, c * model) < 1e-6, f"shape {E2.shape}, scale {c:.6f}, "
              f"max rel diff {rel(E2, c * model):.1e}")

        # T4 R3 compares its model with the very same array (only the model's window differs)
        E3 = engine_data(r2 / "data_dp.hdf5", 0)
        check("T4 R3 (no upsampling) sees the same 256 px data as R2", rel(E3, E2) < 1e-6, f"{rel(E3, E2):.1e}")

        # T5 R4: full engine and its presolve (guard: +-96 mrad = 256 detector px = 1024 fine px). Its truth is the full
        # data cropped 1420 -> 1416 in the engine frame (the bin-4 crop keeps 708 = 177 x 4 px a side, not 710)
        E4f = engine_data(full, 0, presolve=1416)
        E4 = engine_data(r4 / "data_dp.hdf5", 2)
        check("T5 R4 data == engine binning of its 1416 px truth", rel(E4, 16 * binning_2D(E4f, 4)) < 1e-6,
              f"{rel(E4, 16 * binning_2D(E4f, 4)):.1e}")
        E4p = engine_data(r4 / "data_dp.hdf5", 2, presolve=1024)
        check("T5 R4 presolve (1416 -> 1024 fine, then binned) == R2's data", rel(E4p, E2) < 1e-6, f"{rel(E4p, E2):.1e}")

        # T5b the other cameras of the limits wave (2026-10-26): 512 px (2 x 2, k = 1) and 128 px (8 x 8, k = 3), same truth
        Pf = loadmat(r2 / "probe_initial_true.mat")["probe"].astype(np.complex128)   # R2's fine probe (= R1's, T8)
        for rr, B, k in ((r512, 2, 1), (r128, 8, 3)):
            m = loadmat(rr / "sim_meta.mat", simplify_cells=True)["meta"]
            Ek = engine_data(rr / "data_dp.hdf5", k)
            mk = binning_2D(E1, B)
            ck = Ek.sum() / mk.sum()
            Pk = loadmat(rr / f"probe_initial_true_up{k}.mat")["probe"].astype(np.complex128)
            check(f"T5b {1024 // B} px ({B} x {B}, k = {k}): data == engine binning of the R1-frame truth, fine probe = R1's",
                  m["Ndpx"] == 1024 // B and Ek.shape[:2] == (1024 // B,) * 2 and abs(ck - B * B) < 1e-6
                  and rel(Ek, ck * mk) < 1e-6 and overlap(Pk, Pf) > 1 - 1e-9,
                  f"{m['Ndpx']} px, scale {ck:.4f}, max rel diff {rel(Ek, ck * mk):.1e}")

        # T6 the modulus constraint maps each detector pixel's ratio back onto the same fine block the model summed
        rng = np.random.default_rng(3)
        X = rng.random((1024, 1024, 1))
        U = unbinning_2D(binning_2D(X, 4), 4)
        blocks = X.reshape(256, 4, 256, 4).mean(axis=(1, 3))
        check("T6 unbinning_2D(binning_2D(x)) = block mean over the same 4 x 4 blocks (from pixel 0)",
              np.allclose(U[:, :, 0], np.kron(blocks, np.ones((4, 4)))))

        # T7 the test can fail: blocks offset by 1-3 fine pixels, or the engine's centered=true binning
        with h5py.File(full, "r") as f:
            dpf = f["dp"][...].astype(np.float64)
        for off in (1, 2, 3):
            s = N_IN // 2 - 512 - off                      # the 1024 crop, shifted by off
            bad = dpf[:, s:s + 1024, s:s + 1024].reshape(NPOS, 256, 4, 256, 4).sum(axis=(2, 4)).astype(np.float32)
            p = Path(tmp) / f"off{off}.hdf5"
            with h5py.File(p, "w") as f:
                f.create_dataset("dp", data=bad)
            Eb = engine_data(p, 2)
            check(f"T7 control: blocks offset by {off} fine px are caught", rel(Eb, 16 * model) > 0.05,
                  f"max rel diff {rel(Eb, 16 * model):.2f}")
        mc = binning_2D(E1, 4, centered=True)
        check("T7 control: the engine's centered=true binning would NOT match (it is not what the engine calls)",
              rel(E2, 16 * mc) > 0.05, f"{rel(E2, 16 * mc):.2f}")

        # T8 probes
        Pf = loadmat(r2 / "probe_initial_true.mat")["probe"].astype(np.complex128)   # files hold complex64
        P_eng = interpolateFT(P.astype(np.complex128), (1024, 1024))      # the reference presolve's resampling
        check("T8 fine probe (1024) == the engine's interpolateFT of the true probe", overlap(Pf, P_eng) > 1 - 1e-6,
              f"overlap {overlap(Pf, P_eng):.9f}")
        P1 = loadmat(r1 / "probe_initial_true.mat")["probe"].astype(np.complex128)
        check("T8 R1 and R2 get the same probe file content", overlap(P1, Pf) > 1 - 1e-9, f"{overlap(P1, Pf):.12f}")
        P0 = loadmat(r2 / "probe_initial_true_up0.mat")["probe"].astype(np.complex128)
        check("T8 R3 probe = the central 256 px of the fine probe", P0.shape == (256, 256)
              and overlap(P0, Pf[384:640, 384:640]) > 1 - 1e-9)
        old = crop_pad(interpolateFT(Pf, (256, 256)), (1024, 1024))       # the unfixed loader: resample, then zero-pad
        check("T8 why ptycho_model_probe needed the fix: the old path shrinks the probe 4x",
              overlap(old, Pf) < 0.5, f"overlap of the old loader's probe with the true one {overlap(old, Pf):.3f}")

    print(f"\n{'ALL PASSED' if not FAILS else str(len(FAILS)) + ' FAILED: ' + '; '.join(FAILS)}")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
