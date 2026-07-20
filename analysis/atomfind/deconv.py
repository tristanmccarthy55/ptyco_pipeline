#!/usr/bin/env python
"""Richardson-Lucy deconvolution of the reconstructed phase volume.

RL is the Poisson-appropriate, non-negative deconvolution -> the right choice for
electron-counting data, and it CANNOT push the object negative (physical for a phase
that is >= 0 above vacuum). We deconvolve by the measured 3-D system PSF to sharpen the
missing-cone axial blur (the z-crutch): adjacent atoms ~4 layers apart overlap through
the tight kernel into a weakly-modulated streak, and RL / the spike fit pull them apart.

Preprocessing matters (this is what fixed the "solid-colour" figure):
  * per-layer median already subtracted (align.load_phase) -> vacuum ~ 0;
  * CLIP negatives to 0 (do NOT global-min pedestal-shift: that injects a ~0.36 background
    ~3x the atom contrast that RL then redistributes into a flat wash);
  * TRIM to interior layers (cfg.trim_z_A): the entrance/exit "dumping-ground" planes carry
    huge residuals that RL amplifies (seen: interior max ~2 -> 14.5). Deconvolve the interior,
    re-embed into a full-size volume (artifact planes zeroed) so downstream indexing is intact;
  * crop the kernel to its compact in-plane support (it is ~0.1 A / 2 px wide but stored 57 px)
    -> far faster and no ringing from the empty margins.

HONESTY: deconvolution amplifies noise. Iterations are capped (more when noiseless, scaled
down at low dose via cfg.effective_rl_iters), a filter_epsilon floors the divisor, and a
per-iteration DIVERGENCE GUARD stops early if the estimate's peak/median blows past
cfg.rl_blowup_ratio x the input's. The quantitative detector is the model-based spike fit
(find.py / fit.py); RL is the complementary image-space view + a deconvolved volume.
"""
from __future__ import annotations
import numpy as np
from scipy.signal import fftconvolve


def crop_kernel_inplane(psf, frac=0.01, cap=10):
    """Crop a 3-D kernel to its compact in-plane support (keep full z), renormalise."""
    hz = (psf.shape[0]-1)//2; hxy = (psf.shape[1]-1)//2
    proj = psf.max(0)
    yy, xx = np.where(proj >= frac * proj.max())
    rad = int(np.ceil(max(np.abs(yy-hxy).max(), np.abs(xx-hxy).max()))) if len(yy) else 2
    rad = min(max(rad, 2), hxy, cap)
    k = psf[:, hxy-rad:hxy+rad+1, hxy-rad:hxy+rad+1]
    k = np.clip(k, 0, None)
    return k / k.sum()


def _rl_guarded(image, psf, iters, eps, blowup_ratio):
    """3-D Richardson-Lucy with a per-iteration divergence guard.

    Healthy deconvolution CONCENTRATES a blob into a tighter peak (intensity conserved),
    so peak value legitimately grows a few-fold -- a peak/median metric would false-trigger
    on that. We guard on runaway only: stop if the estimate's max exceeds blowup_ratio x the
    input's max (true noise amplification), which sharpening alone never reaches."""
    psf_m = psf[::-1, ::-1, ::-1]
    est = np.full(image.shape, max(image.mean(), eps), dtype=float)
    cap = blowup_ratio * float(image.max())
    done, blew = 0, False
    for i in range(iters):
        conv = np.clip(fftconvolve(est, psf, mode="same"), eps, None)
        est = est * fftconvolve(image / conv, psf_m, mode="same")
        est = np.clip(est, 0, None)
        done = i + 1
        if float(est.max()) > cap:
            blew = True
            break
    return est, done, blew


def mem_3d(V, psf, cfg, iters=200, alpha=0.5, tol=1e-4):
    """Maximum-entropy 3-D deconvolution (Gull-Daniell-type multiplicative-exponential
    iteration), the second routine of the field-standard STEM depth-sectioning workflow
    (Ishizuka et al., Microscopy 70, 241 (2021), where MEM outperformed RL).

        f <- f * exp( alpha * K^T (g - K (*) f) )

    Positivity is enforced by construction; the entropy prior pulls unconstrained voxels
    toward flat. Same interior-trim + re-embed treatment as richardson_lucy_3d. Stops on
    chi^2 stagnation (rel. improvement < tol over 5 checks) or divergence."""
    lo, hi = cfg.trim_z_A
    l0 = max(int(round(lo / cfg.dz)), 0)
    l1 = min(int(round(hi / cfg.dz)), V.shape[0])
    interior = np.clip(V[l0:l1], 0, None)
    scale = np.percentile(interior, 99.9)
    scale = scale if scale > 0 else 1.0
    g = interior / scale
    k = crop_kernel_inplane(psf)
    km = k[::-1, ::-1, ::-1]
    f = np.full(g.shape, max(g.mean(), 1e-6))
    chi_prev = np.inf
    stall = 0
    done = 0
    for i in range(iters):
        r = g - fftconvolve(f, k, mode="same")
        chi = float(np.mean(r * r))
        if chi > chi_prev * 1.5:                          # divergence guard
            break
        if chi_prev - chi < tol * chi_prev:
            stall += 1
            if stall >= 5:
                done = i + 1
                break
        else:
            stall = 0
        chi_prev = chi
        step = np.clip(alpha * fftconvolve(r, km, mode="same"), -2.0, 2.0)
        f = np.clip(f * np.exp(step), 0, None)
        done = i + 1
    out = np.zeros_like(V, dtype=float)
    out[l0:l1] = f * scale
    info = dict(iters_done=done, chi2=chi_prev, trim_layers=[l0, l1], alpha=alpha,
                backend="mem-gull-daniell")
    return out.astype(np.float32), info


def richardson_lucy_3d(V, psf, cfg):
    """Deconvolve the interior of phase volume V by psf. Returns (full-size volume, info)."""
    lo, hi = cfg.trim_z_A
    l0 = max(int(round(lo / cfg.dz)), 0)
    l1 = min(int(round(hi / cfg.dz)), V.shape[0])
    interior = np.clip(V[l0:l1], 0, None)                # vacuum -> 0, keep contrast
    scale = np.percentile(interior, 99.9)
    scale = scale if scale > 0 else 1.0
    img = interior / scale
    k = crop_kernel_inplane(psf)
    iters = cfg.effective_rl_iters()
    dec_int, done, blew = _rl_guarded(img, k, iters, cfg.rl_filter_epsilon, cfg.rl_blowup_ratio)
    dec = np.zeros_like(V, dtype=float)
    dec[l0:l1] = dec_int * scale
    info = dict(iters=iters, iters_done=done, blewup=blew, trim_layers=[l0, l1],
                kernel_inplane=k.shape[1], backend="local-guarded", dose=cfg.dose_e_per_A2)
    return dec.astype(np.float32), info
