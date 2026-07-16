#!/usr/bin/env python
"""3-D point-spread function estimation.

Every atom reconstructs as the system PSF: tight in-plane (~0.5-1 A) but stretched
along the beam (~2 A+, the missing cone). So  recon ~ true_atoms (*) PSF.

Three ways to get the PSF, in order of trustworthiness:
  1. EMPIRICAL from the sim (gold standard, owned by the sims thread): reconstruct a
     single isolated atom through the identical pipeline; the blob IS the PSF. Drops
     in via cfg.single_atom_vol -> empirical_psf(). Nothing else changes.
  2. DATA-DERIVED (default stand-in): average many isolated Pb columns' blobs from THIS
     recon. Real, anisotropic, data-matched. Limitation: the z-tails beyond +-period/2
     are contaminated by the neighbouring in-column Pb (a crystal has no isolated atom
     in z), so we taper them -> the sim PSF (1) is what fixes the true axial tails.
  3. SYNTHETIC (physics stand-in): an anisotropic kernel from the optics (in-plane
     Gaussian, elongated missing-cone z-profile). Controlled + reproducible, assumed.

The user chose "both, compare": we build (2) as default and (3) alongside, and the
pipeline reports detection sensitivity to the PSF choice. All kernels share the SAME
voxel grid so they are drop-in interchangeable.
"""
from __future__ import annotations
import os
import numpy as np
from scipy.ndimage import shift as nd_shift

from . import align as _align


# ---------------------------------------------------------------- grid
def kernel_shape(cfg, dx):
    hz = int(round(cfg.psf_half_z_A / cfg.dz))
    hxy = int(round(cfg.psf_half_xy_A / dx))
    return hz, hxy


def _axes(cfg, dx, hz, hxy):
    z = np.arange(-hz, hz + 1) * cfg.dz
    xy = np.arange(-hxy, hxy + 1) * dx
    return z, xy


# ---------------------------------------------------------------- synthetic (physics)
def synthetic_psf(cfg, dx):
    """Anisotropic Gaussian with an asymmetric missing-cone z-tail. Normalised to sum 1."""
    hz, hxy = kernel_shape(cfg, dx)
    z, xy = _axes(cfg, dx, hz, hxy)
    sxy = cfg.psf_inplane_fwhm_A / 2.3548
    sz = cfg.psf_axial_fwhm_A / 2.3548
    Z, Y, X = np.meshgrid(z, xy, xy, indexing="ij")
    core = np.exp(-0.5 * ((X / sxy) ** 2 + (Y / sxy) ** 2 + (Z / sz) ** 2))
    # missing cone -> extra elongation along z (broad exponential skirt in z only)
    skirt = cfg.psf_missingcone_tail * np.exp(-np.abs(Z) / (1.5 * sz)) * \
        np.exp(-0.5 * ((X / sxy) ** 2 + (Y / sxy) ** 2))
    k = core + skirt
    k /= k.sum()
    return k.astype(float)


# ---------------------------------------------------------------- data-derived
def _isolated_pb_atoms(pos, Z, al, cfg, hz):
    """Individual in-field Pb atoms usable as PSF samples.

    The extraction box (+-psf_half_xy_A in-plane, +-hz layers) is small enough that
    the nearest neighbouring heavy COLUMN (~2.76 A away) never enters it, and the
    in-column Pb neighbours (+-3.9 A) sit at/just outside the z-edge (and are tapered).
    We still apply a light guard: skip a Pb atom if ANOTHER heavy atom lands inside the
    box (covers distorted / domain-wall regions). Returns list of (X, Y, z)."""
    win = _align.in_window(pos, cfg)
    pb = pos[win & (Z == 82)]
    heavy = pos[win & ((Z == 82) | (Z == 22))]
    guard_xy = cfg.psf_half_xy_A + 0.3
    guard_z = cfg.psf_data_zhalf_A          # only reject neighbours inside the TRUSTED
    #                                         central window; farther ones are tapered
    out = []
    for a in pb:
        zr = al.SGN * a[2] + al.OFF
        if zr >= cfg.zmax_show_A or zr <= 0:
            continue
        dxy = np.hypot(heavy[:, 0] - a[0], heavy[:, 1] - a[1])
        dzz = np.abs(heavy[:, 2] - a[2])
        near = (dxy < guard_xy) & (dzz < guard_z) & (dxy > 1e-6)
        if near.any():
            continue
        out.append(a)
    return out


def data_psf(V, dx, al, pos, Z, cfg, recenter=True):
    """Average isolated-Pb blobs into an empirical 3-D kernel on the standard grid."""
    hz, hxy = kernel_shape(cfg, dx)
    nL, Ny, Nx = V.shape
    zw = cfg.psf_data_zhalf_A                       # trusted z half-window (neighbour-free)
    acc = np.zeros((2*hz+1, 2*hxy+1, 2*hxy+1))
    n = 0
    for a in _isolated_pb_atoms(pos, Z, al, cfg, hz):
        rf, cf, lf = al.site_to_index(a[0], a[1], a[2])
        r, c, l = int(round(rf)), int(round(cf)), int(round(lf))
        if not (hz <= l < nL-hz and hxy <= r < Ny-hxy and hxy <= c < Nx-hxy):
            continue
        patch = V[l-hz:l+hz+1, r-hxy:r+hxy+1, c-hxy:c+hxy+1].copy()
        patch -= patch.min()
        if recenter and patch.max() > 0:
            # sub-pixel recentre on the PEAK's local 3x3x3 neighbourhood only (so the
            # tapered z-tails / any residual neighbour do not bias the shift)
            pz, py, px = np.unravel_index(np.argmax(patch), patch.shape)
            z0, z1 = max(pz-1, 0), min(pz+2, patch.shape[0])
            y0, y1 = max(py-1, 0), min(py+2, patch.shape[1])
            x0, x1 = max(px-1, 0), min(px+2, patch.shape[2])
            w = patch[z0:z1, y0:y1, x0:x1]
            gz, gy, gx = np.indices(w.shape)
            cz = z0 + (gz*w).sum()/w.sum(); cy = y0 + (gy*w).sum()/w.sum(); cx = x0 + (gx*w).sum()/w.sum()
            patch = nd_shift(patch, (hz-cz, hxy-cy, hxy-cx), order=1, mode="nearest")
        acc += patch
        n += 1
    if n == 0:
        raise RuntimeError("no isolated Pb blobs found for the data PSF")
    k = acc / n
    # taper the (neighbour-contaminated) z-tails beyond the trusted window
    zax = np.arange(-hz, hz+1) * cfg.dz
    taper = np.ones_like(zax)
    out = np.abs(zax) > zw
    taper[out] = np.clip(1 - (np.abs(zax[out]) - zw) / (cfg.psf_half_z_A - zw + 1e-9), 0, 1)
    k *= taper[:, None, None]
    k = np.clip(k, 0, None)
    k /= k.sum()
    return k.astype(float), n


# ---------------------------------------------------------------- empirical (sim, gold)
def empirical_psf(cfg, dx, vol_path=None):
    """Load a reconstructed single-isolated-atom volume and extract its blob as the PSF.

    vol_path defaults to cfg.single_atom_vol (the gold Pb kernel); pass another path for a
    per-species kernel. The atom is at the scan centre; we take the brightest voxel's
    neighbourhood on the standard grid. This is the DEFENSIBLE, measured system PSF
    (anisotropy + missing cone included). Unlike the data-derived stand-in there is NO
    neighbour to taper (a truly isolated atom), so the axial tails are kept -- that is the
    whole point of the sim kernel. Vacuum is clipped to zero rather than pedestal-shifted
    (per-layer median already centres it)."""
    hz, hxy = kernel_shape(cfg, dx)
    vol = np.load(os.path.expanduser(vol_path or cfg.single_atom_vol))
    V = np.angle(vol).astype(float)
    V -= np.median(V, axis=(1, 2), keepdims=True)
    l, r, c = np.unravel_index(np.argmax(V), V.shape)
    l = int(np.clip(l, hz, V.shape[0]-hz-1))
    r = int(np.clip(r, hxy, V.shape[1]-hxy-1))
    c = int(np.clip(c, hxy, V.shape[2]-hxy-1))
    k = V[l-hz:l+hz+1, r-hxy:r+hxy+1, c-hxy:c+hxy+1].copy()
    k = np.clip(k, 0, None)                       # vacuum -> 0, keep the real blob shape
    k /= k.sum()
    return k.astype(float)


def species_kernels(cfg, dx):
    """Per-species measured kernels: {82: K_Pb, 22: K_Ti} (both normalised, same grid).

    Measured (this data): the Ti single-atom axial response is genuinely broader than Pb's
    (weaker channeling of the lighter atom; 3-D corr 0.87) -- the shape difference is the
    species-classification signal. O's own kernel is below the recon noise floor
    (PSF_SIM_RESPONSE sec.2) so O shares the Ti (light-atom) shape."""
    out = {}
    if cfg.single_atom_vol and os.path.exists(os.path.expanduser(cfg.single_atom_vol)):
        out[82] = empirical_psf(cfg, dx)
    if cfg.ti_kernel_vol and os.path.exists(os.path.expanduser(cfg.ti_kernel_vol)):
        out[22] = empirical_psf(cfg, dx, vol_path=cfg.ti_kernel_vol)
    return out


# ---------------------------------------------------------------- effective 1-D axial kernel
def axial_kernel(k3d, navg=1):
    """The effective 1-D axial response: what a column z-profile (extracted by a (2*navg+1)^2
    in-plane MEAN about the blob centre) sees from a unit atom. This is the operator the
    finder's spike deconvolution must invert -- NOT the centre-voxel line, which for the tight
    gold kernel is a near-delta. Returns (k1d, zoff_layers) with zoff in integer layer units
    relative to the kernel centre; k1d is the mean over the central 3x3 (navg=1) in-plane box.
    Not renormalised: amplitude stays in the same integrated-phase units as fit.py betas."""
    hz = (k3d.shape[0]-1)//2
    hxy = (k3d.shape[1]-1)//2
    k1d = k3d[:, hxy-navg:hxy+navg+1, hxy-navg:hxy+navg+1].mean((1, 2))
    zoff = np.arange(-hz, hz+1)
    return k1d, zoff


# ---------------------------------------------------------------- FWHM report
def measure_fwhm(k, cfg, dx, navg=1):
    """FWHM (A) along z (from the 3x3-integrated axial response -- the meaningful blur the
    finder inverts) and in-plane (through the peak layer). navg sets the axial in-plane box."""
    hz = (k.shape[0]-1)//2; hxy = (k.shape[1]-1)//2
    zline = k[:, hxy-navg:hxy+navg+1, hxy-navg:hxy+navg+1].mean((1, 2))
    xline = k[hz, hxy, :]

    def fwhm(line, step):
        pk = line.max()
        if pk <= 0:
            return np.nan
        above = np.where(line >= pk/2)[0]
        return (above[-1]-above[0]) * step if len(above) > 1 else step
    return fwhm(zline, cfg.dz), fwhm(xline, dx)


# ---------------------------------------------------------------- dispatcher
def build_psfs(V, dx, al, pos, Z, cfg):
    """Return {name: kernel}. 'empirical' present iff cfg.single_atom_vol is set."""
    out = {}
    synth = synthetic_psf(cfg, dx)
    out["synthetic"] = synth
    try:
        data, ndat = data_psf(V, dx, al, pos, Z, cfg)
        out["data"] = data
        out["_data_n"] = ndat
    except RuntimeError as e:
        out["_data_err"] = str(e)
    if cfg.single_atom_vol:
        out["empirical"] = empirical_psf(cfg, dx)
    # default: empirical > data > synthetic (best available measured PSF)
    out["_default"] = "empirical" if "empirical" in out else ("data" if "data" in out else "synthetic")
    return out
