#!/usr/bin/env python
"""@file fit.py
@brief Model-based per-site amplitude fitting -- the (GT-seeded) oxygen detector.

In-plane resolution (~0.1 A) already separates Ti and O columns; oxygen is hard only from
axial OVERLAP (O sits ~1.9 A from Ti along the beam, below the ~2 A axial limit) and
CONTRAST (Z=8). This does not try to RESOLVE O from Ti -- it ATTRIBUTES the signal to known
lattice sites: place the measured PSF at every known site in an in-plane column group and
solve non-negative least squares for a per-site amplitude beta. beta is the detection
statistic; beta_Pb:beta_Ti:beta_O is a calibrated amplitude-vs-Z curve. Per-group solves
stay tiny; a constant-background basis absorbs the vacuum pedestal. Contrast numbers:
RESULTS.md.
"""
from __future__ import annotations
import warnings

import numpy as np
from scipy.optimize import nnls
from scipy.ndimage import shift as nd_shift

from . import align as _align


# ---------------------------------------------------------------- compact kernel
def compact_kernel(psf, dx, cfg, xy_cap_px=4, thresh=0.01):
    """Crop the PSF in-plane to its ACTUAL support (keep full z) for cheap rendering.

    Support = smallest radius holding all in-plane voxels above `thresh` x peak, capped
    at xy_cap_px. Measured from the kernel itself, so the tight (~0.1 A) data PSF gets a
    small support while a broader synthetic one gets a proportionally larger one."""
    hz = (psf.shape[0]-1)//2
    hxy = (psf.shape[1]-1)//2
    inplane = psf.max(axis=0)                       # peak-projected in-plane profile
    yy, xx = np.where(inplane >= thresh * inplane.max())
    rad = int(np.ceil(max(np.abs(yy-hxy).max(), np.abs(xx-hxy).max()))) if len(yy) else 2
    supp = min(max(2, rad), hxy, xy_cap_px)
    k = psf[:, hxy-supp:hxy+supp+1, hxy-supp:hxy+supp+1].copy()
    k = np.clip(k, 0, None)
    k /= k.sum()
    return k, hz, supp


def _nnls_gram(A, b, ridge=0.0):
    """NNLS solved via the normal equations -> reduces an (m x n) problem to (n x n).

    min_{x>=0} ||A x - b||  has a solution depending only on G=A'A and f=A'b. With
    G = R'R (Cholesky), ||Ax-b||^2 = ||R x - R^{-T} f||^2 + const, so nnls(R, R^{-T}f)
    gives the same x in n rows instead of m -- turns a 40k-row solve into a 60-row one."""
    # Floating-point flags are suppressed across this whole solve. Both the Gram products
    # here and scipy's pure-Python NNLS (< 1.15) raise spurious 'overflow'/'invalid value'/
    # 'divide by zero' from matmul on the near-collinear designs this problem generates:
    # ~26k RuntimeWarnings over a full run on numpy 2.0 + scipy 1.13, and NONE on
    # numpy 1.26 + scipy 1.15, for output that agrees to a relative 1e-11 across the whole
    # of report.json with no non-finite value in any export. Suppressing is therefore safe,
    # but it is paired with an explicit finiteness check below so a REAL numerical failure
    # still raises instead of being swallowed by the filter.
    with warnings.catch_warnings(), np.errstate(all="ignore"):
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        G = A.T @ A
        f = A.T @ b
        n = G.shape[0]
        G = G + (ridge + 1e-9 * (np.trace(G) / n + 1e-12)) * np.eye(n)
        R = np.linalg.cholesky(G).T                 # upper-triangular
        y = np.linalg.solve(R.T, f)
        x, _ = nnls(R, y, maxiter=10 * n)
    if not np.all(np.isfinite(x)):
        raise FloatingPointError("NNLS returned a non-finite amplitude vector")
    return x


def _render(box_shape, center, kernel, hz, hxy):
    """Render one unit-amplitude PSF centred at fractional `center` into a zero box."""
    out = np.zeros(box_shape)
    lc, rc, cc = center
    il, ir, ic = int(round(lc)), int(round(rc)), int(round(cc))
    frac = (lc-il, rc-ir, cc-ic)
    k = nd_shift(kernel, frac, order=1, mode="constant", cval=0.0)
    z0, z1 = il-hz, il+hz+1
    y0, y1 = ir-hxy, ir+hxy+1
    x0, x1 = ic-hxy, ic+hxy+1
    kz0 = max(0, -z0); ky0 = max(0, -y0); kx0 = max(0, -x0)
    z0c, y0c, x0c = max(0, z0), max(0, y0), max(0, x0)
    z1c = min(box_shape[0], z1); y1c = min(box_shape[1], y1); x1c = min(box_shape[2], x1)
    if z1c <= z0c or y1c <= y0c or x1c <= x0c:
        return out.ravel()
    out[z0c:z1c, y0c:y1c, x0c:x1c] = k[kz0:kz0+(z1c-z0c), ky0:ky0+(y1c-y0c), kx0:kx0+(x1c-x0c)]
    return out.ravel()


# ---------------------------------------------------------------- the fit
def fit_amplitudes(V, dx, al, pos, Z, psf, cfg, extra_sites=None):
    """Windowed NNLS of the measured PSF at every known site -> per-site amplitudes.

    The field is partitioned into in-plane TILES (robust + reproducible, no clustering
    pathology). Each tile is fit jointly over an EXTENDED box (core + PSF-support margin,
    full depth); each atom's beta is recorded from the tile whose CORE contains it, so
    every atom is fit exactly once. `extra_sites` (Mx3 physical X,Y,z) are added as basis
    functions for the NULL test (their betas are returned tagged Z=-1).

    Returns a record array: idx, Z, X, Y, zGT, row, col, layer, beta, raw, n_basis.
    """
    nL, Ny, Nx = V.shape
    Vnn = V - np.percentile(V, 1.0)               # vacuum baseline -> ~non-negative
    kern, hz, hxy = compact_kernel(psf, dx, cfg)
    margin = hxy + 2
    step = max(8, int(round(cfg.fit_tile_A / dx)))

    win = _align.in_window(pos, cfg)
    idx_all = np.where(win)[0]
    r_all, c_all, l_all = al.site_to_index(pos[idx_all, 0], pos[idx_all, 1], pos[idx_all, 2])
    # optional null / extra basis sites (not recorded as GT atoms)
    if extra_sites is not None and len(extra_sites):
        er, ec, el = al.site_to_index(extra_sites[:, 0], extra_sites[:, 1], extra_sites[:, 2])
    else:
        er = ec = el = np.array([])

    r0f, r1f = int(np.floor(r_all.min())), int(np.ceil(r_all.max()))
    c0f, c1f = int(np.floor(c_all.min())), int(np.ceil(c_all.max()))

    recs = []
    for tr in range(r0f, r1f + 1, step):
        for tc in range(c0f, c1f + 1, step):
            r0, r1 = tr, min(tr + step, r1f + 1)          # tile core
            c0, c1 = tc, min(tc + step, c1f + 1)
            er0, er1 = max(r0 - margin, 0), min(r1 + margin, Ny)
            ec0, ec1 = max(c0 - margin, 0), min(c1 + margin, Nx)

            basis = ((r_all >= er0) & (r_all < er1) & (c_all >= ec0) & (c_all < ec1))
            core = ((r_all >= r0) & (r_all < r1) & (c_all >= c0) & (c_all < c1))
            if not basis.any():
                continue
            bi = np.where(basis)[0]
            # extra (null) sites in this extended box
            xbi = (np.where((er >= er0) & (er < er1) & (ec >= ec0) & (ec < ec1))[0]
                   if len(er) else np.array([], int))

            box = Vnn[:, er0:er1, ec0:ec1]
            bshape = box.shape
            b = box.ravel()
            cols = [_render(bshape, (l_all[i], r_all[i]-er0, c_all[i]-ec0), kern, hz, hxy) for i in bi]
            cols += [_render(bshape, (el[i], er[i]-er0, ec[i]-ec0), kern, hz, hxy) for i in xbi]
            Aatoms = np.column_stack(cols)
            # restrict the solve to voxels within the model support (PSF is tiny in-plane;
            # the vast majority of box voxels are exactly zero for every basis PSF). Keep a
            # sample of background voxels so the constant term is still identifiable.
            supp = Aatoms.any(axis=1)
            ns = int(supp.sum())
            if ns < 8:
                continue
            bg_pool = np.where(~supp)[0]
            if len(bg_pool):
                take = bg_pool[:: max(1, len(bg_pool) // ns)]      # ~ns background voxels
                keep = np.concatenate([np.where(supp)[0], take])
            else:
                keep = np.where(supp)[0]
            A = np.column_stack([Aatoms[keep], np.ones(len(keep))])
            beta = _nnls_gram(A, b[keep], ridge=cfg.fit_ridge)
            nb = len(bi) + len(xbi)
            # record GT atoms whose core-membership is this tile
            for k, i in enumerate(bi):
                if not core[i]:
                    continue
                j = idx_all[i]
                rr, ccl, ll = int(round(r_all[i])), int(round(c_all[i])), int(round(l_all[i]))
                raw = (V[ll, rr-1:rr+2, ccl-1:ccl+2].mean()
                       if 0 <= ll < nL and 1 <= rr < Ny-1 and 1 <= ccl < Nx-1 else np.nan)
                recs.append((j, int(Z[j]), pos[j, 0], pos[j, 1], pos[j, 2],
                             r_all[i], c_all[i], l_all[i], float(beta[k]), raw, nb))
            # record null sites (Z=-1) with core membership in this tile
            for k, i in enumerate(xbi):
                if not (r0 <= er[i] < r1 and c0 <= ec[i] < c1):
                    continue
                recs.append((-1, -1, extra_sites[i, 0], extra_sites[i, 1], extra_sites[i, 2],
                             er[i], ec[i], el[i], float(beta[len(bi)+k]), np.nan, nb))
    dt = np.dtype([("idx", int), ("Z", int), ("X", float), ("Y", float), ("zGT", float),
                   ("row", float), ("col", float), ("layer", float), ("beta", float),
                   ("raw", float), ("n_basis", int)])
    return np.array(recs, dtype=dt)
