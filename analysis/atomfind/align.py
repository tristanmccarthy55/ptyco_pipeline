#!/usr/bin/env python
"""Recon <-> ground-truth alignment (the tested path, factored + generalised).

Reproduces the VALIDATED map from analysis/column_cross_section_overlay.py:
  in-plane : recon (row r, col c) -> GT physical  X = X0 + c*dx ,  Y = Y0 + r*dx
             (recon is the transpose of the GT grid; sub-pixel CAL from blob centroids)
  depth    : z_recon = SGN * z_GT + OFF          (SGN, OFF fitted data-drivenly)

Generalisation over the overlay script: reference columns are picked AUTOMATICALLY
from the brightest ground-truth Pb (A-site) columns instead of two hand-typed pixels,
so this runs unchanged on a bigger/finer volume. On NL70 it reproduces the overlay's
numbers (CAL ~ +1.7 px, depth offset ~ +0.48 A).
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


# ---------------------------------------------------------------- loading
def load_phase(cfg):
    """Load the reconstructed object -> per-layer mean-subtracted phase volume.

    Returns V (nL,Ny,Nx) float, and the effective dx (A/px). dx is derived as
    scan_window/Nx when cfg.dx is None (matches the validated overlay, generalises).
    """
    vol = np.load(cfg.recon_vol)
    V = np.angle(vol).astype(float)
    V -= np.median(V, axis=(1, 2), keepdims=True)     # remove per-layer offset
    dx = cfg.dx if cfg.dx is not None else cfg.scan_window_A / V.shape[2]
    return V, dx


def load_gt(cfg):
    """GT atoms in the RECON physical frame (identical prep to the sim). -> pos, Z."""
    import ase.io, abtem
    a = ase.io.read(cfg.vasp())
    a.rotate(-90, "y", rotate_cell=True)
    a = abtem.orthogonalize_cell(a)
    Lx, Ly, Lz = a.cell.lengths()
    s = max(Lx, Ly)
    a.cell[0, 0] = s
    a.cell[1, 1] = s
    a.center(axis=0)
    a.center(axis=1)
    a.center(axis=2, vacuum=2.0)
    return a.get_positions(), a.get_atomic_numbers()


def in_window(pos, cfg):
    """Boolean mask of atoms inside the scan field of view."""
    cx, cy = cfg.scan_center_xy
    h = cfg.scan_window_A / 2.0
    return ((np.abs(pos[:, 0] - cx) < h) & (np.abs(pos[:, 1] - cy) < h))


# ---------------------------------------------------------------- alignment object
@dataclass
class Alignment:
    dx: float
    dz: float
    X0: float
    Y0: float
    SGN: int
    OFF: float
    CAL_X: float          # in-plane sub-pixel calibration (A), applied to X (col)
    CAL_Y: float          # in-plane sub-pixel calibration (A), applied to Y (row)
    corr_depth: float     # joint depth-registration correlation (diagnostic)

    def site_to_index(self, X, Y, Z_depth):
        """GT physical (X, Y, z) -> fractional recon index (row, col, layer)."""
        col = (np.asarray(X) - self.X0 + self.CAL_X) / self.dx
        row = (np.asarray(Y) - self.Y0 + self.CAL_Y) / self.dx
        layer = (self.SGN * np.asarray(Z_depth) + self.OFF) / self.dz - 0.5
        return row, col, layer

    def index_to_site(self, row, col):
        """Recon (row, col) -> GT physical (X, Y) (inverse in-plane map)."""
        X = self.X0 + np.asarray(col) * self.dx - self.CAL_X
        Y = self.Y0 + np.asarray(row) * self.dx - self.CAL_Y
        return X, Y


# ---------------------------------------------------------------- registration
def _pb_columns(pos, Z, cfg):
    """Unique A-site (Pb) column (X, Y) centres inside the field, with their GT z-lists."""
    win = in_window(pos, cfg)
    pb = pos[win & (Z == 82)]
    # cluster by rounded (X, Y) -> one entry per column
    key = np.round(pb[:, :2] / 0.5) * 0.5
    cols = {}
    for (kx, ky), p in zip(map(tuple, key), pb):
        cols.setdefault((kx, ky), []).append(p[2])
    return [(kx, ky, np.sort(np.array(zs))) for (kx, ky), zs in cols.items()]


def register(V, dx, pos, Z, cfg, n_ref=6):
    """Fit (SGN, OFF, CAL_X, CAL_Y) from the brightest GT Pb columns. Returns Alignment."""
    nL = V.shape[0]
    zrec = (np.arange(nL) + 0.5) * cfg.dz
    dm = V.mean(0)                                    # depth-summed phase (bright at columns)

    # rank Pb columns by recon brightness at their nominal pixel
    cand = _pb_columns(pos, Z, cfg)
    scored = []
    for (X, Y, zs) in cand:
        c = int(round((X - cfg.X0) / dx))
        r = int(round((Y - cfg.Y0) / dx))
        if 1 <= r < V.shape[1] - 1 and 1 <= c < V.shape[2] - 1:
            scored.append((dm[r-1:r+2, c-1:c+2].mean(), X, Y, r, c, zs))
    scored.sort(reverse=True)
    ref = scored[:n_ref]
    if not ref:
        raise RuntimeError("no in-field Pb reference columns found for registration")

    # ---- depth registration: joint comb correlation over the reference columns
    def profile(r, c):
        p = V[:, r-1:r+2, c-1:c+2].mean((1, 2))
        return p - p.mean()

    def comb(zatoms, off, sgn):
        g = np.zeros(nL)
        for za in zatoms:
            g += np.exp(-0.5 * ((zrec - (sgn * za + off)) / 0.7) ** 2)
        return g - g.mean()

    profs = [(profile(r, c), zs) for (_, X, Y, r, c, zs) in ref]
    best = None
    for sgn, lo, hi in cfg.depth_branches:
        for off in np.linspace(lo, hi, 400):
            score = 0.0
            for p, zs in profs:
                g = comb(zs, off, sgn)
                score += float(np.dot(p, g) / (np.linalg.norm(p) * np.linalg.norm(g) + 1e-9))
            if best is None or score > best[0]:
                best = (score, off, sgn)
    corr, OFF, SGN = best
    corr /= len(profs)

    # ---- in-plane sub-pixel CAL: 1-D blob centroids on narrow cross-strips at
    # each Pb plane (matches the tested overlay method: row-averaged strip for the
    # X centroid, column-averaged strip for the Y centroid -> avoids the vertical
    # missing-cone tails biasing a 2-D square-window centroid).
    def centroid_offsets():
        dxs, dys = [], []
        W = int(round(cfg.psf_data_xyhalf_A / dx))
        ax = np.arange(-W, W + 1)
        for (_, X, Y, r, c, zs) in ref:
            for za in zs:
                zr = SGN * za + OFF
                if zr >= cfg.zmax_show_A:
                    continue
                lyr = int(round(zr / cfg.dz - 0.5))
                if not (0 <= lyr < nL) or not (W < r < V.shape[1]-W and W < c < V.shape[2]-W):
                    continue
                xstrip = V[lyr, r-1:r+2, c-W:c+W+1].mean(0)     # row-averaged -> X centroid
                ystrip = V[lyr, r-W:r+W+1, c-1:c+2].mean(1)     # col-averaged -> Y centroid
                wx = np.clip(xstrip - xstrip.min(), 0, None)
                wy = np.clip(ystrip - ystrip.min(), 0, None)
                if wx.sum() > 0:
                    cc = (ax * wx).sum() / wx.sum()
                    dxs.append((c + cc) - (X - cfg.X0) / dx)
                if wy.sum() > 0:
                    rr = (ax * wy).sum() / wy.sum()
                    dys.append((r + rr) - (Y - cfg.Y0) / dx)
        return np.array(dxs), np.array(dys)

    dxs, dys = centroid_offsets()
    CAL_X = float(np.median(dxs) * dx) if len(dxs) else 0.0
    CAL_Y = float(np.median(dys) * dx) if len(dys) else 0.0

    return Alignment(dx=dx, dz=cfg.dz, X0=cfg.X0, Y0=cfg.Y0, SGN=int(SGN),
                     OFF=float(OFF), CAL_X=CAL_X, CAL_Y=CAL_Y, corr_depth=float(corr))


def summarize(al: Alignment):
    return (f"depth: z_recon = {al.SGN:+d}*z_GT + {al.OFF:.2f} A (corr {al.corr_depth:.2f})  |  "
            f"in-plane CAL_X = {al.CAL_X:+.3f} A ({al.CAL_X/al.dx:+.1f} px), "
            f"CAL_Y = {al.CAL_Y:+.3f} A ({al.CAL_Y/al.dx:+.1f} px)  |  dx={al.dx:.4f} A")
