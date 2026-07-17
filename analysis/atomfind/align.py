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
    """In-plane map is AFFINE per axis: nominal index -> measured index
        col_meas = col_nom*(1+mX) + bX ,   col_nom = (X - X0)/dx
    (a constant offset alone cannot register this data: dx = window/N differs from the
    physical object pixel by ~0.6%, which accumulates to ~2 px across the 404-px field).
    CAL_X/CAL_Y report the offset at the FIELD CENTRE (for humans/back-compat)."""
    dx: float
    dz: float
    X0: float
    Y0: float
    SGN: int
    OFF: float
    CAL_X: float          # centre-of-field offset (A), X/col axis (diagnostic)
    CAL_Y: float          # centre-of-field offset (A), Y/row axis (diagnostic)
    corr_depth: float     # joint depth-registration correlation (diagnostic)
    mX: float = 0.0       # affine scale residual, col axis (index units per index)
    bX: float = 0.0       # affine offset, col axis (px)
    mY: float = 0.0
    bY: float = 0.0

    def site_to_index(self, X, Y, Z_depth):
        """GT physical (X, Y, z) -> fractional recon index (row, col, layer)."""
        col_n = (np.asarray(X) - self.X0) / self.dx
        row_n = (np.asarray(Y) - self.Y0) / self.dx
        col = col_n * (1.0 + self.mX) + self.bX
        row = row_n * (1.0 + self.mY) + self.bY
        layer = (self.SGN * np.asarray(Z_depth) + self.OFF) / self.dz - 0.5
        return row, col, layer

    def index_to_site(self, row, col):
        """Recon (row, col) -> GT physical (X, Y) (inverse in-plane map)."""
        col_n = (np.asarray(col) - self.bX) / (1.0 + self.mX)
        row_n = (np.asarray(row) - self.bY) / (1.0 + self.mY)
        return self.X0 + col_n * self.dx, self.Y0 + row_n * self.dx

    def layer_to_z(self, layer):
        """Recon fractional layer -> GT physical z (inverse depth map)."""
        return self.SGN * ((np.asarray(layer) + 0.5) * self.dz - self.OFF)


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

    # ---- in-plane AFFINE calibration from parabolic PEAK offsets across MANY Pb columns.
    # Peak refinement is bias-free (a wide windowed centroid regresses toward zero offset);
    # fitting offset-vs-position as a line per axis captures both the sub-pixel shift AND
    # the ~0.6% scale residual (dx = window/N vs the physical object pixel) that
    # accumulates to ~2 px across the field.
    def _parab(ym, y0, yp):
        d = ym - 2*y0 + yp
        return 0.0 if abs(d) < 1e-12 else float(np.clip(0.5*(ym - yp)/d, -0.5, 0.5))

    def column_peak_offset(X, Y, r, c, zs, search=4):
        """Median (col_off, row_off) in px over this column's planes, or None."""
        offs_c, offs_r = [], []
        for za in zs:
            zr = SGN * za + OFF
            if zr >= cfg.zmax_show_A:
                continue
            lyr = int(round(zr / cfg.dz - 0.5))
            if not (0 <= lyr < nL) or not (search+1 < r < V.shape[1]-search-1
                                           and search+1 < c < V.shape[2]-search-1):
                continue
            box = V[lyr, r-search:r+search+1, c-search:c+search+1]
            pr, pc = np.unravel_index(np.argmax(box), box.shape)
            rr, cc = r - search + pr, c - search + pc
            if not (1 <= rr < V.shape[1]-1 and 1 <= cc < V.shape[2]-1):
                continue
            offs_c.append(cc + _parab(V[lyr, rr, cc-1], V[lyr, rr, cc], V[lyr, rr, cc+1])
                          - (X - cfg.X0) / dx)
            offs_r.append(rr + _parab(V[lyr, rr-1, cc], V[lyr, rr, cc], V[lyr, rr+1, cc])
                          - (Y - cfg.Y0) / dx)
        if len(offs_c) < 4:
            return None
        return float(np.median(offs_c)), float(np.median(offs_r))

    pos_c, off_c, pos_r, off_r = [], [], [], []
    for (_, X, Y, r, c, zs) in scored[:40]:              # many columns, spread over the field
        res = column_peak_offset(X, Y, r, c, zs)
        if res is None:
            continue
        pos_c.append((X - cfg.X0) / dx); off_c.append(res[0])
        pos_r.append((Y - cfg.Y0) / dx); off_r.append(res[1])

    def robust_line(p, o):
        p, o = np.asarray(p), np.asarray(o)
        if len(p) < 3:
            return 0.0, (float(np.median(o)) if len(o) else 0.0)
        m, b = np.polyfit(p, o, 1)
        res = o - (m*p + b)
        keep = np.abs(res - np.median(res)) < 3*np.std(res) + 1e-9
        if keep.sum() >= 3:
            m, b = np.polyfit(p[keep], o[keep], 1)
        return float(m), float(b)

    mX, bX = robust_line(pos_c, off_c)
    mY, bY = robust_line(pos_r, off_r)
    ctr = V.shape[1] / 2.0
    CAL_X = float((mX*ctr + bX) * dx)                    # centre-of-field offset (report)
    CAL_Y = float((mY*ctr + bY) * dx)

    return Alignment(dx=dx, dz=cfg.dz, X0=cfg.X0, Y0=cfg.Y0, SGN=int(SGN),
                     OFF=float(OFF), CAL_X=CAL_X, CAL_Y=CAL_Y, corr_depth=float(corr),
                     mX=mX, bX=bX, mY=mY, bY=bY)


def refine_with_atoms(al: Alignment, found, pos, Z, cfg):
    """Refine the in-plane affine map using matched HEAVY found atoms as fiducials.

    Blob-peak calibration leaves a systematic ~1 px residual because the kernel-fit centre
    (where found atoms -- and hence the honest comparison frame -- live) differs from the
    brightest-voxel peak for asymmetric blobs. This composes a residual affine, fitted on
    bright Pb/Ti matches, into the map. It refines the CALIBRATION only: the finder never
    sees the map, so blindness is untouched. Returns a new Alignment."""
    from dataclasses import replace as _replace
    heavy = np.where((Z == 82) | (Z == 22))[0]
    win = in_window(pos, cfg)
    heavy = heavy[win[heavy]]
    gr, gc, gl = al.site_to_index(pos[heavy, 0], pos[heavy, 1], pos[heavy, 2])
    sel = np.where((found["species"] != 8) & (found["quality"] >= 0.5))[0]
    res_c, nom_c, res_r, nom_r = [], [], [], []
    for i in sel:
        d2 = ((gr - found["row"][i])*al.dx)**2 + ((gc - found["col"][i])*al.dx)**2 \
             + ((gl - found["layer"][i])*al.dz)**2
        j = int(np.argmin(d2))
        if d2[j] > 1.0**2:
            continue
        # residual in px vs the NOMINAL (uncorrected) index of the matched GT atom
        col_n = (pos[heavy[j], 0] - al.X0) / al.dx
        row_n = (pos[heavy[j], 1] - al.Y0) / al.dx
        res_c.append(found["col"][i] - gc[j]); nom_c.append(col_n)
        res_r.append(found["row"][i] - gr[j]); nom_r.append(row_n)
    if len(res_c) < 20:
        return al

    def robust_line(p, o):
        p, o = np.asarray(p), np.asarray(o)
        m, b = np.polyfit(p, o, 1)
        r = o - (m*p + b)
        keep = np.abs(r - np.median(r)) < 3*np.std(r) + 1e-9
        if keep.sum() >= 10:
            m, b = np.polyfit(p[keep], o[keep], 1)
        return float(m), float(b)

    m2c, b2c = robust_line(nom_c, res_c)
    m2r, b2r = robust_line(nom_r, res_r)
    ctr = 0.5 * (np.max(nom_c) + np.min(nom_c))
    new = _replace(al, mX=al.mX + m2c, bX=al.bX + b2c, mY=al.mY + m2r, bY=al.bY + b2r)
    new.CAL_X = float(((new.mX)*ctr + new.bX) * al.dx)
    new.CAL_Y = float(((new.mY)*ctr + new.bY) * al.dx)
    return new


def summarize(al: Alignment):
    return (f"depth: z_recon = {al.SGN:+d}*z_GT + {al.OFF:.2f} A (corr {al.corr_depth:.2f})  |  "
            f"in-plane affine: centre offset X {al.CAL_X:+.3f} A ({al.CAL_X/al.dx:+.1f} px), "
            f"Y {al.CAL_Y:+.3f} A;  scale residual mX {al.mX*100:+.2f}% mY {al.mY*100:+.2f}%"
            f"  |  dx={al.dx:.4f} A")
