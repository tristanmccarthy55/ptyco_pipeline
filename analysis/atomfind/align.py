#!/usr/bin/env python
"""@file align.py
@brief Recon <-> ground-truth alignment (calibration only; the finder never sees this map).

Reproduces the validated map from analysis/column_cross_section_overlay.py:
  in-plane : recon (row r, col c) -> GT physical  X = X0 + c*dx ,  Y = Y0 + r*dx
             (recon is the transpose of the GT grid; sub-pixel CAL from blob centroids)
  depth    : z_recon = SGN * z_GT + OFF          (SGN, OFF fitted data-drivenly)
Reference columns are picked automatically from the brightest GT Pb (A-site) columns, so it
runs unchanged on finer/larger volumes. The map must be affine per axis, not a shift --
see METHODS.md ("Calibration infrastructure").
"""
from __future__ import annotations
from dataclasses import dataclass
import os

import numpy as np


# ---------------------------------------------------------------- loading
def load_recon_mat(path):
    """@brief Read a PtychoShelves reconstruction directly -> complex64 (nL, Ny, Nx).

    `outputs.object_roi` is a MATLAB cell array of per-layer complex arrays written
    column-major, so each layer is transposed on the way in. `object_roi` is used rather than
    `object`, which carries an unreliable probe-overlap halo around the reconstructed region.
    Complex data appears either as a compound ('real','imag') dtype or as a native complex,
    depending on how the .mat was written; both are handled.
    """
    import h5py
    with h5py.File(path, "r") as f:
        refs = f["outputs"]["object_roi"][:, 0]
        layers = []
        for r in refs:
            a = f[r][:]
            if a.dtype.names:                       # compound ('real','imag')
                a = a["real"] + 1j * a["imag"]
            layers.append(a.T)                      # column-major -> row-major
    return np.asarray(layers).astype(np.complex64)


def load_recon_h5(path):
    """@brief Read a newer PtychoShelves `*_recons.h5` -> complex64 (nL, Ny, Nx).

    The current engine writes `reconstruction/object` (NL,1,Ny,Nx) complex instead of the old
    `outputs.object_roi` cell array; there is no separate ROI dataset, so the full object is
    returned (extract_psf.py reads the same field). Kept identical to that reader so the empirical
    PSFs and the labyrinth object come off the same code path.
    """
    import h5py
    with h5py.File(path, "r") as f:
        obj = f["reconstruction"]["object"][:]          # (NL,1,Ny,Nx) or (NL,Ny,Nx) complex
    obj = np.asarray(obj).squeeze()
    if obj.ndim == 2:
        obj = obj[None]
    return obj.astype(np.complex64)


def load_object(path):
    """@brief Load a reconstructed object from whatever the reconstruction actually produced.

    The peer-facing entry point: a PtychoShelves run emits `Niter<N>.mat`, so that is accepted
    directly rather than requiring a pre-converted array. Accepts

      *.mat  an older PtychoShelves reconstruction (outputs.object_roi)
      *.h5   a newer PtychoShelves `*_recons.h5` (reconstruction/object)
      *.npy  a complex (nL, Ny, Nx) array already extracted from one

    @return complex array (nL, Ny, Nx).
    """
    ext = os.path.splitext(path)[1].lower()
    if ext == ".mat":
        return load_recon_mat(path)
    if ext == ".h5":
        return load_recon_h5(path)
    if ext == ".npy":
        return np.load(path)
    raise ValueError(f"unsupported reconstruction format {ext!r} for {path!r}; "
                     f"expected a PtychoShelves .mat / *_recons.h5 or an extracted .npy")


def load_phase(cfg):
    """Load the reconstructed object -> per-layer mean-subtracted phase volume.

    Accepts a PtychoShelves `.mat` or an extracted `.npy` (see load_object). Returns V
    (nL,Ny,Nx) float, and the effective dx (A/px). dx is derived as scan_window/Nx when
    cfg.dx is None (matches the validated overlay, generalises).
    """
    vol = load_object(cfg.recon_vol)
    V = np.angle(vol).astype(float)
    V -= np.median(V, axis=(1, 2), keepdims=True)     # remove per-layer offset
    dx = cfg.dx if cfg.dx is not None else cfg.scan_window_A / V.shape[2]
    return V, dx


GT_CACHE = "gt_prepared.npz"


def _prepare_gt(vasp_path):
    """The sim's exact preparation of the reference structure: rotate so the beam is +z,
    orthogonalize, pad to a square cell, centre, add vacuum. Needs ase + abtem."""
    import ase.io, abtem
    a = ase.io.read(vasp_path)
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


def _prepare_gt_thin(vasp_path, n_cells, z_vacuum=4.0, cell_z=3.905):
    """[thin-ab] Ground truth for the THIN aberration-campaign slab.

    Mirrors sim/simulate_4dstem.py::build_thin_sample(): same rotate/orthogonalize/square-pad
    as _prepare_gt, then crop the beam (z) axis to a central slab n_cells cells thick and
    re-centre with z_vacuum padding each side. This MUST track build_thin_sample exactly --
    the full-box GT (18 cells, 70 A) does not describe a 5-cell 27.5 A sim, and validating
    against the wrong slab silently produces meaningless recall / z-RMS.
    """
    import ase.io, abtem
    a = ase.io.read(vasp_path)
    a.rotate(-90, "y", rotate_cell=True)
    a = abtem.orthogonalize_cell(a)
    Lx, Ly, _ = a.cell.lengths()
    s = max(Lx, Ly)
    a.cell[0, 0] = s
    a.cell[1, 1] = s
    a.center(axis=0)
    a.center(axis=1)
    thick = n_cells * cell_z
    z = a.get_positions()[:, 2]
    zc = 0.5 * (z.min() + z.max())
    a = a[np.abs(z - zc) <= thick / 2.0]
    a.center(axis=2, vacuum=z_vacuum)
    return a.get_positions(), a.get_atomic_numbers()


def load_gt(cfg):
    """@brief Ground-truth atoms in the RECON physical frame (identical prep to the sim).

    Prefers the precomputed `gt_prepared.npz` cache if one is on the data path. The cache
    exists so a reproduction run needs only numpy: preparing the frame from the .vasp calls
    abtem.orthogonalize_cell, and abtem is a heavy dependency to install for one function.
    Build or refresh the cache with `python -m atomfind.make_gt_cache`.

    @return (pos, Z) -- positions (N,3) in angstrom and atomic numbers (N,).
    """
    from . import config as _cfg
    cache = _cfg.data_path(GT_CACHE)
    if os.path.exists(cache):
        d = np.load(cache)
        return d["pos"], d["Z"]
    return _prepare_gt(cfg.vasp())


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
    mZ: float = 0.0       # depth scale residual (layer per layer); OFF carries the shift

    def site_to_index(self, X, Y, Z_depth):
        """GT physical (X, Y, z) -> fractional recon index (row, col, layer)."""
        col_n = (np.asarray(X) - self.X0) / self.dx
        row_n = (np.asarray(Y) - self.Y0) / self.dx
        col = col_n * (1.0 + self.mX) + self.bX
        row = row_n * (1.0 + self.mY) + self.bY
        layer_n = (self.SGN * np.asarray(Z_depth) + self.OFF) / self.dz - 0.5
        layer = layer_n * (1.0 + self.mZ)
        return row, col, layer

    def index_to_site(self, row, col):
        """Recon (row, col) -> GT physical (X, Y) (inverse in-plane map)."""
        col_n = (np.asarray(col) - self.bX) / (1.0 + self.mX)
        row_n = (np.asarray(row) - self.bY) / (1.0 + self.mY)
        return self.X0 + col_n * self.dx, self.Y0 + row_n * self.dx

    def layer_to_z(self, layer):
        """Recon fractional layer -> GT physical z (inverse depth map)."""
        layer_n = np.asarray(layer) / (1.0 + self.mZ)
        return self.SGN * ((layer_n + 0.5) * self.dz - self.OFF)


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


def resolve_origin(V, dx, cfg):
    """@brief Fill in cfg.X0/cfg.Y0 when they are None: the object is CENTRED on the scan centre.

    X0/Y0 place the object's (0,0) pixel in GT coordinates. The historical constants (30, 10)
    are the scan-window corner, which is only the object origin when the object exactly covers
    the scan (dx = scan_window/Nx). A full-field reconstruction is larger than the scan by the
    probe halo -- e.g. a 20 A scan reconstructed onto a 37 A (753 px) or 54.6 A (1109 px) object
    -- so the corner is wrong by (N*dx - scan_window)/2, which is NOT a lattice vector and maps
    every atom onto the wrong site. Deriving it from the object size fixes that for any box.
    Verified against p/illum_sum (the illumination centroid sits at the object centre).
    """
    cx, cy = cfg.scan_center_xy
    if cfg.X0 is None:
        cfg.X0 = cx - 0.5 * V.shape[2] * dx
    if cfg.Y0 is None:
        cfg.Y0 = cy - 0.5 * V.shape[1] * dx
    return cfg


def crop_to_fov(V, dx, cfg):
    """@brief Crop the analysis to the illuminated scan region (cfg.fov_A). Returns the new V.

    A full-field reconstruction is much larger than the scan: the object extends by the probe
    halo (a 20 A scan lands on a 37-55 A object here). Only the scanned area is constrained by
    data -- the halo is reconstructed from partial illumination and is full of junk, so running
    the finder over it costs precision (measured: 11-31%) for detections that were never
    measurable. Standard practice is to analyse the scan field only.

    Crops a cfg.fov_A box centred on cfg.scan_center_xy and shifts X0/Y0 to the new corner, so
    the recon<->GT map still refers to the same physical frame. None => no crop.
    """
    resolve_origin(V, dx, cfg)
    if not cfg.fov_A:
        return V
    cx, cy = cfg.scan_center_xy
    half = 0.5 * cfg.fov_A
    c0 = int(round((cx - half - cfg.X0) / dx)); c1 = int(round((cx + half - cfg.X0) / dx))
    r0 = int(round((cy - half - cfg.Y0) / dx)); r1 = int(round((cy + half - cfg.Y0) / dx))
    c0, r0 = max(c0, 0), max(r0, 0)
    c1, r1 = min(c1, V.shape[2]), min(r1, V.shape[1])
    if c1 - c0 < 8 or r1 - r0 < 8:
        raise RuntimeError(f"fov_A={cfg.fov_A} leaves no usable region in a "
                           f"{V.shape[1]}x{V.shape[2]} object at dx={dx}")
    cfg.X0 += c0 * dx
    cfg.Y0 += r0 * dx
    print(f"[fov] cropped to {cfg.fov_A} A scan field: rows {r0}:{r1}, cols {c0}:{c1} "
          f"({V.shape[1]}x{V.shape[2]} -> {r1-r0}x{c1-c0} px)")
    return V[:, r0:r1, c0:c1]


def _register_depth_atoms(V, dx, pos, Z, cfg, zrec, search_px=6):
    """@brief Depth offset that puts the most GT Pb atoms on reconstructed phase maxima.

    The comb registration correlates a column's depth profile with a comb of its GT depths, which
    fixes OFF only modulo the lattice period c -- and on the thin labyrinth it degenerates
    completely: _pb_columns' 0.5 A rounding splits columns whose atoms wander in-plane (polar
    displacements, median 0.24 A) into ONE-atom fragments, so each comb has one tooth and fits
    equally well on any atom along the column (a90/a100 registered at OFF = -c, with 11% of the
    found atoms then mapped outside the physical slab).

    Here every in-field GT Pb atom is scored at its OWN (x, y, z+OFF): the per-layer maximum of
    the phase within +-search_px of its nominal pixel (tolerant of the pre-refinement in-plane
    offset), linearly interpolated in depth. Atoms mapped outside the atomic band (trim_z_A) score
    zero, so the vacuum-band surface artefacts cannot vote, and a one-cell shift that pushes a
    whole end plane out of the slab loses ~1/N_planes of the score. No column grouping at all.
    @return (OFF, SGN, mean score)
    """
    from scipy.ndimage import maximum_filter
    nL, ny, nx = V.shape
    P = pos[in_window(pos, cfg) & (Z == 82)]
    c = np.round((P[:, 0] - cfg.X0) / dx).astype(int)
    r = np.round((P[:, 1] - cfg.Y0) / dx).astype(int)
    ok = (r >= search_px) & (r < ny - search_px) & (c >= search_px) & (c < nx - search_px)
    P, r, c = P[ok], r[ok], c[ok]
    if len(P) == 0:
        raise RuntimeError("no in-field Pb atoms for depth registration")
    prof = maximum_filter(V, size=(1, 2 * search_px + 1, 2 * search_px + 1))[:, r, c].T   # (n, nL)
    lo_b, hi_b = cfg.trim_z_A
    best = None
    for sgn, lo, hi in cfg.depth_branches:
        for off in np.linspace(lo, hi, 481):
            zm = sgn * P[:, 2] + off                                  # recon-frame depth
            k = np.interp(zm, zrec, np.arange(nL))                    # fractional layer
            k0 = np.clip(np.floor(k).astype(int), 0, nL - 2); w = k - k0
            val = (1 - w) * prof[np.arange(len(P)), k0] + w * prof[np.arange(len(P)), k0 + 1]
            val[(zm < lo_b) | (zm > hi_b)] = 0.0
            score = float(val.mean())
            if best is None or score > best[0]:
                best = (score, off, sgn)
    return best[1], best[2], best[0]


def register(V, dx, pos, Z, cfg, n_ref=6):
    """Fit (SGN, OFF, CAL_X, CAL_Y) from the brightest GT Pb columns. Returns Alignment."""
    resolve_origin(V, dx, cfg)
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

    # ---- depth registration
    if getattr(cfg, "depth_register", "comb") == "atoms":
        OFF, SGN, corr = _register_depth_atoms(V, dx, pos, Z, cfg, zrec)
    else:
        # joint comb correlation over the reference columns
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


def _robust_line(p, o, min_keep=10):
    """Least-squares line with one 3-sigma rejection pass. Returns (slope, intercept)."""
    p, o = np.asarray(p, float), np.asarray(o, float)
    m, b = np.polyfit(p, o, 1)
    r = o - (m*p + b)
    keep = np.abs(r - np.median(r)) < 3*np.std(r) + 1e-9
    if keep.sum() >= min_keep:
        m, b = np.polyfit(p[keep], o[keep], 1)
    return float(m), float(b)


def refine_with_atoms(al: Alignment, found, pos, Z, cfg, iters=4, fit_depth_scale=True):
    """Refine the recon<->GT map (in-plane affine AND depth) on matched HEAVY atoms.

    Two systematic residuals motivate this. (1) In-plane: blob-peak calibration differs
    from the kernel-fit centre by ~1 px for asymmetric blobs, plus a ~0.6% pixel-scale
    residual. (2) DEPTH: the comb-correlation registration can sit ~1 A off, and on a
    lattice whose Ti/O alternate every 1.95 A a 1 A depth error SWAPS SPECIES LABELS
    wholesale (measured: dose1e10 confusion 19% with OFF=1.26 A, 14% at OFF=0.20 A).

    Must ITERATE: the correspondences used to diagnose the offset are themselves computed
    with the (wrong) offset. The z-gate is kept below the heavy-atom column spacing
    (~3.9 A) so heavy->heavy matches cannot jump a lattice site. Refines the CALIBRATION
    only -- the finder never sees the map, so blindness is untouched.

    Returns a new Alignment (OFF absorbs the depth shift so it stays interpretable;
    mZ carries any depth-scale residual)."""
    from dataclasses import replace as _replace
    win = in_window(pos, cfg)
    heavy = np.where(((Z == 82) | (Z == 22)) & win)[0]
    if len(heavy) < 20:
        return al
    # Fiducials chosen by AMPLITUDE, not by species label: on a poorly-registered volume
    # the labels are exactly what is wrong (a mislabelled O sits 1.95 A off and drags the
    # depth fit), so selecting on `species != 8` would feed the error back into itself.
    # The brightest atoms are Pb/Ti whatever they are currently called.
    qual_ok = found["quality"] >= 0.5
    if qual_ok.sum() < 20:
        return al
    amp_cut = np.median(found["amp"][qual_ok])
    sel = np.where(qual_ok & (found["amp"] >= amp_cut))[0]
    if len(sel) < 20:
        return al
    f_row = found["row"][sel]; f_col = found["col"][sel]; f_lay = found["layer"][sel]
    nom_c = (pos[heavy, 0] - al.X0) / al.dx        # nominal (uncorrected) GT indices
    nom_r = (pos[heavy, 1] - al.Y0) / al.dx

    gate_xy = 1.0                                   # A, in-plane: columns are ~2.76 A apart
    gate_z = 1.8                                    # A, < half the 3.9 A heavy spacing
    cur = al
    for _ in range(iters):
        gr, gc, gl = cur.site_to_index(pos[heavy, 0], pos[heavy, 1], pos[heavy, 2])
        mc, mr, ml, mi = [], [], [], []
        for k in range(len(sel)):
            dxy = np.hypot((gr - f_row[k])*cur.dx, (gc - f_col[k])*cur.dx)
            dz = np.abs(gl - f_lay[k]) * cur.dz
            ok = (dxy <= gate_xy) & (dz <= gate_z)
            if not ok.any():
                continue
            cand = np.where(ok)[0]
            j = cand[np.argmin(dxy[cand] + dz[cand])]
            mc.append(f_col[k] - gc[j]); mr.append(f_row[k] - gr[j])
            ml.append(f_lay[k] - gl[j]); mi.append(j)
        if len(mi) < 20:
            break
        mi = np.asarray(mi)
        # ---- depth: shift (into OFF) + optional scale (into mZ) ----
        gl_m = gl[mi]
        if fit_depth_scale and (np.ptp(gl_m) > 10):
            sZ, bZ = _robust_line(gl_m, np.asarray(ml))
        else:
            sZ, bZ = 0.0, float(np.median(ml))
        sZ = float(np.clip(sZ, -0.05, 0.05))        # guard: >5% depth-scale error is a bug
        # layer = layer_n*(1+mZ); absorb the constant into OFF (d layer / d OFF = 1/dz)
        new_mZ = (1.0 + cur.mZ) * (1.0 + sZ) - 1.0
        new_OFF = cur.OFF + bZ * cur.dz / max(1.0 + new_mZ, 1e-6)
        # ---- in-plane affine (composed, as before) ----
        m2c, b2c = _robust_line(nom_c[mi], np.asarray(mc))
        m2r, b2r = _robust_line(nom_r[mi], np.asarray(mr))
        cur = _replace(cur, OFF=float(new_OFF), mZ=float(new_mZ),
                       mX=cur.mX + m2c, bX=cur.bX + b2c,
                       mY=cur.mY + m2r, bY=cur.bY + b2r)
    ctr = 0.5 * (np.max(nom_c) + np.min(nom_c))
    cur.CAL_X = float((cur.mX*ctr + cur.bX) * al.dx)
    cur.CAL_Y = float((cur.mY*ctr + cur.bY) * al.dx)
    return cur


def summarize(al: Alignment):
    return (f"depth: z_recon = {al.SGN:+d}*z_GT + {al.OFF:.2f} A (corr {al.corr_depth:.2f}, "
            f"mZ {al.mZ*100:+.2f}%)  |  "
            f"in-plane affine: centre offset X {al.CAL_X:+.3f} A ({al.CAL_X/al.dx:+.1f} px), "
            f"Y {al.CAL_Y:+.3f} A;  scale residual mX {al.mX*100:+.2f}% mY {al.mY*100:+.2f}%"
            f"  |  dx={al.dx:.4f} A")
