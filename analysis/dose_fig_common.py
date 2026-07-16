#!/usr/bin/env python
"""Shared helpers for the dose-series publication figures.

One place for: the canonical .mat -> phase-volume reader (with a fast .npy cache),
the reference-dose column picker, the ground-truth model load + per-dose depth/in-plane
registration, the kz depth-resolution spectrum, and a clean-publication matplotlib style.

All of this is adapted from the existing analysis scripts so the numbers stay consistent:
  - loader / column pick        : dose_compare.py
  - GT model + registration     : column_cross_section_overlay.py
  - kz spectrum / plane freqs   : depth_resolution.py

Interpreter: ~/hyperspy-bundle/bin/python  (h5py, ase, abtem, scipy).

Self-test:  python dose_fig_common.py --selftest
"""
import os, glob, json, sys
import numpy as np

# ----------------------------------------------------------------- geometry / paths
DATA_ROOT = os.path.expanduser("~/Desktop/dose_series")
FIG_DIR   = os.path.join(DATA_ROOT, "figures")
CACHE_DIR = os.path.join(FIG_DIR, "cache")
VASP      = ("/Users/u2109287/Library/CloudStorage/OneDrive-UniversityofWarwick/"
            "Documents/PtycoShelves/ptychoshelves-clean/sim/PTO6_STO6_18_18_labyrinthPoscar.vasp")

DOSES     = ["1e10", "1e8", "1e6", "1e4"]      # high -> low
DOSE_TeX  = {"1e10": r"$10^{10}$", "1e8": r"$10^{8}$",
             "1e6": r"$10^{6}$", "1e4": r"$10^{4}$"}     # e/Å²
DX_FALLBACK = 0.04916                           # in-plane Å/px (handover)
X0, Y0      = 30.0, 10.0                         # scan-window origin (Å); recon(r,c)->GT
ZMAX_SHOW   = 66.0                               # trim exit-surface artifact (last ~4 Å)

for _d in (FIG_DIR, CACHE_DIR):
    os.makedirs(_d, exist_ok=True)

# per-species marker style for GT overlays (from column_cross_section_overlay.py)
SPECIES = {82: dict(c="#00e5ff", m="o", label="Pb"),   # cyan
           38: dict(c="#ffffff", m="s", label="Sr"),   # white
           22: dict(c="#39ff14", m="x", label="Ti"),   # neon green
           8:  dict(c="#ff35ff", m="D", label="O")}    # magenta


# ----------------------------------------------------------------- publication style
def use_pub_style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "figure.facecolor":  "white",
        "savefig.facecolor": "white",
        "savefig.dpi":       300,
        "savefig.bbox":      "tight",
        "font.family":       "sans-serif",
        "font.sans-serif":   ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size":         10,
        "axes.titlesize":    11,
        "axes.labelsize":    10,
        "axes.linewidth":    0.8,
        "axes.facecolor":    "white",
        "xtick.direction":   "out",
        "ytick.direction":   "out",
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.frameon":    False,
        "legend.fontsize":   8.5,
        "image.interpolation": "nearest",
    })


def savefig(fig, stem, pdf=True):
    """Write <stem>.png (300 dpi) and optionally <stem>.pdf into FIG_DIR."""
    png = os.path.join(FIG_DIR, stem + ".png")
    fig.savefig(png, dpi=300)
    if pdf:
        fig.savefig(os.path.join(FIG_DIR, stem + ".pdf"))
    print("wrote", png + (" (+pdf)" if pdf else ""))
    return png


# ----------------------------------------------------------------- .mat -> phase cache
def _mat_path(dose):
    cand = (glob.glob(os.path.join(DATA_ROOT, f"dose{dose}", "Niter*.mat"))
            or glob.glob(os.path.join(DATA_ROOT, f"*dose{dose}[!0-9]*", "**", "Niter*.mat"),
                         recursive=True))
    if not cand:
        return None
    return sorted(cand, key=lambda p: int("".join(filter(str.isdigit, os.path.basename(p)))))[-1]


def _read_mat(dose):
    """Canonical reader (dose_compare.load_vol): -> phase V[nL,Ny,Nx], dz(Å), dx(Å)."""
    import h5py
    m = _mat_path(dose)
    if m is None:
        raise FileNotFoundError(f"no Niter*.mat for dose {dose} under {DATA_ROOT}")
    with h5py.File(m, "r") as f:
        g = f["outputs"]
        layers = []
        for r in g["object_roi"][:, 0]:            # array of HDF5 object refs
            a = f[r][:]
            a = (a["real"] + 1j * a["imag"]) if a.dtype.names else a
            layers.append(a.T)                     # undo MATLAB column-major
        V = np.angle(np.array(layers)).astype(np.float32)
        V -= np.median(V, (1, 2), keepdims=True)   # per-layer background
        z = g["z_distance"][:, 0]
        dz = float(np.median(z[np.isfinite(z)])) * 1e10
        try:
            dx = float(g["pixel_size"][0, 0]) * 1e10
        except Exception:
            dx = DX_FALLBACK
    return V, dz, dx


_META = os.path.join(CACHE_DIR, "meta.json")


def load_dose(dose, verbose=True):
    """Return (V[nL,Ny,Nx] float32 phase, dz Å, dx Å); cache to CACHE_DIR on first read."""
    npy = os.path.join(CACHE_DIR, f"dose{dose}_phase.npy")
    meta = json.load(open(_META)) if os.path.exists(_META) else {}
    if os.path.exists(npy) and dose in meta:
        V = np.load(npy)
        return V, meta[dose]["dz"], meta[dose]["dx"]
    if verbose:
        print(f"[cache miss] reading dose {dose} .mat (~1.1 GB) ...", flush=True)
    V, dz, dx = _read_mat(dose)
    np.save(npy, V)
    meta[dose] = {"dz": dz, "dx": dx}
    json.dump(meta, open(_META, "w"), indent=0)
    return V, dz, dx


# ----------------------------------------------------------------- column picking
def pick_columns(ref_V, n=1, margin=40, sep=25, pct=97):
    """Strong interior Pb columns on the reference volume (dose_compare.py:49-53).

    Returns a list of (yc, xc), brightest first. n=1 -> just the hero Pb column.
    """
    from scipy.ndimage import maximum_filter
    dm = ref_V.mean(0); dmn = dm - dm.min()
    mask = np.zeros_like(dmn, bool); mask[margin:-margin, margin:-margin] = True
    pk = (dmn == maximum_filter(dmn, sep)) & (dmn > np.percentile(dmn[mask], pct)) & mask
    ys, xs = np.where(pk)
    order = np.argsort(-dmn[ys, xs])
    return [(int(ys[j]), int(xs[j])) for j in order[:n]]


# ----------------------------------------------------------------- ground-truth model
_GT_CACHE = os.path.join(CACHE_DIR, "gt_model.npz")


def load_gt_model():
    """Orientation-corrected GT atoms -> (pos[N,3] Å, Z[N]).  (overlay script:31-34)."""
    if os.path.exists(_GT_CACHE):
        d = np.load(_GT_CACHE)
        return d["pos"], d["Z"]
    import ase.io, abtem
    a = ase.io.read(VASP); a.rotate(-90, "y", rotate_cell=True); a = abtem.orthogonalize_cell(a)
    Lx, Ly, Lz = a.cell.lengths(); s = max(Lx, Ly)
    a.cell[0, 0] = s; a.cell[1, 1] = s
    a.center(axis=0); a.center(axis=1); a.center(axis=2, vacuum=2.0)
    pos = a.get_positions(); Z = a.get_atomic_numbers()
    np.savez(_GT_CACHE, pos=pos, Z=Z)
    return pos, Z


def column_atoms(pos, Z, yc, xc, dx, dxwin=1.3, dywin=0.7):
    """GT atoms sitting in the recon column (yc,xc) -> (Xc, Yc, P[k,3], Zk)."""
    Xc = X0 + xc * dx; Yc = Y0 + yc * dx
    sel = (np.abs(pos[:, 0] - Xc) < dxwin) & (np.abs(pos[:, 1] - Yc) < dywin)
    return Xc, Yc, pos[sel], Z[sel]


# ----------------------------------------------------------------- registration
def inplane_shift(ref, dx, pos, Z, species=(82, 38, 22), zslab=(20, 90)):
    """Residual in-plane shift of the recon vs the GT column map -> (dX0, dY0) in Å.

    Cross-correlates a depth-projected recon slab against a synthetic Gaussian map of
    the GT metal columns at the nominal X0,Y0.  For this sim it is ~0 (the nominal
    origin is right), but it self-corrects any per-dose integer-pixel offset.
    """
    from scipy.ndimage import gaussian_filter
    from skimage.registration import phase_cross_correlation
    Ny, Nx = ref.shape[1:]
    proj = ref[zslab[0]:zslab[1]].mean(0); proj = proj - proj.min()
    sel = np.isin(Z, species)
    ci = np.round((pos[sel, 0] - X0) / dx).astype(int)
    ri = np.round((pos[sel, 1] - Y0) / dx).astype(int)
    ok = (ci >= 0) & (ci < Nx) & (ri >= 0) & (ri < Ny)
    gt = np.zeros((Ny, Nx)); np.add.at(gt, (ri[ok], ci[ok]), 1.0)
    gt = gaussian_filter(gt, 2.0)
    shift, _, _ = phase_cross_correlation(proj, gt, upsample_factor=10, normalization=None)
    return -shift[1] * dx, -shift[0] * dx        # (dX0, dY0)


def slice_atoms(pos, Z, layer, dz, dx, OFF, z_tol=None, dX0=0.0, dY0=0.0):
    """GT atoms in the in-plane slice at recon `layer` -> dict{species: (cols, rows)}.

    recon(row r, col c) -> GT X=X0+c*dx, Y=Y0+r*dx ; z_GT = (layer+0.5)*dz - OFF.
    Returns pixel coords so they overlay directly on imshow of V[layer].
    """
    if z_tol is None:
        z_tol = dz
    z_gt = (layer + 0.5) * dz - OFF
    sel = np.abs(pos[:, 2] - z_gt) < z_tol
    P, Zs = pos[sel], Z[sel]
    out = {}
    for zz in SPECIES:
        m = Zs == zz
        if m.any():
            c = (P[m, 0] - (X0 + dX0)) / dx
            r = (P[m, 1] - (Y0 + dY0)) / dx
            out[zz] = (c, r)
    return out


def register(V, dz, dx, columns, pos, Z, ref_calx=None):
    """Per-dose depth (SGN,OFF) + in-plane sub-pixel (CAL_X) registration.

    Depth: fit one sign+offset for the whole volume by correlating the recon column
    profiles against the GT heavy-atom (Pb/Ti/Sr) comb over physical windows only
    (kills 3.9 Å periodic aliases) — column_cross_section_overlay.py:44-60.
    In-plane: median blob offset of the resolved planes vs GT — overlay:62-75.
    If ref_calx is given it is reused instead of re-fitting (robust at low dose).
    """
    nL = V.shape[0]; zrec = (np.arange(nL) + 0.5) * dz

    def profile(yc, xc):
        p = V[:, yc-1:yc+2, xc-1:xc+2].mean((1, 2)); return p - p.mean()

    heavy, profs = [], []
    for (yc, xc) in columns:
        _, _, P, Zs = column_atoms(pos, Z, yc, xc, dx)
        hz = P[(Zs == 82) | (Zs == 22) | (Zs == 38), 2]
        if hz.size:
            heavy.append(hz); profs.append(profile(yc, xc))

    def comb(zatoms, off, sgn):
        g = np.zeros(nL)
        for za in zatoms:
            g += np.exp(-0.5 * ((zrec - (sgn * za + off)) / 0.7) ** 2)
        return g - g.mean()

    best = None
    for sgn, lo, hi in [(+1, -8, 4), (-1, 66, 78)]:
        for off in np.linspace(lo, hi, 400):
            c = sum(np.dot(p, comb(z, off, sgn)) /
                    (np.linalg.norm(p) * np.linalg.norm(comb(z, off, sgn)) + 1e-9)
                    for p, z in zip(profs, heavy))
            if best is None or c > best[0]:
                best = (c, off, sgn)
    corr, OFF, SGN = best

    if ref_calx is not None:
        return SGN, OFF, ref_calx, corr

    W = 22
    offs = []
    for (yc, xc) in columns:
        cs = V[:, yc-1:yc+2, xc-W:xc+W].mean(1)
        _, _, P, Zs = column_atoms(pos, Z, yc, xc, dx)
        xax = (np.arange(2 * W) - W) * dx
        sel = (Zs == 82) | (Zs == 22)
        for X_at, z_at in zip(P[sel, 0], P[sel, 2]):
            lyr = int(round((SGN * z_at + OFF) / dz - 0.5))
            if 0 <= lyr < nL and (SGN * z_at + OFF) < ZMAX_SHOW:
                w = np.clip(cs[lyr] - cs[lyr].min(), 0, None)
                if w.sum() > 0:
                    offs.append((xax * w).sum() / w.sum() - (X_at - X0 - xc * dx))
    CAL_X = float(np.median(offs)) if offs else 0.0
    return SGN, OFF, CAL_X, corr


# ----------------------------------------------------------------- kz depth spectrum
def hero_depth_offset(V, dz, dx, yc, xc, pos, Z, W=22):
    """Robust depth offset for one column: match the recon plane-blobs (found on a
    drift-following max profile) to the GT Pb comb.  -> (OFF Å, median residual Å).

    More robust than the comb correlation for a single column: the max over the
    in-plane strip follows the column even where it drifts across a domain wall.
    """
    from scipy.signal import find_peaks
    nL = V.shape[0]; zrec = (np.arange(nL) + 0.5) * dz
    strip = V[:, yc-1:yc+2, xc-W:xc+W].mean(1)
    prof = strip.max(1); prof = prof - prof.min()
    pk, _ = find_peaks(prof, distance=max(1, int(2.6 / dz)),
                       prominence=prof.max() * 0.10)
    zblob = zrec[pk]
    _, _, P, Zs = column_atoms(pos, Z, yc, xc, dx)
    zpb = np.sort(P[Zs == 82, 2])
    if zblob.size == 0 or zpb.size == 0:
        return 0.0, np.nan
    offs = np.linspace(-6, 6, 481)
    r = [np.median([np.min(np.abs((z - o) - zpb)) for z in zblob]) for o in offs]
    j = int(np.argmin(r))
    return float(offs[j]), float(r[j])


def kz_spectrum(V, dz):
    """On-column vs vacuum z-power spectrum (depth_resolution.kz_spectrum)."""
    dm = V.mean(0); dmn = dm - dm.min()
    col = dmn > np.percentile(dmn, 95)
    vac = dmn < np.percentile(dmn, 35)
    nL = V.shape[0]
    win = np.hanning(nL)[:, None, None]
    Vc = (V - V.mean(0, keepdims=True)) * win
    P = np.abs(np.fft.rfft(Vc, axis=0)) ** 2
    kz = np.fft.rfftfreq(nL, d=dz)
    return kz, P[:, col].mean(1), P[:, vac].mean(1)


def gt_plane_freqs(pos, Z, cx=40.0, cy=20.0, w=20.0):
    """Dominant along-beam plane spacing for Pb and Ti (depth_resolution.gt_plane_freqs)."""
    win = (np.abs(pos[:, 0] - cx) < w / 2) & (np.abs(pos[:, 1] - cy) < w / 2)
    out = {}
    for nm, zz in [("Pb", 82), ("Ti", 22)]:
        zc = pos[win & (Z == zz), 2]
        zs = np.sort(np.unique(np.round(zc, 1)))
        gaps = np.diff(zs)
        out[nm] = np.median(gaps[gaps > 0.5]) if len(gaps) else np.nan
    return out


def plane_peak_prominence(kz, Pcol, f0, bw=0.03):
    """Depth-resolution scalar: on-column power at the plane freq / local baseline."""
    m = np.abs(kz - f0) < bw
    if not m.any():
        return np.nan
    base = np.median(Pcol[(kz > f0 - 0.12) & (kz < f0 + 0.12)]) + 1e-30
    return float(Pcol[m].max() / base)


# ----------------------------------------------------------------- self-test
def _selftest():
    print("=== dose_fig_common self-test ===")
    pos, Z = load_gt_model()
    print(f"GT atoms: {len(Z)}  species {sorted(set(Z.tolist()))}")
    gaps = gt_plane_freqs(pos, Z)
    f_pb = 1.0 / gaps["Pb"]
    print(f"GT plane spacing  Pb {gaps['Pb']:.2f} Å ({f_pb:.3f} Å⁻¹)  Ti {gaps['Ti']:.2f} Å")

    ref, dz0, dx0 = load_dose("1e10")
    cols = pick_columns(ref, n=3)
    print(f"reference dose 1e10: V{ref.shape} dz={dz0:.3f} dx={dx0:.4f}  columns={cols}")
    _, _, calx_ref, _ = register(ref, dz0, dx0, cols, pos, Z)
    print(f"reference CAL_X = {calx_ref:+.3f} Å")

    for d in DOSES:
        V, dz, dx = load_dose(d)
        SGN, OFF, CALX, corr = register(V, dz, dx, cols, pos, Z, ref_calx=calx_ref)
        kz, Pcol, Pvac = kz_spectrum(V, dz)
        prom = plane_peak_prominence(kz, Pcol, f_pb)
        print(f"  dose {d:>4}: dz={dz:.3f}  reg z={SGN:+d}*zGT+{OFF:5.2f} (corr {corr:.2f})"
              f"  CAL_X={CALX:+.3f}  Pb-plane prominence={prom:5.1f}x")
    print("cache dir:", CACHE_DIR)


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        _selftest()
    else:
        print(__doc__)
