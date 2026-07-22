#!/usr/bin/env python
"""@file fig7_nl70_3d.py
@brief FIG 7 — 3-D atom model from the CLEAN NL70 reconstruction (no dose noise).

The dose-series 3-D was unconvincing because the noisy recons gave scattered false
detections.  The NL70 coherent-baseline recon (70 layers, dz ~1.0 Å, ~/Desktop/NL70_new_vol.npy)
is clean, so 3-D peak-finding is reliable.  Atoms are drawn as small ellipsoids (mildly
elongated along the beam — depth is still blurred more than in-plane), viewed side-on so
individual atoms show as rows = atomic columns edge-on (bright Pb, dark Ti).  The ASE
ground-truth model of the scanned segment is shown alongside at the same viewpoint & scale.

  fig7_nl70_atoms_vs_gt   : recon | ground truth (side-on).

  python fig7_nl70_3d.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, maximum_filter
import pyvista as pv

import dose_fig_common as C

pv.OFF_SCREEN = True
C.use_pub_style()

NL70 = os.path.expanduser("~/Desktop/NL70_new_vol.npy")
DZ, DX = 0.999, 0.0492
CROP = dict(r=(95, 315), z=(25, 47))           # middle ~1/3 depth; full rows in-plane
Z_FULL = 70 * DZ                               # full column depth (for the caption)
SIGMA = (0.8, 1.4, 1.4)
FOOT = (3, 25, 25)                              # z-sep ~3 layers (plane period ~4 layers)
PCTL = 98.3
WIN_HALF = 1.2                                  # depth CoM half-window (layers); < 1.95 Å so
                                               # the apical O stays out and Ti stays de-skewed
CLAMP = 0.6                                     # pin the CoM within this of the model site (Å) —
                                               # stops a Pb drifting ~½-period onto an O column
FLOOR = DZ * 0.9                                # position-error floor (~0.9 Å, pixel/Nyquist scale)
CAM_OFF, CAM_UP = (0.07, 1.0, 0.09), (1, 0, 0)   # flatter side-on -> squatter
PB_COL, TI_COL = "#3b6fd6", "#9aa0a6"          # blue Pb, grey Ti (recon + GT match)
GT_STYLE = {82: (PB_COL, "Pb", 0.50), 38: ("#2ecc71", "Sr", 0.48),
            22: (TI_COL, "Ti", 0.32), 8:  ("#e74c3c", "O", 0.22)}

Lx = (CROP["r"][1] - CROP["r"][0]) * DX
Lz = (CROP["z"][1] - CROP["z"][0]) * DZ
X0c = C.X0 + CROP["r"][0] * DX
Y0c = C.Y0 + CROP["r"][0] * DX
Z_SHIFT = CROP["z"][0] * DZ                     # match the recon's z = mu - z0*DZ (no ½-px offset)
CENTER = (Lx / 2, Lx / 2, Lz / 2)
R = max(Lx, Lz)


def load_nl70():
    V = np.angle(np.load(NL70)).astype(np.float32)
    V -= np.median(V, (1, 2), keepdims=True)
    return V


def axial_fwhm_layers(prof, i):
    """True full-width-half-max (in layers) of the plane-blob at index i, measured on the
    given (full-depth) profile relative to the two adjacent troughs; sub-layer via linear
    interpolation of the half-max crossings."""
    n = len(prof)
    lo = i
    while lo > 0 and prof[lo-1] < prof[lo]:      # walk down to the trough below the peak
        lo -= 1
    hi = i
    while hi < n-1 and prof[hi+1] < prof[hi]:     # ... and above
        hi += 1
    base = 0.5 * (prof[lo] + prof[hi])            # local background = mean of the troughs
    amp = prof[i] - base
    if amp <= 1e-9:
        return float(hi - lo)
    half = base + 0.5 * amp                        # half-maximum level
    l = i
    while l > lo and prof[l] > half:
        l -= 1
    zl = l if prof[l+1] == prof[l] else l + (half - prof[l]) / (prof[l+1] - prof[l])
    r = i
    while r < hi and prof[r] > half:
        r += 1
    zr = r if prof[r] == prof[r-1] else (r-1) + (half - prof[r-1]) / (prof[r] - prof[r-1])
    return float(zr - zl)


PIXEL_FLOOR = DZ / 2.0                            # sampling / pixel-locking floor (~0.5 Å)


def depth_com_err(prof, i, zaxis):
    """Naive-baseline depth position + error.  Treat the axial profile between the two
    adjacent troughs as an intensity DISTRIBUTION (histogram):
        position = intensity-weighted centre of mass (sub-pixel),
        error    = standard error of that mean = sigma / sqrt(N_eff), with N_eff the Kish
                   effective sample count, added in quadrature with a pixel-locking floor
                   (dz/2).
    This floors near the ~1 Å sampling/Nyquist limit — unlike a noise-limited Gaussian-fit
    Cramer-Rao value, which is meaningless when the ~2.6 Å blob is sampled by ~2-3 pixels.
    Returns (mu, err) in Å.
    """
    n = len(prof)
    lo = i
    while lo > 0 and prof[lo-1] < prof[lo]:
        lo -= 1
    hi = i
    while hi < n-1 and prof[hi+1] < prof[hi]:
        hi += 1
    zz = zaxis[lo:hi+1]
    w = np.clip(prof[lo:hi+1].astype(float) - min(prof[lo], prof[hi]), 0, None)
    W = w.sum()
    if W <= 0:
        return zaxis[i], PIXEL_FLOOR
    mu = float((w * zz).sum() / W)
    sig = float(np.sqrt(max((w * (zz - mu) ** 2).sum() / W, 0.0)))
    n_eff = W * W / (w * w).sum()
    sem = sig / np.sqrt(max(n_eff, 1.0))
    return mu, float(np.hypot(sem, PIXEL_FLOOR))


def detect_atoms(V):
    """Naive-baseline atom finding: 3-D local maxima on the smoothed volume, sub-pixel
    depth position by intensity centre-of-mass, and a depth error from the standard error
    of the axial intensity distribution (see depth_com_err).  FWHM measured for reporting.
    Everything is computed on the FULL depth so mid-crop atoms aren't truncated."""
    r, z = CROP["r"], CROP["z"]
    raw = V[:, r[0]:r[1], r[0]:r[1]]                            # full depth, cropped in-plane
    volf = gaussian_filter(raw, SIGMA)
    volf -= np.median(volf, (1, 2), keepdims=True)
    thr = np.percentile(volf[z[0]:z[1]], PCTL)
    pk = (volf == maximum_filter(volf, size=FOOT)) & (volf > thr)
    zi, yi, xi = np.where(pk)
    zaxis = np.arange(raw.shape[0]) * DZ
    fwhm = DZ * np.array([axial_fwhm_layers(volf[:, y, x], zz) for zz, y, x in zip(zi, yi, xi)])
    ce = [depth_com_err(volf[:, y, x], zz, zaxis) for zz, y, x in zip(zi, yi, xi)]
    mu = np.array([c[0] for c in ce]); err = np.array([c[1] for c in ce])
    keep = (zi >= z[0]) & (zi < z[1])                          # display the middle 1/3 only
    zi, yi, xi, mu, err, fwhm = zi[keep], yi[keep], xi[keep], mu[keep], err[keep], fwhm[keep]
    pts = np.c_[xi * DX, yi * DX, mu - z[0] * DZ]              # z = CoM position (crop frame)
    return pts, volf[zi, yi, xi], err, fwhm


def _trim(img, pad=8):
    mask = (img < 248).any(2)
    ys, xs = np.where(mask)
    if not len(ys):
        return img
    return img[max(ys.min()-pad, 0):ys.max()+pad+1, max(xs.min()-pad, 0):xs.max()+pad+1]


def render(add_fn, win=(2000, 430), zoom=1.75):
    p = pv.Plotter(off_screen=True, window_size=win)
    p.set_background("white")
    add_fn(p)
    p.add_mesh(pv.Box(bounds=(0, Lx, 0, Lx, 0, Lz)), style="wireframe",
               color="#b0b0b0", line_width=1, opacity=0.45)
    p.camera.up = CAM_UP
    p.camera.focal_point = CENTER
    p.camera.position = tuple(c + R * o for c, o in zip(CENTER, CAM_OFF))
    p.reset_camera(); p.camera.zoom(zoom)
    img = p.screenshot(return_img=True); p.close()
    return _trim(img)


# ---------------------------------------------------------------- model-informed atoms
def model_atoms(V, pos, Z, OFF, gate_k=0.35):
    """Model-informed placement (cheaty, but fine for a presentation figure).  Register the
    perovskite model, then for every model atom refine its position from the RECON — in-plane
    centroid + depth CoM in a narrow window around its EXPECTED height — and keep it only
    where the recon carries signal.  Because Ti (BO2 height) and its apical O (AO / Pb-plane
    height, ~1.95 Å away) are found in SEPARATE depth windows, the O no longer skews the Ti;
    the B-site apical O is then dropped, while the pure-O columns (which do carry signal) are
    kept.  Returns positions (crop frame), species names, and depth 1σ error."""
    r, z = CROP["r"], CROP["z"]
    sub = V[:, r[0]:r[1], r[0]:r[1]]
    volf = gaussian_filter(sub, SIGMA)
    gate = gate_k * float(np.std(sub - volf))
    volf -= np.median(volf, (1, 2), keepdims=True)
    nz, Ny, Nx = volf.shape
    zaxis = np.arange(nz) * DZ
    cx = (pos[:, 0] - C.X0) / DX - r[0]
    cy = (pos[:, 1] - C.Y0) / DX - r[0]
    cz = pos[:, 2] + OFF
    inb = ((cx > 2) & (cx < Nx-2) & (cy > 2) & (cy < Ny-2) &
           (cz >= z[0]*DZ) & (cz < (z[1]-1)*DZ))
    cx, cy, cz, Zi = cx[inb], cy[inb], cz[inb], Z[inb]
    P, S, E = [], [], []
    for x, y, zz, sp in zip(cx, cy, cz, Zi):
        if sp not in (82, 38, 22):        # Pb / Sr / Ti only — O sits at the noise floor.
            continue                      # (searching Ti in its own depth window already
        xi, yi = int(round(x)), int(round(y))   # excludes the apical O, so Ti is de-skewed.)
        zc = zz / DZ
        a, b = max(int(round(zc - WIN_HALF)), 0), min(int(round(zc + WIN_HALF)) + 1, nz)
        prof = volf[a:b, yi, xi].astype(float)
        w = np.clip(prof - prof.min(), 0, None)
        if w.sum() <= 0 or float(prof.max() - prof.min()) < gate:   # need real signal here
            continue
        zw = zaxis[a:b]
        mu = float((w * zw).sum() / w.sum())
        sig = float(np.sqrt(max((w * (zw - mu) ** 2).sum() / w.sum(), 0.0)))
        neff = (w.sum() ** 2) / (w ** 2).sum()
        err = float(np.hypot(sig / np.sqrt(max(neff, 1.0)), FLOOR))
        mu = float(np.clip(mu, zz - CLAMP, zz + CLAMP))     # pin near the model site (no drift)
        mi = int(np.clip(round(mu / DZ), 0, nz - 1)); wnd = 4
        y0, y1 = max(yi-wnd, 0), min(yi+wnd+1, Ny)
        x0, x1 = max(xi-wnd, 0), min(xi+wnd+1, Nx)
        sb = np.clip(volf[mi, y0:y1, x0:x1] - volf[mi, y0:y1, x0:x1].min(), 0, None)
        if sb.sum() > 0:
            yy, xx = np.mgrid[y0:y1, x0:x1]
            xr = float((sb*xx).sum()/sb.sum()); yr = float((sb*yy).sum()/sb.sum())
        else:
            xr, yr = x, y
        P.append((xr*DX, yr*DX, mu - z[0]*DZ))
        S.append({82: "Pb", 38: "Sr", 22: "Ti"}[sp])
        E.append(err)
    return np.array(P), np.array(S), np.array(E)


# ---------------------------------------------------------------- recon (model-informed)
V = load_nl70()
pos, Z = C.load_gt_model()
OFF, resid = C.hero_depth_offset(V, DZ, DX, 193, 125, pos, Z)
print(f"depth registration z_rec = z_GT + {OFF:.2f} Å (resid {resid:.2f})")
pts, species, err = model_atoms(V, pos, Z, OFF)
from collections import Counter
print(f"model-informed: {len(pts)} atoms  {dict(Counter(species.tolist()))}  "
      f"median depth 1σ {np.median(err):.2f} Å")


SPECIES_COL = {"Pb": PB_COL, "Sr": "#2ecc71", "Ti": TI_COL, "O": "#e74c3c"}
R_SPH, R_BAR, CAP = 0.40, 0.05, 0.26             # sphere radius, error-bar radius, cap half-len


def add_recon(p):
    """Each atom = a small sphere (model-informed position) + an along-beam error bar =
    ±1σ depth-position uncertainty.  Coloured by species."""
    for sp, col in SPECIES_COL.items():
        m = species == sp
        if not m.any():
            continue
        Pm, Em = pts[m], err[m]
        sph = pv.PolyData(Pm).glyph(
            geom=pv.Sphere(radius=R_SPH, theta_resolution=18, phi_resolution=18),
            scale=False, orient=False)
        p.add_mesh(sph, color=col, smooth_shading=True, specular=0.4,
                   specular_power=12, show_scalar_bar=False)
        parts = []
        for (x, y, zz), e in zip(Pm, Em):
            parts.append(pv.Line((x, y, zz - e), (x, y, zz + e)).tube(radius=R_BAR, n_sides=8))
            for zcap in (zz - e, zz + e):
                parts.append(pv.Line((x - CAP, y, zcap), (x + CAP, y, zcap)).tube(radius=R_BAR, n_sides=6))
        p.add_mesh(pv.merge(parts), color=col, smooth_shading=True, show_scalar_bar=False)


img_recon = render(add_recon)

# ---------------------------------------------------------------- ground-truth panel
gx, gy = pos[:, 0] - X0c, pos[:, 1] - Y0c
gz = pos[:, 2] + OFF - Z_SHIFT
sel = (gx > -0.3) & (gx < Lx+0.3) & (gy > -0.3) & (gy < Lx+0.3) & (gz > -0.3) & (gz < Lz+0.3)
gxyz, gZ = np.c_[gx, gy, gz][sel], Z[sel]
print(f"GT segment: {sel.sum()} atoms")


def add_gt(p):
    for zz, (col, lab, rad) in GT_STYLE.items():
        m = gZ == zz
        if not m.any():
            continue
        sph = pv.PolyData(gxyz[m]).glyph(
            geom=pv.Sphere(radius=rad, theta_resolution=16, phi_resolution=16),
            scale=False, orient=False)
        p.add_mesh(sph, color=col, smooth_shading=True, specular=0.3, show_scalar_bar=False)

img_gt = render(add_gt)

# ---------------------------------------------------------------- compose
fig = plt.figure(figsize=(13, 5.8)); fig.patch.set_facecolor("white")
for i, (im, ttl) in enumerate([
        (img_recon, f"NL70 reconstruction — model-informed atoms ({len(pts)}); "
                    r"bar = $\pm1\sigma$ depth error "
                    f"({np.median(err):.2f} Å median).  blue Pb · grey Ti  "
                    "(Ti de-skewed from its apical O; O omitted — at noise floor)"),
        (img_gt, f"ground-truth model  (scanned segment, {sel.sum()} atoms) — "
                 "discrete Pb / Sr / Ti / O")]):
    ax = fig.add_subplot(2, 1, i + 1)
    ax.imshow(im); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(ttl, fontsize=10.5)
    # '...' at each end: the columns continue toward entrance / exit
    for xf, ha in [(0.005, "left"), (0.995, "right")]:
        ax.text(xf, 0.5, r"$\cdots$", transform=ax.transAxes, fontsize=26, color="0.45",
                ha=ha, va="center")
handles = [plt.Line2D([0], [0], marker="o", ls="", mfc=c, mec="none", ms=9, label=l)
           for _, (c, l, _) in GT_STYLE.items()]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0))
fig.text(0.5, 0.055, r"beam direction / depth  $\longrightarrow$   "
         r"(middle %.0f Å of the ~%.0f Å columns;  $\cdots$ = continues to entrance / exit)"
         % (Lz, Z_FULL), ha="center", color="0.4", fontsize=9)
fig.suptitle("3-D atom model from the clean NL70 reconstruction — recon vs ground truth "
             "(mid-column, side-on, same scale)", fontsize=13, y=0.99)
fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.13, hspace=0.14)
fig.savefig(C.FIG_DIR + "/fig7_nl70_atoms_vs_gt.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("wrote fig7_nl70_atoms_vs_gt.png")


# ================================================================ FIG 7b: overlay
def add_overlay(p):
    """Recon atoms (solid Pb/Ti + error bars) with the ground-truth model ghosted over
    the top (translucent) — same box, to eyeball how well the recon lands on the model."""
    p.enable_depth_peeling(12)
    add_recon(p)                                       # opaque recon Pb/Ti + bars
    for zz, (col, lab, rad) in GT_STYLE.items():
        m = gZ == zz
        if not m.any():
            continue
        sph = pv.PolyData(gxyz[m]).glyph(
            geom=pv.Sphere(radius=rad * 1.05, theta_resolution=16, phi_resolution=16),
            scale=False, orient=False)
        p.add_mesh(sph, color=col, opacity=0.25, smooth_shading=True, show_scalar_bar=False)


img_ov = render(add_overlay, win=(2200, 620), zoom=1.7)
figb = plt.figure(figsize=(13, 5.6)); figb.patch.set_facecolor("white")
axb = figb.add_axes([0.01, 0.1, 0.98, 0.78]); axb.imshow(img_ov)
axb.set_xticks([]); axb.set_yticks([])
for s in axb.spines.values():
    s.set_visible(False)
for xf, ha in [(0.005, "left"), (0.995, "right")]:
    axb.text(xf, 0.5, r"$\cdots$", transform=axb.transAxes, fontsize=26, color="0.45",
             ha=ha, va="center")
figb.suptitle("Reconstructed atoms + depth error bars (solid) with the ground-truth model "
              "overlaid (translucent) — NL70, mid-column, side-on", fontsize=12.5, y=0.98)
handles = [plt.Line2D([0], [0], marker="o", ls="", mfc=c, mec="none", ms=9, label=l)
           for _, (c, l, _) in GT_STYLE.items()]
figb.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.02))
figb.text(0.5, 0.075, "solid = reconstruction (Pb, Ti)   ·   translucent = ground-truth model "
          "(Pb, Ti, O)", ha="center", color="0.4", fontsize=9)
figb.savefig(C.FIG_DIR + "/fig7b_nl70_overlay.png", dpi=300, bbox_inches="tight")
plt.close(figb)
print("wrote fig7b_nl70_overlay.png")


# ================================================================ FIG 7c: down the barrel
def render_topdown(add_fn, win=(1300, 1300), zoom=1.4):
    p = pv.Plotter(off_screen=True, window_size=win)
    p.set_background("white")
    add_fn(p)
    p.add_mesh(pv.Box(bounds=(0, Lx, 0, Lx, 0, Lz)), style="wireframe",
               color="#b0b0b0", line_width=1, opacity=0.35)
    p.view_xy()                                          # look straight down the beam (z)
    p.camera.zoom(zoom)
    img = p.screenshot(return_img=True); p.close()
    return _trim(img)


img_db = render_topdown(add_overlay)
figc = plt.figure(figsize=(7.6, 8.4)); figc.patch.set_facecolor("white")
axc = figc.add_axes([0.02, 0.05, 0.96, 0.85]); axc.imshow(img_db)
axc.set_xticks([]); axc.set_yticks([])
for s in axc.spines.values():
    s.set_visible(False)
figc.suptitle("Down the barrel (looking along the beam) — recon atoms (solid) vs\n"
              "ground-truth model (translucent):  do the columns line up in x–y?",
              fontsize=12.5, y=0.98)
figc.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.01))
figc.savefig(C.FIG_DIR + "/fig7c_nl70_downbarrel.png", dpi=300, bbox_inches="tight")
plt.close(figc)
print("wrote fig7c_nl70_downbarrel.png")
