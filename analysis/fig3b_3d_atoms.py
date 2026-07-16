#!/usr/bin/env python
"""FIG 3b — reconstructed 3-D atoms (side-on) vs the ground-truth model.

Viewed from the SIDE (beam/depth runs left-right, not down the barrel) so individual
atoms are visible: the reconstructed atoms sit in rows — each row is an atomic column
seen edge-on, and the bright/dark alternation is the heavy Pb vs lighter Ti columns.

Each reconstructed atom is drawn as a small ellipsoid, mildly elongated along the beam
(schematic — the reconstruction blurs more in depth than in-plane).  The ground-truth
structure (ASE, from the .vasp) for the same scanned segment is shown alongside as a
discrete Pb/Sr/Ti/O ball model, at the same viewpoint and scale.

  fig3b_atoms_vs_gt      : recon | ground truth (side-on).
  fig3b_atoms_dose_grid  : recon atoms at each dose.

  python fig3b_3d_atoms.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, maximum_filter
import pyvista as pv

import dose_fig_common as C

pv.OFF_SCREEN = True
C.use_pub_style()

CROP = dict(r=(95, 315), z=(8, 97))            # wide interior region
SIGMA = (1.0, 1.6, 1.6)
FOOT = (5, 25, 25)
PCTL = 98.5
RZ_HALF, RXY_HALF = 1.8, 0.42                  # ellipsoid half-axes (Å): ~3.6 x 0.84
CAM_OFF, CAM_UP = (0.1, 1.0, 0.18), (1, 0, 0)  # side-on, beam horizontal, slight tilt
GT_STYLE = {82: ("#6a5acd", "Pb", 0.80), 38: ("#2ecc71", "Sr", 0.78),
            22: ("#9aa0a6", "Ti", 0.52), 8:  ("#e74c3c", "O", 0.34)}

dx = C.DX_FALLBACK
_, dz, _ = C.load_dose("1e10")
Lx = (CROP["r"][1] - CROP["r"][0]) * dx
Lz = (CROP["z"][1] - CROP["z"][0]) * dz
X0c = C.X0 + CROP["r"][0] * dx
Y0c = C.Y0 + CROP["r"][0] * dx
Z_SHIFT = (CROP["z"][0] + 0.5) * dz
CENTER = (Lx / 2, Lx / 2, Lz / 2)
R = max(Lx, Lz)


def detect_atoms(V, thr_abs=None):
    r, z = CROP["r"], CROP["z"]
    vol = gaussian_filter(V[z[0]:z[1], r[0]:r[1], r[0]:r[1]].astype(np.float32), SIGMA)
    vol -= np.median(vol, (1, 2), keepdims=True)
    thr = thr_abs if thr_abs is not None else np.percentile(vol, PCTL)
    pk = (vol == maximum_filter(vol, size=FOOT)) & (vol > thr)
    zi, yi, xi = np.where(pk)
    return np.c_[xi * dx, yi * dx, zi * dz], vol[pk], thr


def ellipsoid_mesh(pts, inten):
    blocks = []
    for p, w in zip(pts, inten):
        e = pv.Sphere(radius=1.0, theta_resolution=12, phi_resolution=12)
        e.points *= np.array([RXY_HALF, RXY_HALF, RZ_HALF])
        e.points += p
        e["inten"] = np.full(e.n_points, w)
        blocks.append(e)
    return pv.merge(blocks)


def _trim(img, pad=8):
    mask = (img < 248).any(2)
    ys, xs = np.where(mask)
    if not len(ys):
        return img
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad + 1, img.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad + 1, img.shape[1])
    return img[y0:y1, x0:x1]


def render(add_fn, win=(1900, 470), zoom=1.6):
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


# ---------------------------------------------------------------- recon (1e10)
ref, dz, _ = C.load_dose("1e10")
pts, inten, thr10 = detect_atoms(ref)
ptp = inten.max() - inten.min()
CLIM = (inten.min() - 0.8 * ptp, inten.max())    # brighten: no near-black atoms
print(f"1e10: {len(pts)} atoms  box {Lx:.0f}x{Lx:.0f}x{Lz:.0f} Å")


def add_recon(p, pts=pts, inten=inten):
    p.add_mesh(ellipsoid_mesh(pts, inten), scalars="inten", cmap="copper", clim=CLIM,
               smooth_shading=True, specular=0.4, specular_power=12, show_scalar_bar=False)

img_recon = render(add_recon)

# ---------------------------------------------------------------- ground truth
pos, Z = C.load_gt_model()
OFF, _ = C.hero_depth_offset(ref, dz, dx, 193, 209, pos, Z)
gx, gy = pos[:, 0] - X0c, pos[:, 1] - Y0c
gz = pos[:, 2] + OFF - Z_SHIFT
sel = (gx > -0.3) & (gx < Lx + 0.3) & (gy > -0.3) & (gy < Lx + 0.3) & (gz > -0.3) & (gz < Lz + 0.3)
gxyz, gZ = np.c_[gx, gy, gz][sel], Z[sel]
print(f"GT segment: {sel.sum()} atoms")


def add_gt(p, gxyz=gxyz, gZ=gZ):
    for zz, (col, lab, rad) in GT_STYLE.items():
        m = gZ == zz
        if not m.any():
            continue
        sph = pv.PolyData(gxyz[m]).glyph(
            geom=pv.Sphere(radius=rad, theta_resolution=16, phi_resolution=16),
            scale=False, orient=False)
        p.add_mesh(sph, color=col, smooth_shading=True, specular=0.3, show_scalar_bar=False)

img_gt = render(add_gt)

# ---------------------------------------------------------------- compose (stacked)
fig = plt.figure(figsize=(13.5, 6.6)); fig.patch.set_facecolor("white")
for i, (im, ttl) in enumerate([
        (img_recon, f"reconstruction  (10$^{{10}}$ e/Å², {len(pts)} atoms) — "
                    "rows are atomic columns seen edge-on; bright = Pb, dark = Ti"),
        (img_gt, f"ground-truth model  (scanned segment, {sel.sum()} atoms) — "
                 "discrete Pb / Sr / Ti / O")]):
    ax = fig.add_subplot(2, 1, i + 1)
    ax.imshow(im); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(ttl, fontsize=10.5)
handles = [plt.Line2D([0], [0], marker="o", ls="", mfc=c, mec="none", ms=9, label=l)
           for _, (c, l, _) in GT_STYLE.items()]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0))
fig.text(0.5, 0.055, r"beam direction / depth  (entrance $\rightarrow$ exit,  %.0f Å)  $\longrightarrow$"
         % Lz, ha="center", color="0.4", fontsize=9)
fig.suptitle("What ptychography does to the atoms — reconstructed vs ground truth "
             "(side-on, same scale)", fontsize=13, y=0.99)
fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.12, hspace=0.28)
fig.savefig(C.FIG_DIR + "/fig3b_atoms_vs_gt.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("wrote fig3b_atoms_vs_gt.png")

# ---------------------------------------------------------------- dose grid
counts = {"1e10": len(pts)}
imgs = {"1e10": img_recon}
for d in C.DOSES[1:]:
    V, _, _ = C.load_dose(d)
    P, IN, _ = detect_atoms(V, thr_abs=thr10)
    counts[d] = len(P)
    imgs[d] = render(lambda p, P=P, IN=IN: p.add_mesh(
        ellipsoid_mesh(P, IN), scalars="inten", cmap="copper", clim=CLIM,
        smooth_shading=True, specular=0.4, show_scalar_bar=False), win=(1600, 400))
    del V

fig, axes = plt.subplots(4, 1, figsize=(12, 8.6)); fig.patch.set_facecolor("white")
for axx, d in zip(axes, C.DOSES):
    axx.imshow(imgs[d]); axx.set_xticks([]); axx.set_yticks([])
    for s in axx.spines.values():
        s.set_visible(False)
    axx.set_ylabel(f"{C.DOSE_TeX[d]} e/Å²\n{counts[d]} atoms", fontsize=11, rotation=0,
                   ha="right", va="center", labelpad=28)
fig.suptitle("Reconstructed atoms vs dose (side-on, same detection) — ordered rows give "
             "way to scattered noise", fontsize=13, y=0.965)
fig.subplots_adjust(left=0.1, right=0.995, top=0.93, bottom=0.01, hspace=0.05)
fig.savefig(C.FIG_DIR + "/fig3b_atoms_dose_grid.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("wrote fig3b_atoms_dose_grid.png  counts:", counts)
