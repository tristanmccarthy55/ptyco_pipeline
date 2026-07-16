#!/usr/bin/env python
"""FIG 3 — 3-D reconstructed lattice (low phase made transparent).

The reconstructed phase volume is rendered with an opacity transfer function that clips
the low-phase background to fully transparent, so only the atomic columns/planes remain:
a 3-D "forest" of columns with the resolved atomic planes visible as bright beads down
each string.  Deliverables:

  fig3_3d_hero        : the cinematic single-volume render at 1e10 (framed, colorbar).
  fig3_3d_dose_grid   : the SAME render across the four doses — the 3-D lattice
                        dissolving into noise as dose drops.
  fig3_3d_lattice.html: an interactive, rotatable isosurface (open in any browser).

pyvista renders off-screen (VTK); the renders are composited into clean matplotlib
frames.  The interactive scene is a self-contained plotly HTML (no kaleido needed).

  python fig3_3d_lattice.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.ndimage import gaussian_filter
import pyvista as pv

import dose_fig_common as C

pv.OFF_SCREEN = True
C.use_pub_style()

CMAP = "magma"
CROP = dict(r=(35, 375), c=(35, 375), z=(7, 99))          # interior; trim ent/exit
SIGMA = (0.8, 1.3, 1.3)                                    # denoise (z, y, x)
OPAC = [0.0, 0.03, 0.12, 0.32, 0.62, 1.0]                  # low phase -> transparent
VIEW = dict(vec=(1.0, 0.35, 0.5), az=45, el=20)
WIN = (1500, 1250)


def build_grid(V, dz, dx):
    r, c, z = CROP["r"], CROP["c"], CROP["z"]
    vol = V[z[0]:z[1], r[0]:r[1], c[0]:c[1]].astype(np.float32)
    vol = gaussian_filter(vol, SIGMA)
    vol -= np.median(vol, (1, 2), keepdims=True)
    Nz, Ny, Nx = vol.shape
    grid = pv.ImageData(dimensions=(Nx, Ny, Nz), spacing=(dx, dx, dz))
    grid.point_data["phase"] = vol.ravel(order="C")
    return grid, vol


def render_volume(grid, clim, win=WIN, zoom=1.15):
    p = pv.Plotter(off_screen=True, window_size=win)
    p.set_background("black")
    p.add_volume(grid, scalars="phase", cmap=CMAP, clim=clim, opacity=OPAC,
                 shade=True, ambient=0.35, diffuse=0.9, specular=0.2,
                 show_scalar_bar=False)
    p.view_vector(VIEW["vec"])
    p.camera.azimuth = VIEW["az"]; p.camera.elevation = VIEW["el"]
    p.reset_camera(); p.camera.zoom(zoom)
    img = p.screenshot(return_img=True); p.close()
    return img


# ----------------------------------------------------------------- reference & clim
ref, dz, dx = C.load_dose("1e10")
gref, vref = build_grid(ref, dz, dx)
CLIM = [float(np.percentile(vref, 55)), float(np.percentile(vref, 99.6))]
depth_A = (CROP["z"][1] - CROP["z"][0]) * dz
inplane_A = (CROP["r"][1] - CROP["r"][0]) * dx
print(f"crop {vref.shape}  clim {np.round(CLIM,4)}  box {inplane_A:.0f}x{inplane_A:.0f}x{depth_A:.0f} Å")


def frame(ax, img, title=None):
    ax.imshow(img); ax.set_facecolor("black")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if title:
        ax.set_title(title, color="black", fontsize=12)


# ================================================================ HERO (1e10)
img = render_volume(gref, CLIM, zoom=1.3)
fig = plt.figure(figsize=(8.6, 8.0)); fig.patch.set_facecolor("white")
ax = fig.add_axes([0.02, 0.02, 0.82, 0.9]); frame(ax, img)
ax.annotate("", xy=(0.045, 0.06), xytext=(0.045, 0.92), xycoords="axes fraction",
            arrowprops=dict(arrowstyle="-|>", color="0.35", lw=1.6))
ax.text(0.015, 0.5, r"depth  (entrance $\rightarrow$ exit,  %.0f Å)" % depth_A, rotation=90,
        transform=ax.transAxes, va="center", ha="center", color="0.3", fontsize=9.5)
cax = fig.add_axes([0.87, 0.2, 0.025, 0.55])
sm = ScalarMappable(norm=Normalize(*CLIM), cmap=CMAP); sm.set_array([])
cb = fig.colorbar(sm, cax=cax); cb.set_label("reconstructed phase (rad)")
fig.suptitle("3-D reconstructed lattice at 10$^{10}$ e/Å²\n"
             "low phase made transparent — beads down each column are resolved atomic planes",
             fontsize=12.5, y=0.98)
C.savefig(fig, "fig3_3d_hero", pdf=False)
plt.close(fig)

# ================================================================ DOSE GRID
imgs = {"1e10": img}
for d in C.DOSES[1:]:
    V, dzi, dxi = C.load_dose(d)
    g, _ = build_grid(V, dzi, dxi)
    imgs[d] = render_volume(g, CLIM, win=(1100, 950))
    del V

fig, axes = plt.subplots(1, 4, figsize=(17, 4.7)); fig.patch.set_facecolor("white")
for ax, d in zip(axes, C.DOSES):
    frame(ax, imgs[d], title=f"{C.DOSE_TeX[d]} e/Å²")
fig.subplots_adjust(left=0.005, right=0.9, top=0.84, bottom=0.02, wspace=0.02)
cax = fig.add_axes([0.915, 0.12, 0.012, 0.66])
sm = ScalarMappable(norm=Normalize(*CLIM), cmap=CMAP); sm.set_array([])
cb = fig.colorbar(sm, cax=cax); cb.set_label("phase (rad)")
fig.suptitle("The 3-D atomic lattice vs electron dose — same render, low phase transparent "
             "(depth runs into the page)", fontsize=13, y=0.955)
C.savefig(fig, "fig3_3d_dose_grid", pdf=False)
plt.close(fig)

# ================================================================ INTERACTIVE HTML
import plotly.graph_objects as go

rc = (120, 285)                                            # smaller region, ~8 columns
sub = gaussian_filter(ref[CROP["z"][0]:CROP["z"][1], rc[0]:rc[1], rc[0]:rc[1]], (0.8, 1.2, 1.2))
sub = sub - np.median(sub, (1, 2), keepdims=True)
st = 2                                                     # in-plane stride for a light HTML
sub = sub[:, ::st, ::st]
Nz, Ny, Nx = sub.shape
xg = np.arange(Nx) * dx * st
yg = np.arange(Ny) * dx * st
zg = np.arange(Nz) * dz
Xg, Yg, Zg = np.meshgrid(xg, yg, zg, indexing="ij")        # note: value below matches
val = np.transpose(sub, (2, 1, 0))                         # -> [x,y,z] to match meshgrid
iso = float(np.percentile(sub, 98.5))
figp = go.Figure(go.Isosurface(
    x=Xg.ravel(), y=Yg.ravel(), z=Zg.ravel(), value=val.ravel(),
    isomin=iso, isomax=float(sub.max()), surface_count=3, colorscale="Magma",
    opacity=0.55, caps=dict(x_show=False, y_show=False, z_show=False),
    showscale=True, colorbar=dict(title="phase (rad)")))
figp.update_layout(
    title="Reconstructed 3-D atomic lattice (10¹⁰ e/Å²) — drag to rotate",
    scene=dict(xaxis_title="x (Å)", yaxis_title="y (Å)", zaxis_title="depth z (Å)",
               aspectmode="data", zaxis=dict(autorange="reversed")),
    width=1000, height=850, template="plotly_white")
import os
html = os.path.join(C.FIG_DIR, "fig3_3d_lattice.html")
figp.write_html(html, include_plotlyjs=True, full_html=True)
print("wrote", html)
print("fig3 done.")
