#!/usr/bin/env python
"""FIG 8 — NL70 reconstructed phase vs the ground-truth projected potential.

The reconstructed phase (a phase object) should look like the specimen's projected
potential.  Top row: both summed along the beam (x-y, "down the barrel").  Bottom row:
depth cross-sections (z vs x) — the recon is depth-resolved; the GT potential is sliced
and z-registered to it — so you also see how the atomic planes line up along the beam.

  python fig8_nl70_potential.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.gridspec import GridSpec
from skimage.transform import resize
import ase.io
import abtem

import dose_fig_common as C

C.use_pub_style()
DZ, DX = 0.999, 0.04916
X0, Y0 = 30.0, 10.0
YC, YB = 200, 9                                  # cross-section y-band (rows): centre, half-width
SAMP, STH = 0.10, 1.0                            # GT potential sampling / slice thickness (Å)


def load_recon():
    V = np.angle(np.load(os.path.expanduser("~/Desktop/NL70_new_vol.npy"))).astype(np.float32)
    V -= np.median(V, (1, 2), keepdims=True)
    return V


def gt_potential():
    """3-D projected potential of the scanned model, abtem, same orientation as the sim."""
    a = ase.io.read(C.VASP)
    a.rotate(-90, "y", rotate_cell=True)
    a = abtem.orthogonalize_cell(a)
    s = max(a.cell.lengths()[:2])
    a.cell[0, 0] = s; a.cell[1, 1] = s
    a.center(axis=0); a.center(axis=1); a.center(axis=2, vacuum=2.0)
    pot = abtem.Potential(a, sampling=SAMP, projection="infinite",
                          parametrization="kirkland", slice_thickness=STH)
    arr = np.asarray(pot.build().array)          # [n_slice, gx, gy]  (abtem x,y order)
    samp = s / arr.shape[1]
    sth = float(a.cell[2, 2]) / arr.shape[0]
    return arr, samp, sth


# ---------------------------------------------------------------- data
V = load_recon()
nL, Ny, Nx = V.shape
recon_xy = V.sum(0)                                           # [Y, X]
recon_cs = V[:, YC-YB:YC+YB, :].mean(1)                       # [z, X]

pos, Z = C.load_gt_model()
OFF, _ = C.hero_depth_offset(V, DZ, DX, 193, 125, pos, Z)
arr, samp, sth = gt_potential()
proj = arr.sum(0)                                            # [gx, gy]

# in-plane ROI (X in [30,50], Y in [10,30]) -> transpose to [Y, X], resize to recon grid
xi = slice(int(round(X0/samp)), int(round((X0+Nx*DX)/samp)))
yi = slice(int(round(Y0/samp)), int(round((Y0+Ny*DX)/samp)))
gt_xy = resize(proj[xi, yi].T, (Ny, Nx), order=1, preserve_range=True)

# GT depth cross-section: same y-band + X-ROI, sum over y -> [z_cell, X]; z -> recon frame
yb = slice(int(round((Y0+(YC-YB)*DX)/samp)), int(round((Y0+(YC+YB)*DX)/samp)))
gt_cs_raw = arr[:, xi, yb].sum(2)                            # [n_slice, X]
gt_cs = resize(gt_cs_raw, (nL, Nx), order=1, preserve_range=True)  # onto recon x-grid
zc_axis = (np.arange(arr.shape[0]) + 0.5) * sth + OFF        # GT slice z -> recon z (Å)


def scaled(im, lo=1, hi=99.5):
    return dict(vmin=np.percentile(im, lo), vmax=np.percentile(im, hi))


# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(11.8, 9.2))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1.25, 1.0], hspace=0.22, wspace=0.24,
              left=0.07, right=0.9, top=0.9, bottom=0.07)
ext_xy = [X0, X0 + Nx*DX, Y0 + Ny*DX, Y0]                     # X, Y (Å)
ext_cs = [X0, X0 + Nx*DX, (nL-0.5)*DZ, 0.5*DZ]               # X (Å), z entrance->exit

# --- top: x-y (summed along beam)
a00 = fig.add_subplot(gs[0, 0])
im0 = a00.imshow(recon_xy, cmap="inferno", extent=ext_xy, **scaled(recon_xy, lo=40))
a00.set_title("reconstruction — Σ phase along beam", fontsize=11)
a00.set_ylabel("Y (Å)")
fig.colorbar(im0, ax=a00, fraction=0.046, pad=0.03, label="Σ phase (rad)")

a01 = fig.add_subplot(gs[0, 1])
im1 = a01.imshow(gt_xy, cmap="inferno", extent=ext_xy, **scaled(gt_xy))
a01.set_title("ground truth — projected potential", fontsize=11)
fig.colorbar(im1, ax=a01, fraction=0.046, pad=0.03, label="proj. potential")
for a in (a00, a01):
    a.set_xlabel("X (Å)")

# --- bottom: z vs x cross-sections (same y-band)
a10 = fig.add_subplot(gs[1, 0])
im2 = a10.imshow(recon_cs, cmap="magma", extent=ext_cs, aspect="auto",
                 **scaled(recon_cs, lo=35, hi=99.7))
a10.set_title(f"reconstruction — depth cross-section (Y ≈ {Y0+YC*DX:.0f} Å)", fontsize=11)
a10.set_ylabel("depth  z (Å)  [entrance $\\rightarrow$ exit]"); a10.set_xlabel("X (Å)")
fig.colorbar(im2, ax=a10, fraction=0.046, pad=0.03, label="phase (rad)")

a11 = fig.add_subplot(gs[1, 1])
im3 = a11.imshow(gt_cs, cmap="magma", aspect="auto", **scaled(gt_cs),
                 extent=[X0, X0 + Nx*DX, zc_axis[-1], zc_axis[0]])
a11.set_ylim((nL-0.5)*DZ, 0.5*DZ)                            # match recon z-range
a11.set_title("ground truth — potential cross-section (z-registered)", fontsize=11)
a11.set_xlabel("X (Å)")
fig.colorbar(im3, ax=a11, fraction=0.046, pad=0.03, label="potential")

fig.suptitle("NL70 reconstructed phase vs ground-truth projected potential — "
             "x-y (top) and depth cross-section (bottom)", fontsize=13, y=0.965)
C.savefig(fig, "fig8_nl70_potential", pdf=False)
print(f"wrote fig8  |  xy NCC(recon,GT) = "
      f"{np.corrcoef(recon_xy.ravel(), gt_xy.ravel())[0,1]:.3f}  OFF={OFF:.2f}")
