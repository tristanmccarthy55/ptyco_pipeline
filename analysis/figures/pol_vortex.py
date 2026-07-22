#!/usr/bin/env python
"""Ground-truth polarisation vortex map for the PTO/STO labyrinth model.

Polarisation proxy = B-site off-centering: for each Ti, delta = r_Ti - centroid(6 nearest O)
(the Ti displacement from its oxygen-octahedron charge centre).  Computed in the RAW POSCAR
frame with full PBC, then mapped into the PREPARED sim frame (rotate -90 y, orthogonalize,
pad-to-square, centre) so it lines up with the 4D-STEM scan box.  Beam is along prepared Z, so
the prepared X-Y plane is the (stacking x in-plane) plane where the vortices curl; average each
in-plane column along the beam -> one arrow per column.

Scan area (from sim/simulate_4dstem.py): centre (40,20), 20 A window -> X[30,50] Y[10,30].

  python pol_vortex.py
"""
import os
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ase.io
import abtem

VASP = "/Users/u2109287/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PtycoShelves/PTO6_STO6_18_18_labyrinthPoscar.vasp"
SCAN = dict(cx=40.0, cy=20.0, w=20.0)


# ---------------------------------------------------------------- polarisation (raw frame)
a = ase.io.read(VASP)
Z = a.get_atomic_numbers()
P = a.get_positions()
cell = a.cell.array
ti = P[Z == 22]
ox = P[Z == 8]

# tile O over the 27 periodic images so nearest-neighbour search respects PBC
shifts = np.array([i * cell[0] + j * cell[1] + k * cell[2]
                   for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
ox_t = (ox[None] + shifts[:, None]).reshape(-1, 3)
d6, i6 = cKDTree(ox_t).query(ti, k=6)
o_cen = ox_t[i6].mean(1)                       # oxygen octahedron centre per Ti
delta = ti - o_cen                             # Ti off-centering  ~ local polarisation
mag = np.linalg.norm(delta, axis=1)
print(f"Ti-O6: nn dist {d6.mean():.2f} A | |delta| mean {mag.mean():.3f} "
      f"med {np.median(mag):.3f} max {mag.max():.3f} A")


# ---------------------------------------------------------------- map raw -> prepared frame
# The sim prep is rotate(-90,y) -> orthogonalize -> pad-square -> centre.  orthogonalize only
# redefines the (near-orthogonal) cell, it does NOT move cartesian atoms, so prepared = R @ raw
# + a pure centring translation.  R for rotate(-90,y): (x,y,z) -> (-z, y, x).
b = ase.io.read(VASP)
b.rotate(-90, "y", rotate_cell=True)
b = abtem.orthogonalize_cell(b)
side = float(max(b.cell.lengths()[:2]))
b.cell[0, 0] = side; b.cell[1, 1] = side
b.center(axis=0); b.center(axis=1); b.center(axis=2, vacuum=2.0)
Pp = b.get_positions()                           # prepared cloud (reordered; positions valid)

# the rigid raw->prepared map is a signed axis permutation M (rotation + orthogonalize's axis
# relabel/reflection) plus a translation; find it by matching to the real prepared Ti cloud so
# the vortex HANDEDNESS matches what the sim/recon sees.
from itertools import permutations, product
tip_b = Pp[b.get_atomic_numbers() == 22]
tree_b = cKDTree(tip_b)
best = None
for perm in permutations(range(3)):
    for sg in product((1, -1), repeat=3):
        M = np.zeros((3, 3))
        for i, c in enumerate(perm):
            M[i, c] = sg[i]
        q = (M @ ti.T).T
        tt = tip_b.min(0) - q.min(0)
        med = np.median(tree_b.query(q + tt)[0])
        if best is None or med < best[0]:
            best = (med, M, tt)
_, M, _ = best
q = (M @ ti.T).T
idx = tree_b.query(q + (tip_b.min(0) - q.min(0)))[1]
t = np.median(tip_b[idx] - q, axis=0)            # refined centring translation
ti_p = q + t                                     # prepared Ti positions
dl_p = (M @ delta.T).T                             # prepared polarisation vectors (rigid part)
nn = tree_b.query(ti_p, k=1)[0]
print(f"prepared map: M rows {M.tolist()} | median atom match {np.median(nn):.3f} A | side {side:.2f}")
Xp, Yp, Zp = ti_p[:, 0], ti_p[:, 1], ti_p[:, 2]
dX, dY = dl_p[:, 0], dl_p[:, 1]

# average each in-plane column along the beam (prepared Z)
key = np.round(np.c_[Xp, Yp], 1)
uniq, inv = np.unique(key, axis=0, return_inverse=True)
cnt = np.bincount(inv)
gX = np.bincount(inv, Xp) / cnt; gY = np.bincount(inv, Yp) / cnt
gdX = np.bincount(inv, dX) / cnt; gdY = np.bincount(inv, dY) / cnt
gmag = np.hypot(gdX, gdY)
print(f"{len(uniq)} in-plane columns | in-plane |P| mean {gmag.mean():.3f} A "
      f"| box 0..{side:.0f} A | scan X[{SCAN['cx']-SCAN['w']/2:.0f},{SCAN['cx']+SCAN['w']/2:.0f}]")


# ---------------------------------------------------------------- publication figure
plt.rcParams.update({"figure.facecolor": "white", "font.family": "sans-serif",
                     "font.sans-serif": ["DejaVu Sans"], "font.size": 15})
ang = np.arctan2(gdY, gdX)                        # -pi..pi (arrow direction)
fig, ax = plt.subplots(figsize=(9.6, 9.6))
fig.subplots_adjust(left=0.09, right=0.99, top=0.88, bottom=0.08)
q = ax.quiver(gX, gY, gdX, gdY, ang, cmap="hsv", clim=(-np.pi, np.pi),
              angles="xy", scale_units="xy", scale=0.06, width=0.0042,
              headwidth=4, headlength=4.5, pivot="mid")
# box + scan area
ax.add_patch(plt.Rectangle((0, 0), side, side, fill=False, ec="0.35", lw=1.4))
x0, y0 = SCAN["cx"] - SCAN["w"] / 2, SCAN["cy"] - SCAN["w"] / 2
ax.add_patch(plt.Rectangle((x0, y0), SCAN["w"], SCAN["w"], fill=False, ec="k", lw=2.4))
ax.text(SCAN["cx"], y0 + SCAN["w"] + 0.8, "4D-STEM scan area", ha="center", va="bottom",
        fontsize=14, weight="medium",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.6"))
ax.set_aspect("equal"); ax.set_xlim(-2, side + 2); ax.set_ylim(-2, side + 2)
ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)")
fig.suptitle("Ground-truth polarisation — PbTiO$_3$/SrTiO$_3$ vortices", fontsize=17,
             y=0.965, weight="medium")
ax.set_title("arrow = Ti offset from its O$_6$ charge centre  (projected along the beam)",
             fontsize=13, color="0.35")

# angle colour-wheel inset (bottom-right, over vacuum)
wa = fig.add_axes([0.775, 0.10, 0.15, 0.15], projection="polar")
th = np.linspace(-np.pi, np.pi, 360); rr = np.linspace(0.4, 1, 2)
TH, _ = np.meshgrid(th, rr)
wa.pcolormesh(np.tile(th, (2, 1)), np.tile(rr[:, None], (1, 360)), np.tile(th, (2, 1)),
              cmap="hsv", shading="auto", clim=(-np.pi, np.pi))
wa.set_yticks([]); wa.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2])
wa.set_xticklabels(["→", "↑", "←", "↓"], fontsize=11)
wa.set_title("P direction", fontsize=10, color="0.35", pad=1)

fig.savefig(os.path.expanduser("~/Desktop/dose_series/slides/slide_pol_vortex.png"), dpi=200,
            bbox_inches="tight")
print("wrote slide_pol_vortex.png")
