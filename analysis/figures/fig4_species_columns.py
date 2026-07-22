#!/usr/bin/env python
"""@file fig4_species_columns.py
@brief FIG 4 — Depth cross-sections down Pb, Ti+O and O columns, for every dose.

Three column TYPES (rows) x four doses (columns).  Each panel is a depth cross-section
(z vs in-plane) down a representative column of that species, with the ground-truth atoms
of that column overlaid.  Shows which sublattices are depth-resolved and how each degrades
with dose: the heavy Pb (A-site) planes are strongly resolved and survive to ~1e6; the
lighter Ti+O (B-site) planes are weaker; the pure-oxygen columns are barely modulated even
at the highest dose.  Per-row display scaling (O is ~6x fainter than Pb).

  python fig4_species_columns.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.gridspec import GridSpec
from scipy.ndimage import gaussian_filter

import dose_fig_common as C


def cross(V, yc, xc):
    """Thin depth cross-section strip down a column, mildly smoothed for display
    (sigma < the 5.9-layer plane spacing, so resolved planes are preserved)."""
    cs = V[:, yc-1:yc+2, xc-W:xc+W].mean(1)
    return gaussian_filter(cs, (0.6, 0.9))

C.use_pub_style()
W, CMAP = 22, "magma"

# representative central interior column of each type: label -> (yc, xc)
COLUMNS = [("Pb  (A-site)",       (193, 209)),
           ("Ti + O  (B-site)",   (230, 244)),
           ("O  (oxygen)",        (233, 201))]

ref, dz, dx = C.load_dose("1e10")
pos, Z = C.load_gt_model()
OFF, resid = C.hero_depth_offset(ref, dz, dx, 193, 209, pos, Z)      # global depth reg
nL = ref.shape[0]
zrec = (np.arange(nL) + 0.5) * dz
ext = [-W * dx, W * dx, zrec[-1], zrec[0]]
stroke = [pe.withStroke(linewidth=1.1, foreground="black", alpha=0.5)]
print(f"depth registration z_rec = z_GT + {OFF:.2f} Å (resid {resid:.2f})")

# per-row (per-species) display scale from the clean 1e10 cross-section
row_scale = {}
for _, (yc, xc) in COLUMNS:
    cs = cross(ref, yc, xc)
    shown = cs[zrec <= C.ZMAX_SHOW]
    row_scale[(yc, xc)] = (np.percentile(shown, 4), np.percentile(shown, 99.4))

fig = plt.figure(figsize=(12.4, 10.6))
gs = GridSpec(len(COLUMNS), len(C.DOSES), figure=fig, left=0.11, right=0.985,
              top=0.9, bottom=0.07, wspace=0.1, hspace=0.16)

# preload dose volumes once
vols = {d: C.load_dose(d) for d in C.DOSES}

for r, (label, (yc, xc)) in enumerate(COLUMNS):
    vmin, vmax = row_scale[(yc, xc)]
    Xc, Yc, Patoms, Zatoms = C.column_atoms(pos, Z, yc, xc, dx)
    for c, d in enumerate(C.DOSES):
        V, dzi, dxi = vols[d]
        cs = cross(V, yc, xc)
        ax = fig.add_subplot(gs[r, c])
        ax.imshow(cs, extent=ext, aspect="auto", cmap=CMAP, vmin=vmin, vmax=vmax)
        for zz, st in C.SPECIES.items():
            m = Zatoms == zz
            if not m.any():
                continue
            xin = Patoms[m, 0] - Xc
            zin = Patoms[m, 2] + OFF
            keep = zin <= C.ZMAX_SHOW
            if not keep.any():
                continue
            line = st["m"] in ("x", "+")
            sc = ax.scatter(xin[keep], zin[keep], s=17, marker=st["m"], linewidths=1.0,
                            alpha=0.6, c=st["c"] if line else None,
                            facecolors=(None if line else
                                        ("none" if st["m"] == "o" else st["c"])),
                            edgecolors=(None if line else st["c"]),
                            label=st["label"] if (r == 0 or c == 0) else None)
            sc.set_path_effects(stroke)
        ax.set_ylim(C.ZMAX_SHOW, zrec[0]); ax.set_xlim(-W * dx, W * dx)
        if r == 0:
            ax.set_title(f"{C.DOSE_TeX[d]} e/Å²", fontsize=12, pad=6)
        if r == len(COLUMNS) - 1:
            ax.set_xlabel("in-plane (Å)")
        else:
            ax.set_xticklabels([])
        if c == 0:
            ax.set_ylabel(f"{label}\n\ndepth z (Å)", fontsize=10.5)
        else:
            ax.set_yticklabels([])
        # per-panel species legend on the first column
        if c == 0:
            ax.legend(loc="lower left", framealpha=0.85, facecolor="white",
                      edgecolor="none", handletextpad=0.2, labelspacing=0.2,
                      fontsize=8, markerscale=2.2)

fig.suptitle("Depth cross-sections by column type and dose — which sublattices resolve in depth?\n"
             r"markers = ground-truth atoms (Pb $\circ$, Ti $\times$, O $\diamond$);  "
             r"z: entrance $\rightarrow$ exit", fontsize=12.5, y=0.975)
C.savefig(fig, "fig4_species_columns")
print("wrote fig4  columns:", [(l, c) for l, c in COLUMNS])
