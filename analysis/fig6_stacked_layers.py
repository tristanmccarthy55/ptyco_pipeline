#!/usr/bin/env python
"""FIG 6 — Aesthetic: spaced depth-slices as floating transparent planes (1e10).

A few well-separated in-plane slices are drawn as horizontal planes stacked in depth,
with low phase made transparent so only the atoms remain and you can see THROUGH each
plane to the atoms at other depths.  Because the columns are resolved in depth, different
planes light up at different in-plane positions — the picture makes the depth resolution
tangible.  Purely for looks / talks.

  python fig6_stacked_layers.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.ndimage import gaussian_filter, zoom

import dose_fig_common as C

LAYERS = list(range(47, 54))                  # 7 neighbouring slices, centred on layer 50
CROP = (70, 335)                              # interior in-plane window
DS = 3                                        # downsample for a light mesh
CMAP = cm.magma


def prep(V, layer):
    sl = V[layer, CROP[0]:CROP[1], CROP[0]:CROP[1]].astype(np.float32)
    sl = gaussian_filter(sl, 0.8)
    return sl[::DS, ::DS]


ref, dz, dx = C.load_dose("1e10")
slices = [prep(ref, l) for l in LAYERS]
allv = np.concatenate([s.ravel() for s in slices])
lo, hi = np.percentile(allv, 40), np.percentile(allv, 99.6)

n, m = slices[0].shape
xs = np.arange(m) * dx * DS
ys = np.arange(n) * dx * DS
X, Y = np.meshgrid(xs, ys)

fig = plt.figure(figsize=(9.5, 11.5), facecolor="black")
ax = fig.add_subplot(111, projection="3d")
ax.set_facecolor("black")

GAP = 7.0                                      # visual spacing between displayed planes
for i, (l, sl) in enumerate(zip(LAYERS, slices)):
    norm = np.clip((sl - lo) / (hi - lo), 0, 1)
    rgba = CMAP(norm)
    rgba[..., 3] = 1.0                                    # opaque (stacked solid cards)
    Zc = np.full_like(X, i * GAP)
    ax.plot_surface(X, Y, Zc, facecolors=rgba, rcount=n, ccount=m,
                    shade=False, linewidth=0, antialiased=False)
    # faint sheet frame so each depth reads as a stacked layer
    x0, x1, y0, y1 = xs[0], xs[-1], ys[0], ys[-1]
    ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], [i * GAP] * 5,
            color="0.55", lw=0.8, alpha=0.6)
    ax.text(xs[-1] + 1.2, ys[-1] * 0.5, i * GAP, f"z ≈ {(l+0.5)*dz:.1f} Å",
            color="white", fontsize=10, ha="left", va="center")

ax.set_box_aspect((1, 1, 1.7))
ax.view_init(elev=22, azim=-60)
ax.set_axis_off()
ax.set_zlim(-2, (len(LAYERS) - 1) * GAP + 2)
fig.text(0.5, 0.95, "Seven neighbouring depth-slices centred on z ≈ 34 Å  (10$^{10}$ e/Å²)",
         color="white", ha="center", fontsize=15)
fig.text(0.5, 0.915, "adjacent reconstruction layers, Δz = 0.67 Å (shown exploded) — "
         "how little the atoms change layer to layer",
         color="0.7", ha="center", fontsize=10.5)
fig.savefig(C.FIG_DIR + "/fig6_stacked_layers.png", dpi=300,
            facecolor="black", bbox_inches="tight")
print("wrote fig6_stacked_layers.png  layers", LAYERS, " z(Å)",
      [round((l+0.5)*dz, 1) for l in LAYERS])
