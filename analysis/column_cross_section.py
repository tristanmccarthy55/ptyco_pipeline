#!/usr/bin/env python
"""Depth cross-sections down a Pb and a Ti column (NL70, the clean depth run).

For a chosen column at (yc,xc), cut a thin in-plane strip and stack it over depth
to get a [depth z  x  in-plane] image: the column is the central vertical streak,
and resolved atomic planes show up as periodic bright blobs DOWN its length.
Also a 1-D phase-vs-depth profile with tick marks at the measured 3.9 A plane period.

Usage:  python column_cross_section.py
Output: ~/Desktop/column_cross_section.png
"""
import os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

OUT = os.path.expanduser("~/Desktop")
RUN = "NL70"; DZ = 0.999; DX = 0.0492          # A per layer / per pixel
PERIOD = 3.9                                     # measured plane spacing (A)
PB = (193, 125); TI = (153, 158)                 # (row,col) in object_roi
W  = 22                                          # half in-plane strip (px) ~1.1 A

vol = np.load(os.path.join(OUT, f"{RUN}_new_vol.npy"))
V = np.angle(vol).astype(float); V -= np.median(V, (1, 2), keepdims=True)
nL = V.shape[0]; z = (np.arange(nL) + 0.5) * DZ

def cross(yc, xc):
    strip = V[:, yc-1:yc+2, xc-W:xc+W].mean(1)   # [nL, 2W], avg +-1 row
    prof  = V[:, yc-1:yc+2, xc-1:xc+2].mean((1, 2))  # column-centre phase vs z
    return strip, prof

cs_pb, p_pb = cross(*PB)
cs_ti, p_ti = cross(*TI)
ext = [-W*DX, W*DX, z[-1], z[0]]                 # in-plane (A) ; depth top=entrance

fig = plt.figure(figsize=(12, 9))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.4], wspace=0.35)

for ax, cs, name in [(fig.add_subplot(gs[0]), cs_pb, "Pb (lead) column"),
                     (fig.add_subplot(gs[1]), cs_ti, "Ti column")]:
    vmax = np.percentile(cs, 99.5)
    ax.imshow(cs, extent=ext, aspect="auto", cmap="inferno", vmin=np.percentile(cs, 5), vmax=vmax)
    for zp in np.arange(z[0], z[-1], PERIOD):     # 3.9 A plane ticks
        ax.plot([-W*DX*0.95], [zp], marker="_", color="cyan", ms=10, mew=1.5)
    ax.set_xlabel("in-plane (A)"); ax.set_title(name)
    ax.set_ylabel("depth z (A)  [entrance -> exit]")

axp = fig.add_subplot(gs[2])
axp.plot(p_pb, z, lw=1.6, color="tab:orange", label="Pb (lead)")
axp.plot(p_ti, z, lw=1.6, color="tab:cyan", label="Ti")
for zp in np.arange(z[0], z[-1], PERIOD):
    axp.axhline(zp, color="gray", ls=":", lw=0.5, alpha=0.6)
axp.invert_yaxis(); axp.set_xlabel("phase (mean-sub)"); axp.set_ylabel("depth z (A)")
axp.set_title(f"column-centre phase vs depth\n(dotted = {PERIOD} A plane period)")
axp.legend()

fig.suptitle(f"{RUN}: depth cross-section down a Pb and a Ti column "
             f"(dz={DZ:.2f} A) — cyan ticks = {PERIOD} A planes", fontsize=13)
fig.savefig(os.path.join(OUT, "column_cross_section.png"), dpi=140, bbox_inches="tight")
print("wrote column_cross_section.png   Pb", PB, " Ti", TI)
