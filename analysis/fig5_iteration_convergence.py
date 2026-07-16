#!/usr/bin/env python
"""FIG 5 — When does the reconstruction stop improving?  (1e10, single mid slice)

Reads the SAME mid-depth in-plane slice from every iteration checkpoint
(Niter 25, 50, ... 200) of the highest-dose recon and shows it as a filmstrip, plus
two convergence curves computed on that slice:
  * RMS change between consecutive checkpoints  (how much the image still moves)
  * correlation to the final (Niter 200) image   (how "done" it is)
so you can read off the iteration beyond which refinement is visually worthless.

  python fig5_iteration_convergence.py
"""
import os, glob
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patheffects as pe

import dose_fig_common as C

C.use_pub_style()
DOSE = "1e10"
LAYER = 52                                   # mid depth (~z 35 Å)
CMAP = "inferno"


def read_layer(mat, layer):
    with h5py.File(mat, "r") as f:
        a = f[f["outputs"]["object_roi"][layer, 0]][:]
        a = (a["real"] + 1j * a["imag"]) if a.dtype.names else a
        ph = np.angle(a.T).astype(np.float32)
    return ph - np.median(ph)


folder = os.path.expanduser(f"~/Desktop/dose_series/dose{DOSE}")
mats = sorted(glob.glob(os.path.join(folder, "Niter*.mat")),
              key=lambda p: int("".join(filter(str.isdigit, os.path.basename(p)))))
iters = [int("".join(filter(str.isdigit, os.path.basename(m)))) for m in mats]
slices = [read_layer(m, LAYER) for m in mats]
print("iterations:", iters)

final = slices[-1]
vmin, vmax = np.percentile(final, 3), np.percentile(final, 99.5)
rms_change = [np.nan] + [float(np.sqrt(np.mean((slices[i] - slices[i-1]) ** 2)))
                         for i in range(1, len(slices))]
corr_final = [float(np.corrcoef(s.ravel(), final.ravel())[0, 1]) for s in slices]

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(15.5, 6.4))
gs = GridSpec(2, len(mats), figure=fig, height_ratios=[1.5, 1.0],
              left=0.055, right=0.99, top=0.9, bottom=0.09, wspace=0.05, hspace=0.32)

for i, (it, sl) in enumerate(zip(iters, slices)):
    ax = fig.add_subplot(gs[0, i])
    ax.imshow(sl, cmap=CMAP, vmin=vmin, vmax=vmax, origin="lower")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"Niter {it}", fontsize=11)
    if i == len(mats) - 1:
        ax.text(0.5, -0.08, "(final)", transform=ax.transAxes, ha="center",
                va="top", fontsize=8.5, color="0.4")
ax0 = fig.axes[0]
sb = 5.0 / 0.04916
ax0.plot([12, 12 + sb], [16, 16], "-", color="white", lw=3,
         path_effects=[pe.withStroke(linewidth=5, foreground="black")])
ax0.text(12 + sb / 2, 26, "5 Å", color="white", ha="center", va="bottom", fontsize=8.5,
         path_effects=[pe.withStroke(linewidth=2.5, foreground="black")])

# convergence curves
axc = fig.add_subplot(gs[1, :4])
axc.plot(iters[1:], rms_change[1:], "-o", color="#b2182b", lw=2, ms=6)
axc.set_xlabel("iteration"); axc.set_ylabel("RMS change\nvs previous", color="#b2182b")
axc.tick_params(axis="y", labelcolor="#b2182b")
axc.set_title("How much the slice still changes", fontsize=10.5)
axc.grid(alpha=0.25)

axk = fig.add_subplot(gs[1, 4:])
axk.plot(iters, corr_final, "-s", color="#2166ac", lw=2, ms=6)
axk.axhline(0.99, color="0.6", ls="--", lw=1)
axk.text(iters[1], 0.9905, "0.99", fontsize=8, color="0.4", va="bottom")
axk.set_xlabel("iteration"); axk.set_ylabel("correlation\nto final", color="#2166ac")
axk.tick_params(axis="y", labelcolor="#2166ac")
axk.set_title("How close to the final image", fontsize=10.5)
axk.grid(alpha=0.25)

# annotate the "good enough" iteration (first corr-to-final >= 0.99)
good = next((it for it, cc in zip(iters, corr_final) if cc >= 0.99), iters[-1])
axk.annotate(f"≥0.99 by Niter {good}", (good, 0.99), textcoords="offset points",
             xytext=(8, -18), fontsize=8.5, color="#2166ac",
             arrowprops=dict(arrowstyle="->", color="#2166ac", lw=1))

fig.suptitle(f"Convergence of the {C.DOSE_TeX[DOSE]} e/Å² reconstruction — same mid-depth "
             f"slice (z ≈ {(LAYER+0.5)*0.666:.0f} Å) vs iteration", fontsize=13, y=0.965)
C.savefig(fig, "fig5_iteration_convergence")
print("RMS change:", [round(x, 4) if x == x else None for x in rms_change])
print("corr to final:", [round(x, 4) for x in corr_final], "| good-enough Niter:", good)
