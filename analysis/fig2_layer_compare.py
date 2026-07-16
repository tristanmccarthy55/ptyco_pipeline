#!/usr/bin/env python
"""FIG 2 — In-plane depth slices: dose comparison + depth sectioning.

Produces three figures:
  fig2_layers          : dose (rows) x depth (cols) grid of in-plane phase slices.
                         Down a column = the SAME depth at falling dose; across a row
                         = optical depth-sectioning at fixed dose.
  fig2_depth_montage   : a filmstrip through depth at the best dose (1e10) — the 3-D
                         object really is a stack of resolved planes.
  fig2_slice_overlay   : one mid-depth slice at 1e10 with the ground-truth atoms of
                         that plane overlaid (Pb/Sr/Ti/O) — the atoms are real.

  python fig2_layer_compare.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.gridspec import GridSpec

import dose_fig_common as C

C.use_pub_style()
CMAP = "inferno"

# depths to show (avoid the entrance/exit artifact bands); layers over 105, dz=0.666
ref, dz, dx = C.load_dose("1e10")
nL = ref.shape[0]
DEPTH_LAYERS = [16, 38, 60, 82]                       # z ~ 11, 25, 40, 55 Å
zA = lambda l: (l + 0.5) * dz


def robust(sl, lo=2, hi=99.6):
    return np.percentile(sl, lo), np.percentile(sl, hi)


# ============================================================ FIG 2a: dose x depth grid
# ONE global scale from the clean reference -> exact colorbar & honest dose comparison
VMIN, VMAX = robust(np.stack([ref[l] for l in DEPTH_LAYERS]))

fig = plt.figure(figsize=(11.2, 11.6))
gs = GridSpec(len(C.DOSES), len(DEPTH_LAYERS), figure=fig,
              left=0.075, right=0.88, top=0.92, bottom=0.035, wspace=0.06, hspace=0.06)
im = None
for r, d in enumerate(C.DOSES):
    V, dzi, dxi = C.load_dose(d)
    for c, l in enumerate(DEPTH_LAYERS):
        ax = fig.add_subplot(gs[r, c])
        im = ax.imshow(V[l], cmap=CMAP, vmin=VMIN, vmax=VMAX, origin="lower")
        ax.set_xticks([]); ax.set_yticks([])
        if r == 0:
            ax.set_title(f"z ≈ {zA(l):.0f} Å", fontsize=11)
        if c == 0:
            ax.set_ylabel(f"{C.DOSE_TeX[d]} e/Å²", fontsize=11.5)
    del V

# in-plane scale bar (5 Å) on the top-left panel
ax0 = fig.axes[0]
sb = 5.0 / dx
ax0.plot([12, 12 + sb], [18, 18], "-", color="white", lw=3,
         path_effects=[pe.withStroke(linewidth=5, foreground="black")])
ax0.text(12 + sb / 2, 30, "5 Å", color="white", ha="center", va="bottom", fontsize=9,
         path_effects=[pe.withStroke(linewidth=2.5, foreground="black")])

cax = fig.add_axes([0.895, 0.035, 0.02, 0.885])
cb = fig.colorbar(im, cax=cax)
cb.set_label("reconstructed phase (rad)")
fig.suptitle("In-plane slices vs dose and depth — atoms sharpen with dose, "
             "the structure sections with depth", fontsize=13, y=0.955)
C.savefig(fig, "fig2_layers")
plt.close(fig)

# ============================================================ FIG 2b: depth montage (1e10)
MONT = np.linspace(12, 94, 8).astype(int)             # 8 depths through the crystal
fig = plt.figure(figsize=(15.5, 2.85))
gs = GridSpec(1, len(MONT), figure=fig, left=0.015, right=0.995, top=0.80, bottom=0.02,
              wspace=0.05)
for i, l in enumerate(MONT):
    ax = fig.add_subplot(gs[0, i])
    vmin, vmax = robust(ref[l])
    ax.imshow(ref[l], cmap=CMAP, vmin=vmin, vmax=vmax, origin="lower")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"z ≈ {zA(l):.0f} Å", fontsize=11)
    for s in ax.spines.values():
        s.set_edgecolor("0.6"); s.set_linewidth(0.6)
ax0 = fig.axes[0]
sb = 5.0 / dx
ax0.plot([12, 12 + sb], [16, 16], "-", color="white", lw=3,
         path_effects=[pe.withStroke(linewidth=5, foreground="black")])
ax0.text(12 + sb / 2, 26, "5 Å", color="white", ha="center", va="bottom", fontsize=8.5,
         path_effects=[pe.withStroke(linewidth=2.5, foreground="black")])
fig.suptitle("Optical depth-sectioning at 10$^{10}$ e/Å² — slices stepping "
             r"entrance $\rightarrow$ exit through the crystal", fontsize=13, y=0.97)
C.savefig(fig, "fig2_depth_montage")
plt.close(fig)

# ============================================================ FIG 2c: slice + GT overlay
# Overlay a Pb (A-site) plane: the recon's bright columns are the heavy A-site columns,
# so Pb/Sr markers land on the bright dots (depth res ~3.9 Å can't split the A-B 1.95 Å).
pos, Z = C.load_gt_model()
OFF, resid = C.hero_depth_offset(ref, dz, dx, *C.pick_columns(ref, 1, margin=55)[0], pos, Z)
dX0, dY0 = C.inplane_shift(ref, dx, pos, Z)
npb = lambda l: len(C.slice_atoms(pos, Z, l, dz, dx, OFF, z_tol=0.9,
                                  dX0=dX0, dY0=dY0).get(82, ([],))[0])
L = max(range(35, 70), key=npb)                          # mid-depth Pb-plane layer
atoms = C.slice_atoms(pos, Z, L, dz, dx, OFF, z_tol=1.2, dX0=dX0, dY0=dY0)
stroke = [pe.withStroke(linewidth=2.2, foreground="black")]

fig, axes = plt.subplots(1, 2, figsize=(12.6, 6.6))
vmin, vmax = robust(ref[L])
for ax, overlay in zip(axes, [False, True]):
    ax.imshow(ref[L], cmap=CMAP, vmin=vmin, vmax=vmax, origin="lower")
    ax.set_xticks([]); ax.set_yticks([])
    if overlay:
        for zz, st in C.SPECIES.items():
            if zz not in atoms:
                continue
            cc, rr = atoms[zz]
            m = (cc > 2) & (cc < 402) & (rr > 2) & (rr < 402)
            line = st["m"] in ("x", "+")
            sc = ax.scatter(cc[m], rr[m], s=46, marker=st["m"], linewidths=1.5,
                            c=st["c"] if line else None,
                            facecolors=(None if line else
                                        ("none" if st["m"] == "o" else st["c"])),
                            edgecolors=(None if line else st["c"]),
                            label=st["label"])
            sc.set_path_effects(stroke)
        ax.legend(loc="upper right", ncol=4, columnspacing=0.8, handletextpad=0.2,
                  framealpha=0.9, facecolor="white", edgecolor="none")
        ax.set_title("with ground-truth atoms of this (Pb) plane", fontsize=11)
    else:
        ax.set_title("reconstructed phase", fontsize=11)
sb = 5.0 / dx
axes[0].plot([12, 12 + sb], [16, 16], "-", color="white", lw=3,
             path_effects=[pe.withStroke(linewidth=5, foreground="black")])
axes[0].text(12 + sb / 2, 26, "5 Å", color="white", ha="center", va="bottom", fontsize=9,
             path_effects=[pe.withStroke(linewidth=2.5, foreground="black")])
fig.suptitle(f"A single reconstructed depth-slice at z ≈ {zA(L):.0f} Å (10$^{{10}}$ e/Å²) "
             f"— the bright dots are atomic columns", fontsize=13)
fig.tight_layout()
C.savefig(fig, "fig2_slice_overlay")
plt.close(fig)
print("fig2 done. showcase layer", L, "z≈", round(zA(L), 1), "Å  OFF", round(OFF, 2))
