#!/usr/bin/env python
"""@file fig1_depth_resolution.py
@brief FIG 1 — Depth resolution vs dose (the headline).

Top:    depth cross-sections down the SAME Pb column at each dose, ground-truth Pb
        atomic-plane positions overlaid.  The periodic bright plane-blobs persist at
        1e10/1e8, fade at 1e6, and vanish at 1e4 — depth resolution is the first casualty.
Bottom: the quantitative backing — on-column kz power spectrum (peak at the Pb-plane
        frequency) per dose, and the plane-peak prominence vs dose curve.

Registration is fit once on the cleanest dose (1e10) and reused: all four recons share
the sim geometry, so the GT markers are FIXED and you watch the recon signal fade under
them.  Reuses everything from dose_fig_common.

  python fig1_depth_resolution.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.gridspec import GridSpec

import dose_fig_common as C

C.use_pub_style()

W = 22                                   # half in-plane strip (px) ~1.1 Å
CMAP = "magma"

# ---- reference: pick hero Pb column + registration columns, fit once ----------------
ref, dz, dx = C.load_dose("1e10")
cols = C.pick_columns(ref, n=3, margin=55)
hero = cols[0]
pos, Z = C.load_gt_model()
SGN = 1                                  # dose series is entrance -> exit
OFF, resid = C.hero_depth_offset(ref, dz, dx, *hero, pos, Z)   # robust blob match
CALX = 0.0                               # markers at true GT in-plane positions
print(f"hero Pb column {hero}  registration z_rec=z_GT{OFF:+.2f} Å (resid {resid:.2f} Å)")

f_pb = 1.0 / C.gt_plane_freqs(pos, Z)["Pb"]     # Pb-plane frequency (Å^-1)
Xc, Yc, Patoms, Zatoms = C.column_atoms(pos, Z, *hero, dx)

nL = ref.shape[0]
zrec = (np.arange(nL) + 0.5) * dz
ext = [-W * dx, W * dx, zrec[-1], zrec[0]]
stroke = [pe.withStroke(linewidth=2.2, foreground="black")]

# shared display scale from the reference cross-section (honest amplitude comparison)
_csref = ref[:, hero[0]-1:hero[0]+2, hero[1]-W:hero[1]+W].mean(1)
_shown = _csref[zrec <= C.ZMAX_SHOW]
VMIN, VMAX = np.percentile(_shown, 3), np.percentile(_shown, 99.3)

# ------------------------------------------------------------------ figure skeleton
fig = plt.figure(figsize=(13, 9.2))
gs = GridSpec(2, 4, height_ratios=[1.55, 1.0], hspace=0.28, wspace=0.12,
              left=0.06, right=0.985, top=0.9, bottom=0.08)

kz_all = {}
for i, d in enumerate(C.DOSES):
    V, dzi, dxi = C.load_dose(d)
    cs = V[:, hero[0]-1:hero[0]+2, hero[1]-W:hero[1]+W].mean(1)
    ax = fig.add_subplot(gs[0, i])
    ax.imshow(cs, extent=ext, aspect="auto", cmap=CMAP, vmin=VMIN, vmax=VMAX)

    # GT Pb (and any Ti/O) planes for this column, fixed registration
    for zz, st in C.SPECIES.items():
        m = Zatoms == zz
        if not m.any():
            continue
        xin = Patoms[m, 0] - Xc + CALX
        zin = SGN * Patoms[m, 2] + OFF
        keep = zin <= C.ZMAX_SHOW
        if not keep.any():
            continue
        sc = ax.scatter(xin[keep], zin[keep], s=42,
                        facecolors="none" if st["m"] == "o" else st["c"],
                        edgecolors=st["c"], marker=st["m"], linewidths=1.5,
                        label=st["label"] if i == 0 else None)
        sc.set_path_effects(stroke)

    ax.set_ylim(C.ZMAX_SHOW, zrec[0])
    ax.set_xlim(-W * dx, W * dx)
    ax.set_title(f"{C.DOSE_TeX[d]} e/Å²", fontsize=12, pad=6)
    ax.set_xlabel("in-plane (Å)")
    if i == 0:
        ax.set_ylabel(r"depth  z (Å)   [entrance $\rightarrow$ exit]")
        ax.legend(loc="lower left", markerscale=1.1, handletextpad=0.2,
                  labelspacing=0.25, framealpha=0.85, frameon=True,
                  facecolor="white", edgecolor="none")
    else:
        ax.set_yticklabels([])

    kz, Pcol, Pvac = C.kz_spectrum(V, dzi)
    kz_all[d] = (kz, Pcol, Pvac)

# scale bar on first panel (2 Å in-plane)
ax0 = fig.axes[0]
xb0 = -W * dx + 0.15
ax0.plot([xb0, xb0 + 1.0], [C.ZMAX_SHOW - 2.2, C.ZMAX_SHOW - 2.2], "-", color="white",
         lw=3, solid_capstyle="butt", path_effects=[pe.withStroke(linewidth=5, foreground="black")])
ax0.text(xb0 + 0.5, C.ZMAX_SHOW - 3.6, "1 Å", color="white", ha="center", va="bottom",
         fontsize=8.5, path_effects=[pe.withStroke(linewidth=2.5, foreground="black")])

# ------------------------------------------------------------------ bottom-left: kz spectra
axk = fig.add_subplot(gs[1, :2])
palette = plt.cm.viridis(np.linspace(0.12, 0.86, len(C.DOSES)))
from scipy.ndimage import uniform_filter1d
proms = []
ymax = 1.0
for col, d in zip(palette, C.DOSES):
    kz, Pcol, Pvac = kz_all[d]
    ratio = uniform_filter1d(Pcol / (Pvac + 1e-30), 2)   # on-column excess over vacuum
    axk.plot(kz, ratio, lw=2, color=col, label=f"{C.DOSE_TeX[d]}")
    ymax = max(ymax, ratio[(kz > 0.15) & (kz < 0.45)].max())
    proms.append(C.plane_peak_prominence(kz, Pcol, f_pb))
axk.axhline(1.0, color="0.7", ls="-", lw=0.8)
axk.axvline(f_pb, color="k", ls="--", lw=1.2)
axk.annotate(f"Pb planes\n{1/f_pb:.1f} Å", (f_pb, ymax * 1.03),
             ha="center", va="bottom", fontsize=8.4)
axk.set_xlim(0, 0.55)
axk.set_ylim(0, ymax * 1.18)
axk.set_xlabel(r"axial spatial frequency  $k_z$  (Å$^{-1}$)")
axk.set_ylabel("on-column / vacuum z-power")
axk.set_title("Planes oscillate down the column only at high dose", fontsize=10.5)
axk.legend(title="dose (e/Å²)", ncol=2, loc="upper left")

# ------------------------------------------------------------------ bottom-right: prominence curve
axp = fig.add_subplot(gs[1, 2:])
dvals = [float(x) for x in C.DOSES]
axp.axhline(2.0, color="0.6", ls="--", lw=1)
axp.text(1.3e4, 2.08, "resolution floor (≈2×)", fontsize=8, color="0.4", va="bottom")
axp.plot(dvals, proms, "-o", color="#b2182b", lw=2.2, ms=9, mfc="white", mew=2)
for x, y in zip(dvals, proms):
    axp.annotate(f"{y:.1f}×", (x, y), textcoords="offset points", xytext=(0, 11),
                 ha="center", fontsize=8.5)
axp.set_xscale("log")
axp.invert_xaxis()
axp.set_xlabel(r"dose  (e/Å²)   [high $\rightarrow$ low]")
axp.set_ylabel("Pb-plane peak prominence  (×)")
axp.set_title("Depth resolution collapses as dose drops", fontsize=10.5)
axp.set_ylim(0, max(proms) * 1.18)
axp.grid(True, which="both", axis="y", alpha=0.25)

fig.suptitle("Depth resolution vs electron dose — cross-section down a Pb column "
             "(markers = ground-truth Pb planes)", fontsize=13.5, y=0.965)

C.savefig(fig, "fig1_depth_resolution")
print("Pb-plane prominence by dose:", {d: round(p, 1) for d, p in zip(C.DOSES, proms)})
