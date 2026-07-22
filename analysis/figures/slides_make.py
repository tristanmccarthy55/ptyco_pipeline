#!/usr/bin/env python
"""@file slides_make.py
@brief Presentation (PowerPoint) versions of the NL70 figures -- clean white slide style.

Bigger fonts, simple titles, no dense captions (talk detail is spoken). Reuses the exact
3-D renders, depth registration and projected potential from fig7_nl70_3d and
column_cross_section_overlay, so nothing is recomputed differently.

Outputs -> ~/Desktop/dose_series/slides/: slide_atoms_overlay (7b side + 7c end-on),
slide_depth_columns (Pb & Ti depth cross-sections + GT), slide_potential (recon vs GT).

  python figures/slides_make.py
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import pyvista as pv

# Importing these modules runs their renders once and hands us the trimmed 3-D image arrays
# + the depth registration. fig7 is a sibling in figures/; column_cross_section_overlay
# stays one level up in analysis/, so put that on the path too.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # analysis/
import fig7_nl70_3d as F                 # F.img_recon, F.img_gt, F.img_ov, F.img_db
import column_cross_section_overlay as Q # Q.V, Q.column_atoms, Q.SGN/OFF/CAL_X, ...

OUT = os.path.expanduser("~/Desktop/dose_series/slides")
os.makedirs(OUT, exist_ok=True)

_G = None


def _fig8():
    """Lazy, cached import of the potential script (builds the abtem 3-D potential once)."""
    global _G
    if _G is None:
        import fig8_nl70_potential as m
        _G = m
    return _G

PB, TI, O, SR = F.PB_COL, F.TI_COL, "#e74c3c", "#2ecc71"
MK_PB, MK_TI, MK_O = "#00e5ff", "#39ff14", "#ff21ff"     # markers that pop on inferno


def slide_style():
    plt.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "savefig.dpi": 200,
        "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
        "font.size": 16, "axes.titlesize": 20, "axes.labelsize": 18,
        "xtick.labelsize": 14, "ytick.labelsize": 14, "axes.linewidth": 1.1,
        "legend.fontsize": 17,
    })


def dot(color, label, ms=15):
    return Line2D([0], [0], marker="o", ls="", mfc=color, mec="none", ms=ms, label=label)


def flatten_plane(im, pct=55):
    """Remove a low-order (planar) phase ramp for display.  The NL70 recon phase carries an
    undetermined low-frequency background — a diagonal ~2.3 rad tilt (~20% of the atomic
    signal), a probe/object phase-ambiguity artifact, NOT specimen structure.  Fit a plane to
    the vacuum/background pixels (below `pct` percentile) and subtract so the background is flat
    and matches the GT panel.  (The GT projected potential itself has no such gradient.)"""
    im = im.astype(float); ny, nx = im.shape
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    m = im < np.percentile(im, pct)
    A = np.c_[xx[m], yy[m], np.ones(int(m.sum()))]
    coef, *_ = np.linalg.lstsq(A, im[m], rcond=None)
    return im - (coef[0] * xx + coef[1] * yy + coef[2])


def blank(ax, img, title=None):
    ax.imshow(img)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if title:
        ax.set_title(title, fontsize=20, pad=8)


# ============================================================ SLIDE 1: overlay + barrel
# Slide-only ghost styling (publication fig7b/7c keep the plain 0.25 ghosts).  The ground-truth
# ghosts are kept FAINT (transparent) but given a slightly LIGHTER, distinct hue so they read as
# a halo around the solid recon blob without being confused for it: Pb -> lighter azure blue,
# Ti -> lighter grey.  O is the faintest of all (a pale, barely-there ghost) to hammer home that
# the recon could not locate the oxygen sublattice at all.  Halo radius is enlarged for Pb/Ti.
GHOST_COLOR = {82: "#8fb8ff", 22: "#cfd4d8", 8: "#e74c3c", 38: "#8fe0a0"}
GHOST_OPACITY = {82: 0.30, 22: 0.30, 8: 0.16, 38: 0.30}
GHOST_RADIUS = {82: 1.16, 22: 1.24, 8: 1.0, 38: 1.16}


def _overlay_enhanced(p):
    p.enable_depth_peeling(12)
    F.add_recon(p)                                     # opaque recon Pb/Ti + error bars
    for zz, (col, lab, rad) in F.GT_STYLE.items():
        m = F.gZ == zz
        if not m.any():
            continue
        sph = pv.PolyData(F.gxyz[m]).glyph(
            geom=pv.Sphere(radius=rad * GHOST_RADIUS.get(zz, 1.16),
                           theta_resolution=16, phi_resolution=16),
            scale=False, orient=False)
        p.add_mesh(sph, color=GHOST_COLOR.get(zz, col), opacity=GHOST_OPACITY.get(zz, 0.30),
                   smooth_shading=True, show_scalar_bar=False)


def slide_atoms_overlay():
    slide_style()
    img_ov = F.render(_overlay_enhanced, win=(2200, 620), zoom=1.7)
    img_db = F.render_topdown(_overlay_enhanced)
    fig = plt.figure(figsize=(15, 6.6))
    fig.suptitle("Reconstructed atoms vs. ground-truth model", fontsize=27, y=0.965,
                 weight="medium")
    axL = fig.add_axes([0.010, 0.15, 0.635, 0.70])
    axR = fig.add_axes([0.660, 0.15, 0.330, 0.70])
    blank(axL, img_ov, "Side view — along the columns")
    blank(axR, img_db, "End-on — down the beam")
    # plain scale numbers next to the box (Å): width top-left, depth bottom-right
    axL.text(0.05, 0.93, "11", transform=axL.transAxes, fontsize=17, color="0.3",
             ha="center", va="center")
    axL.text(0.93, 0.07, "22", transform=axL.transAxes, fontsize=17, color="0.3",
             ha="center", va="center")
    # ellipses at the depth ends: this is the cropped middle section of the ~70 Å columns
    for xf, ha in [(0.015, "left"), (0.985, "right")]:
        axL.text(xf, 0.52, r"$\cdots$", transform=axL.transAxes, fontsize=30, color="0.45",
                 ha=ha, va="center")
    fig.legend(handles=[dot(PB, "Pb"), dot(TI, "Ti"), dot(O, "O")],
               loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.045),
               columnspacing=2.6, handletextpad=0.4)
    fig.text(0.5, 0.012, "solid = reconstruction      faint = ground-truth model",
             ha="center", color="0.4", fontsize=15)
    fig.savefig(f"{OUT}/slide_atoms_overlay.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote slide_atoms_overlay.png")


# ============================================================ SLIDE 2: depth columns
def slide_depth_columns():
    slide_style()
    V, W, DX, DZ = Q.V, Q.W, Q.DX, Q.DZ
    zrec, ZMAX, SGN, OFF, CAL = Q.zrec, Q.ZMAX_SHOW, Q.SGN, Q.OFF, Q.CAL_X
    ext = [-W * DX, W * DX, zrec[-1], zrec[0]]
    cols = {82: (MK_PB, "o", "Pb"), 22: (MK_TI, "x", "Ti"), 8: (MK_O, "D", "O")}
    stroke = [pe.withStroke(linewidth=2.6, foreground="black")]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 9.8))
    fig.suptitle("Locating atoms along the beam", fontsize=27, y=0.985, weight="medium")
    for ax, (yc, xc), title in zip(axes, (Q.PB, Q.TI), ("Pb column", "Ti column")):
        cs = V[:, yc-1:yc+2, xc-W:xc+W].mean(1)
        shown = cs[zrec <= ZMAX]
        ax.imshow(cs, extent=ext, aspect="auto", cmap="inferno",
                  vmin=np.percentile(shown, 5), vmax=np.percentile(shown, 99.3))
        Xc, Yc, P, Zs = Q.column_atoms(yc, xc)
        for zz, (col, mk, lab) in cols.items():
            m = Zs == zz
            if not m.any():
                continue
            xin = P[m, 0] - Xc + CAL
            zin = SGN * P[m, 2] + OFF
            keep = zin <= ZMAX
            sc = ax.scatter(xin[keep], zin[keep], s=48,
                            facecolors="none" if mk == "o" else col, edgecolors=col,
                            marker=mk, linewidths=1.5)
            sc.set_path_effects(stroke)
        ax.set_xlim(-W * DX, W * DX); ax.set_ylim(ZMAX, zrec[0])
        ax.set_xlabel("position (Å)")
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.set_title(title, pad=8)
    axes[0].set_ylabel("depth  (Å)")
    axes[1].set_yticklabels([])
    # legend on the white strip below, so the text never fights the dark image
    handles = [Line2D([], [], marker="o", ls="", mfc="none", mec=MK_PB, mew=2.2, ms=13,
                      label="Pb"),
               Line2D([], [], marker="x", ls="", mec=MK_TI, mew=3, ms=13, label="Ti"),
               Line2D([], [], marker="D", ls="", mfc=MK_O, mec=MK_O, ms=12, label="O")]
    for h in handles:
        h.set_path_effects(stroke)
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 0.02), columnspacing=2.6, handletextpad=0.4)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.9, bottom=0.16, wspace=0.06)
    fig.savefig(f"{OUT}/slide_depth_columns.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote slide_depth_columns.png")


# ============================================================ SLIDE 2b: columns vs GT truth
O_COL = (153, 198)                     # pure-oxygen column neighbouring the Ti column (18 O)


def _recon_cs(ax, yc, xc, vmin=None, vmax=None):
    """Reconstructed depth cross-section (z vs in-plane x) down a column — inferno heatmap."""
    V, W, DX, zrec, ZMAX = Q.V, Q.W, Q.DX, Q.zrec, Q.ZMAX_SHOW
    ext = [-W * DX, W * DX, zrec[-1], zrec[0]]
    cs = V[:, yc-1:yc+2, xc-W:xc+W].mean(1)
    shown = cs[zrec <= ZMAX]
    if vmin is None:
        vmin = np.percentile(shown, 5)
    if vmax is None:
        vmax = np.percentile(shown, 99.3)
    ax.imshow(cs, extent=ext, aspect="auto", cmap="inferno", vmin=vmin, vmax=vmax)
    ax.set_xlim(-W * DX, W * DX); ax.set_ylim(ZMAX, zrec[0])
    return vmin, vmax


def _truth_cs(ax, yc, xc, lo=74):
    """Ground-truth column from the abtem PROJECTED POTENTIAL (not model dots): a z-vs-x
    cross-section of the true potential down the same column, sampled from fig8's 3-D potential
    and z-registered into the recon depth frame (same SGN/OFF as the recon-column overlay).  It
    shows the true, sharp atomic planes to compare against the recon's blurred column."""
    G = _fig8()
    arr, samp, sth = G.arr, G.samp, G.sth
    W, DX, zrec, ZMAX, SGN, OFF = Q.W, Q.DX, Q.zrec, Q.ZMAX_SHOW, Q.SGN, Q.OFF
    ax.set_facecolor("black")                               # black behind any uncovered strip
    Xc = Q.X0 + xc * DX; Yc = Q.Y0 + yc * DX                # physical column centre
    gxc, gyc = Xc / samp, Yc / samp                         # -> potential pixels (gx<->X, gy<->Y)
    Wp = int(round(W * DX / samp)); bp = max(int(round(0.7 / samp)), 1)
    x0, x1 = int(round(gxc)) - Wp, int(round(gxc)) + Wp
    y0, y1 = int(round(gyc)) - bp, int(round(gyc)) + bp + 1
    cs = arr[:, x0:x1, y0:y1].mean(2)                       # [n_slice, x]  (z vs X)
    zc = SGN * (np.arange(arr.shape[0]) + 0.5) * sth + OFF  # potential slice depth -> recon frame
    ax.imshow(cs, extent=[-Wp * samp, Wp * samp, zc[-1], zc[0]], aspect="auto", cmap="inferno",
              interpolation="bilinear", vmin=np.percentile(cs, lo), vmax=np.percentile(cs, 99.7))
    ax.set_xlim(-W * DX, W * DX); ax.set_ylim(ZMAX, zrec[0])


def _ti_scale():
    V, W, zrec, ZMAX = Q.V, Q.W, Q.zrec, Q.ZMAX_SHOW
    yc, xc = Q.TI
    cs = V[:, yc-1:yc+2, xc-W:xc+W].mean(1); shown = cs[zrec <= ZMAX]
    return np.percentile(shown, 5), np.percentile(shown, 99.3)


def slide_depth_columns_gt():
    slide_style()
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 9.2))
    fig.suptitle("Locating atoms along the beam — reconstruction vs. ground truth",
                 fontsize=23, y=0.985, weight="medium")
    order = [(Q.PB, "recon", "Pb column"), (Q.PB, "truth", "Pb — ground truth"),
             (Q.TI, "recon", "Ti column"), (Q.TI, "truth", "Ti — ground truth")]
    for ax, (col, kind, title) in zip(axes, order):
        yc, xc = col
        (_recon_cs if kind == "recon" else _truth_cs)(ax, yc, xc)
        ax.set_xlabel("position (Å)"); ax.xaxis.set_major_locator(MaxNLocator(3))
        ax.set_title(title, pad=8, fontsize=18)
    axes[0].set_ylabel("depth  (Å)"); axes[0].yaxis.set_major_locator(MaxNLocator(6))
    for ax in axes[1:]:
        ax.set_yticklabels([])
    fig.subplots_adjust(left=0.06, right=0.99, top=0.9, bottom=0.09, wspace=0.14)
    fig.savefig(f"{OUT}/slide_depth_columns_gt.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote slide_depth_columns_gt.png")


def slide_depth_columns_O():
    slide_style()
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 9.2))
    fig.suptitle("Neighbouring oxygen column — reconstruction vs. ground truth",
                 fontsize=21, y=0.985, weight="medium")
    vmin, vmax = _ti_scale()                            # honest scale: same as the Ti column
    _recon_cs(axes[0], *O_COL, vmin=vmin, vmax=vmax)
    axes[0].set_title("O column — reconstruction", pad=8, fontsize=18)
    _truth_cs(axes[1], *O_COL)
    axes[1].set_title("O column — ground truth", pad=8, fontsize=18)
    axes[0].set_ylabel("depth  (Å)"); axes[0].yaxis.set_major_locator(MaxNLocator(6))
    axes[1].set_yticklabels([])
    for ax in axes:
        ax.set_xlabel("position (Å)"); ax.xaxis.set_major_locator(MaxNLocator(3))
    fig.subplots_adjust(left=0.12, right=0.98, top=0.9, bottom=0.09, wspace=0.1)
    fig.savefig(f"{OUT}/slide_depth_columns_O.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote slide_depth_columns_O.png")


# ============================================================ SLIDE 3: potential vs recon
def slide_potential():
    """Down-the-beam (x-y) only: reconstructed Σ phase vs the GT projected potential.  The
    depth cross-section row was dropped — the recon side view carried a brightness gradient
    across the band (phase-ramp artifact) that distracted from the (convincing) x-y match."""
    slide_style()
    G = _fig8()                           # runs abtem potential build once; gives arrays
    fig = plt.figure(figsize=(12.5, 7.4))
    fig.suptitle("Reconstructed phase vs. ground-truth potential", fontsize=25, y=0.975,
                 weight="medium")
    axl = fig.add_axes([0.015, 0.035, 0.475, 0.80])
    axr = fig.add_axes([0.510, 0.035, 0.475, 0.80])
    recon_flat = flatten_plane(G.recon_xy)             # kill the low-order phase ramp
    axl.imshow(recon_flat, cmap="inferno", **G.scaled(recon_flat, lo=50))
    axr.imshow(G.gt_xy, cmap="inferno", **G.scaled(G.gt_xy))
    for ax in (axl, axr):
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    axl.set_title("Reconstruction", fontsize=22, pad=8)
    axr.set_title("Ground truth", fontsize=22, pad=8)
    fig.text(0.5, 0.90, "summed along the beam  (top view)", ha="center", color="0.4",
             fontsize=15)
    fig.savefig(f"{OUT}/slide_potential.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote slide_potential.png")


if __name__ == "__main__":
    slide_atoms_overlay()
    slide_depth_columns()
    slide_depth_columns_gt()
    slide_depth_columns_O()
    slide_potential()
    # the recon-vs-GT ball slide was dropped (superseded by the overlay slide)
    old = f"{OUT}/slide_recon_vs_gt.png"
    if os.path.exists(old):
        os.remove(old); print("removed slide_recon_vs_gt.png")
    print("slides ->", OUT)
