#!/usr/bin/env python
"""@file make_figs.py
@brief Generate the figures for sim_recon_methods.tex.

Writes vector PDFs (+ PNG previews) into paper/figs/:
  sample          the PbTiO3/SrTiO3 model: ferroelectric domains + atomic columns & probe
  technique       4D-STEM acquisition (side view) + the ptychographic overlap principle
  debye_waller    coherent fraction e^{-2M} vs angle, 4 species (real RT B), diffuse band
  ptycho_inverse  the LSQ-ML inverse loop (probe+slices -> multislice -> compare -> update)

  ~/hyperspy-bundle/bin/python analysis/atomfind/paper/make_figs.py

The sample figure loads the real structure via the simulation's own load_and_prepare_atoms
(needs abtem+ase -> use the hyperspy-bundle interpreter). Debye-Waller is computed from
tabulated room-temperature B factors.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Polygon, FancyArrowPatch, Circle

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
FIGS = os.path.join(HERE, "figs"); os.makedirs(FIGS, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 8.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.edgecolor": "#333333", "text.color": "#1a1a1a",
    "axes.labelcolor": "#1a1a1a", "xtick.color": "#333333", "ytick.color": "#333333",
})
C_Pb, C_O, C_Ti, C_Sr = "#D55E00", "#0072B2", "#009E73", "#E69F00"   # Okabe-Ito (CVD-safe)
INK, MUTED, PANEL, ACCENT = "#1a1a1a", "#666666", "#ececec", "#0072B2"
LAMBDA = 0.01969


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), bbox_inches="tight", dpi=300, facecolor="white")
    plt.close(fig); print(f"wrote figs/{name}.pdf (+ .png)")


# ============================================================ Fig: the sample
def _polarisation():
    """Ti off-centring vs its O cage = the true local polarisation (the ground-truth
    notebook method), in the prepared beam=z frame. Returns vec(N,3), loc(N,3), pos, syms, S."""
    sys.path.insert(0, os.path.join(REPO, "sim"))
    import simulate_4dstem as S
    from scipy.spatial import KDTree
    atoms, _ = S.load_and_prepare_atoms()
    pos = atoms.get_positions(); syms = np.array(atoms.get_chemical_symbols())
    Ti = pos[syms == "Ti"]; O = pos[syms == "O"]
    tree = KDTree(O); vec = []; loc = []
    for t in Ti:
        d, ind = tree.query(t, k=6)
        if d.mean() > 2.8:
            continue
        vec.append(t - O[ind].mean(0)); loc.append(t)
    return np.array(vec), np.array(loc), pos, syms, S


def fig_sample():
    from scipy.interpolate import griddata
    from scipy.ndimage import gaussian_filter
    vec, loc, pos, syms, S = _polarisation()
    cx, cy, win = S.SCAN_CENTER_X_A, S.SCAN_CENTER_Y_A, S.SCAN_WINDOW_A
    ewr = (74.0 + 20.0) * np.tan(0.1); halo = win / 2 + ewr        # exit-wave broadening radius

    fig, (a, b) = plt.subplots(1, 2, figsize=(6.8, 3.5),
                               gridspec_kw=dict(width_ratios=[1, 1], wspace=0.26))

    # (a) ONE mean-direction arrow per atomic column (unit arrows -> clean grid, not a hairball),
    #     over a transparent filled scan window + exit-wave halo (fills read better than outlines).
    zc = pos[:, 2].mean()
    key = np.round(loc[:, :2]).astype(int)                    # group Ti into columns (~3.9 Å apart)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    U = np.zeros(len(uniq)); V = np.zeros(len(uniq)); C = np.zeros(len(uniq))
    np.add.at(U, inv, vec[:, 0]); np.add.at(V, inv, vec[:, 1]); np.add.at(C, inv, 1)
    U /= C; V /= C; mg = np.hypot(U, V) + 1e-9
    m = (syms == "Pb") & (pos[:, 2] < zc); a.scatter(pos[m, 0], pos[m, 1], s=6, color="#dcd4ca", ec="none", zorder=1)
    a.add_patch(Circle((cx, cy), halo, fc="#e0954a", ec="none", alpha=0.14, zorder=0))
    a.add_patch(Rectangle((cx - win / 2, cy - win / 2), win, win, fc="#3b7dd8", ec="none", alpha=0.20, zorder=2))
    a.quiver(uniq[:, 0], uniq[:, 1], U / mg, V / mg, scale=30, width=0.005, headwidth=4,
             pivot="mid", color="#111111", zorder=3)
    _bb = dict(fc="white", ec="none", alpha=0.8, pad=1)
    a.annotate("scan window", xy=(cx, cy - win / 2), xytext=(29, 3.5), fontsize=6.5, color="#1a4a8a",
               ha="center", bbox=_bb, arrowprops=dict(arrowstyle="-|>", color="#1a4a8a", lw=0.9))
    a.annotate("beam exit halo", xy=(cx + halo * 0.72, cy - halo * 0.72), xytext=(55, 4), fontsize=6.5,
               color="#a05a17", ha="center", bbox=_bb, arrowprops=dict(arrowstyle="-|>", color="#a05a17", lw=0.9))
    a.set_xlim(3, 67); a.set_ylim(3, 67); a.set_aspect("equal")
    a.set_xlabel("x (Å)"); a.set_ylabel("y (Å)")
    a.set_title("(a) polarisation map + scan footprint", fontsize=7.6, loc="left", pad=3)

    # (b) atomic columns (beam projection) in the scan region + the probe -> Pb polar dumbbells
    zoom = 11
    m0 = (np.abs(pos[:, 0] - cx) < zoom + 1) & (np.abs(pos[:, 1] - cy) < zoom + 1) & (pos[:, 2] < zc)
    for el, cc, ss in [("Pb", C_Pb, 34), ("Sr", C_Sr, 26), ("Ti", C_Ti, 16), ("O", C_O, 5)]:
        mm = m0 & (syms == el); b.scatter(pos[mm, 0], pos[mm, 1], s=ss, color=cc, alpha=0.55, ec="none", label=el)
    b.add_patch(Rectangle((cx - win / 2, cy - win / 2), win, win, fc="none", ec=INK, lw=1.0, ls="--"))
    b.add_patch(Circle((cx, cy), 2.0, fc=ACCENT, alpha=0.20, ec=ACCENT, lw=1.4))
    b.annotate("probe\n(~4 Å)", xy=(cx + 1.4, cy + 1.4), xytext=(cx + 5, cy + 6), fontsize=6.8,
               color="#0a4e77", ha="center", arrowprops=dict(arrowstyle="-|>", color="#0a4e77", lw=0.8))
    b.set_aspect("equal"); b.set_xlim(cx - zoom, cx + zoom); b.set_ylim(cy - zoom, cy + zoom)
    b.set_xlabel("x (Å)"); b.set_title("(b) atomic columns + probe", fontsize=8.0, loc="left", pad=3)
    b.legend(loc="upper right", fontsize=6.4, handletextpad=0.2, borderpad=0.3, framealpha=0.85, markerscale=1.3)
    save(fig, "sample")


# ============================================================ Fig: the tornado (3-D depth-evolving polarisation)
def fig_tornado():
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers the 3d projection)
    vec, loc, pos, syms, S = _polarisation()
    cx, cy, win = S.SCAN_CENTER_X_A, S.SCAN_CENTER_Y_A, S.SCAN_WINDOW_A
    mag = np.linalg.norm(vec, axis=1)
    # the SCANNED sub-volume (x,y within the scan window; full depth): the polarisation twists
    # into a vortex tube through the beam direction -> exactly what a depth-resolved recon recovers.
    sub = (np.abs(loc[:, 0] - cx) < win / 2 - 1) & (np.abs(loc[:, 1] - cy) < win / 2 - 1) & (mag < 0.7) & (mag > 0.05)
    s = np.where(sub)[0]; s = s[::max(1, len(s) // 320)]
    zc = loc[s, 2]; cols = plt.cm.viridis((zc - zc.min()) / max(np.ptp(zc), 1e-9))

    fig = plt.figure(figsize=(3.9, 3.9)); ax = fig.add_subplot(111, projection="3d")
    ax.quiver(loc[s, 0], loc[s, 1], loc[s, 2], vec[s, 0], vec[s, 1], vec[s, 2],
              length=2.6, normalize=True, colors=cols, lw=0.7)
    ax.set_xlim(29, 51); ax.set_ylim(9, 31); ax.set_zlim(5, 69)      # tight -> less empty box
    ax.set_xticks([30, 40, 50]); ax.set_yticks([10, 20, 30]); ax.set_zticks([10, 30, 50])
    ax.set_xlabel("x (Å)", fontsize=7); ax.set_ylabel("y (Å)", fontsize=7)
    ax.set_zlabel("z = depth (Å)", fontsize=7); ax.tick_params(labelsize=6)
    try:
        ax.set_box_aspect((1, 1, 1.3), zoom=1.12)     # zoom trims the 3-D whitespace; leaves room for labels
    except TypeError:
        ax.set_box_aspect((1, 1, 1.3))
    ax.view_init(elev=18, azim=-68)
    # depth is the z-axis AND the arrow colour, so no colourbar is needed (avoids the z-axis clash)
    ax.set_title("polar vortex twisting through depth\n(scanned volume; colour = depth)", fontsize=8.2)
    for _ext in ("pdf", "png"):                        # small pad so mplot3d axis labels aren't cropped
        fig.savefig(os.path.join(FIGS, f"tornado.{_ext}"), bbox_inches="tight", pad_inches=0.28,
                    dpi=300, facecolor="white")
    plt.close(fig); print("wrote figs/tornado.pdf (+ .png)")


# ============================================================ Fig: technique (4D-STEM + depth parallax)
def _ripple(ax, x, y, col, r0=0.011, rings=(0.024, 0.038)):
    ax.add_patch(Circle((x, y), r0, color=col, zorder=4))
    for rr in rings:
        ax.add_patch(Circle((x, y), rr, fc="none", ec=col, lw=0.6, alpha=0.85, zorder=4))


def fig_technique():
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.9, 3.35),
                               gridspec_kw=dict(width_ratios=[0.92, 1.12], wspace=0.02))
    for ax in (a, b):
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_aspect("equal")                         # square panels -> square detector/pixels
        for s in ax.spines.values(): s.set_visible(False)

    # (a) 4D-STEM: ONE symmetric probe cone (same angle above/below crossover), square detector
    a.text(0.02, 0.985, "(a) 4D-STEM acquisition", fontsize=8.2, color=INK)
    xc, ycross, slope = 0.36, 0.80, 0.52
    ytop, yent, yexit, ydet = 0.95, 0.64, 0.42, 0.34
    hw = lambda y: slope * abs(ycross - y)               # half-width = angle x distance (constant angle)
    for s_ in (-1, 1):                                    # cone edges: straight lines through the crossover
        a.plot([xc + s_ * hw(ytop), xc, xc + s_ * hw(ydet)], [ytop, ycross, ydet], color=ACCENT, lw=1.1)
    a.add_patch(Polygon([(xc - hw(ytop), ytop), (xc + hw(ytop), ytop), (xc, ycross)], color=ACCENT, alpha=0.12, ec="none"))
    a.add_patch(Polygon([(xc, ycross), (xc - hw(ydet), ydet), (xc + hw(ydet), ydet)], color=ACCENT, alpha=0.10, ec="none"))
    a.add_patch(Rectangle((0.05, yexit), 0.60, yent - yexit, fc="#dfe6ee", ec="#9a9a9a", lw=0.7))   # sample slab
    a.text(0.67, (yent + yexit) / 2, "sample\n~70 Å", fontsize=6.8, va="center", color=INK)
    a.annotate("", xy=(xc - 0.135, ycross), xytext=(xc - 0.135, yent),
               arrowprops=dict(arrowstyle="<->", color="#0a4e77", lw=0.8))
    a.text(xc - 0.155, (ycross + yent) / 2, "overfocus\n20 Å", fontsize=6.3, color="#0a4e77", ha="right", va="center")
    a.annotate("", xy=(xc + hw(ytop), 0.90), xytext=(xc - hw(ytop), 0.90),
               arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.8))
    a.text(xc + hw(ytop) + 0.02, 0.90, "scan", fontsize=6.4, color=INK, ha="left", va="center")
    # pixelated detector: a square grid of pixels showing a (synthetic) diffraction pattern
    ds = 0.26; dx0, dy0 = xc - ds / 2, 0.04; npx = 11; pw = ds / npx
    yy, xx = np.mgrid[0:npx, 0:npx]; rr = np.hypot(xx - (npx - 1) / 2, yy - (npx - 1) / 2)
    dp = np.exp(-(rr / 1.5) ** 2) + 0.35 * np.exp(-((rr - 3.7) / 0.9) ** 2); dp /= dp.max()
    for i in range(npx):
        for j in range(npx):
            a.add_patch(Rectangle((dx0 + i * pw, dy0 + j * pw), pw * 0.86, pw * 0.86,
                                  fc=plt.cm.inferno(dp[j, i]), ec="none"))          # gaps -> pixels
    a.annotate("", xy=(xc, dy0 + ds + 0.004), xytext=(xc, ydet - 0.004), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.7))
    a.text(dx0 + ds + 0.015, dy0 + ds / 2, "pixelated\ndetector:\nfull DP\nper position", fontsize=6.0, va="center", color=INK)

    # (b) DEPTH BY PARALLAX: two atoms at different z shift at different rates between adjacent DPs.
    # A Ronchigram is a shadow image with the viewpoint at the crossover, so the atom NEARER the
    # surface (nearer the crossover) sweeps FASTER as the probe steps (lit: parallax depth sectioning).
    b.text(0.02, 0.985, "(b) depth by parallax", fontsize=8.2, color=INK)
    yS, yShal, yDeep, yDet = 0.86, 0.66, 0.50, 0.22
    S1, S2 = 0.43, 0.57                                   # two probe-crossover positions (step delta)
    shadow = lambda xs, ya: xs + (0.5 - xs) * (yS - yDet) / (yS - ya)   # atom = pivot -> shadow on detector
    b.plot([0.5, 0.5], [yDet, yS], color=MUTED, lw=0.6, ls=":")         # optic axis
    for xs, dsh in [(S1, False), (S2, True)]:                            # rays: each source through each atom
        for ya, col in [(yShal, C_Pb), (yDeep, C_O)]:
            xsh = shadow(xs, ya)
            b.plot([xs, 0.5, xsh], [yS, ya, yDet], color=col, lw=0.7,
                   ls="--" if dsh else "-", alpha=0.9, zorder=3)
            b.add_patch(Circle((xsh, yDet), 0.010, color=col, zorder=6))
    b.add_patch(Rectangle((0.13, yDet - 0.028), 0.74, 0.028, fc="#0f1216", ec="#9a9a9a", lw=0.5, zorder=4))
    b.text(0.875, yDet - 0.014, "detector", fontsize=6.0, va="center", ha="left", color=MUTED)
    b.add_patch(Circle((S1, yS), 0.013, color=INK, zorder=7))            # crossover pos 1 (solid)
    b.add_patch(Circle((S2, yS), 0.013, fc="white", ec=INK, lw=1.0, zorder=7))   # pos 2 (open)
    b.annotate("", xy=(S2, yS + 0.055), xytext=(S1, yS + 0.055), arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.7))
    b.text(0.5, yS + 0.085, r"probe crossover (step $\delta$)", fontsize=6.3, ha="center", color=INK)
    b.add_patch(Circle((0.5, yShal), 0.014, color=C_Pb, zorder=7))
    b.text(0.5 + 0.028, yShal, "near-surface", fontsize=6.0, va="center", color=C_Pb)
    b.add_patch(Circle((0.5, yDeep), 0.014, color=C_O, zorder=7))
    b.text(0.5 + 0.028, yDeep, "deep", fontsize=6.0, va="center", color=C_O)
    for ya, col, yb in [(yShal, C_Pb, yDet - 0.070), (yDeep, C_O, yDet - 0.115)]:   # shadow-shift spans
        x1, x2 = shadow(S1, ya), shadow(S2, ya)
        b.annotate("", xy=(max(x1, x2), yb), xytext=(min(x1, x2), yb), arrowprops=dict(arrowstyle="<->", color=col, lw=1.1))
    b.text(0.5, 0.045, "near-surface shadow shifts more", fontsize=6.2, color=C_Pb, ha="center")
    b.text(0.5, 0.008, r"$\Rightarrow$ shift rate encodes depth", fontsize=6.8, color=INK, ha="center")
    save(fig, "technique")


# ============================================================ Fig: Debye-Waller
def fig_debye_waller():
    B = {"Pb": 0.90, "O": 0.80, "Sr": 0.55, "Ti": 0.45}
    col = {"Pb": C_Pb, "O": C_O, "Sr": C_Sr, "Ti": C_Ti}
    theta = np.linspace(0, 120, 400); s = (theta * 1e-3) / (2 * LAMBDA)
    coh = {el: np.exp(-2 * B[el] * s ** 2) for el in B}
    fig, ax = plt.subplots(figsize=(3.5, 2.9))
    ax.fill_between(theta, 0, coh["O"], color=C_O, alpha=0.22, lw=0)
    ax.fill_between(theta, coh["O"], 1, color="#b8b8b8", alpha=0.30, lw=0)
    ax.text(88, 0.62, "diffuse\n(TDS)", color=MUTED, fontsize=7.2, ha="center", va="center")
    ax.text(20, 0.12, "O coherent", color=C_O, fontsize=7.2, ha="center")
    for el in ["Ti", "Sr", "Pb", "O"]:
        ax.plot(theta, coh[el], color=col[el], lw=1.4)
    for el, tt in [("Ti", 78), ("Sr", 62), ("Pb", 40), ("O", 46)]:
        yy = np.interp(tt, theta, coh[el])
        ax.annotate(el, (tt, yy), color=col[el], fontsize=7.6, fontweight="bold",
                    xytext=(3, 3), textcoords="offset points")
    ax.set_xlim(0, 120); ax.set_ylim(0, 1.0)
    ax.set_xlabel("scattering angle (mrad)"); ax.set_ylabel(r"coherent fraction $e^{-2M}$")
    ax.grid(True, lw=0.4, color="#e6e6e6"); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    save(fig, "debye_waller")


# ============================================================ Fig: ptycho inverse loop
def _box(ax, x, y, w, h, text, fc=PANEL, ec="#9a9a9a", fs=7.6, tc=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.02",
                                fc=fc, ec=ec, lw=0.8))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=tc)


def _arrow(ax, p0, p1, style="-|>", color=None, lw=1.1, rad=0.0):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=10,
                                 color=color or INK, lw=lw, connectionstyle=f"arc3,rad={rad}"))


def fig_ptycho_inverse():
    fig, ax = plt.subplots(figsize=(6.6, 2.75)); ax.set_xlim(0, 1); ax.set_ylim(-0.2, 1.0)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    _box(ax, 0.02, 0.60, 0.17, 0.24, r"probe $P$")
    _box(ax, 0.02, 0.16, 0.17, 0.24, "")
    for xx in np.linspace(0.055, 0.15, 5):
        ax.add_patch(Rectangle((xx, 0.20), 0.006, 0.16, color="#b9c7d6", ec="#8fa3b6", lw=0.3))
    ax.text(0.105, 0.11, r"object slices $\{q_l\}$", ha="center", fontsize=7.4, color=INK)
    _box(ax, 0.33, 0.38, 0.19, 0.28, "multislice\n" + r"$\mathcal{M}$", fc="#eef3ef", ec="#9ab0a2")
    _box(ax, 0.60, 0.38, 0.17, 0.28, "model\n" + r"$|\mathcal{F}\psi_j|^2$")
    _box(ax, 0.82, 0.38, 0.16, 0.28, "measured\namplitude", fc="#111318", ec="#333", tc="white")
    _arrow(ax, (0.19, 0.70), (0.35, 0.60)); _arrow(ax, (0.19, 0.30), (0.35, 0.44))
    _arrow(ax, (0.52, 0.52), (0.60, 0.52)); _arrow(ax, (0.77, 0.52), (0.82, 0.52))
    ax.text(0.795, 0.57, "compare", fontsize=6.6, color=MUTED, ha="center")
    _arrow(ax, (0.90, 0.38), (0.195, 0.26), color=C_O, lw=1.1, rad=-0.20)
    ax.text(0.52, -0.12, r"residual $\;\rightarrow\;$ least-squares update of $\{q_l\}$ and $P$",
            ha="center", fontsize=7.4, color=C_O)
    ax.text(0.5, 0.93, r"minimise $\;\mathcal{D}=\sum_{j,\mathbf{k}}\left(\sqrt{I^{\mathrm{model}}}-\sqrt{I^{\mathrm{meas}}}\right)^2$",
            ha="center", fontsize=8.0, color=INK)
    save(fig, "ptycho_inverse")


if __name__ == "__main__":
    fig_sample()
    fig_tornado()
    fig_technique()
    fig_debye_waller()
    fig_ptycho_inverse()
    print("done ->", FIGS)
