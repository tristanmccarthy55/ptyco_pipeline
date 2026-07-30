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
def fig_sample():
    sys.path.insert(0, os.path.join(REPO, "sim"))
    import simulate_4dstem as S
    atoms, _ = S.load_and_prepare_atoms()
    p = atoms.get_positions(); sym = np.array(atoms.get_chemical_symbols())
    cx, cy, win = S.SCAN_CENTER_X_A, S.SCAN_CENTER_Y_A, S.SCAN_WINDOW_A

    fig, (a, b) = plt.subplots(1, 2, figsize=(6.8, 3.25),
                               gridspec_kw=dict(width_ratios=[1, 1], wspace=0.30))

    # (a) ferroelectric polarisation: B-site (Ti) off-centring vs the A-site cage, mid-depth slab
    from scipy.spatial import cKDTree
    Asite = p[(sym == "Pb") | (sym == "Sr")]; Ti = p[sym == "Ti"]
    _, idx = cKDTree(Asite).query(Ti, k=8)
    disp = Ti - Asite[idx].mean(1); mag = np.linalg.norm(disp, axis=1)
    zc = p[:, 2].mean(); sl = (np.abs(Ti[:, 2] - zc) < 8) & (mag < 0.45) & (mag > 0.02)
    ang = np.arctan2(disp[sl, 1], disp[sl, 0])
    a.quiver(Ti[sl, 0], Ti[sl, 1], disp[sl, 0], disp[sl, 1], ang, cmap="twilight",
             scale=6, width=0.005, pivot="mid", clim=(-np.pi, np.pi))
    a.add_patch(Rectangle((cx - win / 2, cy - win / 2), win, win, fc="none", ec=INK, lw=1.3, ls="--"))
    a.text(cx, cy - win / 2 - 2.2, "scan window", ha="center", fontsize=6.8, color=INK)
    a.set_aspect("equal"); a.set_xlim(2, 68); a.set_ylim(2, 68)
    a.set_xlabel("x (Å)"); a.set_ylabel("y (Å)")
    a.set_title("(a) ferroelectric polarisation", fontsize=8.2, loc="left", pad=3)

    # (b) atomic columns (beam projection) in the scan region + the probe
    style = [("Pb", C_Pb, 34), ("Sr", C_Sr, 26), ("Ti", C_Ti, 16), ("O", C_O, 5)]
    zoom = 11
    m0 = (np.abs(p[:, 0] - cx) < zoom + 1) & (np.abs(p[:, 1] - cy) < zoom + 1)
    for el, cc, ss in style:
        mm = m0 & (sym == el)
        b.scatter(p[mm, 0], p[mm, 1], s=ss, color=cc, alpha=0.55, ec="none", label=el)
    b.add_patch(Rectangle((cx - win / 2, cy - win / 2), win, win, fc="none", ec=INK, lw=1.0, ls="--"))
    # probe: overfocus 20 A x conv 100 mrad -> ~2 A geometric radius at the entrance
    b.add_patch(Circle((cx, cy), 2.0, fc=ACCENT, alpha=0.20, ec=ACCENT, lw=1.4))
    b.annotate("probe\n(~4 Å)", xy=(cx + 1.4, cy + 1.4), xytext=(cx + 5, cy + 6),
               fontsize=6.8, color="#0a4e77", ha="center",
               arrowprops=dict(arrowstyle="-|>", color="#0a4e77", lw=0.8))
    b.set_aspect("equal"); b.set_xlim(cx - zoom, cx + zoom); b.set_ylim(cy - zoom, cy + zoom)
    b.set_xlabel("x (Å)")
    b.set_title("(b) atomic columns + probe", fontsize=8.2, loc="left", pad=3)
    b.legend(loc="upper right", fontsize=6.4, handletextpad=0.2, borderpad=0.3,
             framealpha=0.85, markerscale=1.3)
    save(fig, "sample")


# ============================================================ Fig: technique (4D-STEM + ptychography)
def fig_technique():
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.8, 3.15),
                               gridspec_kw=dict(width_ratios=[1.05, 1.0], wspace=0.28))
    for ax in (a, b):
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        for s in ax.spines.values(): s.set_visible(False)

    # (a) 4D-STEM side view: overfocused probe -> slab -> pixelated detector records a DP
    a.text(0.02, 0.98, "(a) 4D-STEM acquisition", fontsize=8.2, color=INK)
    foc = (0.32, 0.74)                                   # crossover 2 nm ABOVE the entrance (overfocus)
    a.plot([0.14, foc[0]], [0.92, foc[1]], color=ACCENT, lw=1.0)
    a.plot([0.50, foc[0]], [0.92, foc[1]], color=ACCENT, lw=1.0)
    a.add_patch(Polygon([(0.14, 0.92), (0.50, 0.92), foc], color=ACCENT, alpha=0.12, ec="none"))
    a.add_patch(Polygon([foc, (0.16, 0.40), (0.48, 0.40)], color=ACCENT, alpha=0.10, ec="none"))
    a.plot([foc[0], 0.16], [foc[1], 0.40], color=ACCENT, lw=1.0)
    a.plot([foc[0], 0.48], [foc[1], 0.40], color=ACCENT, lw=1.0)
    a.add_patch(Rectangle((0.06, 0.40), 0.56, 0.16, fc="#dfe6ee", ec="#9a9a9a", lw=0.7))  # slab
    a.text(0.64, 0.48, "sample\n~70 Å", fontsize=6.8, color=INK, va="center")
    a.annotate("overfocus\n20 Å", xy=foc, xytext=(0.66, 0.76), fontsize=6.6, color="#0a4e77",
               ha="center", arrowprops=dict(arrowstyle="-|>", color="#0a4e77", lw=0.7))
    a.annotate("", xy=(0.50, 0.90), xytext=(0.14, 0.90),
               arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.8))            # scan
    a.text(0.53, 0.90, "scan", fontsize=6.6, color=INK, ha="left", va="center")
    a.add_patch(Rectangle((0.12, 0.06), 0.40, 0.20, fc="#111318", ec="#9a9a9a", lw=0.5))  # detector
    a.add_patch(Circle((0.32, 0.16), 0.028, color="#f2c14e"))
    a.add_patch(Circle((0.32, 0.16), 0.055, fc="none", ec="#7a6a3a", lw=0.4))
    for dx in (0.16, 0.24):
        a.add_patch(Circle((0.32, 0.16), dx * 0.6 + 0.05, fc="none", ec="#5a4a2a", lw=0.25))
    a.annotate("", xy=(0.32, 0.28), xytext=(0.32, 0.39), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.7))
    a.text(0.55, 0.16, "pixelated\ndetector:\nfull DP per\nposition", fontsize=6.6, color=INK, va="center")

    # (b) ptychography principle: overlapping probes -> interfering diffraction discs
    b.text(0.02, 0.98, "(b) ptychographic overlap", fontsize=8.2, color=INK)
    # real space: two overlapping probe discs over a row of atoms
    for x0 in (0.30, 0.44):
        b.add_patch(Circle((x0, 0.74), 0.13, fc=ACCENT, alpha=0.18, ec=ACCENT, lw=1.0))
    for ax_ in np.linspace(0.16, 0.60, 8):
        b.add_patch(Circle((ax_, 0.74), 0.012, color=INK))
    b.annotate("", xy=(0.44, 0.90), xytext=(0.30, 0.90), arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.7))
    b.text(0.37, 0.925, r"shift $\ll$ probe", fontsize=6.4, color=INK, ha="center")
    b.text(0.68, 0.74, "overlapping\nillumination", fontsize=6.6, color=INK, va="center")
    # reciprocal space: two overlapping CBED discs with interference fringes in the overlap
    c1, c2, R = (0.34, 0.30), (0.50, 0.30), 0.15
    for cc in (c1, c2):
        b.add_patch(Circle(cc, R, fc="#cfe0ee", ec=ACCENT, lw=1.0, alpha=0.7))
    # fringes in the lens-shaped overlap
    xm = 0.42
    for yy in np.linspace(0.19, 0.41, 7):
        hw = np.sqrt(max(R**2 - (yy - 0.30)**2, 0)) - (xm - c1[0])
        if hw > 0:
            b.plot([xm - hw, xm + hw], [yy, yy], color="#0a4e77", lw=0.8)
    b.annotate("interference fringes\n$\\Rightarrow$ object phase", xy=(0.42, 0.30), xytext=(0.60, 0.14),
               fontsize=6.6, color="#0a4e77", ha="center",
               arrowprops=dict(arrowstyle="-|>", color="#0a4e77", lw=0.7))
    b.text(0.16, 0.10, "diffraction plane", fontsize=6.6, color=MUTED)
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
    fig_technique()
    fig_debye_waller()
    fig_ptycho_inverse()
    print("done ->", FIGS)
