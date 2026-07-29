#!/usr/bin/env python
"""@file make_figs.py
@brief Generate the four figures for sim_recon_methods.tex.

Writes vector PDFs (for LaTeX) + PNG previews into paper/figs/:
  acquisition     multislice 4D-STEM forward operator (schematic)
  debye_waller    coherent -> diffuse redistribution with angle (computed from real B)
  ptycho_inverse  the ptychographic inverse loop (schematic)
  missing_cone    reciprocal-space null space + the REAL measured Pb PSF

  python analysis/atomfind/paper/make_figs.py

The missing-cone panel embeds the real kernel ~/Desktop/psf_Pb_NL70_vol.npy; if it is
unreadable a synthetic anisotropic Gaussian stands in (a warning is printed).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Polygon, FancyArrowPatch, Circle, Wedge
from matplotlib.collections import PatchCollection

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figs"); os.makedirs(FIGS, exist_ok=True)

# ---- paper-consistent style (LaTeX-like serif + cm math) --------------------
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 8.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.edgecolor": "#333333", "text.color": "#1a1a1a",
    "axes.labelcolor": "#1a1a1a", "xtick.color": "#333333", "ytick.color": "#333333",
})
# Okabe-Ito (CVD-safe by construction) — species identities
C_Pb, C_O, C_Ti, C_Sr = "#D55E00", "#0072B2", "#009E73", "#E69F00"
INK, MUTED, PANEL, ACCENT = "#1a1a1a", "#666666", "#ececec", "#0072B2"
LAMBDA = 0.01969   # electron wavelength at 300 keV [A]


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), bbox_inches="tight",
                    dpi=300, facecolor="white")
    plt.close(fig)
    print(f"wrote figs/{name}.pdf (+ .png)")


# ============================================================ Fig 1: acquisition
def fig_acquisition():
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.0, 2.35),
                                  gridspec_kw=dict(width_ratios=[1.25, 0.9, 1.05], wspace=0.32))
    for ax in (a, b, c):
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)

    # (a) multislice + convergent probe -> exit wave
    a.set_xlim(0, 1); a.set_ylim(0, 1)
    xs = np.linspace(0.42, 0.72, 6)
    for x in xs:
        a.add_patch(Rectangle((x, 0.18), 0.012, 0.64, color="#c7c7c7", ec="#9a9a9a", lw=0.4))
    # probe cone converging to a focus just above the entrance, then diverging
    foc = (0.40, 0.5)
    a.add_patch(Polygon([(0.02, 0.86), (0.02, 0.14), foc], closed=True,
                        color=ACCENT, alpha=0.16, ec="none"))
    a.add_patch(Polygon([foc, (0.72, 0.80), (0.72, 0.20)], closed=True,
                        color=ACCENT, alpha=0.10, ec="none"))
    for y0, y1 in [(0.86, 0.80), (0.14, 0.20)]:
        a.plot([0.02, foc[0], 0.72], [y0, foc[1], y1], color=ACCENT, lw=1.0)
    a.annotate("", xy=(0.90, 0.5), xytext=(0.74, 0.5),
               arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.1))
    a.text(0.05, 0.06, "convergent\nprobe", ha="left", va="center", fontsize=7.2, color=INK)
    a.text(0.57, 0.92, "potential slices", ha="center", fontsize=7.2, color=INK)
    a.text(0.91, 0.60, "exit\nwave", ha="left", va="center", fontsize=7.2, color=INK)
    a.annotate(r"$t_n=e^{\,i\sigma v_n}$", xy=(xs[3], 0.2), xytext=(0.68, 0.05),
               fontsize=7.0, color=MUTED, ha="center",
               arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5))
    a.annotate(r"$\mathcal{P}_{\Delta z}$", xy=(0.5, 0.74), fontsize=7.2, color=MUTED, ha="center")
    a.set_title("(a) multislice propagation", fontsize=8.0, loc="left", pad=2)

    # (b) scan raster
    b.set_xlim(0, 1); b.set_ylim(0, 1)
    b.add_patch(Rectangle((0.14, 0.16), 0.72, 0.68, fc="none", ec="#9a9a9a", lw=0.7))
    gx = np.linspace(0.20, 0.80, 7); gy = np.linspace(0.22, 0.78, 7)
    for j, y in enumerate(gy):
        row = gx if j % 2 == 0 else gx[::-1]
        b.plot(row, [y] * len(row), color="#cfcfcf", lw=0.6, zorder=1)
    XX, YY = np.meshgrid(gx, gy)
    b.scatter(XX, YY, s=6, color=ACCENT, zorder=2)
    b.text(0.5, 0.05, r"scan positions $\mathbf{r}_p$  (0.1–0.15 Å step)",
           ha="center", fontsize=7.2, color=INK)
    b.set_title("(b) raster scan", fontsize=8.0, loc="left", pad=2)

    # (c) 4D data: a diffraction pattern at every position
    c.set_xlim(0, 1); c.set_ylim(0, 1)
    cx = np.linspace(0.14, 0.70, 3); cy = np.linspace(0.62, 0.18, 3)
    for x in cx:
        for y in cy:
            c.add_patch(Rectangle((x, y), 0.20, 0.20, fc="#111318", ec="#9a9a9a", lw=0.4))
            c.add_patch(Circle((x + 0.10, y + 0.10), 0.028, color="#f2c14e"))
            c.add_patch(Circle((x + 0.10, y + 0.10), 0.055, fc="none", ec="#7a6a3a", lw=0.4))
    c.annotate("", xy=(0.92, 0.16), xytext=(0.92, 0.82),
               arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.8))
    c.text(0.95, 0.5, r"$\mathbf{r}_p$", fontsize=7.5, color=MUTED, va="center")
    c.text(0.44, 0.045, r"$I(\mathbf{r}_p,\mathbf{k})$: 4D dataset", ha="center",
           fontsize=7.2, color=INK)
    c.set_title("(c) a pattern per position", fontsize=8.0, loc="left", pad=2)
    save(fig, "acquisition")


# ============================================================ Fig 2: Debye-Waller
def fig_debye_waller():
    B = {"Pb": 0.90, "O": 0.80, "Sr": 0.55, "Ti": 0.45}     # RT isotropic B [A^2]
    col = {"Pb": C_Pb, "O": C_O, "Sr": C_Sr, "Ti": C_Ti}
    theta = np.linspace(0, 120, 400)                          # scattering angle [mrad]
    s = (theta * 1e-3) / (2 * LAMBDA)                         # s = sin(theta/2)/lambda [1/A]
    coh = {el: np.exp(-2 * B[el] * s ** 2) for el in B}       # coherent surviving fraction e^{-2M}

    fig, ax = plt.subplots(figsize=(3.5, 2.9))
    # O as a stacked area: coherent (blue) shrinking, diffuse (grey) growing
    ax.fill_between(theta, 0, coh["O"], color=C_O, alpha=0.22, lw=0)
    ax.fill_between(theta, coh["O"], 1, color="#b8b8b8", alpha=0.30, lw=0)
    ax.text(88, 0.62, "diffuse\n(TDS)", color=MUTED, fontsize=7.2, ha="center", va="center")
    ax.text(20, 0.12, "O coherent", color=C_O, fontsize=7.2, ha="center")
    # species comparison: coherent fraction as thin lines
    for el in ["Ti", "Sr", "Pb", "O"]:
        ax.plot(theta, coh[el], color=col[el], lw=1.4)
    # direct labels at a readable angle
    for el, tt in [("Ti", 78), ("Sr", 62), ("Pb", 40), ("O", 46)]:
        yy = np.interp(tt, theta, coh[el])
        ax.annotate(el, (tt, yy), color=col[el], fontsize=7.6, fontweight="bold",
                    xytext=(3, 3), textcoords="offset points")
    ax.set_xlim(0, 120); ax.set_ylim(0, 1.0)
    ax.set_xlabel("scattering angle (mrad)"); ax.set_ylabel(r"coherent fraction $e^{-2M}$")
    ax.grid(True, lw=0.4, color="#e6e6e6"); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    save(fig, "debye_waller")


# ============================================================ Fig 3: ptycho loop
def _box(ax, x, y, w, h, text, fc=PANEL, ec="#9a9a9a", fs=7.6, tc=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.02",
                                fc=fc, ec=ec, lw=0.8))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=tc)


def _arrow(ax, p0, p1, style="-|>", color=None, lw=1.1, rad=0.0):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=10,
                                 color=color or INK, lw=lw,
                                 connectionstyle=f"arc3,rad={rad}"))


def fig_ptycho_inverse():
    fig, ax = plt.subplots(figsize=(6.6, 2.75)); ax.set_xlim(0, 1); ax.set_ylim(-0.2, 1.0)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

    _box(ax, 0.02, 0.60, 0.17, 0.24, r"probe $P$")
    # object box with a mini slice-stack glyph
    _box(ax, 0.02, 0.16, 0.17, 0.24, "")
    for i, xx in enumerate(np.linspace(0.055, 0.15, 5)):
        ax.add_patch(Rectangle((xx, 0.20), 0.006, 0.16, color="#b9c7d6", ec="#8fa3b6", lw=0.3))
    ax.text(0.105, 0.11, r"object slices $\{q_l\}$", ha="center", fontsize=7.4, color=INK)
    _box(ax, 0.33, 0.38, 0.19, 0.28, r"multislice" + "\n" + r"$\mathcal{M}$", fc="#eef3ef", ec="#9ab0a2")
    _box(ax, 0.60, 0.38, 0.17, 0.28, r"model" + "\n" + r"$|\mathcal{F}\psi_j|^2$")
    _box(ax, 0.82, 0.38, 0.16, 0.28, "measured\namplitude", fc="#111318", ec="#333", tc="white")

    _arrow(ax, (0.19, 0.70), (0.35, 0.60))     # probe -> M
    _arrow(ax, (0.19, 0.30), (0.35, 0.44))     # object -> M
    _arrow(ax, (0.52, 0.52), (0.60, 0.52))     # M -> model
    _arrow(ax, (0.77, 0.52), (0.82, 0.52))     # model -> measured (compare)
    ax.text(0.795, 0.57, "compare", fontsize=6.6, color=MUTED, ha="center")
    # feedback: residual updates probe & object
    _arrow(ax, (0.90, 0.38), (0.195, 0.26), color=C_O, lw=1.1, rad=-0.20, style="-|>")
    ax.text(0.52, -0.12, r"residual $\;\rightarrow\;$ least-squares update of $\{q_l\}$ and $P$",
            ha="center", fontsize=7.4, color=C_O)
    ax.text(0.5, 0.93, r"minimise $\;\mathcal{D}=\sum_{j,\mathbf{k}}\left(\sqrt{I^{\mathrm{model}}}-\sqrt{I^{\mathrm{meas}}}\right)^2$",
            ha="center", fontsize=8.0, color=INK)
    save(fig, "ptycho_inverse")


# ============================================================ Fig 4: missing cone + PSF
def _load_psf():
    # prefer a repo-local copy (Desktop is often OS-sandbox-blocked); then Desktop
    for p in (os.path.join(FIGS, "psf_Pb_NL70_vol.npy"),
              os.path.join(HERE, "psf_Pb_NL70_vol.npy"),
              os.path.expanduser("~/Desktop/psf_Pb_NL70_vol.npy")):
        try:
            v = np.load(p)
            V = np.angle(v).astype(float); V -= np.median(V, (1, 2), keepdims=True)
            print(f"  loaded real PSF: {p}")
            return V, 0.0492, 0.999, True
        except Exception:
            continue
    print("  ! no readable PSF (Desktop OS-blocked?); synthetic fallback — drop "
          "psf_Pb_NL70_vol.npy into figs/ and re-run to bake in the real kernel")
    if True:
        nL, N = 46, 60; z = (np.arange(nL) - nL // 2) * 0.999
        y = (np.arange(N) - N // 2) * 0.0492
        Z, Y, X = np.meshgrid(z, y, y, indexing="ij")
        V = np.exp(-0.5 * ((X / 0.35) ** 2 + (Y / 0.35) ** 2 + (Z / 1.1) ** 2))
        return V, 0.0492, 0.999, False


def _fwhm(line, step):
    pk = line.max()
    above = np.where(line >= pk / 2)[0]
    return (above[-1] - above[0]) * step if len(above) > 1 else step


def fig_missing_cone():
    V, dx, dz, real = _load_psf()
    l, r, c = np.unravel_index(np.argmax(V), V.shape)
    inplane = np.abs(V).max(0)
    axial = np.abs(V[:, r, :])
    fwhm_xy = _fwhm(np.abs(V[l, r, :]), dx)
    fwhm_z = _fwhm(np.abs(V[:, r, c]), dz)

    fig = plt.figure(figsize=(6.6, 2.95))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1.0], height_ratios=[1, 1],
                          wspace=0.28, hspace=0.40)
    axd = fig.add_subplot(gs[:, 0]); axi = fig.add_subplot(gs[0, 1]); axa = fig.add_subplot(gs[1, 1])

    # --- reciprocal-space schematic: accessible slab + missing cone ---
    axd.set_xlim(-1.15, 1.15); axd.set_ylim(-1.18, 1.42); axd.set_aspect("equal")
    axd.set_xticks([]); axd.set_yticks([])
    # missing double cone about K_z (grey, hatched)
    for sgn in (1, -1):
        axd.add_patch(Polygon([(0, 0), (-0.5, sgn * 1.0), (0.5, sgn * 1.0)], closed=True,
                              fc="#d9d9d9", ec="#b0b0b0", lw=0.5, hatch="////", alpha=0.9))
    # accessible thin lens near K_z = 0 (blue)
    th = np.linspace(0, 2 * np.pi, 200)
    axd.add_patch(Polygon(np.c_[0.96 * np.cos(th), 0.30 * np.sin(th)], closed=True,
                          fc=ACCENT, alpha=0.28, ec=ACCENT, lw=1.0))
    axd.annotate("", xy=(1.08, 0), xytext=(-1.08, 0), arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.8))
    axd.annotate("", xy=(0, 1.10), xytext=(0, -1.10), arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.8))
    axd.text(1.05, -0.16, r"$K_\parallel$", fontsize=8.5, color=INK)
    axd.text(0.10, 1.13, r"$K_z$", fontsize=8.5, color=INK)
    bb = dict(fc="white", ec="none", alpha=0.75, pad=1)
    axd.text(0.0, 0.03, "transferred", color="#0a4e77", fontsize=7.2, ha="center", va="bottom", bbox=bb)
    axd.text(0.0, 0.62, "missing cone\n(null space)", color=MUTED, fontsize=7.0, ha="center", va="center", bbox=bb)
    axd.text(0.58, 0.40, r"$\delta z\!\approx\!\lambda/\alpha^2$", color="#0a4e77", fontsize=8.0, ha="left")
    axd.text(-1.12, 1.27, "reciprocal space (single orientation)", fontsize=8.0, ha="left", color=INK)

    # --- the real measured PSF ---
    ext_i = [-inplane.shape[1] / 2 * dx, inplane.shape[1] / 2 * dx] * 2
    axi.imshow(inplane, cmap="inferno", extent=ext_i, aspect="equal")
    axi.set_title(f"measured PSF — in-plane  (FWHM {fwhm_xy:.2f} Å)", fontsize=7.2, pad=2)
    axi.set_xticks([]); axi.set_yticks([])
    ext_a = [-axial.shape[1] / 2 * dx, axial.shape[1] / 2 * dx, axial.shape[0] / 2 * dz, -axial.shape[0] / 2 * dz]
    axa.imshow(axial, cmap="inferno", extent=ext_a, aspect="auto")
    axa.set_title(f"axial ($z$–$x$)  (FWHM {fwhm_z:.1f} Å)", fontsize=7.2, pad=2)
    axa.set_xlabel("x (Å)", fontsize=7.0); axa.set_ylabel("z (Å)", fontsize=7.0)
    axa.tick_params(labelsize=6.5)
    if not real:
        axi.text(0.5, 0.5, "SYNTHETIC", transform=axi.transAxes, color="white",
                 ha="center", fontsize=8, alpha=0.6)
    save(fig, "missing_cone")


if __name__ == "__main__":
    fig_acquisition()
    fig_debye_waller()
    fig_ptycho_inverse()
    fig_missing_cone()
    print("done ->", FIGS)
