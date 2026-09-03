#!/usr/bin/env python
"""@file make_eels_fig.py
@brief Figure 5: orientation dichroism of the PbTiO3 apical-oxygen O K edge.

Replaces the cropped slide PNG with vector art in the manuscript's style. Panel (a) is a
schematic of the two limiting orientations; panel (b) is the four measured OptaDOS core-loss
spectra, read straight from the CASTEP/OptaDOS output rather than re-plotted from a slide.

  polar cell    tet_Pz_Oap.{qc,qperp}     -- q || c and q perp c
  zero-P ref    scan_0.00_Oap.{qc,qperp}  -- s = 0, same tetragonal strain (the backbone)

All four share one OptaDOS normalisation, which is what makes max|D|/max S meaningful, so
the axis is a common scale rather than four separately normalised curves.

Run:  python report/make_eels_fig.py            (cwd ptychoshelves-clean/)
Writes report/figs/eels_dichroism.pdf.
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrow

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(os.path.dirname(HERE), "eels", "runs")
FIGDIR = os.path.join(HERE, "figs")

plt.rcParams.update({
    "font.size": 8.5, "font.family": "serif", "mathtext.fontset": "cm",
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "axes.labelsize": 9, "axes.titlesize": 9, "legend.fontsize": 7.2,
    "figure.dpi": 150, "savefig.dpi": 300,
})

PAR, PERP, REF = "#2166ac", "#b2182b", "0.55"     # q||c, q_|_c, zero-P backbone
PB, TI, OX = "#4d4d4d", "#8bbdd9", "#c0392b"
EMIN, EMAX = -1.0, 32.0                            # eV relative to the O K edge


def load(seed):
    """(energy relative to the edge, intensity) for the excited atom's block."""
    col, e, y, grab = 2, [], [], False
    with open(os.path.join(RUNS, f"{seed}_core_edge.dat")) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("#"):
                grab = ":exc" in s                 # select the core-hole atom, never sum
                continue
            p = s.split()
            if grab and len(p) > col:
                e.append(float(p[0])); y.append(float(p[col]))
    e, y = np.asarray(e), np.asarray(y)
    m = (e >= EMIN) & (e <= EMAX)
    return e[m], y[m] * 1e5                        # 1e-5 folded into the axis label


def cell(ax, x0, y0, s, polar_along_beam, label):
    """One projected perovskite cell: Pb corners, Ti centre, O at the edge midpoints."""
    for dx in (0, s):
        for dy in (0, s):
            ax.add_patch(Circle((x0 + dx, y0 + dy), 0.115 * s, fc=PB, ec="none", zorder=3))
    # oxygen at the four edge midpoints; the Ti displacement shortens one bond
    d = 0.13 * s
    ti = (x0 + s / 2, y0 + s / 2 - d) if polar_along_beam else (x0 + s / 2 - d, y0 + s / 2)
    for ox, oy in [(x0 + s / 2, y0), (x0 + s / 2, y0 + s),
                   (x0, y0 + s / 2), (x0 + s, y0 + s / 2)]:
        ax.plot([ti[0], ox], [ti[1], oy], color="0.75", lw=0.7, zorder=1)
        ax.add_patch(Circle((ox, oy), 0.085 * s, fc=OX, ec="none", zorder=3))
    ax.add_patch(Circle(ti, 0.105 * s, fc=TI, ec="none", zorder=4))
    # the polarisation arrow, along the Ti off-centring
    col = PAR if polar_along_beam else PERP
    if polar_along_beam:
        ax.add_patch(FancyArrow(ti[0], ti[1] + 0.34 * s, 0, -0.30 * s, width=0.018 * s,
                                head_width=0.08 * s, head_length=0.09 * s, fc=col, ec="none",
                                length_includes_head=True, zorder=5))
        ax.text(ti[0] + 0.10 * s, ti[1] + 0.16 * s, "$P$", color=col, fontsize=8, zorder=6)
    else:
        ax.add_patch(FancyArrow(ti[0] + 0.34 * s, ti[1], -0.30 * s, 0, width=0.018 * s,
                                head_width=0.08 * s, head_length=0.09 * s, fc=col, ec="none",
                                length_includes_head=True, zorder=5))
        ax.text(ti[0] + 0.02 * s, ti[1] + 0.13 * s, "$P$", color=col, fontsize=8, zorder=6)
    # the electron beam, always vertical
    xb = x0 + 1.40 * s
    ax.annotate("", xy=(xb, y0 - 0.18 * s), xytext=(xb, y0 + 1.18 * s),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=0.9,
                                shrinkA=0, shrinkB=0, mutation_scale=9))
    ax.text(x0 + 0.62 * s, y0 + 1.30 * s, label, ha="center", va="bottom",
            color=col, fontsize=8)


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    e_par, y_par = load("tet_Pz_Oap.qc")
    e_perp, y_perp = load("tet_Pz_Oap.qperp")
    e_rp, y_rp = load("scan_0.00_Oap.qc")
    e_rt, y_rt = load("scan_0.00_Oap.qperp")

    fig = plt.figure(figsize=(6.9, 2.85), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 2.35])

    # ---------------------------------------------------------------- (a) schematic
    ax = fig.add_subplot(gs[0, 0]); ax.set_axis_off()
    ax.set_xlim(-0.30, 2.20); ax.set_ylim(-0.60, 4.45); ax.set_aspect("equal")
    cell(ax, 0.15, 2.45, 1.0, True, r"polar axis $\parallel$ beam")
    cell(ax, 0.15, 0.40, 1.0, False, r"polar axis $\perp$ beam")
    for i, (c, lab) in enumerate([(PB, "Pb"), (TI, "Ti"), (OX, "O")]):
        ax.add_patch(Circle((0.18 + 0.62 * i, -0.42), 0.078, fc=c, ec="none"))
        ax.text(0.30 + 0.62 * i, -0.42, lab, va="center", fontsize=7)
    ax.text(-0.30, 4.42, "(a)", fontsize=8.5, va="top", ha="left")

    # ---------------------------------------------------------------- (b) spectra
    ax = fig.add_subplot(gs[0, 1])
    ax.fill_between(e_par, y_par, y_perp, color="#9ecae1", alpha=0.45, lw=0,
                    label=r"$|\Delta|$, polar cell")
    ax.plot(e_rp, y_rp, color=REF, lw=0.9, ls="--",
            label=r"zero-P, $\mathbf{q}\parallel c$")
    ax.plot(e_rt, y_rt, color=REF, lw=0.9, ls=":",
            label=r"zero-P, $\mathbf{q}\perp c$")
    ax.plot(e_perp, y_perp, color=PERP, lw=1.3, label=r"polar, $\mathbf{q}\perp c$")
    ax.plot(e_par, y_par, color=PAR, lw=1.3, label=r"polar, $\mathbf{q}\parallel c$")

    ax.annotate(r"$\pi^*$", xy=(8.64, y_perp[np.argmin(abs(e_perp - 8.64))]),
                xytext=(6.0, 9.4), color=PERP, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=PERP, lw=0.6))
    ax.annotate(r"$\sigma^*$", xy=(18.67, y_par[np.argmin(abs(e_par - 18.67))]),
                xytext=(21.0, 9.4), color=PAR, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=PAR, lw=0.6))

    ax.set_xlim(EMIN, EMAX); ax.set_ylim(0, 10.6)
    ax.set_xlabel("energy loss relative to the O K edge (eV)")
    ax.set_ylabel(r"core-loss intensity ($\times10^{-5}$)")
    ax.legend(loc="upper right", framealpha=1.0, edgecolor="0.8",
              handlelength=1.5, handletextpad=0.5, borderpad=0.35, labelspacing=0.32)
    ax.text(-0.115, 1.115, "(b)", transform=ax.transAxes, fontsize=8.5, va="top")

    out = os.path.join(FIGDIR, "eels_dichroism.pdf")
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.replace(".pdf", ".png"), bbox_inches="tight")
    plt.close(fig)

    num = np.max(np.abs(y_par - y_perp)) / max(y_par.max(), y_perp.max())
    print(f"wrote {out}\n  max|D|/max S over the plotted window = {100*num:.1f}% "
          f"(RESULTS.md: 78.2% over the full edge)")


if __name__ == "__main__":
    main()
