#!/usr/bin/env python
"""@file make_figures.py
@brief Figures for the hollow-detector fusion study, from the real simulations, reconstructions and spectra.

  1  the headline       -- the specimen, the detector split into its two channels, and the sign:
                           what the model test picks against what the reconstructed atoms support
  2  the EELS axis      -- the CASTEP ladder, the domain spectra, and the inversion to |delta_z|
  3  sign vs hole size  -- the blind likelihood-ratio sign test at 50-95 mrad: evidence per recorded
                           and per incident electron, and the electrons per pattern that implies
  4  depth sections     -- the blind reconstruction and the reference volume through the Ti-O and
                           equatorial-O columns: atomfind picks (95 % depth intervals) against truth
  5  Ti vs its oxygens  -- z(Ti) - mean z(equatorial O) per located Ti, with atomfind's 95 % intervals

Inputs are the Blythe pulls under ~/Desktop/fusion_recons/round2-5 (fusion/runs is gitignored).

    ~/hyperspy-bundle/bin/python make_figures.py [--fig 1 2 3 4 5]
"""
from __future__ import annotations

import argparse
import functools
import glob
import json
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "analysis"))
import analyze_fusion as AF     # noqa: E402
import eels_forward as F        # noqa: E402
import atomfind_sign as AS      # noqa: E402

DESK = "/Users/u2109287/Desktop/fusion_recons"
R5 = os.path.join(DESK, "round5", "fusion", "runs")
S015_BUDGET = os.path.join(DESK, "round4", "fusion", "runs", "fusion_s015", "hollow_budget.json")
NL70 = "/Users/u2109287/Desktop/NL70_new_vol.npy"

COL = {"A": "#c0392b", "B": "#2980b9", "C": "#e67e22", "D": "#16a085"}
SPEC = {"Pb": "#6D4AA6", "TiO": "#16786A", "Oeq": "#B3701A"}
EELS_C, PTY_C, HOLE_C = "#6D4AA6", "#16786A", "#58D6C2"
RECON_C = "#B3701A"
OVER = {82: "#C9A7FF", 22: "#3FE0C0", 8: "#FFB347"}          # overlays on the bone colormap
AF_DIR = os.path.join(DESK, "round5", "atomfind")
AF_PREFIX = os.path.join(AF_DIR, "atomfind_sign_015")
IDEAL_NOISE = 0.1                                             # x the reconstruction background
AF_LABEL = {"recon": "blind reconstruction, 0.15 Å step",
            "ideal": "reference: true atoms × measured kernels (not a reconstruction)"}
AF_SHORT = {"recon": "blind reconstruction", "ideal": "reference (true atoms × kernels)"}
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "white"})


# ---------------------------------------------------------------- inputs
def blind015_path():
    return glob.glob(os.path.join(R5, "recon_hsa0_NL26", "01", "*step02*", "Niter200.mat"))[0]


@functools.lru_cache(maxsize=None)
def af_context():
    """@brief The blind 0.15 A reconstruction, its lattice frame and the atomfind set-up."""
    return AS.context(blind015_path(), S015_BUDGET, os.path.join(HERE, "sample", "toy_truth.npz"),
                      os.path.join(HERE, "sample", "toy_membrane.vasp"))


@functools.lru_cache(maxsize=None)
def af_volume(tag):
    """@brief 'recon' = the blind reconstruction; 'ideal' = the reference volume atomfind ran on."""
    ctx = af_context()
    return ctx["V"] if tag == "recon" else AS.ideal_volume(ctx, IDEAL_NOISE, seed=0)[0]


@functools.lru_cache(maxsize=None)
def af_result(tag):
    """@brief (found atoms, Ti-equatorial-O offsets with atomfind's 95 % intervals), from the cache."""
    ctx = af_context()
    name = "recon" if tag == "recon" else "ideal_noise%g" % IDEAL_NOISE
    path = f"{AF_PREFIX}_{name}_found.npy"
    found = AS.find_or_load(ctx, None if os.path.exists(path) else af_volume(tag), path, reuse=True)
    return found, AS.equatorial_offsets(found, ctx)


def dp_samples():
    return np.load(os.path.join(R5, "fig_dp_step0p3.npz"), allow_pickle=True)


def names_of(truth):
    return [str(x) for x in truth["names"]]


# ---------------------------------------------------------------- figure 1
INK = "#14131A"


def _vector3d(ax, delta, dmax):
    """@brief One domain's polar displacement in 3-D: the vector, the beam axis, and its projection.

    The dashed shadow on the image plane is everything projected ptychography measures (A and B cast
    the same shadow); the vertical drop from it to the tip is |delta_z|, what EELS measures.
    """
    v = 0.95 * np.asarray(delta, float) / dmax
    plane = "#D9D5E0"
    sq = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1], [-1, -1]], float)
    ax.plot(sq[:, 0], sq[:, 1], np.zeros(5), color=plane, lw=0.9)
    ax.plot([-1, 1], [0, 0], [0, 0], color=plane, lw=0.7)
    ax.plot([0, 0], [-1, 1], [0, 0], color=plane, lw=0.7)
    ax.text(1.12, 0, 0, "x", color="0.55", fontsize=8, ha="center", va="center")
    ax.text(0, 1.12, 0, "y", color="0.55", fontsize=8, ha="center", va="center")
    ax.plot([0, 0], [0, 0], [-1.1, 1.1], color="#8F8A99", lw=1.1)
    ax.text(0.12, 0, 1.05, "z", color="0.4", fontsize=8.5, ha="left", va="center")
    ax.plot([0, v[0]], [0, v[1]], [0, 0], color=PTY_C, lw=2.2, ls=(0, (3, 1.6)))   # projection
    ax.plot([v[0], v[0]], [v[1], v[1]], [0, v[2]], color=EELS_C, lw=1.8)          # |delta_z|
    ax.quiver(0, 0, 0, v[0], v[1], v[2], color=INK, linewidth=2.4, arrow_length_ratio=0.22)
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
    ax.set_proj_type("ortho")
    ax.view_init(elev=24, azim=-58)
    ax.set_box_aspect((1, 1, 1), zoom=1.9)
    ax.set_axis_off()


def figure_headline(truth, out, hole=75.0):
    """@brief The whole argument in three panels: the specimen, the split detector, and the fusion."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Wedge
    names = names_of(truth)
    dmax = float(truth["delta_Ti_A"])

    with plt.rc_context({"font.size": 10, "axes.titlesize": 10.5, "axes.labelsize": 10}):
        fig = plt.figure(figsize=(13.2, 4.9), constrained_layout=True)
        gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 0.95, 1.25])

        # (a) the polar displacement in each domain, in 3-D, beam along z
        sub = gs[0].subgridspec(4, 2, height_ratios=[0.09, 1, 1, 0.20], hspace=0.02, wspace=0.02)
        tax = fig.add_subplot(sub[0, :])
        tax.axis("off")
        tax.text(0.0, 0.5, "(a) four domains: same |δ|, different direction", transform=tax.transAxes,
                 ha="left", va="center", fontsize=10.5)
        for k, n in enumerate(names):
            d = np.asarray(truth["deltas"][k], float)
            ax3 = fig.add_subplot(sub[1 + k // 2, k % 2], projection="3d")
            _vector3d(ax3, d, dmax)
            ax3.text2D(0.0, 1.0, f"{n} · {'up' if d[2] > 0 else 'down'}", transform=ax3.transAxes,
                       fontsize=11, weight="bold", color=INK, va="top")
        fax = fig.add_subplot(sub[3, :])
        fax.axis("off")
        fax.text(0.0, 0.95, f"arrow: polar displacement δ, |δ| = {dmax:.3f} Å; z is the beam axis\n"
                 "teal dashed: in-plane projection (ptychography)   purple: |δz| (EELS)",
                 transform=fax.transAxes, ha="left", va="top", fontsize=8.5, color="0.35",
                 linespacing=1.6)

        # (b) the detector: the measured mean pattern, split into the two channels
        ax = fig.add_subplot(gs[1])
        dps = dp_samples()
        th, M = dps["theta_mrad"], dps["mean"]
        alpha = float(dps["convergence_mrad"])
        R = th.shape[0] * float(dps["d_alpha_mrad"]) / 2
        eels = M[th < hole].sum() / M.sum()
        norm = LogNorm(vmin=M.max() * 1e-5, vmax=M.max())
        ax.imshow(np.maximum(M, norm.vmin), cmap="gray", norm=norm, extent=[-R, R, R, -R])
        ax.add_patch(Wedge((0, 0), R, 0, 360, width=R - hole, fc=PTY_C, alpha=0.18, ec="none"))
        ax.add_patch(Circle((0, 0), hole, fc=EELS_C, alpha=0.55, ec="none"))
        ax.add_patch(Circle((0, 0), alpha, fill=False, ec="white", lw=1.0, ls=(0, (4, 3))))
        ax.text(0, 0, f"hole → EELS\n{100 * eels:.0f} % of beam", ha="center", va="center",
                fontsize=10, weight="bold", color="white")
        ax.text(0, -(hole + R) / 2 - 8, f"annulus → ptychography  {100 * (1 - eels):.0f} %",
                ha="center", va="center", fontsize=9.5, weight="bold", color=PTY_C,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.92))
        ax.text(alpha * 0.71 + 6, alpha * 0.71 + 6, "α", color="white", fontsize=10)
        ax.set_xlim(-R, R); ax.set_ylim(R, -R)
        ax.set_xticks([-200, -100, 0, 100, 200]); ax.set_yticks([-200, -100, 0, 100, 200])
        ax.set_xlabel("scattering angle (mrad)")
        ax.set_title(f"(b) one detector, two channels: {hole:.0f} mrad hole", loc="left")

        # (c) the sign: what the model test picks, and what the reconstructed atoms support
        ax = fig.add_subplot(gs[2])
        theta = AF.eels_theta_from_contrast(truth, alpha, hole)
        sweep = {r["hole"]: r for r in sign_sweep()}[hole]
        rec = AS.domain_summary(af_result("recon")[1], truth)
        ys = {n: len(names) - 1 - k for k, n in enumerate(names)}
        for k, n in enumerate(names):
            y = ys[n]
            d = np.asarray(truth["deltas"][k], float)
            mag = float(np.linalg.norm(d[:2]) / np.tan(np.radians(theta[n])))
            s = 1.0 if sweep["dom"][n]["decided"] == "up" else -1.0
            ax.plot([-mag, mag], [y, y], color=EELS_C, lw=1.2, alpha=0.45, zorder=1)
            ax.scatter([-mag, mag], [y, y], s=130, facecolors="white", edgecolors=EELS_C,
                       linewidths=2.0, zorder=3)
            ax.annotate("", xy=(s * mag * 0.93, y + 0.32), xytext=(0, y + 0.32),
                        arrowprops=dict(arrowstyle="-|>", color=PTY_C, lw=2.2, mutation_scale=14))
            ax.plot([d[2], d[2]], [y - 0.46, y + 0.46], color="0.55", lw=1.6, zorder=2)
            if n in rec:
                r = rec[n]
                ax.errorbar([r["mean"]], [y - 0.32], xerr=[r["hw_mean"]], fmt="D", ms=6.5,
                            color=RECON_C, ecolor=RECON_C, elinewidth=2.0, capsize=4, zorder=4)
                print(f"[fig1] {n}: EELS |dz| {mag:.4f}; model test {sweep['dom'][n]['decided']}; "
                      f"reconstructed atoms {r['mean']:+.3f} +- {r['hw_mean']:.3f} A (95 %, n={r['n']})")
        ax.axvline(0, color="0.75", lw=1.0, zorder=0)
        ax.set_yticks([ys[n] for n in names])
        ax.set_yticklabels(names, fontsize=12, weight="bold")
        ax.set_ylim(-0.75, len(names) - 0.3)
        ax.set_xlim(-0.9, 0.75)
        ax.set_xticks([-0.6, -0.3, 0, 0.3, 0.6])
        ax.set_xlabel("δz (Å)")
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        handles = [Line2D([], [], ls="none", marker="o", ms=9, mfc="white", mec=EELS_C, mew=2,
                          label="EELS: |δz|, either sign"),
                   Line2D([], [], color=PTY_C, lw=2.2, marker=">", ms=7,
                          label="model test on simulated patterns"),
                   Line2D([], [], color=RECON_C, lw=2.0, marker="D", ms=6,
                          label="from reconstructed atoms, 95 %"),
                   Line2D([], [], color="0.55", lw=1.6, label="ground truth")]
        ax.legend(handles=handles, frameon=False, fontsize=8.5, ncol=2, loc="upper center",
                  bbox_to_anchor=(0.45, -0.14))
        ax.set_title("(c) the sign: in the data, not yet in the reconstruction", loc="left")

        fig.savefig(out, dpi=170)
    print(f"wrote {out}")


# ---------------------------------------------------------------- figure 2
def figure_eels(truth, out, alpha=100.0, beta=75.0):
    """@brief The CASTEP ladder, the domain spectra, and the inversion to |delta_z|."""
    lad = F.load_ladder()
    s_grid, e, S_par, S_perp = lad
    w = F.aperture_weights(alpha, beta)
    dmax = float(truth["delta_Ti_A"])
    names = names_of(truth)
    m = (e >= -3) & (e <= 30)

    fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.1), constrained_layout=True)

    # (a) the displacement ladder, q perp c -- the channel a 100 mrad probe lands on
    ramp = plt.get_cmap("viridis")(np.linspace(0.05, 0.78, len(s_grid)))   # skip the pale yellow
    for i, sv in enumerate(s_grid):
        ax[0].plot(e[m], S_perp[i][m] / S_perp.max(), color=ramp[i], lw=1.7,
                   label=f"s = {sv:.2f}")
    ax[0].set_xlabel("energy − O-K onset (eV)"); ax[0].set_ylabel("intensity (norm.)")
    ax[0].set_title("(a) CASTEP ladder, q ⊥ c", loc="left")
    ax[0].legend(frameon=False, fontsize=7.5, title="along-chain distortion", title_fontsize=7.5,
                 loc="upper right")
    ax[0].annotate("π* grows with the\npolar displacement", (8.9, 1.0), xytext=(7.2, 0.97),
                   textcoords="data", fontsize=8, color="0.4", ha="right", va="top",
                   arrowprops=dict(arrowstyle="-", color="0.55", lw=0.8))

    # (b) the four domains, and the difference
    spec = {n: F.column_spectrum(truth["deltas"][k], dmax, w, lad) for k, n in enumerate(names)}
    ref = spec[names[0]]
    for n in names[1:]:
        ax[1].plot(e[m], (spec[n] - ref)[m] / ref[m].max() * 100, color=COL[n],
                   lw=2.6 if n == "C" else 1.5,
                   label=f"{n} − A   |δz| = {abs(truth['deltas'][names.index(n)][2]):.3f} Å")
    ax[1].annotate("C − A and D − A coincide", (0.97, 0.10), xycoords="axes fraction",
                   ha="right", fontsize=8, color="0.4")
    ax[1].axhline(0, color="0.4", lw=0.8)
    ax[1].set_xlabel("energy − O-K onset (eV)"); ax[1].set_ylabel("difference (% of edge max)")
    ax[1].set_title("(b) domain spectra, differenced", loc="left")
    ax[1].legend(frameon=False, fontsize=7.5)
    ax[1].annotate("B − A ≡ 0 : sign-blind", (0.04, 0.16), xycoords="axes fraction",
                   fontsize=8.5, color=COL["B"], weight="bold")

    # (c) the inversion: contrast against |delta_z|
    th = np.linspace(0, 90, 46)
    lib = [F.column_spectrum(dmax * np.array([np.sin(np.radians(t)), 0, np.cos(np.radians(t))]),
                             dmax, w, lad) for t in th]
    base = lib[0]
    c = np.array([F.contrast(x, base, e) * 100 for x in lib])
    dz = dmax * np.cos(np.radians(th))
    ax[2].plot(dz, c, color="#6D4AA6", lw=2)
    for k, n in enumerate(names):
        dzk = abs(truth["deltas"][k][2])
        ck = np.interp(-dzk, -dz, c)
        ax[2].scatter([dzk], [ck], s=48, color=COL[n], zorder=3, edgecolors="white", linewidths=1.2)
        ax[2].annotate(n, (dzk, ck), xytext=(7 if n in ("A", "C") else -13, 5),
                       textcoords="offset points", fontsize=9, color=COL[n], weight="bold")
    ax[2].set_xlabel("|δz| (Å)"); ax[2].set_ylabel("O-K contrast vs fully along-beam (%)")
    ax[2].set_title("(c) EELS is monotonic in |δz|", loc="left")
    for lbl, k in (("A, B coincide — sign-blind", 0), ("C, D coincide", 2)):
        dzk = abs(truth["deltas"][k][2])
        ax[2].annotate(lbl, (dzk, np.interp(-dzk, -dz, c)), xytext=(-30, -24) if k == 0 else (-14, -26),
                       textcoords="offset points", fontsize=8, color="0.4", ha="right",
                       arrowprops=dict(arrowstyle="-", color="0.55", lw=0.8))

    fig.suptitle("The EELS axis: O-K near-edge structure carries |δz|, and nothing about its sign",
                 fontsize=10)
    fig.savefig(out, dpi=170)
    print(f"wrote {out}")


# ---------------------------------------------------------------- figure 3
_LLR = re.compile(r"  ([ABCD]): truth\s+(\w+), decided\s+(\w+)\s+(?:OK|WRONG)\s+\|\s+LLR ([-+0-9.e]+) "
                  r"\(([-+0-9.e]+)/pattern\), (\d+)% of individual")


def sign_sweep():
    """@brief The four sign-test logs -> per hole: fraction recorded, and per-domain LLR results.

    `sign_test.log_likelihood_ratio` normalises every pattern, so LLR/pattern is evidence per
    RECORDED electron (~ the KL divergence between the hypotheses). Multiplying by the fraction the
    annulus keeps -- measured on the same patterns -- gives evidence per INCIDENT electron.
    """
    d = dp_samples()
    M, th = d["mean"], d["theta_mrad"]
    rows = []
    for f in glob.glob(os.path.join(R5, "sign_test_12744*.out")):
        s = open(f).read()
        hole = float(re.search(r"hole (\d+) mrad", s).group(1))
        dom = {m.group(1): dict(truth=m.group(2), decided=m.group(3), llr=float(m.group(4)),
                                per=float(m.group(5)), frac=int(m.group(6)))
               for m in _LLR.finditer(s)}
        rows.append(dict(hole=hole, kept=float(M[th >= hole].sum() / M.sum()), dom=dom))
    return sorted(rows, key=lambda r: r["hole"])


def figure_sign(truth, out):
    """@brief The blind sign test against hole size: the evidence barely notices the hole."""
    rows = sign_sweep()
    names = names_of(truth)
    x = np.array([100 * (1 - r["kept"]) for r in rows])
    kept = np.array([r["kept"] for r in rows])
    n_ok = sum(r["dom"][n]["decided"] == r["dom"][n]["truth"] for r in rows for n in names)
    worst = min(r["dom"][n]["frac"] for r in rows for n in names)

    fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.3), constrained_layout=True)
    ax[1].sharey(ax[0])
    for k, n in enumerate(names):
        dz = float(truth["deltas"][k][2])
        per = np.array([abs(r["dom"][n]["per"]) for r in rows])
        inc = per * kept
        style = dict(color=COL[n], marker="o" if dz > 0 else "s", ms=6.5, lw=1.8,
                     mfc=COL[n] if dz > 0 else "white", mew=1.6)
        ax[0].semilogy(x, per, **style, label=f"{n}   δz {dz:+.3f} Å")
        ax[1].semilogy(x, inc, **style)
        ax[2].semilogy(x, np.log(100) / inc, **style)
        print(f"[fig2] {n}: per recorded e- {per.round(6)}, per incident e- {inc.round(6)} "
              f"({100 * (inc[-1] / inc[0] - 1):+.0f} %), e-/pattern for 100:1 "
              f"{(np.log(100) / inc).round(0)}")
        if n in ("A", "C"):
            ax[1].annotate(f"{100 * (inc[-1] / inc[0] - 1):+.0f} %", (x[-1], inc[-1]),
                           xytext=(6, -3 if n == "A" else 6), textcoords="offset points",
                           fontsize=8.5, color=COL[n], weight="bold")
    ticks = [f"{xx:.0f} %\n({r['hole']:.0f})" for xx, r in zip(x, rows)]
    for a_ in ax:
        a_.set_xticks(x)
        a_.set_xticklabels(ticks, fontsize=7.5)
        a_.set_xlim(15, 106)
        a_.set_xlabel("share of beam sent to the spectrometer\n(hole semi-angle, mrad)")
        a_.grid(axis="y", which="major", color="0.92", lw=0.8)
    ax[0].set_ylabel("|log-likelihood ratio| per electron (nats)")
    ax[0].set_title("(a) sign evidence per recorded electron", loc="left")
    ax[0].legend(frameon=False, fontsize=8, loc="upper left", title="filled = up, open = down",
                 title_fontsize=7.5)
    ax[1].set_title("(b) per incident electron: nearly flat", loc="left")
    ax[1].tick_params(labelleft=False)
    ax[2].set_ylabel("incident electrons per pattern")
    ax[2].set_title("(c) dose for 100 : 1 odds from a single pattern", loc="left")
    ax[2].set_yticks([2e3, 5e3, 1e4, 2e4, 4e4])
    ax[2].set_yticklabels(["2×10³", "5×10³", "10⁴", "2×10⁴", "4×10⁴"])
    ax[2].minorticks_off()
    ax[2].text(0.03, 0.45, "expected LLR = N · KL, so N = ln 100 / KL\n(no Poisson draw about it)",
               transform=ax[2].transAxes, fontsize=7.5, color="0.45")
    pats = "100 %" if worst == 100 else f"≥ {worst} %"
    fig.suptitle(f"Blind sign test, hollow data only, no depth reconstruction: {n_ok}/{len(rows) * 4} "
                 f"domain decisions correct at 4 hole sizes; {pats} of individual patterns correct "
                 f"at every hole", fontsize=10)
    fig.savefig(out, dpi=170)
    print(f"wrote {out}")


# ---------------------------------------------------------------- figure 4
def figure_sections(truth, out, half_width_A=0.3, near_A=0.5):
    """@brief Depth sections through the Ti-O and equatorial-O columns: atomfind picks against truth.

    Each section runs along x at y = (j + 1/2) a, the plane holding every Ti-O column and the
    equatorial O either side of it, and crosses a domain wall so an up and a down domain sit side by
    side. Squares are the true atoms, circles atomfind's picks with their calibrated 95 % depth
    intervals.
    """
    from matplotlib.lines import Line2D
    ctx = af_context()
    a = float(truth["a"])
    x0, y0, ox, oy, dx = ctx["frame"]
    T = ctx["T"]
    wall = 0.5 * float(truth["box_A"])
    names = names_of(truth)
    up = {nm: truth["deltas"][k][2] > 0 for k, nm in enumerate(names)}
    rows = ((4.5 * a, "A", "B"), (7.5 * a, "C", "D"))

    fig = plt.figure(figsize=(13.2, 13.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    for col, tag in enumerate(("recon", "ideal")):
        V = af_volume(tag)
        _, E = af_result(tag)
        xyz, spec, hw = E["xyz_all"], E["spec_all"], E["hw_all"]["z"]
        X0 = x0 - ox
        X1 = X0 + V.shape[2] * dx
        for rr, (ys, left, right) in enumerate(rows):
            ax = fig.add_subplot(gs[rr, col])
            r = int(round((ys + oy - y0) / dx))
            h = max(1, int(round(half_width_A / dx)))
            img = V[:, r - h:r + h + 1, :].mean(1)
            lo, hi = np.percentile(img, [2, 99.7])
            ax.imshow(img, cmap="bone", vmin=lo, vmax=hi, extent=[X0, X1, T, 0], aspect="equal")
            mt = np.abs(ctx["pos_t"][:, 1] - ys) < near_A
            for sym in ("Pb", "Ti", "O"):
                m = mt & (ctx["sym_t"] == sym)
                ax.scatter(ctx["pos_t"][m, 0], ctx["pos_t"][m, 2], marker="s", s=56,
                           facecolors="none", edgecolors=OVER[AS.Z_OF[sym]], linewidths=1.3, zorder=3)
            mf = np.abs(xyz[:, 1] - ys) < near_A
            for zc in (82, 22, 8):
                m = mf & (spec == zc)
                if m.any():
                    ax.errorbar(xyz[m, 0], xyz[m, 2], yerr=hw[m], fmt="o", ms=5, mfc="none",
                                mec=OVER[zc], mew=1.3, ecolor=OVER[zc], elinewidth=1.0, capsize=2,
                                zorder=4)
            ax.axvline(wall, color="white", lw=1.0, ls=(0, (5, 3)), alpha=0.8)
            for nm, xc in ((left, 0.5 * (X0 + AS.EDGE_A + wall)), (right, 0.5 * (wall + X1 - AS.EDGE_A))):
                ax.text(xc, 0.8, f"{nm} · {'up' if up[nm] else 'down'}", ha="center", va="top",
                        color="white", fontsize=10, weight="bold",
                        bbox=dict(boxstyle="round,pad=0.25", fc="black", alpha=0.55, ec="none"))
            ax.set_xlim(X0 + AS.EDGE_A, X1 - AS.EDGE_A)
            ax.set_ylim(T, 0)
            if rr == 1:
                ax.set_xlabel("x (Å)")
            if col == 0:
                ax.set_ylabel(f"depth z (Å) · beam ↓   (section at y = {ys:.1f} Å)")
            if rr == 0:
                zr = E["z_rms"]
                ax.set_title(f"({'ab'[col]}) {AF_LABEL[tag]}\ndepth error RMS  Pb {zr.get('Pb', np.nan):.2f} · "
                             f"Ti {zr.get('Ti', np.nan):.2f} · O {zr.get('O', np.nan):.2f} Å", loc="left")
        print(f"[fig4] {tag}: depth RMS {E['z_rms']}")
    handles = [Line2D([], [], ls="none", marker="s", ms=8, mfc="none", mec="0.25", mew=1.3, label="ground truth"),
               Line2D([], [], ls="none", marker="o", ms=7, mfc="none", mec="0.25", mew=1.3,
                      label="atomfind pick, 95 % depth interval"),
               Line2D([], [], ls="none", marker="s", ms=8, mfc=OVER[22], mec="0.3", label="Ti"),
               Line2D([], [], ls="none", marker="s", ms=8, mfc=OVER[8], mec="0.3", label="O"),
               Line2D([], [], ls="none", marker="s", ms=8, mfc=OVER[82], mec="0.3", label="Pb")]
    fig.legend(handles=handles, loc="outside lower center", ncol=5, frameon=False, fontsize=9)
    fig.suptitle("Depth sections through the Ti–O and equatorial-O columns (±0.3 Å slab), each crossing a "
                 "domain wall: where atomfind puts the atoms, against where they are", fontsize=10.5)
    fig.savefig(out, dpi=170)
    print(f"wrote {out}")


# ---------------------------------------------------------------- figure 5
def figure_offsets(truth, out):
    """@brief z(Ti) - mean z(equatorial O) per located Ti, with atomfind's propagated 95 % intervals."""
    from matplotlib.lines import Line2D
    names = names_of(truth)
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.6), constrained_layout=True, sharey=True)
    for ax, tag, letter in zip(axes, ("recon", "ideal"), "ab"):
        _, E = af_result(tag)
        S = AS.domain_summary(E, truth)
        n_clear = 0
        for k, nm in enumerate(names):
            m = np.where(E["domain"] == nm)[0]
            if not len(m):
                continue
            m = m[np.argsort(E["dz"][m])]
            v, h = E["dz"][m], E["hw"][m]
            sgn = S[nm]["truth"]
            xs = k - 0.36 + 0.56 * (np.arange(len(m)) + 0.5) / len(m)
            clear = np.abs(v) > h
            for xi, vi, hi, ok_, right in zip(xs, v, h, clear, np.sign(v) == sgn):
                c = (PTY_C if right else "#A63D3D") if ok_ else "#8F8A99"
                ax.plot([xi, xi], [vi - hi, vi + hi], color=c, lw=0.9, alpha=0.8)
                ax.plot([xi], [vi], "o", ms=3.2, color=c)
            n_clear += int(clear.sum())
            t = float(np.nanmedian(E["truth_dz"][m])) if np.isfinite(E["truth_dz"][m]).any() \
                else float(truth["deltas"][k][2])
            ax.plot([k - 0.42, k + 0.42], [t, t], color=INK, lw=1.3, ls=(0, (4, 2)), zorder=3)
            ax.errorbar([k + 0.32], [S[nm]["mean"]], yerr=[S[nm]["hw_mean"]], fmt="D", ms=7,
                        color=RECON_C, ecolor=RECON_C, elinewidth=2.4, capsize=4, zorder=5)
            ax.text(k, 2.12, f"n = {len(m)}", ha="center", va="center", fontsize=8.5, color="0.3",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9), zorder=6)
            print(f"[fig5] {tag} {nm}: n {len(m)}, clear of zero {int(clear.sum())}, mean "
                  f"{S[nm]['mean']:+.3f} +- {S[nm]['hw_mean']:.3f} (true {t:+.3f})")
        ax.axhline(0, color="0.35", lw=1.0, zorder=0)
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels([f"{nm}\n{'up' if truth['deltas'][k][2] > 0 else 'down'}"
                            for k, nm in enumerate(names)])
        ax.set_xlim(-0.55, len(names) - 0.45)
        zr = E["z_rms"]
        ax.set_title(f"({letter}) {AF_SHORT[tag]}: {n_clear} of {len(E['dz'])} Ti intervals clear zero\n"
                     f"depth RMS  Pb {zr.get('Pb', np.nan):.2f} · Ti {zr.get('Ti', np.nan):.2f} · "
                     f"O {zr.get('O', np.nan):.2f} Å  ·  95 % coverage {100 * E['coverage']['z']:.0f} %",
                     loc="left")
    axes[0].set_ylabel("z(Ti) − mean z(equatorial O)  (Å)")
    axes[0].set_ylim(-2.3, 2.3)
    handles = [Line2D([], [], color="#8F8A99", lw=1.2, marker="o", ms=4, label="Ti, 95 % interval spans zero"),
               Line2D([], [], color=PTY_C, lw=1.2, marker="o", ms=4, label="clears zero, right sign"),
               Line2D([], [], color="#A63D3D", lw=1.2, marker="o", ms=4, label="clears zero, wrong sign"),
               Line2D([], [], color=INK, lw=1.3, ls=(0, (4, 2)), label="true offset"),
               Line2D([], [], color=RECON_C, lw=2.2, marker="D", ms=6, label="domain mean, 95 %")]
    fig.legend(handles=handles, loc="outside lower center", ncol=5, frameon=False, fontsize=9)
    fig.suptitle("Ti against its equatorial oxygens along the beam, every located Ti with ≥ 2 of them, "
                 "atomfind's calibrated 95 % intervals propagated", fontsize=10.5)
    fig.savefig(out, dpi=170)
    print(f"wrote {out}")


FIGURES = {"1": (figure_headline, "fig1_headline.png"),
           "2": (figure_eels, "fig2_eels_axis.png"),
           "3": (figure_sign, "fig3_sign_vs_hole.png"),
           "4": (figure_sections, "fig4_sections_atomfind.png"),
           "5": (figure_offsets, "fig5_ti_oxygen_offsets.png")}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default=os.path.join(HERE, "figs"))
    ap.add_argument("--truth", default=os.path.join(HERE, "sample", "toy_truth.npz"))
    ap.add_argument("--fig", nargs="*", default=list(FIGURES))
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    truth = np.load(args.truth, allow_pickle=True)
    for k in args.fig:
        fn, name = FIGURES[k]
        fn(truth, os.path.join(args.out_dir, name))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
