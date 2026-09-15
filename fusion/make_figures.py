#!/usr/bin/env python
"""@file make_figures.py
@brief Figures for the hollow-detector fusion study, from the real simulations, reconstructions and spectra.

  1  the headline       -- the specimen, the detector split into its two channels, and the fusion
                           of the EELS |delta_z| with the sign the hollow data decide
  2  the EELS axis      -- the CASTEP ladder, the domain spectra, and the inversion to |delta_z|
  3  sign vs hole size  -- the blind likelihood-ratio sign test at 50-95 mrad: evidence per recorded
                           and per incident electron, and the electrons per pattern that implies
  4  the phase volume   -- projected phase at the true column positions, depth sections through Pb
                           and Ti-O columns against the true atom depths, and the k_z spectra
  5  atoms in depth     -- fitted atom depths and the per-cell sign readout, blind against known start

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

DESK = "/Users/u2109287/Desktop/fusion_recons"
R5 = os.path.join(DESK, "round5", "fusion", "runs")
S015_BUDGET = os.path.join(DESK, "round4", "fusion", "runs", "fusion_s015", "hollow_budget.json")
NL70 = "/Users/u2109287/Desktop/NL70_new_vol.npy"

COL = {"A": "#c0392b", "B": "#2980b9", "C": "#e67e22", "D": "#16a085"}
SPEC = {"Pb": "#6D4AA6", "TiO": "#16786A", "Oeq": "#B3701A"}
EELS_C, PTY_C, HOLE_C = "#6D4AA6", "#16786A", "#58D6C2"
RUN_C = {"blind03": "#9C8F7A", "blind015": "#B3701A", "known": "#16786A"}
LABEL = {"blind03": "blind · 0.3 Å step", "blind015": "blind · 0.15 Å step",
         "known": "known start · diagnostic"}
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "white"})


# ---------------------------------------------------------------- inputs
def pixel_size(path):
    import h5py
    with h5py.File(path, "r") as f:
        return float(np.asarray(f["outputs"]["pixel_size"]).ravel()[0]) * 1e10


@functools.lru_cache(maxsize=None)
def recon(name):
    """@brief (volume, budget) for a named reconstruction; the object pixel always from the file."""
    if name == "blind015":
        p = glob.glob(os.path.join(R5, "recon_hsa0_NL26", "01", "*step02*", "Niter200.mat"))[0]
        b = json.load(open(S015_BUDGET))
    else:
        p = (os.path.join(DESK, "hsa0_Niter200.mat") if name == "blind03" else
             glob.glob(os.path.join(DESK, "round3", "**", "recon_hsa0_NL16_known", "01", "*",
                                    "Niter200.mat"), recursive=True)[0])
        b = json.load(open(os.path.join(DESK, "round2", "hollow_budget_production.json")))
    b["dx_object_A"] = pixel_size(p)
    return AF.load_volume(p), b


def dp_samples():
    return np.load(os.path.join(R5, "fig_dp_step0p3.npz"), allow_pickle=True)


def names_of(truth):
    return [str(x) for x in truth["names"]]


def kz_spectrum(P, dz):
    R = P.reshape(P.shape[0], -1)
    R = R - R.mean(0)
    t = np.polynomial.polynomial.polyvander(np.arange(P.shape[0]) * dz, 3)
    R = R - t @ np.linalg.lstsq(t, R, rcond=None)[0]
    S = (np.abs(np.fft.rfft(R, axis=0)) ** 2).mean(1)
    f = np.fft.rfftfreq(P.shape[0], d=dz)
    return f, S / np.median(S[f > 0.12])


def eels_dz_prior(truth):
    """@brief |delta_z| per domain as the EELS channel reports it (theta inverted, delta_xy known)."""
    theta = AF.eels_theta_from_contrast(truth, 100.0, 75.0)
    return {n: float(np.linalg.norm(truth["deltas"][k][:2]) / np.tan(np.radians(theta[n])))
            for k, n in enumerate(names_of(truth))}


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
            ax3.text2D(0.0, 1.0, n, transform=ax3.transAxes, fontsize=13, weight="bold",
                       color=INK, va="top")
            ax3.text2D(0.13, 0.985, "up" if d[2] > 0 else "down", transform=ax3.transAxes,
                       fontsize=9.5, color="0.35", va="top")
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

        # (c) the fusion, per domain
        ax = fig.add_subplot(gs[2])
        theta = AF.eels_theta_from_contrast(truth, alpha, hole)
        sweep = {r["hole"]: r for r in sign_sweep()}[hole]
        ys = {n: len(names) - 1 - k for k, n in enumerate(names)}
        for k, n in enumerate(names):
            y = ys[n]
            d = np.asarray(truth["deltas"][k], float)
            mag = float(np.linalg.norm(d[:2]) / np.tan(np.radians(theta[n])))
            res = sweep["dom"][n]
            s = 1.0 if res["decided"] == "up" else -1.0
            ax.plot([-mag, mag], [y, y], color=EELS_C, lw=1.2, alpha=0.45, zorder=1)
            ax.annotate("", xy=(s * mag * 0.93, y + 0.30), xytext=(0, y + 0.30),
                        arrowprops=dict(arrowstyle="-|>", color=PTY_C, lw=2.2, mutation_scale=14))
            ax.scatter([-mag, mag], [y, y], s=130, facecolors="white", edgecolors=EELS_C,
                       linewidths=2.0, zorder=3)
            ax.plot([d[2], d[2]], [y - 0.22, y + 0.18], color="0.55", lw=1.6, zorder=2)
            ax.scatter([s * mag], [y], s=42, color=INK, zorder=4)
            ax.text(0.47, y, f"{res['frac']} %", ha="right", va="center", fontsize=9.5,
                    color=PTY_C, weight="bold")
            print(f"[fig1] {n}: EELS |dz| {mag:.4f} (truth {abs(d[2]):.4f}), sign {res['decided']} "
                  f"({res['frac']} % of patterns), fused {s * mag:+.4f} vs truth {d[2]:+.4f}")
        ax.axvline(0, color="0.75", lw=1.0, zorder=0)
        ax.set_yticks([ys[n] for n in names])
        ax.set_yticklabels(names, fontsize=12, weight="bold")
        ax.set_ylim(-0.7, len(names) - 0.15)
        ax.set_xlim(-0.42, 0.50)
        ax.set_xticks([-0.3, -0.15, 0, 0.15, 0.3])
        ax.set_xlabel("δz (Å)")
        ax.text(0.47, len(names) - 0.38, "patterns correct", ha="right", va="center",
                fontsize=8.5, color=PTY_C)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        handles = [Line2D([], [], ls="none", marker="o", ms=9, mfc="white", mec=EELS_C, mew=2,
                          label="EELS: |δz|, either sign"),
                   Line2D([], [], color=PTY_C, lw=2.2, marker=">", ms=7, label="hollow data: the sign"),
                   Line2D([], [], ls="none", marker="o", ms=6, color=INK, label="fused δz"),
                   Line2D([], [], color="0.55", lw=1.6, label="ground truth")]
        ax.legend(handles=handles, frameon=False, fontsize=8.5, ncol=2, loc="upper center",
                  bbox_to_anchor=(0.45, -0.14))
        ax.set_title("(c) fusion: EELS leaves ±|δz|, the recorded annulus picks one", loc="left")

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
def _cell(truth, X, Y):
    a, n = float(truth["a"]), int(truth["n_lat"])
    return int(np.floor(X / a + 1e-6)) % n, int(np.floor(Y / a + 1e-6)) % n


def column_line(truth, frame, shape, name, n_max=5, margin_A=1.2):
    """@brief The longest [110] run of alternating Pb and Ti-O columns inside one domain and the ROI."""
    a, n = float(truth["a"]), int(truth["n_lat"])
    x0, y0, ox, oy, dx = frame
    _, Ny, Nx = shape
    m = margin_A / dx

    def ok(X, Y):
        px, py = (X + ox - x0) / dx, (Y + oy - y0) / dx
        return (m < px < Nx - 1 - m and m < py < Ny - 1 - m
                and truth["domain_grid"][_cell(truth, X, Y)] == name)

    best = []
    for i in range(n):
        for j in range(n):
            for start in ((i * a, j * a), ((i + .5) * a, (j + .5) * a)):
                run = []
                while ok(start[0] + len(run) * a / 2, start[1] + len(run) * a / 2):
                    run.append((start[0] + len(run) * a / 2, start[1] + len(run) * a / 2))
                if len(run) > len(best):
                    best = run
    k = max(0, (len(best) - n_max) // 2)
    return best[k:k + n_max]


def section(V, frame, run, half_width_A=0.3, pad_A=1.0):
    """@brief Depth section (layers x distance) along a column run, averaged across a thin slab."""
    from scipy.ndimage import map_coordinates
    x0, y0, ox, oy, dx = frame
    (Xa, Ya), (Xb, Yb) = run[0], run[-1]
    L = np.hypot(Xb - Xa, Yb - Ya)
    u = np.array([Xb - Xa, Yb - Ya]) / L
    s = np.arange(-pad_A, L + pad_A + 1e-9, dx)
    S, W = np.meshgrid(s, np.linspace(-half_width_A, half_width_A, 7), indexing="ij")
    px = (Xa + S * u[0] - W * u[1] + ox - x0) / dx
    py = (Ya + S * u[1] + W * u[0] + oy - y0) / dx
    sec = np.stack([map_coordinates(V[l], [py.ravel(), px.ravel()], order=1)
                    .reshape(S.shape).mean(1) for l in range(V.shape[0])])
    return s, sec


def figure_volume(truth, out):
    """@brief The reconstructed phase volume against the true atoms, and its depth spectrum."""
    a_lat, c_lat = float(truth["a"]), float(truth["c"])
    fig = plt.figure(figsize=(13.2, 9.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.55], width_ratios=[1, 1.75])

    # (a) projected phase of the blind full-resolution object, with the three column types
    from scipy.ndimage import gaussian_filter
    V, b = recon("blind015")
    x0, y0, ox, oy, dx = AF.lattice_frame(V, b, truth)
    ax = fig.add_subplot(gs[0, 0])
    P = gaussian_filter(V.sum(0), 0.08 / dx)
    n = int(round(3.2 * a_lat / dx))
    o = P.shape[0] // 2
    sub = P[o:o + n, o:o + n]
    lo, hi = np.percentile(sub, [1, 99.7])
    ax.imshow(sub, cmap="bone", vmin=lo, vmax=hi, extent=[0, n * dx, n * dx, 0])
    # markers at the TRUE column positions: registered site + the species' in-plane polar shift,
    # referred to the Pb sublattice the registration anchors on
    cols = {k: AF.inplane_offsets(V, b, truth, k) for k in AF.KINDS}
    pb_mean = cols["Pb"][2].mean(0)
    for lab, key, mk in (("Pb", "Pb", "o"), ("Ti + O$_{ap}$", "TiO", "s"), ("O$_{eq}$", "Oeq", "^")):
        site, _, pred, _ = cols[key]
        X = ((site[:, 0] + ox - x0) / dx - o + 0.5) * dx + pred[:, 0] - pb_mean[0]
        Y = ((site[:, 1] + oy - y0) / dx - o + 0.5) * dx + pred[:, 1] - pb_mean[1]
        k = (X > .3) & (X < n * dx - .3) & (Y > .3) & (Y < n * dx - .3)
        ax.scatter(X[k], Y[k], s=30, marker=mk, facecolors="none", edgecolors=SPEC[key],
                   linewidths=1.4, label=lab)
    meas, pred = cols["Pb"][1], cols["Pb"][2]
    meas, pred = meas - meas.mean(0), pred - pred.mean(0)
    before = np.sqrt((meas ** 2).sum(1).mean())
    after = np.sqrt(((meas - pred) ** 2).sum(1).mean())
    print(f"[fig3] Pb in-plane: {before:.3f} A RMS about the undistorted lattice, {after:.3f} A about "
          f"the true polar shift ({len(meas)} columns)")
    ax.text(0.02, 0.02, f"Pb peaks: {after:.2f} Å RMS from their true in-plane positions\n"
            f"({before:.2f} Å from the undistorted lattice)", transform=ax.transAxes, fontsize=7.5,
            color="white", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", fc="black", alpha=0.55, ec="none"))
    ax.set_xlim(0, n * dx); ax.set_ylim(n * dx, 0)
    ax.legend(fontsize=7.5, loc="upper right", labelcolor="white", facecolor="black",
              framealpha=0.55, edgecolor="none")
    ax.set_title("(a) projected phase · blind 0.15 Å, true column positions", loc="left")
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)")

    # (b) the depth spectrum
    ax = fig.add_subplot(gs[0, 1])
    for name, dz in (("blind03", None), ("blind015", None), ("known", None), ("NL70", 0.999)):
        if name == "NL70":
            W = np.angle(np.load(NL70)).astype(float)
            lbl, col = "70 Å labyrinth reference (NL70)", "#6D4AA6"
        else:
            W, bb = recon(name)
            dz = float(bb["beam_thickness_A"]) / W.shape[0]
            lbl, col = f"{LABEL[name]}  (NL{W.shape[0]})", RUN_C[name]
        W = W - np.median(W, axis=(1, 2), keepdims=True)
        f, S = kz_spectrum(W, dz)
        m = f > 0.10
        ax.semilogy(f[m], S[m], color=col, lw=2.0 if name == "blind015" else 1.5, label=lbl)
        ax.axvline(1 / (2 * dz), color=col, lw=0.9, ls=":", alpha=0.8)
    ax.axvline(1 / c_lat, color="0.35", lw=1.2, ls="--")
    ax.annotate(f"PbTiO₃ c = {c_lat:.2f} Å", (1 / c_lat, ax.get_ylim()[1]), xytext=(4, -12),
                textcoords="offset points", fontsize=8, color="0.35")
    ax.set_xlabel("depth spatial frequency $k_z$ (Å$^{-1}$)")
    ax.set_ylabel("power / median band power")
    ax.set_title("(b) the depth axis: 0.15 Å step moves the peak off the layer Nyquist and onto "
                 "the lattice", loc="left")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.text(0.01, 0.02, "dotted = each run's layer-grid Nyquist", transform=ax.transAxes,
            fontsize=7.5, color="0.45")

    # (c) depth sections along [110] through alternating Pb and Ti-O columns, truth overlaid
    z_off = 0.5 * (float(b["beam_thickness_A"]) - float(truth["slab_thickness_A"]))
    groups = gs[1, :].subgridspec(1, 3)
    first = None
    for g, name in enumerate(("blind03", "blind015", "known")):
        V, bb = recon(name)
        frame = AF.lattice_frame(V, bb, truth)
        T = float(bb["beam_thickness_A"])
        pair = groups[g].subgridspec(1, 2, wspace=0.04)
        for p, dom in enumerate(("A", "B")):
            run = column_line(truth, frame, V.shape, dom)
            s, sec = section(V, frame, run)
            ax = fig.add_subplot(pair[p], sharey=first) if first else fig.add_subplot(pair[p])
            first = first or ax
            lo, hi = np.percentile(sec, [2, 99.7])
            ax.imshow(sec, cmap="bone", vmin=lo, vmax=hi, extent=[s[0], s[-1], T, 0],
                      aspect="auto")
            Xa, Ya = run[0]
            for X, Y in run:
                sk = np.hypot(X - Xa, Y - Ya)
                d = truth["delta_grid"][_cell(truth, X, Y)]
                is_pb = abs(X / a_lat - round(X / a_lat)) < 0.25
                zs, wt = AF.column_sites(d, truth, "Pb" if is_pb else "TiO")
                zs = zs + z_off
                heavy = np.isclose(wt, wt.max())
                if is_pb:
                    ax.scatter(np.full(heavy.sum(), sk), zs[heavy], s=34, marker="o",
                               facecolors="none", edgecolors=SPEC["Pb"], linewidths=1.3)
                else:
                    ax.scatter(np.full(heavy.sum(), sk), zs[heavy], s=30, marker="s",
                               facecolors="none", edgecolors=SPEC["TiO"], linewidths=1.3)
                    ax.scatter(np.full((~heavy).sum(), sk), zs[~heavy], s=16, marker="^",
                               facecolors="none", edgecolors=SPEC["Oeq"], linewidths=1.1)
            head = f"{LABEL[name]}\n" if p == 0 else "\n"
            ax.set_title(f"{head}domain {dom} · δz {'> 0' if dom == 'A' else '< 0'}", loc="left",
                         fontsize=8.5, color=RUN_C[name] if p == 0 else "0.2")
            ax.set_xlabel("along [110] (Å)")
            if g == 0 and p == 0:
                ax.set_ylabel("depth z (Å) · beam ↓")
            else:
                ax.tick_params(labelleft=False)
    fig.supxlabel("(c) depth sections through Pb (○) and Ti–O columns: true Ti □ and apical O △ "
                  "overlaid. A and B differ only in the stacking order of Ti against Pb.",
                  fontsize=8.5, x=0.01, ha="left")
    fig.savefig(out, dpi=170)
    print(f"wrote {out}")


# ---------------------------------------------------------------- figure 5
def fit_z(profile, z, z0, half):
    """@brief Sub-layer z centroid of the peak nearest z0, by parabolic fit on the maximum."""
    m = np.abs(z - z0) <= half
    if not m.any():
        return np.nan
    idx = np.where(m)[0]
    k = idx[np.argmax(profile[idx])]
    if k == 0 or k == len(z) - 1:
        return float(z[k])
    y0, y1, y2 = profile[k - 1], profile[k], profile[k + 1]
    den = y0 - 2 * y1 + y2
    d = 0.5 * (y0 - y2) / den if abs(den) > 1e-12 else 0.0
    return float(z[k] + np.clip(d, -1, 1) * (z[1] - z[0]))


def atom_fits(name, truth):
    """@brief Fitted against true depth for every Pb and Ti atom -> (kind, z_true, z_fit)."""
    V, b = recon(name)
    a_lat, c_lat, n = float(truth["a"]), float(truth["c"]), int(truth["n_lat"])
    z_off = 0.5 * (float(b["beam_thickness_A"]) - float(truth["slab_thickness_A"]))
    cells, xy, z = AF.column_maps(V, b, truth)
    gi = np.floor(xy[:, 0] / a_lat).astype(int) % n
    gj = np.floor(xy[:, 1] / a_lat).astype(int) % n
    rows = []
    for i in range(len(cells)):
        d = truth["delta_grid"][gi[i], gj[i]]
        for kind in ("Pb", "TiO"):
            zs, wt = AF.column_sites(d, truth, kind)
            for zt in zs[np.isclose(wt, wt.max())] + z_off:       # Pb, or the Ti of the Ti-O column
                if z.min() + 1.2 <= zt <= z.max() - 1.2:
                    zf = fit_z(cells[i][kind], z, zt, 0.45 * c_lat)
                    if np.isfinite(zf):
                        rows.append((kind, zt, zf))
    return (np.array([r[0] for r in rows]), np.array([r[1] for r in rows]),
            np.array([r[2] for r in rows]))


def figure_atoms(truth, out):
    """@brief Atoms placed in depth, and the per-cell sign, blind against a known start."""
    c_lat = float(truth["c"])
    names = names_of(truth)
    dzp = eels_dz_prior(truth)
    fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.4), constrained_layout=True)
    guess = 0.45 * c_lat / np.sqrt(3)

    for name in ("blind015", "known"):
        kind, zt, zf = atom_fits(name, truth)
        res = zf - zt
        for kd, mk in (("Pb", "o"), ("TiO", "s")):
            m = kind == kd
            ax[0].scatter(zt[m], zf[m], s=15, marker=mk, alpha=0.6, color=RUN_C[name],
                          edgecolors="none",
                          label=f"{LABEL[name]} · {'Pb' if kd == 'Pb' else 'Ti'}")
        rms = {kd: float(np.sqrt(np.mean(res[kind == kd] ** 2))) for kd in ("Pb", "TiO")}
        ax[1].hist(res, bins=np.linspace(-2, 2, 33), alpha=0.6, color=RUN_C[name],
                   label=f"{LABEL[name]}\nPb {rms['Pb']:.2f} Å · Ti {rms['TiO']:.2f} Å RMS")
        print(f"[fig4] {name}: {len(res)} atoms, RMS Pb {rms['Pb']:.3f} Ti {rms['TiO']:.3f} "
              f"all {np.sqrt(np.mean(res ** 2)):.3f} A (uniform guess {guess:.2f})")
    V0, _ = recon("known")
    z = (np.arange(V0.shape[0]) + 0.5) * float(recon("known")[1]["beam_thickness_A"]) / V0.shape[0]
    ax[0].plot([z.min(), z.max()], [z.min(), z.max()], color="0.4", lw=1, ls="--", zorder=0)
    ax[0].set_xlabel("true atom depth (Å)"); ax[0].set_ylabel("fitted depth (Å)")
    ax[0].set_title("(a) atom depths: fitted against truth", loc="left")
    ax[0].legend(frameon=False, fontsize=7.5, loc="upper left", markerscale=1.6)
    ax[1].axvline(0, color="0.4", lw=1, ls="--")
    ax[1].set_xlabel("fitted − true depth (Å)"); ax[1].set_ylabel("atoms")
    ax[1].set_title("(b) depth residual", loc="left")
    ax[1].legend(frameon=False, fontsize=7.5, loc="upper left")
    ax[1].text(0.98, 0.97, f"fit window ±{0.45 * c_lat:.2f} Å\nuniform guess: RMS {guess:.2f} Å",
               transform=ax[1].transAxes, ha="right", va="top", fontsize=7.5, color="0.4")

    # (c) the per-cell matched-filter readout, each run scaled by its own RMS (only the sign is used)
    for r, name in enumerate(("blind015", "known")):
        V, b = recon(name)
        _, obs, dom, _ = AF.read_sign(V, b, truth, dz_prior=dzp)
        obs = obs / np.sqrt(np.mean(obs ** 2))
        right = []
        for k, dn in enumerate(names):
            m = dom == dn
            if not m.any():
                continue
            want = np.sign(truth["deltas"][k][2])
            right.append(np.sign(obs[m]) == want)
            xk = k + (-0.19 if r == 0 else 0.19)
            rng = np.random.default_rng(10 * k + r)
            ax[2].scatter(xk + rng.normal(0, .045, m.sum()), obs[m], s=22, alpha=0.8,
                          color=RUN_C[name], edgecolors="none")
            ax[2].plot([xk - .13, xk + .13], [np.median(obs[m])] * 2, color=RUN_C[name], lw=2.4)
        acc = 100 * np.mean(np.concatenate(right))
        ax[2].scatter([], [], color=RUN_C[name], s=22, label=f"{LABEL[name]}: {acc:.0f} % of cells")
        print(f"[fig4] {name}: sign readout {acc:.0f} % of cells")
    ax[2].axhline(0, color="0.35", lw=1.1)
    ax[2].set_xticks(range(4))
    ax[2].set_xticklabels([f"{n}\n{'δz > 0' if truth['deltas'][i][2] > 0 else 'δz < 0'}"
                           for i, n in enumerate(names)])
    ax[2].set_ylabel("matched-filter discriminant / RMS")
    ax[2].set_title("(c) per-cell sign from the depth profile", loc="left")
    ax[2].legend(frameon=False, fontsize=7.5, loc="lower left")
    ax[2].text(0.98, 0.97, "correct: > 0 for A, C; < 0 for B, D", transform=ax[2].transAxes,
               ha="right", va="top", fontsize=7.5, color="0.4")
    fig.suptitle("Placing atoms in depth: the known-start object does it and reads every cell's "
                 "sign; the blind object, even at 0.15 Å step, does neither yet", fontsize=10)
    fig.savefig(out, dpi=170)
    print(f"wrote {out}")


FIGURES = {"1": (figure_headline, "fig1_headline.png"),
           "2": (figure_eels, "fig2_eels_axis.png"),
           "3": (figure_sign, "fig3_sign_vs_hole.png"),
           "4": (figure_volume, "fig4_phase_volume.png"),
           "5": (figure_atoms, "fig5_atoms_depth.png")}


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
