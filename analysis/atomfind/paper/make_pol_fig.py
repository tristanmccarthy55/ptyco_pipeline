#!/usr/bin/env python
"""@file make_pol_fig.py
@brief Figure for the polarisation result: Ti-O6 off-centring recovered from the located atoms.

  fig_polarisation.pdf -- (a) in-plane off-centring map (located atoms vs ground truth);
                          (b) recovered vs true off-centring, in-plane and along the beam;
                          (d) signal-to-noise per component (true spread / propagated sigma,
                              log axis, decision line at 1) -- the blind statement that the
                              along-beam component is not measured.

Run:  ~/hyperspy-bundle/bin/python atomfind/paper/make_pol_fig.py   (cwd analysis/)
Needs <out_dir>/polarisation.npz from atomfind/polarisation.py.
"""
from __future__ import annotations
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
from atomfind import config

FIGDIR = os.path.join(HERE, "figs")
os.makedirs(FIGDIR, exist_ok=True)
plt.rcParams.update({
    "font.size": 8.5, "font.family": "serif", "mathtext.fontset": "cm",
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "axes.labelsize": 9, "axes.titlesize": 9, "legend.fontsize": 7.5,
    "figure.dpi": 150, "savefig.dpi": 300,
})
REC, GT = "#2166ac", "#b0b4b8"


def main():
    cfg = config.preset("NL70_coherent")
    d = np.load(os.path.join(cfg.out_dir, "polarisation.npz"))
    ti, A, B, sig, tail = d["ti"], d["delta"], d["delta_gt"], d["sigma"], d["tail"]

    fig = plt.figure(figsize=(7.4, 2.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 4, width_ratios=[1.3, 1, 1, 1.15])

    # (a) in-plane off-centring map, averaged per atomic column so the texture is legible
    ax = fig.add_subplot(gs[0, 0])
    key = np.round(ti[:, :2] / 2.0).astype(int)
    _, inv = np.unique(key, axis=0, return_inverse=True)
    P = np.array([[ti[inv == u, 0].mean(), ti[inv == u, 1].mean(),
                   A[inv == u, 0].mean(), A[inv == u, 1].mean(),
                   B[inv == u, 0].mean(), B[inv == u, 1].mean()]
                  for u in range(inv.max() + 1)])
    ax.quiver(P[:, 0], P[:, 1], P[:, 4], P[:, 5], color=GT, angles="xy",
              scale_units="xy", scale=0.10, width=0.014, label="true")
    ax.quiver(P[:, 0], P[:, 1], P[:, 2], P[:, 3], color=REC, angles="xy",
              scale_units="xy", scale=0.10, width=0.006, label="located atoms")
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)")
    ax.set_title("(a) in-plane off-centring", fontsize=8.5)
    ax.set_aspect("equal")
    ax.legend(loc="lower left", framealpha=0.95, handlelength=1.0,
              handletextpad=0.4, borderpad=0.25, fontsize=6.8)

    # (b,c) recovered vs true, in-plane and along the beam
    for n, (k, lab, ttl) in enumerate([(0, r"$\delta_x$", "(b) in-plane"),
                                       (2, r"$\delta_z$", "(c) along beam")]):
        ax = fig.add_subplot(gs[0, 1 + n])
        lo, hi = -0.55, 0.55
        ax.plot([lo, hi], [lo, hi], color="0.6", lw=0.8, zorder=0)
        ax.errorbar(B[:, k], A[:, k], yerr=1.96 * sig[:, k], fmt="none",
                    ecolor="#c8d4e3", elinewidth=0.7, zorder=1)
        ax.scatter(B[~tail, k], A[~tail, k], s=5, color=REC, lw=0, zorder=2)
        ax.scatter(B[tail, k], A[tail, k], s=11, facecolor="none", edgecolor="#d6604d",
                   lw=0.7, zorder=3, label="cage error")
        r = np.corrcoef(A[:, k], B[:, k])[0, 1]
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
        ax.set_xlabel(f"true {lab} (Å)", labelpad=1)
        ax.set_ylabel(f"recovered {lab} (Å)", labelpad=1)
        ax.set_title(f"{ttl},  $r={r:.2f}$", fontsize=8.5)
        if n == 1:
            ax.legend(loc="lower right", framealpha=0.95, markerscale=1.1,
                      handletextpad=0.3, borderpad=0.3, fontsize=7)

    # (d) is each component actually measured?  Signal-to-noise = the true spread of that
    # component divided by the propagated sigma on it, on a log axis with the decision line
    # at 1, so the reader is asked for no arithmetic: a bar below the line is a component
    # whose error bar exceeds the entire signal it is meant to resolve.
    ax = fig.add_subplot(gs[0, 3])
    x = np.arange(3)
    snr = np.array([B[:, k].std() / np.median(sig[:, k]) for k in range(3)])
    ax.bar(x, snr, 0.55, color=[REC if v >= 1 else "#d6604d" for v in snr], zorder=2)
    ax.axhline(1.0, color="#d6604d", lw=0.9, ls="--", zorder=1)
    ax.set_yscale("log")
    ax.set_ylim(0.15, 90); ax.set_xlim(-0.6, 2.6)
    ax.set_xticks(x); ax.set_xticklabels([r"$\delta_x$", r"$\delta_y$", r"$\delta_z$"])
    ax.set_ylabel("signal-to-noise", labelpad=2)
    ax.set_title("(d) is it measured?", fontsize=8.5)
    for xi, v in zip(x, snr):
        if v >= 1:                       # value sits above the bar
            ax.text(xi, v * 1.25, f"{v:.0f}", fontsize=7, ha="center",
                    va="bottom", color=REC)
        else:                            # short bar: value inside, verdict above
            ax.text(xi, v * 0.78, f"{v:.2f}", fontsize=7, ha="center", va="top",
                    color="white")
            ax.text(xi, v * 1.22, "not measured", fontsize=6.5, ha="center",
                    va="bottom", color="#d6604d")

    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"fig_polarisation.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("wrote fig_polarisation.pdf/.png")


if __name__ == "__main__":
    main()
