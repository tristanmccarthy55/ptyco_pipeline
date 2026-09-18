#!/usr/bin/env python
"""@file fig1_probe.py
@brief Figure 1 -- what opening the aperture costs the probe.

ONE CLAIM: past ~70 mrad the Ronchigram is no longer flat and the probe grows, so conventional
imaging is finished in exactly the regime the rest of the paper works in.

(a-c) simulated Ronchigrams (probe through a thin amorphous film) at 50 / 70 / 100 mrad, each
      scaled to its own aperture so the shadow structure is comparable;
(d)   Ronchigram flatness -- peak-to-valley of the aberration phase with its defocus part removed,
      against the lambda/4 criterion -- and the 90%-enclosed probe diameter, both vs alpha.

Probe parameters and sizes come from campaign/round_sweep.tsv (the planner's own output, so the
figure cannot drift from what was simulated); the flatness metric comes from plan_probe.py.

    ~/hyperspy-bundle/bin/python analysis/paper/fig1_probe.py
"""
from __future__ import annotations
import argparse, importlib.util, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pubstyle as ps                                                    # noqa: E402


def _mod(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alphas", type=int, nargs=3, default=[50, 70, 100],
                    help="the three alphas shown as Ronchigrams")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    ps.use()
    import matplotlib.pyplot as plt
    R = _mod(os.path.join(os.path.dirname(HERE), "make_ronchigram_fig.py"), "ronchi")
    plan = {al: (c3, c1) for al, c3, c1 in R.PTS}

    # sizes + flatness straight from the planner's table
    tsv = os.path.join(os.path.dirname(os.path.dirname(HERE)), "campaign", "round_sweep.tsv")
    rows = []
    for line in open(tsv):
        if line.lstrip("#").startswith("a") and "\t" in line:
            f = line.lstrip("#").rstrip("\n").split("\t")
            if f[0].startswith("a") and f[1].isdigit():
                rows.append((int(f[1]), float(f[10]), float(f[13])))     # alpha, d90, flat_pv
    rows.sort()
    al = np.array([r[0] for r in rows]); d90 = np.array([r[1] for r in rows])
    pv = np.array([r[2] for r in rows])

    fig = plt.figure(figsize=(ps.COL2, 2.7))
    # column 3 is a spacer: panel (d) carries a y-label on each side, which would otherwise
    # collide with the last Ronchigram
    gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 1, 0.34, 1.85], wspace=0.20)

    rng = np.random.default_rng(7)
    for j, alpha in enumerate(a.alphas):
        C3, C1 = plan[alpha]
        P = R.probe(alpha, C3, C1)
        img, lim = R.ronchigram(P, R.L_BOX, alpha, rng)
        ax = fig.add_subplot(gs[0, j])
        ps.imshow_clean(ax, img, cmap="gray", extent=[-lim, lim, -lim, lim],
                        vmax=np.percentile(img, 99.5))
        th = np.linspace(0, 2 * np.pi, 200)
        ax.plot(alpha * np.cos(th), alpha * np.sin(th), color=ps.YELLOW, lw=0.7, ls=(0, (3, 2)))
        ax.set_title(f"{alpha} mrad", pad=3)
        ps.panel(ax, "abc"[j], dx=-0.06, dy=1.13)
        flat = pv[al == alpha][0]
        ax.set_xlabel(("flat" if flat <= R.PLAN.FLAT_TOL else "not flat") + f"  ({flat:.1f} rad)",
                      labelpad=2, color=(ps.GREEN if flat <= R.PLAN.FLAT_TOL else ps.GREY))

    axd = fig.add_subplot(gs[0, 4])
    axd.semilogy(al, pv, "-o", color=ps.PURPLE, label="Ronchigram flatness")
    axd.axhline(R.PLAN.FLAT_TOL, color=ps.PURPLE, lw=0.7, ls=":")
    axd.text(27, R.PLAN.FLAT_TOL * 0.52, "λ/4", color=ps.PURPLE, ha="left", fontsize=7)
    axd.set_xlabel("convergence semi-angle α (mrad)")
    axd.set_ylabel("non-defocus P–V (rad)", color=ps.PURPLE)
    axd.tick_params(axis="y", colors=ps.PURPLE)
    axd.set_ylim(0.1, 400); axd.set_xlim(25, 125)
    ps.panel(axd, "d", dx=-0.22)

    axr = axd.twinx()
    axr.semilogy(al, d90, "-s", color=ps.ORANGE, label="probe diameter")
    axr.axhline(4.0, color=ps.ORANGE, lw=0.7, ls=":")
    axr.text(123, 4.5, "4 Å target", color=ps.ORANGE, fontsize=7, ha="right")
    axr.set_ylabel("probe d90 (Å)", color=ps.ORANGE)
    axr.tick_params(axis="y", colors=ps.ORANGE)
    axr.set_ylim(0.5, 200)
    axr.spines["right"].set_visible(True); axr.spines["right"].set_linewidth(0.6)
    axr.spines["top"].set_visible(False)

    h1, l1 = axd.get_legend_handles_labels(); h2, l2 = axr.get_legend_handles_labels()
    axd.legend(h1 + h2, l1 + l2, loc="upper left", bbox_to_anchor=(0.02, 1.0))

    ps.save(fig, "fig1_probe", a.out)


if __name__ == "__main__":
    main()
