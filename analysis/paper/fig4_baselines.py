#!/usr/bin/env python
"""@file fig4_baselines.py
@brief Figure 4 -- where the matched kernel actually earns its place.

ONE CLAIM: on the heavy columns any reasonable detector works, and peak picking is marginally
better; the matched single-atom kernel wins on the LIGHT atoms and on precision. Deconvolving
first -- the standard workflow -- makes oxygen worse, not better.

(a) atoms recovered per species; (b) precision, i.e. what fraction of the detections are real.
All methods run on the SAME reconstruction, from one report.json, so they see identical data.
Baselines are 3-D local-maxima peak picking on the raw phase volume and on its Richardson-Lucy and
maximum-entropy deconvolutions (the Ishizuka-style workflow).

    ~/hyperspy-bundle/bin/python analysis/paper/fig4_baselines.py --run ~/Desktop/thin_ab_af_final/out/atomfind_a90
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pubstyle as ps                                                    # noqa: E402

METHODS = [("peaks3d_raw", "peak-pick (raw)", "0.82"),
           ("peaks3d_rl", "peak-pick (RL)", "0.62"),
           ("peaks3d_mem", "peak-pick (MEM)", "0.42"),
           ("v3", "matched kernel (this work)", ps.PURPLE)]
SPECIES = ["Pb", "Ti", "O"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="~/Desktop/thin_ab_af_final/out/atomfind_a90")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    ps.use()
    import matplotlib.pyplot as plt
    run = os.path.expanduser(a.run)
    fin = json.load(open(os.path.join(run, "report.json")))["finder"]
    alpha = "".join(c for c in os.path.basename(run.rstrip("/")) if c.isdigit())
    use = [(k, lab, col) for k, lab, col in METHODS if k in fin]

    fig, (axr, axp) = plt.subplots(1, 2, figsize=(ps.COL15, 2.45),
                                   gridspec_kw=dict(width_ratios=[2.1, 1], wspace=0.42))

    w = 0.2
    x = np.arange(len(SPECIES), dtype=float)
    for i, (key, lab, col) in enumerate(use):
        vals = [100 * fin[key][s]["recall_bulk"] for s in SPECIES]
        axr.bar(x + (i - (len(use) - 1) / 2) * w, vals, w * 0.9, label=lab,
                color=col, edgecolor="0.25", linewidth=0.4)
        print(f"  {lab:28s} precision {fin[key]['precision']:.2f}   " +
              "  ".join(f"{s} {v:.0f}%" for s, v in zip(SPECIES, vals)))
    axr.set_xticks(x); axr.set_xticklabels(SPECIES)
    axr.set_ylabel("atoms recovered (%)"); axr.set_ylim(0, 108)
    # above the axes: inside, it sits on top of the bars
    axr.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2, fontsize=6.2,
               handlelength=1.0, columnspacing=0.8, handletextpad=0.4, borderpad=0.2)
    ps.panel(axr, "a", dx=-0.13)

    xp = np.arange(len(use), dtype=float)
    axp.bar(xp, [100 * fin[k]["precision"] for k, _, _ in use], 0.62,
            color=[c for _, _, c in use], edgecolor="0.25", linewidth=0.4)
    axp.set_xticks(xp)
    axp.set_xticklabels(["raw", "RL", "MEM", "this\nwork"][:len(use)])
    axp.set_ylabel("precision (%)"); axp.set_ylim(0, 108)
    axp.set_title(f"α = {alpha} mrad", pad=4)
    ps.panel(axp, "b", dx=-0.30)

    ps.save(fig, "fig4_baselines", a.out)


if __name__ == "__main__":
    main()
