#!/usr/bin/env python
"""@file fig2_depth.py
@brief Figure 2 -- depth recovery improves as the aperture opens.

ONE CLAIM: atom recovery and depth accuracy both improve with alpha, through the regime where the
Ronchigram has stopped being flat (Fig. 1).

(a) fraction of atoms recovered, per species, over the interior planes;
(b) depth error (RMS against the known structure) with the diffraction depth resolution
    delta_z = lambda/alpha^2 for reference -- the finder fits atom centres, so it sits far below it.

Open symbols mark alphas where the species labels are not reliable (atomfind's own species-confusion
health check, > 5%): there the split between species is label swapping, not physics.

    ~/hyperspy-bundle/bin/python analysis/paper/fig2_depth.py --glob '~/Desktop/thin_ab_af_final/out/atomfind_a*'
"""
from __future__ import annotations
import argparse, glob, importlib.util, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pubstyle as ps                                                    # noqa: E402

LAM = 0.0196877
CONF_MAX = 0.05


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--glob", default="~/Desktop/thin_ab_af_final/out/atomfind_a*")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    ps.use()
    import matplotlib.pyplot as plt
    spec = importlib.util.spec_from_file_location(
        "collate", os.path.join(os.path.dirname(HERE), "collate_atomfind_depth.py"))
    C = importlib.util.module_from_spec(spec); spec.loader.exec_module(C)

    runs = sorted(glob.glob(os.path.expanduser(a.glob)))
    rows = sorted((C.load_run(p) for p in runs), key=lambda r: r["alpha"])
    if not rows:
        raise SystemExit(f"no runs matched {a.glob}")
    al = np.array([r["alpha"] for r in rows], float)
    bad = np.array([(r.get("confusion") or 0) > CONF_MAX for r in rows])
    print("alphas:", al.astype(int), " labels unreliable at:", al[bad].astype(int))

    fig, (axr, axz) = plt.subplots(2, 1, figsize=(ps.COL1, 3.9), sharex=True,
                                   gridspec_kw=dict(hspace=0.12))

    for sp in ("Pb", "Ti", "O"):
        col = ps.SPECIES_NAME[sp]
        y = np.array([100 * (r[f"{sp}_recall_bulk"] or np.nan) for r in rows])
        axr.plot(al, y, "-", color=col, label=sp)
        axr.plot(al[~bad], y[~bad], "o", color=col, mfc=col)
        axr.plot(al[bad], y[bad], "o", color=col, mfc="white", mew=1.0)
    axr.set_ylabel("atoms recovered (%)")
    axr.set_ylim(30, 104)
    axr.legend(loc="lower right", ncol=3, columnspacing=0.9, handlelength=1.4)
    ps.panel(axr, "a")

    zr = np.array([r["z_rms"] for r in rows])
    axz.semilogy(al, zr, "-o", color=ps.GREY, mfc=ps.GREY, label="depth error (RMS)")
    aa = np.linspace(al.min(), al.max(), 100)
    axz.semilogy(aa, LAM / (aa * 1e-3) ** 2, ":", color="0.55", lw=1.0,
                 label="δz = λ/α²")
    axz.set_xlabel("convergence semi-angle α (mrad)")
    axz.set_ylabel("depth (Å)")
    axz.set_ylim(0.2, 12)
    axz.set_yticks([0.3, 1, 3, 10])
    axz.get_yaxis().set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter("%g"))
    axz.legend(loc="upper right")
    ps.panel(axz, "b")

    if bad.any():
        axr.text(0.03, 0.06, "open: species labels unreliable", transform=axr.transAxes,
                 fontsize=6.5, color=ps.GREY)

    ps.save(fig, "fig2_depth", a.out)


if __name__ == "__main__":
    main()
