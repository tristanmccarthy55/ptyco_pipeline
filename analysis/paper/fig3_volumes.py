#!/usr/bin/env python
"""@file fig3_volumes.py
@brief Figure 3 -- the reconstructions behind Fig. 2, with the answer drawn on.

ONE CLAIM: at low alpha a column is one unbroken streak; by 90 mrad it separates into individual
atoms at the right depths, and the blind finder puts them there.

Depth (x-z) sections through a row of B-site columns -- Ti alternating with apical O every 1.95 A
along the same column, the hardest depth target in the structure -- at two alphas. Rings are
atomfind's blind detections, dots the known atom positions, so a ring with no dot is a spurious
atom and a dot with no ring a miss. Cropped to the atomic slab; the z-vacuum the full-box
reconstruction uses to absorb surface artefacts is outside the frame.

    ~/hyperspy-bundle/bin/python analysis/paper/fig3_volumes.py \
        --recons ~/Desktop/thin_ab_af_final --atomfind ~/Desktop/thin_ab_af_final/out \
        --gt ~/Desktop/thin_ab_af/gtdata
"""
from __future__ import annotations
import argparse, glob, importlib.util, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pubstyle as ps                                                    # noqa: E402
sys.path.insert(0, os.path.dirname(HERE))
from atomfind import align                                               # noqa: E402

BOXZ, ZVAC = 27.525, 4.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recons", required=True)
    ap.add_argument("--atomfind", required=True)
    ap.add_argument("--gt", required=True)
    ap.add_argument("--alphas", type=int, nargs=2, default=[50, 90])
    ap.add_argument("--half-A", type=float, default=0.25, help="half-thickness of the x-z slab (Å)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    ps.use()
    import matplotlib.pyplot as plt
    M = importlib.util.spec_from_file_location(
        "mep", os.path.join(os.path.dirname(HERE), "make_mep_volumes_fig.py"))
    mep = importlib.util.module_from_spec(M); M.loader.exec_module(mep)

    fig, axes = plt.subplots(1, 2, figsize=(ps.COL2, 3.0), sharey=True,
                            gridspec_kw=dict(wspace=0.06))

    for j, alpha in enumerate(a.alphas):
        d = glob.glob(os.path.join(os.path.expanduser(a.recons), f"recon_af_a{alpha}_lab_NL*"))
        if not d:
            raise SystemExit(f"no recon dir for a{alpha}")
        NL = int(re.search(r"_NL(\d+)$", d[0].rstrip("/")).group(1))
        h5 = sorted(glob.glob(os.path.join(d[0], "**", "*_recons.h5"), recursive=True),
                    key=os.path.getmtime)[-1]
        V, dx, cfg, pos, Z = mep.load(h5, NL, os.path.expanduser(a.gt))
        V = mep.deplane(V); dz = cfg.dz; nL, ny, nx = V.shape

        run = os.path.join(os.path.expanduser(a.atomfind), f"atomfind_a{alpha}")
        found = np.load(os.path.join(run, "found_atoms.npy"))
        al = align.refine_with_atoms(align.register(V, dx, pos, Z, cfg), found, pos, Z, cfg)
        gr, gc, gl = al.site_to_index(pos[:, 0], pos[:, 1], pos[:, 2])

        row = mep.pick_row(found, 22, ny, dx)              # a row of B-site (Ti + apical O) columns
        half = a.half_A / dx
        sec = mep.section(V, row, half)

        ax = axes[j]
        ax.imshow(sec, cmap="magma", aspect="equal", origin="upper",
                  extent=[0, nx * dx, (nL - 0.5) * dz, -0.5 * dz],
                  vmin=np.percentile(sec, 1), vmax=np.percentile(sec, 99.7))
        m = (np.abs(gr - row) <= 2 * half) & (gc >= 0) & (gc < nx)
        fm = np.abs(found["row"] - row) <= 2 * half
        for sp, (nm, col) in ps.SPECIES.items():
            q = m & (Z == sp)
            if q.any():
                ax.scatter(gc[q] * dx, gl[q] * dz, s=7, marker="o", c=col, alpha=0.95,
                           linewidths=0, zorder=2)
            f = found[fm & (found["species"] == sp)]
            if len(f):
                ax.scatter(f["col"] * dx, f["z_A"], s=26, marker="o", facecolors="none",
                           edgecolors=col, linewidths=0.9, zorder=3,
                           label=nm if j == 0 else None)
        ax.set_ylim(BOXZ - ZVAC, ZVAC)                    # the atomic slab only
        ax.set_xlim(0, nx * dx)
        ax.set_xlabel("x (Å)")
        ax.set_title(f"{alpha} mrad", pad=3)
        ps.panel(ax, "ab"[j], dx=-0.02 if j else -0.14, dy=1.04)
        if j == 0:
            ax.set_ylabel("depth z (Å)")
            leg = ax.legend(loc="lower left", ncol=3, handlelength=1.0, columnspacing=0.8,
                            labelcolor="white", bbox_to_anchor=(0.0, -0.02))
            for t in leg.get_texts():
                t.set_fontsize(6.5)

    fig.text(0.5, -0.04, "rings: atoms found blind    ·    dots: true positions",
             ha="center", fontsize=7, color=ps.GREY)
    ps.save(fig, "fig3_volumes", a.out)


if __name__ == "__main__":
    main()
