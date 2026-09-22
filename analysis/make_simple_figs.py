#!/usr/bin/env python
"""
@file make_simple_figs.py
@brief The plain-pictures figure set: show the reconstruction and what the atom finder did to it,
       rather than summary statistics about them.

Written 2026-09-22 after the round-two page was judged too abstract. Every figure here is an image
of real output -- a Ronchigram, a phase volume, a kernel -- with the ground truth drawn on. The only
numbers are recall percentages in the panel titles. The reader is meant to be able to look at a
depth section, see the peaks, see which ones carry a detection ring, and judge for themselves
whether the reconstruction or the finder is at fault. That question cannot be answered by a chart
of recall against aberration, which is why it kept being asked.

  fig_probe.png      the probe: Ronchigram the operator would see, and the probe at the specimen
  fig_nonround.png   depth sections as six-fold astigmatism grows, ground truth + found atoms
  fig_dose.png       the same as the electron dose falls
  fig_psf.png        one matched kernel, in-plane and in depth
  fig_focus.png      the in-solver focus fit, one panel

    ~/hyperspy-bundle/bin/python analysis/make_simple_figs.py            # all of them
    ~/hyperspy-bundle/bin/python analysis/make_simple_figs.py --only nonround dose
"""
from __future__ import annotations

import argparse
import csv
import datetime
import glob
import importlib.util
import json
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from atomfind import config, align                                  # noqa: E402

BOXZ, ZVAC = 27.525, 4.0
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
# ground truth vs found: one colour per species, ground truth a filled dot, a detection a ring
SP = {82: ("Pb", "#2a78d6"), 22: ("Ti", "#eb6834"), 8: ("O", "#1baf7a")}


def _load_mod(name, path):
    s = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m)
    return m


def week_dir():
    d = os.path.join(REPO, "aberration_experiment", "figs",
                     datetime.date.today().strftime("%G-W%V"), "simple")
    os.makedirs(d, exist_ok=True)
    return d


def style():
    plt.rcParams.update({
        "figure.facecolor": "white", "savefig.facecolor": "white", "savefig.dpi": 190,
        "savefig.bbox": "tight", "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 10.5,
        "axes.titlesize": 11, "axes.labelsize": 10, "axes.titleweight": "semibold",
        "axes.labelcolor": INK2, "axes.edgecolor": "#c3c2b7", "axes.grid": False,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelsize": 9, "ytick.labelsize": 9,
        "text.color": INK, "legend.frameon": False, "legend.fontsize": 9,
    })


def save(fig, name):
    p = os.path.join(week_dir(), name)
    fig.savefig(p); plt.close(fig)
    print(f"wrote {p}")
    return p


# ------------------------------------------------------------------ loading a leg
def load_leg(h5, nl, gt_dir, af_dir=None):
    """The volume in atomfind's own frame, plus the ground truth and the blind detections mapped
    into it -- so what is drawn is exactly what the finder was looking at."""
    config.set_data_dir(os.path.expanduser(gt_dir))
    cfg = config.preset("thin")
    cfg.recon_vol, cfg.dz = h5, BOXZ / nl
    cfg = cfg.resolve()
    V, dx = align.load_phase(cfg)
    V = align.crop_to_fov(V, dx, cfg)
    pos, Z = align.load_gt(cfg)
    mep = _load_mod("mep", os.path.join(HERE, "make_mep_volumes_fig.py"))
    V = mep.deplane(V)
    found = gtidx = None
    if af_dir and os.path.exists(os.path.join(af_dir, "found_atoms.npy")):
        found = np.load(os.path.join(af_dir, "found_atoms.npy"))
        al = align.refine_with_atoms(align.register(V, dx, pos, Z, cfg), found, pos, Z, cfg)
        gtidx = al.site_to_index(pos[:, 0], pos[:, 1], pos[:, 2])
    return dict(V=V, dx=dx, dz=cfg.dz, pos=pos, Z=Z, found=found, gt=gtidx, cfg=cfg)


def xz_panel(ax, leg, row, half_A=0.25, xlim=None, show_y=True):
    """One x-z depth section with the ground truth drawn on it and every blind detection ringed."""
    V, dx, dz = leg["V"], leg["dx"], leg["dz"]
    nL, ny, nx = V.shape
    hp = half_A / dx
    r0, r1 = int(round(row - hp)), int(round(row + hp)) + 1
    s = V[:, max(r0, 0):min(r1, ny), :].mean(1)
    ax.imshow(s, cmap="magma", aspect="equal", origin="upper",
              extent=[0, nx * dx, (nL - 0.5) * dz, -0.5 * dz],
              vmin=np.percentile(s, 1), vmax=np.percentile(s, 99.7), interpolation="nearest")
    # ground truth: every atom whose column falls inside this slab
    if leg["gt"] is not None:
        gr, gc, gl = leg["gt"]
        m = np.abs(gr - row) <= hp
        seen = set()
        for z_, col, sp in zip(gl[m], gc[m], leg["Z"][m]):
            if sp in SP:
                ax.plot(col * dx, (z_ + 0.5) * dz, ".", color=SP[sp][1], ms=5.0, alpha=.95, zorder=3)
                seen.add(int(sp))
        leg.setdefault("species_seen", set()).update(seen)
    # blind detections: a ring, so a ring with no dot under it is spurious and a dot with no ring a miss
    if leg["found"] is not None:
        f = leg["found"]
        m = np.abs(f["row"] - row) <= hp * 2.2
        for a_ in f[m]:
            c = SP.get(int(a_["species"]), (None, "white"))[1]
            ax.plot(a_["col"] * dx, (a_["layer"] + 0.5) * dz, "o", mfc="none", mec=c, mew=1.2,
                    ms=7.5, zorder=4)
    for zb in (ZVAC, BOXZ - ZVAC):                      # the vacuum bands: artefacts belong here
        ax.axhline(zb, color="#7fd4ff", lw=0.8, ls=(0, (4, 3)), alpha=.85)
    if xlim:
        ax.set_xlim(*xlim)
    ax.set_ylim(BOXZ, 0)
    ax.set_xlabel("x (Å)")
    if show_y:
        ax.set_ylabel("depth z (Å)")
    else:
        ax.set_yticklabels([])


def legend_axes(fig, y=0.005, species=None):
    keep = SP if species is None else {k: v for k, v in SP.items() if k in species}
    h = [plt.Line2D([], [], ls="none", marker=".", color=c, ms=10, label=f"{n} in the structure")
         for _, (n, c) in keep.items()]
    h += [plt.Line2D([], [], ls="none", marker="o", mfc="none", mec=INK2, mew=1.5, ms=9,
                     label="atom the finder reported")]
    h += [plt.Line2D([], [], color="#7fd4ff", ls=(0, (4, 3)), lw=1.2, label="edge of the vacuum")]
    fig.legend(handles=h, loc="lower center", ncol=5, bbox_to_anchor=(0.5, y), fontsize=9.5)


def recall(af_dir):
    p = os.path.join(af_dir, "report.json")
    if not os.path.exists(p):
        return None
    r = json.load(open(p))["finder"]["v3"]
    return {k: 100 * r[k]["recall_bulk"] for k in ("Pb", "Ti", "O")}


# ------------------------------------------------------------------ the figures
NR_ROOT = "~/Desktop/relax_0922/atomfind_results_nr0p1_C56_0p1w-nr0p2_C56_0p2w-nr0p3_C56_0p3w-nr0p45_C56_0p45w_20260921_225614"
NR_AF = "~/Desktop/relax_0922/analysis_0922"
DOSE_ROOT = "~/Desktop/relax_0917/atomfind_results_a70-90_dose1e7-1e6-1e5-1e4_20260917_181903"
DOSE_AF = "~/Desktop/relax_0917/dose_analysis"
BASE_RECON = "~/Desktop/thin_ab_af_final/recon_af_a70_lab_NL14"
BASE_AF = "~/Desktop/thin_ab_af_final/out/atomfind_a70"


def _h5(d):
    g = glob.glob(os.path.join(os.path.expanduser(d), "analysis", "**", "*_recons.h5"), recursive=True)
    return sorted(g, key=os.path.getmtime)[-1] if g else None


def _series(rows, gt, half_A, title, name, note):
    """One row of depth sections, same scene every panel, only the condition changing."""
    legs = []
    for lab, recon, af_dir in rows:
        h5 = _h5(recon)
        if not h5:
            print(f"  skip {lab}: no h5 under {recon}")
            continue
        af_dir = os.path.expanduser(af_dir) if af_dir else None
        legs.append((lab, load_leg(h5, 14 if "NL14" in recon else 23, gt, af_dir),
                     recall(af_dir) if af_dir else None))
    if not legs:
        raise SystemExit("nothing to draw")
    mep = _load_mod("mep", os.path.join(HERE, "make_mep_volumes_fig.py"))
    # ONE row of columns for every panel: the same physical slice of the same object throughout,
    # chosen on the baseline, so any difference between panels is the condition and nothing else.
    base = legs[0][1]
    row = mep.pick_row(base["found"], 22, base["V"].shape[1], base["dx"]) if base["found"] is not None \
        else base["V"].shape[1] / 2
    nx, dx = base["V"].shape[2], base["dx"]
    xlim = (max(0, row * dx - 9), min(nx * dx, row * dx + 9))
    xlim = (nx * dx / 2 - 9, nx * dx / 2 + 9)

    fig, axes = plt.subplots(1, len(legs), figsize=(3.55 * len(legs), 5.1), squeeze=False)
    for i, (lab, leg, rc) in enumerate(legs):
        xz_panel(axes[0, i], leg, row, half_A=half_A, xlim=xlim, show_y=(i == 0))
        if rc is None:
            t = f"{lab}\nno kernel could be measured,\nso the finder never ran"
            axes[0, i].set_title(t, fontsize=10.5, pad=7, color="#c2490a")
        else:
            axes[0, i].set_title(f"{lab}\nPb {rc['Pb']:.0f}%   Ti {rc['Ti']:.0f}%   O {rc['O']:.0f}%",
                                 fontsize=10.5, pad=7)
    fig.suptitle(title, x=0.008, ha="left", fontsize=12.5, fontweight="semibold", y=1.0)
    fig.text(0.008, -0.055, note, ha="left", va="top", fontsize=9.5, color=INK2, wrap=True)
    fig.tight_layout(rect=[0, 0.055, 1, 0.965])
    seen = set()
    for _, leg, _ in legs:
        seen |= leg.get("species_seen", set())
    legend_axes(fig, y=-0.005, species=seen)
    return save(fig, name)


def fig_nonround(a):
    rows = [("no six-fold\n(the baseline)", BASE_RECON, BASE_AF)]
    for lab, tag in [("0.10 waves", "0p1"), ("0.20 waves", "0p2"), ("0.45 waves", "0p45")]:
        rows.append((f"six-fold\n{lab}",
                     f"{NR_ROOT}/recon_af_nr{tag}_C56_{tag}w_lab_NL14",
                     f"{NR_AF}/atomfind_nr{tag}_C56_{tag}w"))
    return _series(rows, a.gt, a.half, "Six-fold astigmatism: the same columns, the same slice, four probes",
                   "fig_nonround.png",
                   "The same row of B-site columns in every panel, at 70 mrad. Titanium alternates with apical "
                   "oxygen every 1.95 Å down each column. Dots are where the atoms really are; a ring is an atom "
                   "the finder reported. A dot with no ring is a miss, a ring with no dot is spurious.")


def fig_dose(a):
    rows = [("noiseless\n(the baseline)", BASE_RECON, BASE_AF)]
    for d in ("1e7", "1e6", "1e5", "1e4"):
        rows.append((f"{d.replace('1e','10^')} e/Å²",
                     f"{DOSE_ROOT}/recon_af_a70_lab_dose{d}_NL14",
                     f"{DOSE_AF}/atomfind_a70_dose{d}"))
    return _series(rows, a.gt, a.half, "Electron dose: the same columns as the counting noise rises",
                   "fig_dose.png",
                   "The same row of B-site columns at 70 mrad, reconstructed from Poisson-sampled data. "
                   "Dots are the real atoms, rings are the finder's detections.")


def fig_probe(a):
    """What the aberration does to the probe: the Ronchigram an operator would tune on, and the
    probe that lands on the specimen. One identical patch of amorphous film in every panel."""
    mrf = _load_mod("make_ronchigram_fig", os.path.join(HERE, "make_ronchigram_fig.py"))
    mmf = _load_mod("make_meeting_figs", os.path.join(HERE, "make_meeting_figs.py"))
    aw = _load_mod("aberration_waves", os.path.join(REPO, "campaign", "aberration_waves.py"))
    rows = {}
    with open(os.path.join(REPO, "campaign", "nonround_sweep.tsv")) as f:
        for r in csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"):
            rows[r["label"]] = r
    show = ["nr0_round", "nr0p1_C56_0p1w", "nr0p2_C56_0p2w", "nr0p45_C56_0p45w"]
    fig, axes = plt.subplots(2, len(show), figsize=(3.3 * len(show), 7.0))
    for i, lab in enumerate(show):
        r = rows[lab]
        ab = ({"C30": float(r["c3"]), "C50": float(r["c5"])} if r["aber_json"] == "-"
              else {k: float(v) for k, v in json.loads(r["aber_json"]).items()})
        w = max([abs(ab[t]) / aw.per_wave(t, 70.0) for t in ab if t in aw.ORDER and t[-1] != "0"] or [0.0])
        P = mmf.build_probe(70, ab, float(r["c1"]))
        d90 = mrf.enclosed(P, 140.0, (0.9,))[0]
        R, hr = mrf.ronchigram(P, 140.0, 70)
        ax = axes[0, i]
        ax.imshow(R, cmap="gray", extent=[-hr, hr, -hr, hr], interpolation="bilinear")
        ax.add_artist(plt.Circle((0, 0), 70, fill=False, color="#ffd24a", lw=1.2, ls=":"))
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("round probe" if w == 0 else f"six-fold {w:.2f} waves", pad=7, fontsize=11)
        if i == 0:
            ax.set_ylabel("Ronchigram\nwhat the operator sees", fontsize=10, color=INK2)
        I, hI = mmf.probe_crop(P, 140.0, 8.0)
        ax = axes[1, i]
        ax.imshow(I ** 0.3, cmap="magma", extent=[-hI, hI, -hI, hI], interpolation="bilinear")
        ax.set_xticks([]); ax.set_yticks([])
        ax.plot([-hI + 0.8, -hI + 2.8], [-hI + 1.0] * 2, color="white", lw=2.6, solid_capstyle="butt")
        ax.set_xlabel(f"probe diameter {d90:.1f} Å", fontsize=10,
                      color=(INK if w == 0 else "#c2490a"), labelpad=5)
        if i == 0:
            ax.text(-hI + 0.8, -hI + 1.7, "2 Å", color="white", fontsize=9, va="bottom")
            ax.set_ylabel("the probe at the specimen", fontsize=10, color=INK2)
    fig.suptitle("The probe, at 70 mrad, as six-fold astigmatism is added",
                 x=0.008, ha="left", fontsize=12.5, fontweight="semibold")
    fig.text(0.008, 0.012,
             "Same patch of amorphous film in every Ronchigram, so only the probe differs. The quoted "
             "diameter is the circle holding 90 % of the probe's intensity: it barely moves, while the "
             "six-pointed tails that wreck the reconstruction grow.",
             ha="left", va="bottom", fontsize=9.5, color=INK2)
    fig.tight_layout(rect=[0, 0.055, 1, 0.965])
    return save(fig, "fig_probe.png")


def fig_psf(a):
    """One matched kernel: the response of the whole instrument-plus-solver to a single atom."""
    vol = os.path.expanduser(a.psf)
    K = np.load(vol)
    K = np.angle(K) if np.iscomplexobj(K) else K
    l, r, c = np.unravel_index(np.argmax(K), K.shape)
    dz = BOXZ / 14.0
    dx = 0.0492
    half = 14
    xy = K[l, r - half:r + half + 1, c - half:c + half + 1]
    xz = K[:, r, c - half:c + half + 1]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.3))
    e = half * dx
    axes[0].imshow(xy, cmap="magma", extent=[-e, e, e, -e], interpolation="nearest")
    axes[0].set_title("across the beam", fontsize=11, pad=6)
    axes[0].set_xlabel("x (Å)"); axes[0].set_ylabel("y (Å)")
    axes[1].imshow(xz, cmap="magma", aspect="auto", extent=[-e, e, (K.shape[0] - 0.5) * dz, -0.5 * dz],
                   interpolation="nearest")
    axes[1].set_title("along the beam", fontsize=11, pad=6)
    axes[1].set_xlabel("x (Å)"); axes[1].set_ylabel("depth z (Å)")
    fig.suptitle("One lead atom, as this instrument and this solver render it",
                 x=0.008, ha="left", fontsize=12.5, fontweight="semibold")
    fig.text(0.008, 0.01,
             "Measured, not modelled: a grid of lead atoms is simulated and reconstructed through the same "
             "probe, box and solver settings as the data, and the atoms are averaged. This is the shape the "
             "finder looks for. Peak "
             f"{K[l, r, c]:.2f} rad, about {K.shape[0]} depth slices of {dz:.2f} Å.",
             ha="left", va="bottom", fontsize=9.5, color=INK2)
    fig.tight_layout(rect=[0, 0.06, 1, 0.94])
    return save(fig, "fig_psf.png")


def fig_focus(a):
    """The in-solver focus fit, in one panel: where it ended against where it started."""
    mrf = _load_mod("make_relaxation_figs", os.path.join(HERE, "make_relaxation_figs.py"))
    runs = [r for r in mrf.dfo_runs(mrf.RELAX) if not r["crashed"]]
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    lim = 38
    ax.plot([-lim, lim], [-lim, lim], "--", color=MUTED, lw=1.3)
    ax.axhline(0, color="#1baf7a", lw=2.0)
    for alpha, col, mk in ((70, "#2a78d6", "o"), (90, "#eb6834", "s")):
        t = mrf.truth_c1(alpha)
        rs = [r for r in runs if r["alpha"] == alpha]
        ax.plot([r["c1"] - t for r in rs], [r["c1"] + r["shift"] - t for r in rs], mk,
                color=col, ms=9, ls="none", label=f"α = {alpha} mrad", alpha=.9)
    ax.annotate("the fit did nothing", (lim * 0.55, lim * 0.55), rotation=39, fontsize=9,
                color=INK2, ha="center", va="bottom")
    ax.annotate("a working fit would land on this line", (-lim * 0.95, 1.6), fontsize=9.5, color="#1baf7a")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("focus error the solver was given (Å)")
    ax.set_ylabel("focus error it ended with (Å)")
    ax.legend(loc="lower right")
    ax.set_title("Letting the solver fit focus: 22 runs, none of them recover it",
                 loc="left", fontsize=12, pad=9)
    fig.tight_layout()
    return save(fig, "fig_focus.png")


def fig_didfinderfail(a):
    """The question the summary charts could not answer: when recall falls, is the signal gone, or is
    the finder failing to pick up signal that is plainly there? Same slice twice -- once bare, once
    with the answer drawn on -- so the reader judges the peaks before seeing the detections."""
    mep = _load_mod("mep", os.path.join(HERE, "make_mep_volumes_fig.py"))
    tag, lab = a.fail_leg, a.fail_label
    recon = f"{NR_ROOT}/recon_af_nr{tag}_C56_{tag}w_lab_NL14"
    af_dir = os.path.expanduser(f"{NR_AF}/atomfind_nr{tag}_C56_{tag}w")
    base = load_leg(_h5(BASE_RECON), 14, a.gt, os.path.expanduser(BASE_AF))
    leg = load_leg(_h5(recon), 14, a.gt, af_dir)
    row = mep.pick_row(base["found"], 22, base["V"].shape[1], base["dx"])
    nx, dx = leg["V"].shape[2], leg["dx"]
    xlim = (nx * dx / 2 - 9, nx * dx / 2 + 9)
    rc = recall(af_dir)

    fig, axes = plt.subplots(1, 3, figsize=(12.4, 5.4))
    bare = dict(leg); bare["found"] = None; bare["gt"] = None
    xz_panel(axes[0], bare, row, half_A=a.half, xlim=xlim, show_y=True)
    axes[0].set_title("the reconstruction, as it comes out", fontsize=10.5, pad=7)
    gtonly = dict(leg); gtonly["found"] = None
    xz_panel(axes[1], gtonly, row, half_A=a.half, xlim=xlim, show_y=False)
    axes[1].set_title("where the atoms actually are", fontsize=10.5, pad=7)
    xz_panel(axes[2], leg, row, half_A=a.half, xlim=xlim, show_y=False)
    axes[2].set_title(f"what the finder reported\nPb {rc['Pb']:.0f}%   Ti {rc['Ti']:.0f}%   O {rc['O']:.0f}%",
                      fontsize=10.5, pad=7)
    fig.suptitle(f"Did the reconstruction fail, or did the atom finder?   ({lab}, 70 mrad)",
                 x=0.008, ha="left", fontsize=12.5, fontweight="semibold")
    fig.tight_layout(rect=[0, 0.075, 1, 0.95])
    legend_axes(fig, y=0.012, species=leg.get("species_seen"))
    fig.text(0.008, -0.005,
             "Left: no markers, so the columns can be judged on their own. Middle: the true positions. "
             "Right: the blind detections on top. Columns that are still plainly visible but carry no ring "
             "are the finder's failure; columns that have dissolved into the background are the "
             "reconstruction's.",
             ha="left", va="top", fontsize=9.5, color=INK2)
    return save(fig, f"fig_didfinderfail_{tag}.png")


FIGS = {"probe": fig_probe, "nonround": fig_nonround, "dose": fig_dose, "psf": fig_psf,
        "focus": fig_focus, "blame": fig_didfinderfail}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gt", default="~/Desktop/thin_ab_af/gtdata")
    ap.add_argument("--psf", default="~/Desktop/thin_ab_af3/psf/psf_Pb_a70_vol.npy")
    ap.add_argument("--half", type=float, default=0.25, help="half-thickness of the x-z slab (Å)")
    ap.add_argument("--fail-leg", default="0p2", help="which six-fold rung the blame figure uses")
    ap.add_argument("--fail-label", default="six-fold 0.20 waves")
    ap.add_argument("--only", nargs="+", default=list(FIGS))
    a = ap.parse_args()
    style()
    for k in a.only:
        FIGS[k](a)


if __name__ == "__main__":
    main()
