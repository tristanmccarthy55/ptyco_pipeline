#!/usr/bin/env python
"""
@file make_relaxation_figs.py
@brief Figures for the 2026-09-21 relaxation results: focus fitted inside the solver, the completed
       step-1 legs, frozen phonons, and non-round aberrations. Every value is read from run output.

Figures written to aberration_experiment/figs/<ISO-week>/meeting/:
  figA_dfo.png        the one-parameter (focus-only) probe fit, against the landscape it should descend
  figB_step1.png      fitted focus vs known probe: what the focus search costs in the final numbers
  figC_phonons.png    thermal vibration of the atoms, against the same baseline
  figD_nonround.png   six-fold astigmatism escalated, the corrector residual that round knobs cannot fix

Sources:
  results/<week>/c1_objective.csv                     the fixed-probe focus landscape (step 1)
  <relax>/c1dfo_*/recon_c1dfo_*/analysis/**           focus-only runs: error trace + engine shift header
  <relax>/analysis_0921/atomfind_<tag>/report.json    the finder's numbers for each new leg
  <af_final>/out/atomfind_a<A>/report.json            the known-probe baseline
  campaign/nonround_sweep.tsv                          the aberration strength of each non-round leg

    ~/hyperspy-bundle/bin/python analysis/make_relaxation_figs.py
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from make_meeting_figs import (S1, S2, S3, INK, INK2, MUTED, GRIDC, SPECIES,   # noqa: E402
                               style, tag, save, read_csv_rows, read_reports, results_dir,
                               _load_mod, build_probe, aperture_phase, probe_crop)

RELAX = ["~/Desktop/relax_0921", "~/Desktop/relax_0922"]
AF_FINAL = "~/Desktop/thin_ab_af_final"
BAD = "#c1502a"
_AW = None


def _aw():
    """campaign/aberration_waves.py -- one definition of "waves at the edge" for every order,
    shared by the sweep, the figures and the page so they cannot disagree."""
    global _AW
    if _AW is None:
        _AW = _load_mod("aberration_waves", os.path.join(REPO, "campaign", "aberration_waves.py"))
    return _AW

NR_PROBE_HALF = 8.0     # probe panels are cropped to +-8 A, covering the worst d90 (7.4 A)


# The defocus-only runs come in three variants, told apart by the tag after "dfo" in the dir name
# (campaign/run_c1_search.sh): bare = object free and step uncapped (2026-09-21); "oinf" = the object
# FROZEN in the full engine (experiment A); "c<cap>" = the focus step capped (experiment A2).
DFO_VARIANTS = {"": "object free", "oinf": "object frozen", "c0.05": "step capped"}
DFO_RE = re.compile(r"recon_c1dfo_a(?P<alpha>\d+)_lab_df(?P<c1>-?\d+)_ps\d+(?:x\d+)?dfo(?P<tag>[^_]*)_n")


def dfo_runs(relax):
    """Focus-only runs across every root: alpha, start C1, variant, the engine's cumulative
    shift, the final error, and whether the run crashed."""
    out = []
    for root in ([relax] if isinstance(relax, str) else relax):
        for d in sorted(glob.glob(os.path.expanduser(os.path.join(root, "c1dfo_*", "recon_c1dfo_*")))):
            m = DFO_RE.search(os.path.basename(d))
            if not m:
                continue
            rec = dict(alpha=int(m["alpha"]), c1=float(m["c1"]),
                       variant=DFO_VARIANTS.get(m["tag"], m["tag"] or "object free"), dir=d)
            tr = glob.glob(os.path.join(d, "analysis", "**", "*_error_trace.csv"), recursive=True)
            h5 = glob.glob(os.path.join(d, "analysis", "**", "*_recons.h5"), recursive=True)
            if not tr or not h5:
                out.append(dict(rec, shift=np.nan, err=np.nan, crashed=True))
                continue
            shift = np.nan
            for line in open(tr[-1]):
                if line.startswith("# probe_defocus_shift_A="):
                    shift = float(line.split("=")[1])
            rows = read_csv_rows(tr[-1])
            out.append(dict(rec, shift=shift, err=float(rows[-1]["fourier_error"]), crashed=False))
    return out


def dfo_trajectory(run_dir):
    """The focus fit iteration by iteration, from the solver's own log: the cumulative C1 shift
    printed every iteration, plus the index where the presolve hands over to the full engine.

    The shift is the only per-iteration record there is -- no error history is written per
    iteration, and the sample's depth centroid is NOT logged in this batch (the 2026-09-22 handoff
    says it is; it is not). Returns (shift array, handover index or None)."""
    logs = sorted(glob.glob(os.path.join(run_dir, "slurm_*.out")))
    if not logs:
        return np.array([]), None
    shifts, handover = [], None
    for line in open(logs[-1], errors="replace"):
        if "defocus-only probe update: cumulative C1 shift" in line:
            shifts.append(float(line.rsplit("shift", 1)[1].split("A")[0]))
        elif "Finished MLs solver" in line and handover is None and shifts:
            handover = len(shifts)
    return np.array(shifts), handover


def finder(path):
    return json.load(open(os.path.join(path, "report.json")))["finder"]["v3"]


def af(relax, tag_):
    """The finder's report for a leg, from whichever root holds it (each batch writes its own
    analysis_<date> dir, so a fixed one silently returned None for the newer legs)."""
    for root in ([relax] if isinstance(relax, str) else relax):
        for d in sorted(glob.glob(os.path.expanduser(os.path.join(root, "analysis_*", f"atomfind_{tag_}")))):
            if os.path.isfile(os.path.join(d, "report.json")):
                return finder(d)
    return None


def truth_c1(alpha):
    for line in open(os.path.join(REPO, "campaign", "round_sweep.tsv")):
        if re.match(rf"#?a0*{alpha}\t", line):
            return float(line.lstrip("#").split("\t")[4])
    raise SystemExit(f"no a{alpha} row in round_sweep.tsv")


# ----------------------------------------------------------------- A: focus fitted in the solver
VAR_STYLE = {"object free": ("o", "full"), "object frozen": ("s", "open"), "step capped": ("^", "full")}


def _dfo_panel(ax, runs, alpha, lim, letter, title):
    """Focus error after fitting against the focus error handed in, one series per variant.
    On the diagonal the fit did nothing; on y = 0 it recovered the truth."""
    t = truth_c1(alpha)
    ax.plot([-lim, lim], [-lim, lim], "--", color=MUTED, lw=1.3, zorder=1)
    ax.axhline(0, color=S3, lw=1.8, zorder=1)
    for var, (mk, fill) in VAR_STYLE.items():
        rs = [r for r in runs if r["alpha"] == alpha and r["variant"] == var and not r["crashed"]]
        if not rs:
            continue
        x = [r["c1"] - t for r in rs]
        y = [r["c1"] + r["shift"] - t for r in rs]
        kw = dict(mfc="none", mew=1.9, ms=13) if fill == "open" else dict(ms=8)
        ax.plot(x, y, mk, color=INK if fill == "open" else S1, ls="none", label=var, zorder=3, **kw)
    crashed = [r["c1"] - t for r in runs if r["alpha"] == alpha and r["crashed"]]
    if crashed:
        ax.plot(crashed, [-lim * 0.9] * len(crashed), "x", color=BAD, ms=12, mew=2.4, zorder=3,
                label="diverged (no reconstruction)")
    ax.annotate("focus never moved", (lim * 0.58, lim * 0.58), rotation=39, fontsize=8.5,
                color=INK2, ha="center", va="bottom")
    free = [r for r in runs if r["alpha"] == alpha and r["variant"] == "object free"
            and not r["crashed"] and abs(r["shift"]) > 1.0]
    if len(free) >= 3:
        ymean = float(np.mean([r["c1"] + r["shift"] - t for r in free]))
        ax.annotate(f"with the object free, every fit that moves\nends near {ymean:+.0f} Å, whatever it started from",
                    (-lim * 0.96, ymean - lim * 0.10), fontsize=8.5, color=S1)
    ax.annotate("recovering the truth would land here", (-lim * 0.96, lim * 0.06), fontsize=8.5, color=S3)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("focus error handed to the solver (Å)")
    ax.set_ylabel("focus error after fitting (Å)")
    ax.legend(loc="lower right", fontsize=9)
    tag(ax, letter, title)


def figA_dfo(a):
    runs = dfo_runs(a.relax)
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 4.6))
    _dfo_panel(axes[0], runs, 70, 36, "a", "α = 70 mrad: the fit moves, but not to the truth")
    _dfo_panel(axes[1], runs, 90, 40, "b", "α = 90 mrad: the fit does not move at all")

    # (c) the a70 moves drawn on the landscape they should be descending. Each arrow starts at the
    # landscape value for the probe the solver was handed and ends where the fit actually finished.
    ax = axes[2]
    rows = read_csv_rows(os.path.join(results_dir(), "c1_objective.csv"))
    pts = sorted([(float(r["dc1"]), float(r["final_error"])) for r in rows
                  if int(r["alpha"]) == 70 and int(r["niter"]) == 200])
    xs, ys = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
    ax.plot(xs, ys, "-", color=S1, lw=1.8, alpha=.55, label="probe held fixed (the landscape)")
    t = truth_c1(70)
    for r in (r for r in runs if r["alpha"] == 70 and not r["crashed"] and r["variant"] == "object free"):
        x0, x1, y1 = r["c1"] - t, r["c1"] + r["shift"] - t, r["err"]
        if abs(x1 - x0) < 1.0:
            ax.plot([x1], [y1], "o", color=INK2, ms=7, mec="white", mew=1.1, zorder=5)
            continue
        y0 = float(np.interp(x0, xs, ys))
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=BAD, lw=1.6, shrinkA=3, shrinkB=3), zorder=4)
        ax.plot([x1], [y1], "o", color=BAD, ms=7, zorder=5)
    ax.plot([], [], "o", color=BAD, ms=7, label="where the fit ended")
    ax.annotate("Object-free runs only: the frozen-object runs never\nrefine the object at full resolution, so their error is\nnot comparable with this landscape.",
                (0.02, 0.985), xycoords="axes fraction", fontsize=8.2, color=INK2, va="top")
    ax.annotate("wherever it moves, it ends\nworse than it started", (0.30, 0.50),
                xycoords="axes fraction", fontsize=8.8, color=BAD, va="top")
    ax.set_yscale("log"); ax.set_yticks([5, 10, 20, 50, 100, 200, 400]); ax.minorticks_off()
    ax.set_yticklabels(["5", "10", "20", "50", "100", "200", "400"])
    ax.set_xlim(-40, 40); ax.set_ylim(4, 900)
    ax.set_xlabel("focus error of the probe (Å)")
    ax.set_ylabel("mismatch to the measured data")
    ax.legend(loc="lower center", fontsize=9)
    tag(ax, "c", "Measured against the landscape it should descend")
    fig.tight_layout()
    return save(fig, "figA_dfo.png")


# ----------------------------------------------------------------- F: the fit, iteration by iteration
def figF_dfo_traj(a):
    """What the focus actually did while the solver ran. The only per-iteration record the engine
    writes is its own cumulative shift, so that is what is plotted; the dashed line marks the
    handover from the presolve to the full-resolution engine."""
    runs = [r for r in dfo_runs(a.relax) if not r["crashed"]]
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.6))
    for ax, alpha in zip(axes, (70, 90)):
        t = truth_c1(alpha)
        sel = [r for r in runs if r["alpha"] == alpha]
        starts = sorted({r["c1"] - t for r in sel})
        cmap = plt.get_cmap("coolwarm")
        col = {d: cmap(0.5 + 0.5 * d / (max(abs(min(starts)), abs(max(starts))) or 1)) for d in starts}
        drawn = set()
        hands = []
        for r in sorted(sel, key=lambda r: r["c1"]):
            sh, hand = dfo_trajectory(r["dir"])
            if sh.size == 0:
                continue
            d0 = r["c1"] - t
            ls = {"object free": "-", "object frozen": "--", "step capped": ":"}[r["variant"]]
            lab = None
            if (d0, "c") not in drawn:
                lab = f"start {d0:+.0f} Å"; drawn.add((d0, "c"))
            ax.plot(np.arange(sh.size), d0 + sh, ls, color=col[d0], lw=1.7, label=lab)
            if hand:
                hands.append(hand)
        if hands:
            h = int(np.median(hands))
            ax.axvline(h, color=MUTED, lw=1.1, ls=(0, (2, 3)))
            ax.annotate("full-resolution engine starts here", (h, 0.30), xycoords=("data", "axes fraction"),
                        xytext=(-5, 0), textcoords="offset points", fontsize=8.2, color=INK2,
                        rotation=90, ha="right", va="center")
        ax.axhline(0, color=S3, lw=1.8)
        ax.annotate("the true focus", (2, 0), xytext=(0, 4), textcoords="offset points",
                    fontsize=8.5, color=S3)
        lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + 0.42 * (hi - lo))     # headroom for the legend
        ax.set_xlabel("solver iteration (presolve, then the full-resolution engine)")
        ax.set_ylabel("focus error of the probe (Å)")
        ax.legend(loc="upper left", fontsize=8.5, ncol=3)
        tag(ax, "ab"[alpha == 90], f"α = {alpha} mrad")
    axes[0].annotate("solid = object free · dashed = object frozen · dotted = step capped",
                     (0.02, 0.085), xycoords="axes fraction", fontsize=8.5, color=INK2)
    axes[0].annotate("in the presolve every run that starts at or above the true focus\n"
                     "is pulled to the same wrong place, about 12 Å below it",
                     (0.30, 0.60), xycoords="axes fraction", fontsize=8.5, color=INK2)
    axes[1].annotate("no start moves by as much as 1 Å:\nat this aperture the one-parameter\ngradient carries no usable signal",
                     (0.36, 0.50), xycoords="axes fraction", fontsize=8.8, color=BAD)
    fig.tight_layout()
    return save(fig, "figF_dfo_traj.png")


def _species_panel(ax, groups, title, letter):
    """Grouped bars: one group per configuration, one bar per species (found %)."""
    n = len(groups); w = 0.26
    xi = np.arange(n)
    for j, (sp, col) in enumerate(SPECIES.items()):
        ax.bar(xi + (j - 1) * w, [100 * g[1][sp]["recall_bulk"] for g in groups], w, color=col, label=sp)
    ax.set_xticks(xi); ax.set_xticklabels([g[0] for g in groups])
    ax.set_ylim(0, 126); ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_ylabel("atoms found (%)")
    ax.legend(loc="upper left", ncol=3)
    tag(ax, letter, title)


def _error_panel(ax, groups, title, letter, extra=None):
    xi = np.arange(len(groups)); w = 0.34
    ax.bar(xi - w / 2, [g[1]["z_rms_A"] for g in groups], w, color=S1, label="depth error")
    ax.bar(xi + w / 2, [g[1]["xy_rms_A"] for g in groups], w, color=S3, label="in-plane error")
    for i, g in enumerate(groups):
        ax.annotate(f"{g[1]['z_rms_A']:.2f}", (i - w / 2, g[1]["z_rms_A"]), xytext=(0, 3),
                    textcoords="offset points", ha="center", fontsize=8.5, color=INK)
    if extra:
        ax.axhline(extra[0], color=MUTED, ls="--", lw=1.4)
        ax.annotate(extra[1], (len(groups) - 0.5, extra[0]), xytext=(0, 4), textcoords="offset points",
                    fontsize=8.5, color=INK2, ha="right")
    ax.set_xticks(xi); ax.set_xticklabels([g[0] for g in groups])
    ax.set_ylabel("error (Å)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    tag(ax, letter, title)


# ----------------------------------------------------------------- B: step 1 complete
def figB_step1(a):
    base = read_reports(a.af)
    groups = []
    for alpha in (70, 90):
        groups.append((f"α {alpha}\nknown probe", base[alpha]))
        r = af(a.relax, f"a{alpha}_fitC1")
        if r:
            groups.append((f"α {alpha}\nfocus fitted", r))
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.3))
    _species_panel(axes[0], groups, "Atoms found, by species", "a")
    _error_panel(axes[1], groups, "Position error", "b")
    fig.tight_layout()
    return save(fig, "figB_step1.png")


# ----------------------------------------------------------------- C: phonons
def figC_phonons(a):
    base = read_reports(a.af)
    groups = []
    for alpha in (70, 90):
        groups.append((f"α {alpha}\nstatic atoms", base[alpha]))
        r = af(a.relax, f"a{alpha}_ph16")
        if r:
            groups.append((f"α {alpha}\nvibrating", r))
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.3))
    _species_panel(axes[0], groups, "Atoms found, by species", "a")
    # room-temperature 1-D vibration amplitude of Pb, the floor any position measurement inherits
    sig = float(np.sqrt(0.90 / (8 * np.pi ** 2)))
    _error_panel(axes[1], groups, "Position error", "b",
                 extra=(sig, f"Pb vibration amplitude, {sig:.2f} Å"))
    fig.tight_layout()
    return save(fig, "figC_phonons.png")


# ----------------------------------------------------------------- D: non-round aberrations
def nr_kernel(relax, label, el):
    """The matched kernel this rung's own grid legs produced, from whichever batch ran it."""
    for root in ([relax] if isinstance(relax, str) else relax):
        for p in sorted(glob.glob(os.path.expanduser(
                os.path.join(root, "analysis_*", "psf", f"psf_{el}_{label}_vol.npy")))):
            return p
    return None


def nonround_terms(ab):
    """The non-round terms of an aberration dict: C_nm with m > 0, phi angles excluded."""
    return sorted(k for k in ab
                  if len(k) == 3 and k[0] == "C" and k[1:].isdigit() and k[2] != "0")


def nr_legs(a, term="C56"):
    """Every rung carrying ONLY this non-round term, as (waves, label, report), round baseline
    first. Restricted to a single term on purpose: the sweep also holds rungs for other orders
    (two-fold, three-fold, three-lobe) and a combined six-fold-plus-four-fold level, and plotting
    any of those on an axis labelled "six-fold astigmatism" would be a different aberration drawn
    as if it were this one. The wave figure comes from the sweep TSV the sims were run from, so
    the axis cannot drift from what was actually simulated."""
    plan = read_nonround_sweep()
    legs = []
    base = read_reports(a.af).get(70)
    if base:
        legs.append((0.0, "nr0_round", base))
    for label, (ab, c1, w) in plan.items():
        if nonround_terms(ab) != [term]:
            continue
        r = af(a.relax, label)
        if r:
            legs.append((w, label, r))
    return sorted(legs)


def figD_nonround(a):
    legs = nr_legs(a)
    combo = af(a.relax, "nr4_C56_C34")
    x = [l[0] for l in legs]
    fig, axes = plt.subplots(1, 3, figsize=(16.4, 4.5))

    ax = axes[0]
    for sp, col in SPECIES.items():
        ax.plot(x, [100 * l[2][sp]["recall_bulk"] for l in legs], "-o", color=col, label=sp)
    ax.axvline(x[1], color=BAD, lw=1.3, ls=(0, (3, 3)))
    ax.annotate(f"the smallest residual tested is {x[1]:.2f} waves,\nand it has already cost a third of the oxygen\nand failed the species-labelling health check.\nThe usable tolerance is below this line.",
                (0.30, 0.90), xycoords="axes fraction", fontsize=8.5, color=BAD, va="top")
    ax.set_ylim(0, 116); ax.set_xlim(-0.08, max(x) * 1.05)
    ax.set_xlabel("six-fold astigmatism (waves at the aperture edge)")
    ax.set_ylabel("atoms found (%)"); ax.legend(loc="upper right", ncol=3)
    tag(ax, "a", "Atoms found as the corrector's own residual grows")

    ax = axes[1]
    ax.plot(x, [l[2]["z_rms_A"] for l in legs], "-o", color=S1, label="depth error")
    ax.plot(x, [l[2]["xy_rms_A"] for l in legs], "-o", color=S3, label="in-plane error")
    if combo:
        xc = max(x) * 1.14
        ax.plot([xc], [combo["z_rms_A"]], "s", color=S1, mfc="white", mew=1.8, ms=8)
        ax.plot([xc], [combo["xy_rms_A"]], "s", color=S3, mfc="white", mew=1.8, ms=8)
        ax.annotate("+ four-fold", (xc, combo["z_rms_A"]), xytext=(0, 9), textcoords="offset points",
                    fontsize=8.5, color=INK2, ha="center")
    ax.set_ylim(0, None); ax.set_xlim(-0.08, max(x) * 1.25)
    ax.annotate("past the first rung these are the errors of the\nfew atoms still being found, so they flatten out",
                (0.26, 0.16), xycoords="axes fraction", fontsize=8.5, color=INK2)
    ax.set_xlabel("six-fold astigmatism (waves at the aperture edge)")
    ax.set_ylabel("error (Å)"); ax.legend(loc="center right")
    tag(ax, "b", "Position error")

    # (c) the kernel is measured from the SAME reconstruction pipeline, so it says whether the
    # reconstruction or the finder is failing.
    ax = axes[2]
    mmf = _load_mod("make_meeting_figs", os.path.join(HERE, "make_meeting_figs.py"))
    for el, col, ref in (("Pb", S1, a.round_pb), ("Ti", S2, a.round_ti)):
        xs, ys = [], []
        for w, label, _ in legs:
            path = ref if label == "nr0_round" else nr_kernel(a.relax, label, el)
            if path and os.path.exists(os.path.expanduser(path)):
                xs.append(w); ys.append(mmf.kernel_quality(os.path.expanduser(path))[1])
        if xs:
            ax.plot(xs, ys, "-o", color=col, label=f"{el} kernel")
    ax.set_yscale("log"); ax.set_yticks([10, 20, 50, 100, 200, 500]); ax.minorticks_off()
    ax.set_yticklabels(["10", "20", "50", "100", "200", "500"])
    ax.set_xlim(-0.08, max(x) * 1.05)
    ax.set_xlabel("six-fold astigmatism (waves at the aperture edge)")
    ax.set_ylabel("matched-kernel peak / background")
    ax.legend(loc="upper right")
    ax.annotate("the kernel is measured from the same\nreconstruction as the data, so its collapse\nsays the reconstruction is failing,\nnot just the atom finder",
                (0.33, 0.70), xycoords="axes fraction", fontsize=8.5, color=INK2, va="top")
    tag(ax, "c", "Quality of the matched kernel")
    fig.tight_layout()
    return save(fig, "figD_nonround.png")


# --------------------------------------------------- fig E: what non-round aberration looks like
# The tolerance turned out to sit BELOW the smallest rung, so the panels show the threshold region
# (0.1-0.45 waves) and keep one catastrophic level for scale, rather than four catastrophic ones.
NR_SHOW = ["nr0_round", "nr0p1_C56_0p1w", "nr0p2_C56_0p2w", "nr0p45_C56_0p45w", "nr2_C56_1p2w"]


def read_nonround_sweep():
    """The sweep the sims were actually run from: label -> (aberration dict, C1, waves at the edge)."""
    rows = {}
    with open(os.path.join(REPO, "campaign", "nonround_sweep.tsv")) as f:
        rd = csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t")
        for r in rd:
            ab = {"C30": float(r["c3"]), "C50": float(r["c5"])}
            if r["aber_json"] != "-":
                ab = json.loads(r["aber_json"])
            # waves at the aperture edge of the row's LARGEST non-round term, by order, via the
            # campaign helper -- so a two-fold or three-lobe row reports its own waves and not a
            # C56 figure that would be zero for it.
            w = 0.0
            for t in nonround_terms(ab):
                w = max(w, abs(float(ab[t])) / _aw().per_wave(t, float(r["alpha"])))
            rows[r["label"]] = (ab, float(r["c1"]), w)
    return rows


def figE_nonround_probe(a):
    """The six-fold residual, seen three ways: the Ronchigram the operator would judge focus on,
    the wavefront across the aperture, and the probe that reaches the sample."""
    mrf = _load_mod("make_ronchigram_fig", os.path.join(HERE, "make_ronchigram_fig.py"))
    plan = read_nonround_sweep()
    show = [l for l in NR_SHOW if l in plan]
    recall = {l: af(a.relax, l) for l in show}
    base = read_reports(a.af).get(70)
    if base:
        recall[show[0]] = base

    fig, axes = plt.subplots(3, len(show), figsize=(13.0, 9.4))
    for i, lab in enumerate(show):
        ab, c1, w = plan[lab]
        P = build_probe(70, ab, c1)
        d90 = mrf.enclosed(P, 140.0, (0.9,))[0]

        R, hr = mrf.ronchigram(P, 140.0, 70)
        ax = axes[0, i]
        ax.imshow(R, cmap="gray", extent=[-hr, hr, -hr, hr], interpolation="bilinear")
        ax.add_artist(plt.Circle((0, 0), 70, fill=False, color="#ffd24a", lw=1.2, ls=":"))
        ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
        ax.set_title("round probe\n(no six-fold)" if w == 0 else
                     f"six-fold\n{w:.2f} waves" if w < 1 else f"six-fold\n{w:.1f} waves",
                     pad=6, color=INK, fontsize=10.5)
        if i == 0:
            ax.set_ylabel("Ronchigram\n(what the operator sees)", fontsize=10, color=INK2)

        ph, hp = aperture_phase(P, 70)
        ax = axes[1, i]
        ax.imshow(ph, cmap="twilight_shifted", vmin=-np.pi, vmax=np.pi,
                  extent=[-hp, hp, -hp, hp], interpolation="bilinear")
        ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
        if i == 0:
            ax.set_ylabel("aberration across\nthe aperture", fontsize=10, color=INK2)

        I, hI = probe_crop(P, 140.0, NR_PROBE_HALF)
        ax = axes[2, i]
        ax.imshow(I ** 0.3, cmap="magma", extent=[-hI, hI, -hI, hI], interpolation="bilinear")
        ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
        ax.plot([-hI + 0.8, -hI + 2.8], [-hI + 1.0] * 2, color="white", lw=2.4, solid_capstyle="butt")
        if i == 0:
            ax.text(-hI + 0.8, -hI + 1.7, "2 Å", color="white", fontsize=8.5, va="bottom")
            ax.set_ylabel("probe at the sample\n(square root of intensity)", fontsize=10, color=INK2)
        r = recall.get(lab)
        lines = [f"probe diameter {d90:.1f} Å"]
        if r:
            lines.append(f"Pb found {100 * r['Pb']['recall_bulk']:.0f}%")
        ax.set_xlabel("\n".join(lines), fontsize=9,
                      color=(S3 if r and r["Pb"]["recall_bulk"] > 0.9 else INK2), labelpad=4)

    fig.suptitle("a   Six-fold astigmatism: the residual a hexapole corrector cannot remove",
                 x=0.012, ha="left", fontsize=11.5, fontweight="semibold", color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    return save(fig, "figE_nonround_probe.png")


FIGS = {"A": figA_dfo, "B": figB_step1, "C": figC_phonons, "D": figD_nonround,
        "E": figE_nonround_probe, "F": figF_dfo_traj}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--relax", nargs="+", default=RELAX); ap.add_argument("--af", default=AF_FINAL)
    ap.add_argument("--round-pb", default="~/Desktop/thin_ab_af3/psf/psf_Pb_a70_vol.npy",
                    help="the round a70 kernels, the reference the non-round rungs are measured against")
    ap.add_argument("--round-ti", default="~/Desktop/thin_ab_af3/psf/psf_Ti_a70_vol.npy")
    ap.add_argument("--only", nargs="+", default=list(FIGS))
    a = ap.parse_args()
    style()
    for k in a.only:
        FIGS[k.upper()](a)


if __name__ == "__main__":
    main()
