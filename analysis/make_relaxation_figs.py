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

RELAX = "~/Desktop/relax_0921"
AF_FINAL = "~/Desktop/thin_ab_af_final"
BAD = "#c1502a"
C56_WAVE = 1.004e6      # C56 [A] giving one wave of six-fold astigmatism at the 70 mrad edge
NR_PROBE_HALF = 8.0     # probe panels are cropped to +-8 A, covering the worst d90 (7.4 A)


def dfo_runs(relax):
    """Focus-only runs: (alpha, start C1, engine's cumulative shift, final error, crashed)."""
    out = []
    for d in sorted(glob.glob(os.path.expanduser(os.path.join(relax, "c1dfo_*", "recon_c1dfo_*")))):
        m = re.search(r"recon_c1dfo_a(\d+)_lab_df(-?\d+)_", os.path.basename(d))
        if not m:
            continue
        alpha, c1 = int(m.group(1)), float(m.group(2))
        tr = glob.glob(os.path.join(d, "analysis", "**", "*_error_trace.csv"), recursive=True)
        h5 = glob.glob(os.path.join(d, "analysis", "**", "*_recons.h5"), recursive=True)
        if not tr or not h5:
            out.append(dict(alpha=alpha, c1=c1, shift=np.nan, err=np.nan, crashed=True))
            continue
        shift = np.nan
        for line in open(tr[-1]):
            if line.startswith("# probe_defocus_shift_A="):
                shift = float(line.split("=")[1])
        rows = read_csv_rows(tr[-1])
        out.append(dict(alpha=alpha, c1=c1, shift=shift, err=float(rows[-1]["fourier_error"]), crashed=False))
    return out


def finder(path):
    return json.load(open(os.path.join(path, "report.json")))["finder"]["v3"]


def af(relax, tag_):
    p = os.path.expanduser(os.path.join(relax, "analysis_0921", f"atomfind_{tag_}"))
    return finder(p) if os.path.isfile(os.path.join(p, "report.json")) else None


def truth_c1(alpha):
    for line in open(os.path.join(REPO, "campaign", "round_sweep.tsv")):
        if re.match(rf"#?a0*{alpha}\t", line):
            return float(line.lstrip("#").split("\t")[4])
    raise SystemExit(f"no a{alpha} row in round_sweep.tsv")


# ----------------------------------------------------------------- A: focus fitted in the solver
def figA_dfo(a):
    rows = read_csv_rows(os.path.join(results_dir(), "c1_objective.csv"))
    runs = dfo_runs(a.relax)
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.4))

    # left: what the solver did with the focus, against where it started
    ax = axes[0]
    for alpha, col in ((70, S1), (90, S2)):
        t = truth_c1(alpha)
        good = [(r["c1"] - t, r["c1"] + r["shift"] - t) for r in runs if r["alpha"] == alpha and not r["crashed"]]
        if good:
            ax.plot([g[0] for g in good], [g[1] for g in good], "o", color=col, ms=10, label=f"α = {alpha} mrad")
        for r in (r for r in runs if r["alpha"] == alpha and r["crashed"]):
            ax.plot([r["c1"] - t], [0], "x", color=BAD, ms=12, mew=2.4)
    lim = 36
    ax.plot([-lim, lim], [-lim, lim], "--", color=MUTED, lw=1.4)
    ax.annotate("focus never moved", (lim * 0.60, lim * 0.60), rotation=39, fontsize=9,
                color=INK2, ha="center", va="bottom")
    ax.axhline(0, color=S3, lw=1.8)
    ax.annotate("what success would look like", (-lim * 0.95, 1.8), fontsize=9, color=S3)
    ax.plot([], [], "x", color=BAD, ms=10, mew=2.2, label="run crashed")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("focus error handed to the solver (Å)")
    ax.set_ylabel("focus error after fitting (Å)")
    ax.legend(loc="lower right")
    tag(ax, "a", "Fitting focus inside the solver")

    # right: those moves drawn on the landscape they should be descending. Each arrow starts at the
    # landscape value for the probe the solver was handed and ends where the fit actually finished.
    ax = axes[1]
    curves = {}
    for alpha, col in ((70, S1), (90, S2)):
        pts = sorted([(float(r["dc1"]), float(r["final_error"])) for r in rows
                      if int(r["alpha"]) == alpha and int(r["niter"]) == 200])
        if not pts:
            continue
        xs, ys = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
        curves[alpha] = (xs, ys)
        ax.plot(xs, ys, "-", color=col, lw=1.8, alpha=.5, label=f"α = {alpha}: probe held fixed")
    for alpha, col in ((70, S1), (90, S2)):
        t = truth_c1(alpha)
        for r in (r for r in runs if r["alpha"] == alpha and not r["crashed"]):
            x0, x1, y1 = r["c1"] - t, r["c1"] + r["shift"] - t, r["err"]
            if abs(x1 - x0) < 1.0:
                ax.plot([x1], [y1], "o", color=col, ms=8, mec="white", mew=1.2, zorder=5)
                continue
            y0 = float(np.interp(x0, *curves[alpha])) if alpha in curves else y1
            ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                        arrowprops=dict(arrowstyle="-|>", color=BAD, lw=1.7, shrinkA=3, shrinkB=3), zorder=4)
            ax.plot([x1], [y1], "o", color=BAD, ms=7, zorder=5)
    ax.plot([], [], "o", color=INK2, ms=8, label="where the fit ended")
    ax.annotate("at 70 mrad the fit walks the probe\nuphill, to a far worse reconstruction",
                (0.03, 0.97), xycoords="axes fraction", fontsize=9, color=BAD, va="top")
    ax.set_yscale("log"); ax.set_yticks([5, 10, 20, 50, 100, 200, 400]); ax.minorticks_off()
    ax.set_yticklabels(["5", "10", "20", "50", "100", "200", "400"])
    ax.set_xlim(-40, 40); ax.set_ylim(4, 700)
    ax.set_xlabel("focus error of the probe (Å)")
    ax.set_ylabel("mismatch to the measured data")
    ax.legend(loc="lower center", fontsize=9)
    tag(ax, "b", "Measured against the landscape")
    fig.tight_layout()
    return save(fig, "figA_dfo.png")


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
def figD_nonround(a):
    base = read_reports(a.af)[70]
    legs = [("0", base)]
    waves = {"nr1_C56_0p6w": "0.6", "nr2_C56_1p2w": "1.2", "nr3_C56_2p5w": "2.5"}
    for k, w in waves.items():
        r = af(a.relax, k)
        if r:
            legs.append((w, r))
    combo = af(a.relax, "nr4_C56_C34")
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.3))

    ax = axes[0]
    x = [float(l[0]) for l in legs]
    for sp, col in SPECIES.items():
        ax.plot(x, [100 * l[1][sp]["recall_bulk"] for l in legs], "-o", color=col, label=sp)
    if combo:
        for sp, col in SPECIES.items():
            ax.plot([2.9], [100 * combo[sp]["recall_bulk"]], "s", color=col, mfc="white", mew=1.8, ms=8)
        ax.annotate("+ four-fold\nastigmatism", (2.9, 12), fontsize=8.5, color=INK2, ha="center")
    ax.set_ylim(0, 108); ax.set_xlim(-0.15, 3.25)
    ax.set_xlabel("six-fold astigmatism (waves at the aperture edge)")
    ax.set_ylabel("atoms found (%)"); ax.legend(loc="upper right", ncol=3)
    tag(ax, "a", "Atoms found as the corrector's own residual grows")

    ax = axes[1]
    ax.plot(x, [l[1]["z_rms_A"] for l in legs], "-o", color=S1, label="depth error")
    ax.plot(x, [l[1]["xy_rms_A"] for l in legs], "-o", color=S3, label="in-plane error")
    if combo:
        ax.plot([2.9], [combo["z_rms_A"]], "s", color=S1, mfc="white", mew=1.8, ms=8)
        ax.plot([2.9], [combo["xy_rms_A"]], "s", color=S3, mfc="white", mew=1.8, ms=8)
        ax.annotate("+ four-fold", (2.9, combo["z_rms_A"]), xytext=(0, 9), textcoords="offset points",
                    fontsize=8.5, color=INK2, ha="center")
    ax.set_ylim(0, None); ax.set_xlim(-0.15, 3.25)
    ax.annotate("past 0.6 waves these are the errors of\nthe few atoms still being found",
                (0.5, 0.06), xycoords="axes fraction", fontsize=8.5, color=INK2)
    ax.set_xlabel("six-fold astigmatism (waves at the aperture edge)")
    ax.set_ylabel("error (Å)"); ax.legend(loc="center left")
    tag(ax, "b", "Position error, α = 70 mrad")
    fig.tight_layout()
    return save(fig, "figD_nonround.png")




# --------------------------------------------------- fig E: what non-round aberration looks like
NR_SHOW = ["nr0_round", "nr0p3_C56_0p3w", "nr1_C56_0p6w", "nr2_C56_1p2w", "nr3_C56_2p5w"]


def read_nonround_sweep():
    """The sweep the sims were actually run from: label -> (aberration dict, C1, waves at the edge)."""
    rows = {}
    with open(os.path.join(REPO, "campaign", "nonround_sweep.tsv")) as f:
        rd = csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t")
        for r in rd:
            ab = {"C30": float(r["c3"]), "C50": float(r["c5"])}
            if r["aber_json"] != "-":
                ab = json.loads(r["aber_json"])
            rows[r["label"]] = (ab, float(r["c1"]), ab.get("C56", 0.0) / C56_WAVE)
    return rows


def figE_nonround_probe(a):
    """The six-fold residual, seen three ways: the Ronchigram the operator would judge focus on,
    the wavefront across the aperture, and the probe that reaches the sample."""
    mrf = _load_mod("make_ronchigram_fig", os.path.join(HERE, "make_ronchigram_fig.py"))
    plan = read_nonround_sweep()
    rng = np.random.default_rng(3)
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

        R, hr = mrf.ronchigram(P, 140.0, 70, rng)
        ax = axes[0, i]
        ax.imshow(R, cmap="gray", extent=[-hr, hr, -hr, hr], interpolation="bilinear")
        ax.add_artist(plt.Circle((0, 0), 70, fill=False, color="#ffd24a", lw=1.2, ls=":"))
        ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
        ax.set_title("round probe\n(no six-fold)" if w == 0 else f"six-fold\n{w:.1f} waves",
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
        "E": figE_nonround_probe}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--relax", default=RELAX); ap.add_argument("--af", default=AF_FINAL)
    ap.add_argument("--only", nargs="+", default=list(FIGS))
    a = ap.parse_args()
    style()
    for k in a.only:
        FIGS[k.upper()](a)


if __name__ == "__main__":
    main()
