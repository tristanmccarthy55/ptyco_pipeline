#!/usr/bin/env python
"""
@file c1_objective.py
@brief The C1 (defocus) search objective: final Fourier error vs trial C1, its noise floor, the
       flat-bottom estimate of C1, and an emulated blind search from an offset start.

Reads the trial dirs written by campaign/run_c1_search.sh (recon_c1_a<A>_<mode>_df<C1>[_rN]_n<NITER>_NL<NL>)
from one or more extracted tarballs:
  - analysis/**/*_error_trace.csv   (ptycho/run_synthetic_recon_ML.m) -> objective = error at the last iteration
  - analysis/**/*_layer_stats.csv   -> depth profile of the phase std -> slab position in the recon box
  - 01/probe_initial.json           (sim/make_probe.py) -> the sim's true C1, probe d90
Truth (c1_sim_A) is used only to plot dC1 and to score the estimate, never inside the estimator.

Why a flat bottom and not an argmin. Defocus is degenerate with the object's depth in this model: a
probe shifted in focus by dC1 and an object shifted in depth by dC1 give the same far-field
intensities, and the full-box recon has 4 A of vacuum on each side to absorb that shift. So the
objective is expected to be flat over roughly C1_true +- 4 A and to rise outside. The estimate is the
interval of C1 whose error lies within K x sigma of the minimum, its midpoint the fitted C1 and its
half-width the uncertainty. sigma comes from repeated trials at the same C1 (the GPU engine reseeds its
RNG from the clock, so repeats differ in batch order and layer starts); without repeats it falls back to
the scatter of second differences on the finest grid, and the report says which. Panel (c) tests the
degeneracy directly: inside the bottom the slab's depth centroid should move one-for-one with dC1.

    python analysis/c1_objective.py --root ~/Desktop/c1_grid_a70
    python analysis/c1_objective.py --root ~/Desktop/c1_grid_a70 ~/Desktop/c1_n200_a70 --blind-start -40 -20 20 40
"""
from __future__ import annotations

import argparse
import csv
import datetime
import glob
import json
import os
import re
from collections import defaultdict

import numpy as np

NAME_RE = re.compile(r"recon_(?P<camp>[A-Za-z0-9]+)_a(?P<alpha>\d+)_(?P<mode>lab|Pb|Ti)_df(?P<c1>-?\d+(?:\.\d+)?)"
                     r"(?:_r(?P<rep>\d+))?(?P<ps>_ps\d+(?:x\d+)?)?_n(?P<niter>\d+)_NL(?P<nl>\d+)$")
LAMBDA_A = 0.0196877
# chart ink + series (dataviz reference palette, light surface)
SURF, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]               # categorical slots 1-3 (NITER)
DIVERGE = {-2: "#184f95", -1: "#5598e7", 0: "#52514e", 1: "#e34948", 2: "#9e2f2e"}   # dC1 sign x magnitude


def week_dir(sub):
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    d = os.path.join(repo, "aberration_experiment", sub, datetime.date.today().strftime("%G-W%V"))
    os.makedirs(d, exist_ok=True)
    return d


def read_csv_rows(path):
    with open(path) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def _profile_edges(z, w):
    """(centroid, entrance, exit) of a non-negative depth profile w(z); edges at half maximum."""
    if not np.isfinite(w).all() or w.sum() <= 0:
        return np.nan, np.nan, np.nan
    cen = float((z * w).sum() / w.sum())
    half = 0.5 * w.max(); above = np.where(w >= half)[0]
    def cross(i0, i1):                    # linear crossing of half between layers i0 (below) and i1 (above)
        if i0 < 0 or i0 >= len(w):
            return float(z[i1])
        return float(z[i0] + (half - w[i0]) / (w[i1] - w[i0] + 1e-30) * (z[i1] - z[i0]))
    return cen, cross(above[0] - 1, above[0]), cross(above[-1] + 1, above[-1])


def slab_profile(stats_path):
    """Slab position from the sidecar's RAW per-layer phase std (fallback when no h5 was packed).

    Caution, measured on the a70 smoke test: the raw std is dominated by the phase-ramp gauge (vacuum layer
    2 raw 0.095 vs 0.025 deramped) and the ramp is not spread evenly (layer 1 carries none), so this
    centroid is biased. Prefer h5_slab_profile."""
    rows = read_csv_rows(stats_path)
    z = np.array([float(r["z_A"]) for r in rows]); s = np.array([float(r["phase_std"]) for r in rows])
    return _profile_edges(z, np.clip(s - s.min(), 0, None))


def h5_slab_profile(h5_path, box_z=None, hp_sigma_A=1.0):
    """Slab position from the recon h5: per layer, the std over the illuminated field (illum_sum > 0.15 max)
    of the phase after removing a fitted plane (the ramp gauge) and a Gaussian low-pass of hp_sigma_A
    (keeps the atomic contrast). The profile minus its minimum is the weight. On the a70 smoke test this
    puts the slab at 14.39 A against 14.19 A for the 200-iteration known-probe recon (slab centre 13.76 A)."""
    import h5py
    from scipy.ndimage import gaussian_filter
    with h5py.File(h5_path, "r") as f:
        r = f["reconstruction"]
        obj = np.asarray(r["object"]).squeeze()
        il = np.asarray(r["p"]["illum_sum"]["illum_sum_0"])
        dx = float(np.ravel(r["p"]["dx_spec"])[0]) * 1e10
    if obj.ndim == 2:
        obj = obj[None]
    m = il > 0.15 * il.max()
    if m.shape != obj.shape[1:]:
        m = m.T
    yy, xx = np.mgrid[:obj.shape[1], :obj.shape[2]]
    A = np.c_[xx[m], yy[m], np.ones(m.sum())]
    hp = []
    for l in range(obj.shape[0]):
        v = np.angle(obj[l]); c, *_ = np.linalg.lstsq(A, v[m], rcond=None)
        d = v - (c[0] * xx + c[1] * yy + c[2])
        hp.append((d - gaussian_filter(d, hp_sigma_A / dx))[m].std())
    hp = np.array(hp); nl = len(hp)
    z = (np.arange(nl) + 0.5) * ((box_z or 27.525) / nl)
    return _profile_edges(z, np.clip(hp - hp.min(), 0, None))


def load_trials(roots, camps):
    """Trials under `roots`, for each campaign prefix in `camps` (a string or a list)."""
    trials = []; seen = set()
    camps = [camps] if isinstance(camps, str) else camps
    for root in roots:
        for camp in camps:
            for d in sorted(glob.glob(os.path.join(os.path.expanduser(root), "**", f"recon_{camp}_a*"), recursive=True)):
                if d in seen:
                    continue
                seen.add(d)
                m = NAME_RE.match(os.path.basename(d.rstrip("/")))
                if not m or not os.path.isdir(d):
                    continue
                t = dict(dir=d, alpha=int(m["alpha"]), mode=m["mode"], c1=float(m["c1"]), rep=int(m["rep"] or 1),
                         niter=int(m["niter"]), nl=int(m["nl"]))
                tr = glob.glob(os.path.join(d, "analysis", "**", "*_error_trace.csv"), recursive=True)
                if not tr:
                    print(f"  MISSING trace: {d} (job failed or sidecar skipped -- check its slurm log)")
                    continue
                rows = read_csv_rows(sorted(tr)[-1])
                it = np.array([int(r["iteration"]) for r in rows]); err = np.array([float(r["fourier_error"]) for r in rows])
                t.update(trace_it=it, trace_err=err, final_error=float(err[np.argmax(it)]), last_iter=int(it.max()),
                         n_rows=len(it), complete=bool(it.max() == t["niter"]), finite=bool(np.all(np.isfinite(err))))
                pj = os.path.join(d, "01", "probe_initial.json")
                if os.path.isfile(pj):
                    info = json.load(open(pj))
                    t.update(c1_sim=float(info["c1_sim_A"]), d90=info.get("d90_A"), self_check=info.get("self_check"))
                    if abs(float(info["c1_A"]) - t["c1"]) > 1e-6:
                        print(f"  WARNING {d}: probe_initial.json C1 {info['c1_A']} != dir-name C1 {t['c1']}")
                else:
                    t.update(c1_sim=np.nan, d90=None, self_check=None)
                    print(f"  NOTE {d}: no probe_initial.json -- truth unknown for this trial")
                h5 = glob.glob(os.path.join(d, "analysis", "**", "*_recons.h5"), recursive=True)
                ls = glob.glob(os.path.join(d, "analysis", "**", "*_layer_stats.csv"), recursive=True)
                if h5:
                    t["z_centroid"], t["z_entrance"], t["z_exit"] = h5_slab_profile(sorted(h5)[-1]); t["z_source"] = "h5"
                elif ls:
                    t["z_centroid"], t["z_entrance"], t["z_exit"] = slab_profile(sorted(ls)[-1]); t["z_source"] = "sidecar(raw)"
                else:
                    t["z_centroid"], t["z_entrance"], t["z_exit"] = (np.nan,) * 3; t["z_source"] = None
                trials.append(t)
    return trials


MIN_REPEAT_DOF = 4     # below this a repeat-based sigma is too uncertain to set the bottom on its own


def _sigma_d2(groups):
    """Robust sigma from 2nd differences of single-trial values on the finest evenly spaced run (var = 6 sigma^2).
    Includes any curvature, so it is an upper bound where the objective bends."""
    c1 = np.array(sorted(groups)); f = np.array([groups[c][0] for c in c1])
    if len(c1) < 5:
        return np.nan, 0, None
    h = np.diff(c1); fine = h.min()
    idx = [i for i in range(1, len(c1) - 1) if abs(h[i - 1] - fine) < 1e-6 and abs(h[i] - fine) < 1e-6]
    if len(idx) < 3:
        return np.nan, 0, fine
    d2 = np.array([f[i - 1] - 2 * f[i] + f[i + 1] for i in idx])
    return float(1.4826 * np.median(np.abs(d2 - np.median(d2))) / np.sqrt(6)), len(idx), fine


def noise_sigma(groups):
    """sigma of ONE trial's final error. groups: {c1: [values]}. Pooled repeats when they carry >= MIN_REPEAT_DOF
    degrees of freedom; otherwise the larger of the repeat estimate and the 2nd-difference estimate (one pair of
    repeats gave a sigma 20x too small on synthetic data)."""
    reps = [np.asarray(v) for v in groups.values() if len(v) >= 2]
    dof = sum(len(v) - 1 for v in reps)
    s_rep = float(np.sqrt(sum(((v - v.mean()) ** 2).sum() for v in reps) / dof)) if dof else np.nan
    if dof >= MIN_REPEAT_DOF:
        return s_rep, f"repeats ({dof} dof)"
    s_d2, n_d2, fine = _sigma_d2(groups)
    d2_lab = f"2nd differences ({n_d2} on the {fine:g} A grid; upper bound)" if fine is not None else "2nd differences"
    cands = [(x, lab) for x, lab in ((s_rep, f"repeats ({dof} dof)"), (s_d2, d2_lab)) if np.isfinite(x)]
    if not cands:
        return np.nan, "none (no repeats and no evenly spaced fine run)"
    x, lab = max(cands)
    other = [f"{l} {v:.3g}" for v, l in cands if l != lab]
    return x, lab + (f"; larger than {other[0]}" if other else "") + f" -- fewer than {MIN_REPEAT_DOF} repeat dof"


def _run_below(f, i0, thr):
    lo = hi = i0
    while lo - 1 >= 0 and f[lo - 1] <= thr:
        lo -= 1
    while hi + 1 < len(f) and f[hi + 1] <= thr:
        hi += 1
    return lo, hi


def flat_bottom(c1, f, sigma, k):
    """Contiguous run of grid C1 around the minimum whose error is within K sigma of the bottom's level.

    Two passes, because the minimum of several noisy points is biased low: pass 1 takes f <= min + K sigma sqrt(2)
    (a difference of two noisy values), pass 2 re-thresholds at mean(pass-1 run) + K sigma. C1_fit = the run's
    midpoint; half-width = half the run + half the spacing to the first point outside it (the edge lies between).
    Monte Carlo on the planned grid (repeats at 0 x3 and +-10 x2, truth off-grid, K=2, 400 draws): flat half-width
    4 A -> rms error 1.4-1.9 A, coverage 0.82-0.85; half-width 2 A -> 0.8-1.1 A, 0.92-0.93; a sharp minimum ->
    0.5-0.6 A, 0.98-0.99; the one-pass min + K sigma rule covered only 0.67-0.70 at 4 A.
    @return dict(c1_fit, halfwidth, lo, hi, threshold, n_in_bottom, bottom_open, parabola_vertex, c1_argmin)"""
    o = np.argsort(c1); c1, f = np.asarray(c1, float)[o], np.asarray(f, float)[o]
    i0 = int(np.argmin(f)); ks = k * sigma if np.isfinite(sigma) else 0.0
    lo, hi = _run_below(f, i0, f[i0] + ks * np.sqrt(2))
    thr = float(f[lo:hi + 1].mean() + ks)
    lo, hi = _run_below(f, i0, thr)
    gap_lo = (c1[lo] - c1[lo - 1]) / 2 if lo > 0 else np.nan
    gap_hi = (c1[hi + 1] - c1[hi]) / 2 if hi + 1 < len(c1) else np.nan
    edge = np.nanmax([gap_lo, gap_hi]) if np.isfinite([gap_lo, gap_hi]).any() else 0.0
    fit = 0.5 * (c1[lo] + c1[hi]); hw = 0.5 * (c1[hi] - c1[lo]) + edge
    a, b = max(lo - 1, 0), min(hi + 1, len(c1) - 1)
    vertex = np.nan
    if b - a + 1 >= 3:
        p = np.polyfit(c1[a:b + 1], f[a:b + 1], 2)
        if p[0] > 0:
            vertex = float(-p[1] / (2 * p[0]))
    return dict(c1_fit=float(fit), halfwidth=float(hw), lo=float(c1[lo]), hi=float(c1[hi]), threshold=thr,
                n_in_bottom=int(hi - lo + 1), bottom_open=bool(lo == 0 or hi == len(c1) - 1), parabola_vertex=vertex,
                c1_argmin=float(c1[i0]))


def blind_search(means, start, coarse, fine, sigma, k):
    """Operator emulation on measured points: walk downhill in `coarse` steps from `start` until both neighbours
    are higher, then take the `fine` scan inside that bracket and apply flat_bottom. Returns (result, needed C1s)."""
    have = lambda c: min(means, key=lambda x: abs(x - c)) if means and min(abs(x - c) for x in means) < 1e-6 else None
    needed = []; c = start; path = [start]
    for _ in range(25):
        nb = [c - coarse, c, c + coarse]; miss = [x for x in nb if have(x) is None]
        if miss:
            return None, sorted(set(miss)), path
        fl, fc, fr = (means[have(x)] for x in nb)
        if fc <= fl and fc <= fr:
            break
        c = c - coarse if fl < fr else c + coarse; path.append(c)
    else:
        return None, [], path
    grid = np.arange(c - coarse, c + coarse + fine / 2, fine)
    needed = [float(x) for x in grid if have(x) is None]
    if needed:
        return None, needed, path
    pts = np.array([have(x) for x in grid]); vals = np.array([means[x] for x in pts])
    return flat_bottom(pts, vals, sigma, k), [], path


def rank_agreement(a, b):
    """Spearman rank correlation of two {c1: value} dicts over their common C1."""
    common = sorted(set(a) & set(b))
    if len(common) < 3:
        return np.nan, len(common)
    ra = np.argsort(np.argsort([a[c] for c in common])); rb = np.argsort(np.argsort([b[c] for c in common]))
    return float(np.corrcoef(ra, rb)[0, 1]), len(common)


def golden_proposal(c1, f):
    """Next golden-section point from the lowest measured triple (for a bottom narrower than the grid)."""
    o = np.argsort(c1); c1, f = np.asarray(c1)[o], np.asarray(f)[o]; i = int(np.clip(np.argmin(f), 1, len(f) - 2))
    a, b, c = c1[i - 1], c1[i], c1[i + 1]; g = 0.381966
    return float(b + g * (c - b)) if (c - b) > (b - a) else float(b - g * (b - a))


def make_figure(res, trials, out_png, logy):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter, LogLocator, NullFormatter
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.edgecolor": "#c3c2b7",
                         "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
                         "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.grid": True, "grid.color": GRID,
                         "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
                         "legend.frameon": False})
    keys = sorted(res)
    fig, axes = plt.subplots(len(keys), 3, figsize=(13, 3.9 * len(keys)), squeeze=False)
    for row, key in zip(axes, keys):
        alpha, mode = key; R = res[key]; truth = R["c1_sim"]
        axA, axB, axC = row
        # (a) final error vs dC1, one series per NITER; bottom interval of the primary NITER shaded
        for j, n in enumerate(sorted(R["by_niter"])):
            S = R["by_niter"][n]; col = SERIES[j % len(SERIES)]
            ts = [t for t in trials if (t["alpha"], t["mode"], t["niter"]) == (alpha, mode, n)]
            axA.plot([t["c1"] - truth for t in ts], [t["final_error"] for t in ts], "o", ms=3.5, color=col,
                     alpha=0.55, mec="none")
            c = np.array(sorted(S["means"])); axA.plot(c - truth, [S["means"][x] for x in c], "-", lw=2, color=col,
                                                       label=f"NITER {n}")
            if S.get("estimate") and n == R["primary_niter"]:
                E = S["estimate"]
                axA.axvspan(E["lo"] - truth, E["hi"] - truth, color=col, alpha=0.12, lw=0)
                axA.axhline(E["threshold"], color=col, lw=1, ls=":")
                axA.errorbar([E["c1_fit"] - truth], [min(S["means"].values())], xerr=[E["halfwidth"]], fmt="D",
                             ms=7, color=INK, mfc=SURF, mew=1.5, capsize=3, zorder=5,
                             label=f"fit {E['c1_fit'] - truth:+.1f} ± {E['halfwidth']:.1f} Å")
        axA.axvline(0, color=INK2, lw=1, ls="--")
        axA.annotate("true C1", (0, 0.5), xycoords=("data", "axes fraction"), xytext=(4, 0),
                     textcoords="offset points", color=INK2, fontsize=8, va="center")
        if logy:
            axA.set_yscale("log")
        axA.set_xlabel("trial C1 − true C1 (Å)"); axA.set_ylabel("final Fourier error")
        axA.set_title(f"α = {alpha} mrad ({mode}): objective, true C1 = {truth:g} Å", loc="left", fontsize=10)
        axA.legend(loc="upper center", fontsize=8)
        # (b) traces vs iteration at the largest NITER, diverging colour by dC1
        n = max(R["by_niter"]); ts = [t for t in trials if (t["alpha"], t["mode"], t["niter"], t["rep"]) == (alpha, mode, n, 1)]
        d = sorted({t["c1"] - truth for t in ts})
        pick = sorted({min(d, key=lambda x: abs(x - v)) for v in (-20, -6, 0, 6, 20)})
        neg = [x for x in pick if x < 0]; pos = [x for x in pick if x > 0]
        for t in sorted(ts, key=lambda t: t["c1"]):
            dc = t["c1"] - truth
            if dc not in pick:
                continue
            slot = 0 if dc == 0 else (-(len(neg) - neg.index(dc)) if dc < 0 else pos.index(dc) + 1)
            col = DIVERGE[int(np.clip(slot, -2, 2))]
            axB.plot(t["trace_it"], t["trace_err"], "-", lw=2, color=col, label=f"ΔC1 {dc:+g} Å")
        axB.set_yscale("log"); axB.set_xlabel("iteration (full-resolution engine)"); axB.set_ylabel("Fourier error")
        axB.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
        axB.yaxis.set_major_formatter(FormatStrFormatter("%.3g")); axB.yaxis.set_minor_formatter(NullFormatter())
        axB.set_title(f"convergence, NITER {n}", loc="left", fontsize=10); axB.legend(fontsize=8)
        # (c) slab depth centroid vs dC1 (degeneracy test), with slope-1 references through the dC1 = 0 point
        for j, nn in enumerate(sorted(R["by_niter"])):
            ts = [t for t in trials if (t["alpha"], t["mode"], t["niter"]) == (alpha, mode, nn) and np.isfinite(t["z_centroid"])]
            if not ts:
                continue
            x = np.array([t["c1"] - truth for t in ts]); y = np.array([t["z_centroid"] for t in ts])
            axC.plot(x, y, "o", ms=4, color=SERIES[j % len(SERIES)], mec="none", label=f"NITER {nn}")
        z0 = R.get("z0")
        if z0 is not None and np.isfinite(z0):
            # references only across the bottom (+-4 A beyond), so they cannot stretch the depth axis
            lims = axC.get_xlim(), axC.get_ylim()
            w0, w1 = R["z_window"][0] - truth - 4, R["z_window"][1] - truth + 4; xx = np.array([w0, w1])
            for sgn, ls in ((1, "--"), (-1, ":")):
                axC.plot(xx, z0 + sgn * xx, ls, color=MUTED, lw=1, label=f"slope {sgn:+d}")
            axC.set_xlim(lims[0]); axC.set_ylim(lims[1][0] - 2, lims[1][1] + 2)
        if R.get("z_slope") is not None:
            axC.set_title(f"slab depth centroid: slope {R['z_slope']:+.2f} inside the bottom", loc="left", fontsize=10)
        axC.set_xlabel("trial C1 − true C1 (Å)"); axC.set_ylabel("phase-std centroid (Å, recon box)")
        axC.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    print(f"wrote {out_png}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", nargs="+", required=True, help="extracted tarball dir(s) holding recon_c1_* dirs")
    ap.add_argument("--camp", nargs="+", default=["c1"],
                    help="dir prefixes recon_<camp>_ to combine, e.g. --camp c1 c1fine (default c1)")
    ap.add_argument("--k", type=float, default=2.0, help="bottom = error <= min + K*sigma (default 2)")
    ap.add_argument("--niter", type=int, default=None, help="NITER that defines the estimate (default: smallest with the most C1)")
    ap.add_argument("--blind-start", type=float, nargs="*", default=[], help="start offsets from the true C1 to emulate [A]")
    ap.add_argument("--coarse", type=float, default=10.0); ap.add_argument("--fine", type=float, default=2.0)
    ap.add_argument("--logy", action="store_true", help="log y on the objective panel")
    ap.add_argument("--out-csv", default=None); ap.add_argument("--out-png", default=None); ap.add_argument("--out-json", default=None)
    ap.add_argument("--no-fig", action="store_true")
    a = ap.parse_args()

    trials = load_trials(a.root, a.camp)
    if not trials:
        raise SystemExit(f"no recon_{{{','.join(a.camp)}}}_* trials with an error trace under {a.root}")
    print(f"\n{len(trials)} trials")
    for t in trials:
        flags = [] if t["complete"] and t["finite"] else (["INCOMPLETE"] if not t["complete"] else []) + (["NaN"] if not t["finite"] else [])
        if t.get("self_check") == "FAIL":
            flags.append("PROBE SELF-CHECK FAIL")
        if flags:
            print(f"  {'/'.join(flags)}: {os.path.basename(t['dir'])} (last iter {t['last_iter']}, {t['n_rows']} rows)")
    good = [t for t in trials if t["complete"] and t["finite"] and t.get("self_check") != "FAIL"]

    res = {}
    for key in sorted({(t["alpha"], t["mode"]) for t in good}):
        ts = [t for t in good if (t["alpha"], t["mode"]) == key]
        truths = {t["c1_sim"] for t in ts if np.isfinite(t["c1_sim"])}
        if len(truths) > 1:
            raise SystemExit(f"a{key[0]} {key[1]}: trials disagree on the sim's C1 {truths} -- mixed sims in one analysis")
        R = dict(c1_sim=truths.pop() if truths else np.nan, by_niter={})
        for n in sorted({t["niter"] for t in ts}):
            groups = defaultdict(list)
            for t in ts:
                if t["niter"] == n:
                    groups[t["c1"]].append(t["final_error"])
            means = {c: float(np.mean(v)) for c, v in groups.items()}
            sigma, how = noise_sigma(groups)
            S = dict(means=means, n_c1=len(means), n_trials=sum(len(v) for v in groups.values()), sigma=sigma, sigma_from=how)
            if len(means) >= 3 and np.isfinite(sigma):
                S["estimate"] = flat_bottom(list(means), list(means.values()), sigma, a.k)
            R["by_niter"][n] = S
        ns = sorted(R["by_niter"], key=lambda n: (-R["by_niter"][n]["n_c1"], n))
        R["primary_niter"] = a.niter if a.niter in R["by_niter"] else ns[0]
        if len(R["by_niter"]) >= 2:
            lo, hi = min(R["by_niter"]), max(R["by_niter"])
            R["rank_rho"], R["rank_n"] = rank_agreement(R["by_niter"][lo]["means"], R["by_niter"][hi]["means"])
        # degeneracy test on the primary NITER: centroid vs trial C1 INSIDE the estimated bottom, where the recon can
        # still absorb the focus error as a depth shift (outside it the vacuum band is used up)
        Ep = R["by_niter"][R["primary_niter"]].get("estimate")
        R["z_window"] = (Ep["lo"], Ep["hi"]) if Ep else (R["c1_sim"] - 4, R["c1_sim"] + 4)
        tp = [t for t in ts if t["niter"] == R["primary_niter"] and np.isfinite(t["z_centroid"])
              and R["z_window"][0] - 1e-6 <= t["c1"] <= R["z_window"][1] + 1e-6]
        if len(tp) >= 3:
            x = np.array([t["c1"] - R["c1_sim"] for t in tp]); y = np.array([t["z_centroid"] for t in tp])
            p = np.polyfit(x, y, 1); R["z_slope"], R["z0"] = float(p[0]), float(p[1])
        res[key] = R

    # ---- report ----
    summary = {}
    for (alpha, mode), R in res.items():
        truth = R["c1_sim"]; print(f"\n=== a{alpha} {mode}: true C1 {truth:g} A")
        for n, S in sorted(R["by_niter"].items()):
            print(f"  NITER {n}: {S['n_trials']} trials at {S['n_c1']} C1 values; sigma {S['sigma']:.4g} from {S['sigma_from']}")
            E = S.get("estimate")
            if not E:
                print("    no estimate (needs >= 3 C1 values and a noise sigma)")
                continue
            err = E["c1_fit"] - truth; ok = abs(err) <= E["halfwidth"]
            print(f"    bottom [{E['lo']:g}, {E['hi']:g}] A ({E['n_in_bottom']} points <= bottom level + {a.k:g} sigma"
                  f"{'; OPEN at the grid edge -- widen the grid' if E['bottom_open'] else ''})")
            print(f"    C1_fit {E['c1_fit']:g} +- {E['halfwidth']:g} A  (true {truth:g}: error {err:+g} A -> "
                  f"{'within' if ok else 'OUTSIDE'} the half-width); argmin {E['c1_argmin']:g}, parabola vertex "
                  f"{E['parabola_vertex']:.1f}")
            if E["n_in_bottom"] <= 1:
                c1s = list(S["means"]); print(f"    bottom narrower than the grid: golden-section proposal C1 = "
                                              f"{golden_proposal(c1s, [S['means'][c] for c in c1s]):.1f} A")
            summary[f"a{alpha}_{mode}_n{n}"] = dict(alpha=alpha, mode=mode, niter=n, c1_true=truth, sigma=S["sigma"],
                                                   sigma_from=S["sigma_from"], within=bool(ok), **E)
        if "rank_rho" in R:
            print(f"  NITER {min(R['by_niter'])} vs {max(R['by_niter'])}: Spearman rho {R['rank_rho']:.3f} over "
                  f"{R['rank_n']} common C1 (1 = the short run ranks trials as the long one does)")
        if "z_slope" in R:
            print(f"  slab centroid vs trial C1 over [{R['z_window'][0]:g}, {R['z_window'][1]:g}] A: slope {R['z_slope']:+.3f} "
                  f"(degeneracy predicts magnitude 1)")
        Sp = R["by_niter"][R["primary_niter"]]
        for s0 in a.blind_start:
            est, need, path = blind_search(Sp["means"], truth + s0, a.coarse, a.fine, Sp["sigma"], a.k)
            walk = " -> ".join(f"{p - truth:+g}" for p in path)
            if est is None:
                print(f"  blind start {s0:+g} A (walk {walk}): needs trials at C1 = {', '.join(f'{x:g}' for x in need) or '?'}")
            else:
                err = est["c1_fit"] - truth
                print(f"  blind start {s0:+g} A (walk {walk}): C1_fit {est['c1_fit']:g} +- {est['halfwidth']:g} "
                      f"(error {err:+g} A, {'within' if abs(err) <= est['halfwidth'] else 'OUTSIDE'})")
                summary[f"a{alpha}_{mode}_blind{s0:+g}"] = dict(start_offset=s0, c1_true=truth, **est)
        if a.blind_start:
            print("  (emulation caveat: the grid is centred on the true C1, so a bracket always lands on it; a real "
                  "search from an unknown offset puts its fine scan wherever its bracket falls)")

    rows = sorted(trials, key=lambda t: (t["alpha"], t["mode"], t["niter"], t["c1"], t["rep"]))
    cols = ["alpha", "mode", "niter", "c1", "dc1", "rep", "final_error", "last_iter", "n_rows", "complete", "finite",
            "self_check", "d90", "z_centroid", "z_entrance", "z_exit", "z_source", "dir"]
    out_csv = a.out_csv or os.path.join(week_dir("results"), "c1_objective.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for t in rows:
            w.writerow({c: (t["c1"] - t["c1_sim"] if c == "dc1" else t.get(c)) for c in cols})
    print(f"\nwrote {out_csv}")
    out_json = a.out_json or os.path.join(week_dir("results"), "c1_objective_summary.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(f"wrote {out_json}")
    if not a.no_fig and res:
        make_figure(res, good, a.out_png or os.path.join(week_dir("figs"), "c1_objective.png"), a.logy)


if __name__ == "__main__":
    main()
