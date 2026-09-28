#!/usr/bin/env python
"""@file make_ceos_figs.py
@brief The CEOS-approx sweep figure set: an older hexapole-corrected column opened from 40 to 80 mrad, probe known.

Built up in the order the page reads, every number read from the run outputs or computed from the sweep table --
nothing typed in here:
  fig1_aberration.png   the aberration across the aperture in waves (all of it, then only what the corrector cannot
                        null) and the Ronchigram it makes, per aperture
  fig2_probes.png       the probe that lands on the specimen: intensity and phase, CEOS against the round control
  fig3_growth.png       the non-round terms in waves as the aperture opens, and the probe diameter against the scan
  fig4_scan.png         the change that made the large probes work: the scan field and step, and the fit it gives
  fig5_recons.png       the reconstructions: depth sections with the real atoms (dots) and the found ones (rings)
  fig6_numbers.png      recall per species and depth error against aperture, CEOS against round
  fig7_stability.png    why the round 75/80 controls needed a different probe: the presolve diverging, per slice count

Reuses make_meeting_figs (style, probe, aperture phase), make_ronchigram_fig (one shared amorphous film, d90),
make_simple_figs (depth sections with atoms drawn on) and campaign/aberration_waves (waves per term).

    ~/hyperspy-bundle/bin/python analysis/make_ceos_figs.py [--figs 1 2 3 ...]
"""
from __future__ import annotations

import argparse
import csv
import datetime
import glob
import importlib.util
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)


def _load_mod(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


mmf = _load_mod("make_meeting_figs", os.path.join(HERE, "make_meeting_figs.py"))
aw = _load_mod("aberration_waves", os.path.join(REPO, "campaign", "aberration_waves.py"))
LAM = mmf.LAM
INK, INK2, MUTED, GRIDC = mmf.INK, mmf.INK2, mmf.MUTED, mmf.GRIDC
SPECIES = mmf.SPECIES                       # Pb blue, Ti orange, O green: the colours every earlier figure used
TSV = os.path.join(REPO, "campaign", "ceos_sweep.tsv")
ALPHAS = (40, 50, 60, 65, 70, 75, 80)
SCAN_GOLD = "#c99700"                       # the 20 A scan, in every figure (never a species colour)
ROUND_TERMS = ("C10", "C30", "C50")         # what a round (Cs-type) corrector adjusts and leaves; the rest is non-round
# the round control used at each aperture in the results (the planner's 4 A row, or the stable balance where the
# planner's probe diverged at Nyquist slicing -- fig7)
ROUND_CONTROL = {40: "round_a040", 50: "round_a050", 60: "round_a060", 65: "round_a065", 70: "round_a070",
                 75: "round_a075_b7", 80: "round_a080_b9"}


def week_dir():
    d = os.path.join(REPO, "aberration_experiment", "figs", datetime.date.today().strftime("%G-W%V"), "ceos")
    os.makedirs(d, exist_ok=True)
    return d


def clean(ax):
    """An image panel: no ticks, grid or frame (the shared style keeps left/bottom spines for charts)."""
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def save(fig, name):
    p = os.path.join(week_dir(), name)
    fig.savefig(p)
    plt.close(fig)
    print(f"wrote {p}")
    return p


# ------------------------------------------------------------------------------------------------ the sweep table
def rows():
    out = {}
    with open(TSV) as f:
        for r in csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"):
            out[r["label"]] = r
    return out


def aberrations(r):
    """abTEM aberration dict of a sweep row, exactly as the simulator received it (aber_json, else round c3/c5)."""
    if r["aber_json"] not in ("-", ""):
        return {k: float(v) for k, v in json.loads(r["aber_json"]).items()}
    return {"C30": float(r["c3"]), "C50": float(r["c5"])}


def chi_waves(ab, c1, alpha, n=301, keep=1.0, only_nonround=False, part=None):
    """The aberration function across the aperture in WAVES, from abTEM's polar convention:
    chi = 2 pi / lambda * sum C_nm theta^(n+1) / (n+1) cos(m (phi - phi_nm)), defocus = -C10.
    Returned in abTEM's array orientation (first axis = x, as every probe and Ronchigram panel here is drawn), in
    mrad, NaN outside the aperture. Checked against abTEM's own probe by check_chi (the axis order was the catch:
    the round terms agree either way, the non-round ones only with x first)."""
    t = np.linspace(-keep * alpha, keep * alpha, n) * 1e-3
    tx, ty = np.meshgrid(t, t)
    th, ph = np.hypot(tx, ty), np.arctan2(ty, tx)
    terms = dict(ab); terms["C10"] = -float(c1)
    w = np.zeros_like(th)
    for k, v in terms.items():
        part_ = part or ("nonround" if only_nonround else "all")
        if not k.startswith("C") or (part_ == "nonround" and k in ROUND_TERMS) or \
                (part_ == "round" and k not in ROUND_TERMS):
            continue
        nn, m = int(k[1]), int(k[2])
        w += float(v) * th ** (nn + 1) / (nn + 1) * np.cos(m * (ph - float(terms.get("phi" + k[1:], 0.0))))
    w = w / LAM                                  # in waves (chi / 2 pi)
    w[th > alpha * 1e-3] = np.nan
    return w.T, t[-1] * 1e3


def check_chi(ab, c1, alpha, P, L):
    """Median disagreement (rad, modulo 2 pi, after removing a constant) between the analytic chi and the phase of
    the probe abTEM built. Guards the sign, angle and axis conventions the figure relies on (median: the read-back
    phase is itself coarsely sampled where chi is steep, at the edge of the 80 mrad aperture)."""
    ph, hp = mmf.aperture_phase(P, alpha, L, keep=1.0)
    n = ph.shape[0]
    w, _ = chi_waves(ab, c1, alpha, n=n, keep=hp / alpha)
    d = np.angle(np.exp(1j * (ph + 2 * np.pi * w)))          # probe = A exp(-i chi)
    ok = np.isfinite(d)
    d = np.angle(np.exp(1j * (d - np.angle(np.nanmean(np.exp(1j * d[ok]))))))
    return float(np.median(np.abs(d[ok])))


# ------------------------------------------------------------------------------------------------ fig 1
def wrapped_panel(ax, w, h, al):
    """A wavefront in waves, wrapped: one full colour cycle per wave, so the fringes crowd where the phase changes
    fast -- the same places the Ronchigram turns to speckle."""
    ax.imshow(np.mod(w, 1.0), cmap="twilight", vmin=0, vmax=1, extent=[-h, h, -h, h], interpolation="nearest")
    ax.add_artist(plt.Circle((0, 0), al, fill=False, color=INK2, lw=0.8))
    clean(ax)


def fig1_aberration(a):
    import make_ronchigram_fig as mrf           # noqa: E402  (one shared film for every Ronchigram)
    R = rows()
    fig, axes = plt.subplots(2, len(ALPHAS), figsize=(2.35 * len(ALPHAS), 5.4))
    for i, al in enumerate(ALPHAS):
        r = R[f"ceosopt_a{al:03d}"]; ab = aberrations(r); c1 = float(r["c1"])
        L, N = 160.0, 2560
        P = mmf.build_probe(al, ab, c1, L=L, N=N)
        err = check_chi(ab, c1, al, P, L)
        if err > 0.15:
            raise SystemExit(f"analytic chi disagrees with abTEM's probe by {err:.2f} rad at {al} mrad")
        w, h = chi_waves(ab, c1, al, n=601)
        wrapped_panel(axes[0, i], w, h, al)
        axes[0, i].set_title(f"{al} mrad\n{np.nanmax(w) - np.nanmin(w):.0f} waves peak to valley", fontsize=10.5, pad=5)
        Rg, hr = mrf.ronchigram(P, L, al)
        ax = axes[1, i]
        ax.imshow(Rg, cmap="gray", extent=[-hr, hr, -hr, hr], interpolation="bilinear")
        ax.add_artist(plt.Circle((0, 0), al, fill=False, color="#ffd24a", lw=1.0, ls=":"))
        clean(ax)
    axes[0, 0].set_ylabel("every aberration\n(one colour cycle = 1 wave)", fontsize=9.5, color=INK2)
    axes[1, 0].set_ylabel("Ronchigram", fontsize=9.5, color=INK2)
    fig.suptitle("The CEOS-approx column opened from 40 to 80 mrad: the whole wavefront, and the Ronchigram it makes",
                 x=0.01, ha="left", fontsize=12, fontweight="semibold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    return save(fig, "fig1_aberration.png")


# ------------------------------------------------------------------------------------------------ fig 2
def complex_rgb(P, gamma=0.45):
    """A complex probe as colour: hue = phase, brightness = amplitude^gamma (normalised). The usual picture of a
    complex wave -- the phase is only drawn where the probe has intensity to carry it."""
    from matplotlib.colors import hsv_to_rgb
    A = np.abs(P); A = (A / A.max()) ** gamma
    H = (np.angle(P) + np.pi) / (2 * np.pi)
    return hsv_to_rgb(np.dstack([H, np.full_like(H, 0.85), A]))


def window_for(al, P, L):
    """The reconstruction window the planner's rule gives this probe (plan_probe.region_geometry, from d90/d99)."""
    import make_ronchigram_fig as mrf
    sys.path.insert(0, os.path.join(REPO, "campaign"))
    from plan_probe import region_geometry
    d90, d99 = mrf.enclosed(P, L, (0.9, 0.99))
    return region_geometry(d90, d99)[1], d90, d99


def probe_panel(ax, img, half, text, cmap=None):
    ax.imshow(img, cmap=cmap, extent=[-half, half, -half, half], interpolation="bilinear")
    clean(ax)
    bar = 0.2 * 2 * half
    ax.plot([-half * 0.9, -half * 0.9 + bar], [-half * 0.85] * 2, color="white", lw=2.2, solid_capstyle="butt")
    ax.text(-half * 0.9, -half * 0.8, f"{bar:g} Å", color="white", fontsize=7.5, va="bottom")
    ax.text(half * 0.93, half * 0.9, text, color="white", fontsize=8, ha="right", va="top")


CORE_A = 6.0                                     # half-width of the probe-core phase panels (A)


def fig2_probes(a):
    import make_ronchigram_fig as mrf
    R = rows()
    L, N = 220.0, 3584                               # holds the 80 mrad probe's d99 (76 A) with room
    px = L / N; c = N // 2
    fig, axes = plt.subplots(3, len(ALPHAS), figsize=(2.35 * len(ALPHAS), 7.7))
    for i, al in enumerate(ALPHAS):
        r = R[f"ceosopt_a{al:03d}"]
        P = mmf.build_probe(al, aberrations(r), float(r["c1"]), L=L, N=N)
        win, d90, d99 = window_for(al, P, L)
        half = win / 2; h = int(round(half / px))
        Pc = P[c - h:c + h, c - h:c + h]
        cy, cx = np.unravel_index((np.abs(P) ** 2).argmax(), P.shape)
        probe_panel(axes[0, i], (np.abs(Pc) ** 2) ** 0.4, half, f"d90 {d90:.1f} Å", "magma")
        axes[0, i].add_artist(plt.Circle(((cx - c) * px, (cy - c) * px), d90 / 2, fill=False, color="white",
                                         lw=0.9, ls=(0, (3, 2))))
        axes[0, i].set_title(f"{al} mrad\nwindow {win:g} Å", fontsize=10.5, pad=5)
        # the core only, at one fixed zoom: across the halo the phase wraps faster than a panel can show
        zh = int(round(CORE_A / px)); cy, cx = np.unravel_index((np.abs(P) ** 2).argmax(), P.shape)
        probe_panel(axes[1, i], complex_rgb(P[cy - zh:cy + zh, cx - zh:cx + zh]), CORE_A, f"core ±{CORE_A:g} Å")
        rr = R[ROUND_CONTROL[al]]
        Pr = mmf.build_probe(al, aberrations(rr), float(rr["c1"]), L=L, N=N)
        d90r = mrf.enclosed(Pr, L, (0.9,))[0]
        probe_panel(axes[2, i], (np.abs(Pr[c - h:c + h, c - h:c + h]) ** 2) ** 0.4, half, f"d90 {d90r:.1f} Å", "magma")
        cy, cx = np.unravel_index((np.abs(Pr) ** 2).argmax(), Pr.shape)
        axes[2, i].add_artist(plt.Circle(((cx - c) * px, (cy - c) * px), d90r / 2, fill=False, color="white",
                                         lw=0.9, ls=(0, (3, 2))))
    for row, lab in enumerate(("CEOS probe\nintensity", "CEOS probe core\ncolour = phase", "round control\nintensity")):
        axes[row, 0].set_ylabel(lab, fontsize=9.5, color=INK2)
    fig.suptitle("The probe on the specimen, each drawn in the reconstruction window it needs",
                 x=0.01, ha="left", fontsize=12, fontweight="semibold")
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    return save(fig, "fig2_probes.png")


# ------------------------------------------------------------------------------------------------ fig 3
def fig3_growth(a):
    import make_ronchigram_fig as mrf
    R = rows()
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(13.0, 5.0), gridspec_kw=dict(width_ratios=[1.25, 1]))
    terms = {}
    for al in ALPHAS:
        for t, w in aw.tableau_waves(aberrations(R[f"ceosopt_a{al:03d}"]), al).items():
            terms.setdefault(t, []).append((al, abs(w)))
    ends = []
    for t, pts in sorted(terms.items()):
        x, y = np.array(pts).T
        hw = t in aw.CEOS_HARDWARE
        col = INK if hw else MUTED
        lw = 2.4 if t == "C56" else (1.6 if hw else 1.1)
        ax.plot(x, y, "-o", color=col, lw=lw, ms=3.5 if not hw else 4.5, alpha=1 if hw else 0.8)
        ends.append([np.log10(y[-1]), f"{t} {aw.NAME.get(t, t)}", col])
    ends.sort(key=lambda e: e[0]); gap = 0.12
    for k in range(1, len(ends)):
        if ends[k][0] - ends[k - 1][0] < gap:
            ends[k][0] = ends[k - 1][0] + gap
    for y, lab, col in ends:
        ax.annotate(lab, (ALPHAS[-1], 10 ** y), xytext=(6, 0), textcoords="offset points", fontsize=8,
                    color=col, va="center")
    ax.set_yscale("log"); ax.set_xlim(ALPHAS[0] - 2, ALPHAS[-1] + 22); ax.set_xticks(ALPHAS)
    ax.set_yticks([0.01, 0.1, 1, 10]); ax.set_yticklabels(["0.01", "0.1", "1", "10"])
    ax.legend(handles=[Line2D([], [], color=INK, lw=2.4, label="six-fold, 1 mm: no knob (CEOS)"),
                       Line2D([], [], color=INK, lw=1.6, label="other hardware terms: no knob"),
                       Line2D([], [], color=MUTED, lw=1.1, label="retuned by the operator at each aperture")],
              loc="lower right", fontsize=8.5)
    ax.set_xlabel("aperture semi-angle (mrad)"); ax.set_ylabel("waves at the aperture edge")
    ax.set_title("a   What is left at each aperture after the operator retunes", loc="left", pad=8)
    # b: probe diameter against the scan fields
    for fam, col, ls, lab in (("ceosopt_a{:03d}", INK, "-", "CEOS probe"),
                              (None, MUTED, "--", "round control")):
        d = []
        for al in ALPHAS:
            rr = R[fam.format(al)] if fam else R[ROUND_CONTROL[al]]
            P = mmf.build_probe(al, aberrations(rr), float(rr["c1"]), L=220.0, N=3584)
            d.append(mrf.enclosed(P, 220.0, (0.9,))[0])
        bx.plot(ALPHAS, d, ls=ls, marker="o", color=col, label=lab)
    new = float(R["ceosopt_a070_f20"]["win"])
    bx.axhline(new, color=SCAN_GOLD, lw=1.6, ls=(0, (4, 3)), label=f"scanned field, {new:g} Å")
    bx.set_yscale("log"); bx.set_xlabel("aperture semi-angle (mrad)"); bx.set_ylabel("length (Å)")
    bx.set_yticks([4, 10, 20, 50, 100]); bx.set_yticklabels(["4", "10", "20", "50", "100"]); bx.set_xticks(ALPHAS)
    bx.minorticks_off()
    bx.set_title("b   Probe diameter against the scanned field", loc="left", pad=8)
    bx.legend(loc="upper left")
    fig.tight_layout()
    return save(fig, "fig3_growth.png")


# ------------------------------------------------------------------------------------------------ results
D = os.path.expanduser("~/Desktop")
AN_1712 = f"{D}/ceos80_sweeps/analysis_local_fixed"          # first sweep, re-scored locally from the h5s (09-28)
AN_W35 = glob.glob(f"{D}/ceos_0928/analysis_*/analysis_*")   # the 35 A rerun + CEOS 75, Blythe, fixed extraction
AN_NEW = glob.glob(f"{D}/ceos_0928c/analysis_ceosopt_a070_f20-ceosopt_a070_nl10*/analysis_*")   # 09-28 legs
# which run stands for each point, and where its numbers live -- paths only, every number is read from the file
RESULTS = [
    ("CEOS", 50, "ceosopt_a050", f"{AN_1712}/ceos"), ("CEOS", 60, "ceosopt_a060", f"{AN_1712}/ceos6x"),
    ("CEOS", 65, "ceosopt_a065", f"{AN_1712}/ceos6x"), ("CEOS", 70, "ceosopt_a070_f20", AN_NEW),
    ("CEOS", 75, "ceosopt_a075", AN_NEW), ("CEOS", 80, "ceosopt_a080_f20", AN_NEW),
    ("round", 50, "round_a050", f"{AN_1712}/round"), ("round", 60, "round_a060", f"{AN_1712}/round"),
    ("round", 65, "round_a065", AN_W35), ("round", 70, "round_a070", f"{AN_1712}/round"),
    ("round", 75, "round_a075_b8", AN_NEW), ("round", 80, "round_a080_b9", AN_NEW),
]
# the first sweep's large-probe legs on their original scans, which did not reconstruct (triage: SATURATED)
FAILED_FIRST = [("CEOS", 70, "ceosopt_a070", f"{D}/ceos_phase_ceos/analysis_*/analysis_*"),
                ("CEOS", 80, "ceosopt_a080", f"{D}/ceos_phase_ceos/analysis_*/analysis_*")]
CONF_MAX = 0.05                                             # atomfind's own health threshold on species confusion


def summary_row(label, where):
    ws = where if isinstance(where, list) else sorted(glob.glob(where))
    for w in ws:
        p = os.path.join(w, "summary.csv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p)):
                if r["label"] == label:
                    return r, p
    raise SystemExit(f"no summary row for {label} under {where}")


def results_table():
    """One row per plotted point, with the file it came from; also written as a CSV beside the figures."""
    out = []
    for inst, al, lab, where in RESULTS:
        r, src = summary_row(lab, where)
        tsv = rows()[lab]
        out.append(dict(instrument=inst, alpha=al, label=lab, status=r["status"],
                        Pb=float(r["Pb_recall_bulk"]), Ti=float(r["Ti_recall_bulk"]), O=float(r["O_recall_bulk"]),
                        z_rms=float(r["z_rms"]), precision=float(r["precision"]), confusion=float(r["confusion"]),
                        NL=int(r["NL"]), scan=f"{tsv['win']} A at {tsv['step']} A", source=src))
    for inst, al, lab, where in FAILED_FIRST:
        r, src = summary_row(lab, where)
        tsv = rows()[lab]
        out.append(dict(instrument=inst, alpha=al, label=lab, status=r["status"], Pb=np.nan, Ti=np.nan, O=np.nan,
                        z_rms=np.nan, precision=np.nan, confusion=np.nan, NL=0,
                        scan=f"{tsv['win']} A at {tsv['step']} A", source=src))
    p = os.path.join(REPO, "aberration_experiment", "results", datetime.date.today().strftime("%G-W%V"))
    os.makedirs(p, exist_ok=True)
    with open(os.path.join(p, "ceos_sweep_results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    results_fragment(out, os.path.join(p, "ceos_sweep_results_table.html"))
    return out


def results_fragment(T, path):
    """The results table as an HTML fragment for the page (page/build_page.py inserts it), from the same rows."""
    import html as H
    def pct(v):
        return "–" if not np.isfinite(v) else f"{100 * v:.0f}"
    rows_ = []
    for t in sorted(T, key=lambda t: (t["alpha"], t["instrument"] != "CEOS", t["status"] != "ok")):
        ok = t["status"] == "ok"
        conf = t["confusion"]
        flag = ('<span class="chip c-stop">did not reconstruct</span>' if not ok else
                ('<span class="chip c-warn">labels unreliable</span>' if conf > CONF_MAX else ""))
        rows_.append("<tr>" + "".join(f"<td>{c}</td>" for c in (
            f"{t['alpha']}", t["instrument"], f"<code>{H.escape(t['label'])}</code>", H.escape(t["scan"]).replace(" A", " Å"),
            f"{t['NL']}" if ok else "–", f"{pct(t['Pb'])} / {pct(t['Ti'])} / {pct(t['O'])}" if ok else "–",
            f"{t['z_rms']:.2f}" if ok else "–", f"{t['precision']:.2f}" if ok else "–",
            f"{100 * conf:.0f} %" if ok else "–", flag)) + "</tr>")
    head = "".join(f"<th>{h}</th>" for h in ("mrad", "column", "run", "scan", "slices", "Pb / Ti / O %",
                                             "depth err. Å", "precision", "confusion", ""))
    open(path, "w").write(f'<table>\n<thead><tr>{head}</tr></thead>\n<tbody>\n' + "\n".join(rows_) +
                          "\n</tbody>\n</table>\n")


def fig6_numbers(a):
    T = results_table()
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 4.3), sharex=True)
    for k, (sp, ax) in enumerate(zip(("Pb", "Ti", "O", "z_rms"), axes)):
        for inst, col, ls, mk in (("CEOS", SPECIES.get(sp, INK), "-", "o"), ("round", MUTED, "--", "s")):
            pts = sorted((t for t in T if t["instrument"] == inst and t["status"] == "ok"), key=lambda t: t["alpha"])
            x = [t["alpha"] for t in pts]
            y = [100 * t[sp] if sp != "z_rms" else t[sp] for t in pts]
            ax.plot(x, y, ls=ls, color=col, lw=2.0 if inst == "CEOS" else 1.6, zorder=2)
            for t, yy in zip(pts, y):
                bad = t["confusion"] > CONF_MAX
                ax.plot(t["alpha"], yy, marker=mk, ms=7.5, mew=1.6, mec=col, mfc="white" if bad else col, zorder=3)
        if sp == "z_rms":
            ax.set_ylabel("depth error, z-RMS (Å)"); ax.set_title("d   depth error", loc="left")
            ax.set_ylim(0, 1.25)
        else:
            ax.set_ylim(0, 105); ax.set_ylabel("bulk recall (%)" if k == 0 else "")
            ax.set_title(f"{'abc'[k]}   {sp}", loc="left")
        ax.set_xlabel("aperture semi-angle (mrad)"); ax.set_xticks(sorted({t['alpha'] for t in T}))
    h = [Line2D([], [], color=INK, ls="-", marker="o", ms=7, label="CEOS column (species colour)"),
         Line2D([], [], color=MUTED, ls="--", marker="s", ms=7, label="round control"),
         Line2D([], [], color=INK2, ls="none", marker="o", mfc="white", ms=7, label="species labels unreliable (confusion > 5 %)")]
    fig.legend(handles=h, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    return save(fig, "fig6_numbers.png")


# ------------------------------------------------------------------------------------------------ fig 7
def presolve_trace(slurm_glob):
    """Presolve (engine 1) error trace from a recon log, as eps = residual x N / 2e5 with N its pattern width."""
    import re
    for f in sorted(glob.glob(slurm_glob), key=os.path.getmtime, reverse=True):
        t = open(f, errors="ignore").read()
        m = re.search(r"engine 1 error trace \(iteration:error\): ([0-9: .e+\-]+)", t)
        n = re.search(r"Ndp(\d+)_step01", t)
        if m and n:
            pts = [tuple(map(float, p.split(":"))) for p in m.group(1).split()]
            it, v = np.array(pts).T
            return it, v * int(n.group(1)) / 2e5
    return None


def full_trace(recon_dir):
    """Full-engine error trace (sidecar CSV), as eps with N from the file name."""
    import re
    f = glob.glob(os.path.join(recon_dir, "analysis", "*", "*", "*_error_trace.csv"))
    if not f:
        return None
    N = int(re.search(r"_(\d+)x\d+_", os.path.basename(f[0])).group(1))
    d = np.array([l.strip().split(",") for l in open(f[0]) if l[0].isdigit()], float)
    return d[:, 0], d[:, 1] * N / 2e5


def change_per_slice(al, ab, c1, nl, box=27.525, ext=48.0, n=768):
    """How much the probe's intensity pattern changes from one recon slice to the next, through the box:
    1 - mean Pearson r of |P(z)|^2 against |P(z + dz)|^2, dz = box / NL (free-space propagation)."""
    P = mmf.build_probe(al, ab, c1, L=ext, N=n)
    k = np.fft.fftfreq(n, ext / n); k2 = k[:, None] ** 2 + k[None, :] ** 2
    F = np.fft.fft2(P); dz = box / nl; cors = []
    for z in np.arange(0, box - dz + 1e-6, 0.5):
        u = np.abs(np.fft.ifft2(F * np.exp(-1j * np.pi * LAM * z * k2))) ** 2
        v = np.abs(np.fft.ifft2(F * np.exp(-1j * np.pi * LAM * (z + dz) * k2))) ** 2
        cors.append(np.corrcoef(u.ravel(), v.ravel())[0, 1])
    return 1 - float(np.mean(cors))


HIA = f"{D}/ceos_hia/atomfind_results_hia_a080_nl14-hia_a080_nl23-hia_a080_k90_20260927_203444"
B1 = glob.glob(f"{D}/ceos_0928c/atomfind_results_ceosopt_a070_f20-*")
NL10 = glob.glob(f"{D}/ceos_0928c/atomfind_results_ceosopt_a070_nl10_*")
# the compact-probe legs with a known outcome: (label, NL, run dir or None, outcome); the trace is read if logged
STAB = [("round_a070", 14, None, True), ("round_a070", 28, None, False), ("round_a075", 16, None, False),
        ("round_a080", 14, f"{HIA}/recon_af_hia_a080_nl14_lab_NL14", True),
        ("round_a080", 18, None, False),
        ("round_a080", 23, f"{HIA}/recon_af_hia_a080_nl23_lab_NL23", False),
        ("hia_a080_k90", 18, f"{HIA}/recon_af_hia_a080_k90_lab_NL18", True),
        ("round_a075_b8", 16, "B1:recon_af_round_a075_b8_lab_NL16", False),
        ("round_a080_b9", 18, "B1:recon_af_round_a080_b9_lab_NL18", True)]


def _rd(d):
    return os.path.join(B1[0], d[3:]) if d and d.startswith("B1:") else d


def probe_name(r):
    return f"{r['alpha']} mrad: C3 {float(r['c3']) / 1e4:+.0f} µm, C1 {float(r['c1']):+.0f} Å"


PENDING_STAB = [("round_a075_b7", 16)]            # predicted from the rule BEFORE running (drawn as a prediction)
BAD = INK2                               # failed: grey, dashed / hollow (orange is titanium's colour)


def fig7_stability(a):
    R = rows()
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(16.0, 5.2), gridspec_kw=dict(width_ratios=[1.0, 0.9]))
    seen = set()
    for lab, nl, d, ok in STAB:
        if not d:
            continue
        tr = presolve_trace(os.path.join(_rd(d), "slurm_*.out"))
        if tr is None:
            continue
        it, e = tr
        r = R[lab]
        ax.plot(it, e, color=INK if ok else BAD, lw=1.8, ls="-" if ok else (0, (4, 2)))
        key = (r["alpha"], r["c3"], r["c1"], nl)
        if key in seen:                                     # the same probe run twice: one label
            continue
        seen.add(key)
        ax.annotate(f"{r['alpha']} mrad, C3/C1 {float(r['c3']) / 1e4:+.0f} µm/{float(r['c1']):+.0f} Å, {nl} sl.",
                    (it[-1], e[-1]), xytext=(5, 0), textcoords="offset points",
                    fontsize=8, va="center", color=INK if ok else BAD)
    ax.set_yscale("log"); ax.set_xlim(0, 380)
    ax.set_yticks([0.05, 0.1, 0.2, 0.5]); ax.set_yticklabels(["0.05", "0.1", "0.2", "0.5"]); ax.minorticks_off()
    ax.set_xlabel("presolve iteration"); ax.set_ylabel("misfit to the data, ε")
    ax.set_title("a   The presolve: settling, or running away", loc="left")
    ax.legend(handles=[Line2D([], [], color=INK, lw=1.8, label="final reconstruction good"),
                       Line2D([], [], color=BAD, lw=1.8, ls=(0, (4, 2)), label="final reconstruction failed")],
              loc="upper left", bbox_to_anchor=(0.1, 0.8))
    pts, labels = [], []
    for lab, nl, d, ok in STAB:
        r = R[lab]
        pts.append((100 * change_per_slice(float(r["alpha"]), aberrations(r), float(r["c1"]), nl), ok))
        labels.append(f"{probe_name(r)}, {nl} slices")
    for lab, nl in PENDING_STAB:
        r = R[lab]
        pts.append((100 * change_per_slice(float(r["alpha"]), aberrations(r), float(r["c1"]), nl), None))
        labels.append(f"{probe_name(r)}, {nl} slices")
    lo = max(m for m, ok in pts if ok is False); hi = min(m for m, ok in pts if ok is True)
    bx.axvspan(lo, hi, color="#e1e0d9", zorder=0)
    for i, (m, ok) in enumerate(pts):
        if ok is None:
            bx.plot(m, i, marker="D", ms=8, mew=1.6, mec=INK2, mfc="white")
        else:
            bx.plot(m, i, marker="o", ms=8, mew=1.6, mec=INK if ok else BAD, mfc=INK if ok else "white")
    bx.set_yticks(range(len(labels))); bx.set_yticklabels(labels, fontsize=8.5); bx.invert_yaxis()
    bx.set_ylim(len(labels) - 0.5, -0.5)
    bx.set_xlabel("probe intensity change per recon slice (%)")
    bx.set_title("b   A candidate rule: the probe must change enough per slice", loc="left")
    bx.legend(handles=[Line2D([], [], color=INK, ls="none", marker="o", ms=8, label="good"),
                       Line2D([], [], color=BAD, ls="none", marker="o", mfc="white", ms=8, label="failed"),
                       Line2D([], [], color=INK2, ls="none", marker="D", mfc="white", ms=8, label="predicted, running"),
                       plt.Rectangle((0, 0), 1, 1, color="#e1e0d9", label="threshold lies in here")],
              loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.subplots_adjust(left=0.05, right=0.9, wspace=0.62, bottom=0.12, top=0.92)
    return save(fig, "fig7_stability.png")


# ------------------------------------------------------------------------------------------------ legs as images
GT_DIR = f"{D}/ceos_region_gt"                     # make_gt_cache --thin-cells 5 --z-vacuum 4 --region-side 210
AN_FD = sorted(glob.glob(f"{D}/ceos_figdata/analysis_*/analysis_*"))   # Blythe figure data (figdata.py), when pulled
FIRST_ROUND = glob.glob(f"{D}/ceos80_sweeps/atomfind_results_round_a040-*/")
FIRST_CEOS = glob.glob(f"{D}/ceos80_sweeps/atomfind_results_ceosopt_a040-*/")
FIRST_HEAVY = glob.glob(f"{D}/ceos80_sweeps/atomfind_results_ceosopt_a070-*/")
# where each leg's volume and found atoms live locally: (root holding recon_af_<label>_lab_NL*, atomfind dir root)
LOCAL = {"round_a050": (FIRST_ROUND, f"{AN_1712}/round"), "round_a060": (FIRST_ROUND, f"{AN_1712}/round"),
         "round_a070": (FIRST_ROUND, f"{AN_1712}/round"), "ceosopt_a050": (FIRST_CEOS, f"{AN_1712}/ceos"),
         "ceosopt_a060": (FIRST_CEOS, f"{AN_1712}/ceos6x"), "ceosopt_a065": (FIRST_CEOS, f"{AN_1712}/ceos6x")}


def leg(label):
    """A leg as the finder saw it: volume, ground truth and found atoms in one frame (as make_simple_figs.load_leg),
    from a local h5 when there is one, else from the Blythe figure data. None if neither is here yet."""
    import figdata
    from atomfind import align
    from analyse_sweep import geometry, recon_h5
    V = h5 = af = None
    if label in LOCAL:
        roots, af_root = LOCAL[label]
        d = glob.glob(os.path.join(roots[0], f"recon_af_{label}_lab_NL*"))
        if d:
            h5 = recon_h5(d[0]); g = geometry(d[0], h5, 4.0)
            V, dx, _ = figdata.frame(h5, GT_DIR, g["dz_A"], g["dx_A"], g["scan_centre"])
            dz, centre = g["dz_A"], g["scan_centre"]
            af = os.path.join(af_root, f"atomfind_{label}")
    if V is None:
        for an in AN_FD:
            fd = os.path.join(an, "figdata", label)
            if os.path.exists(os.path.join(fd, "phase_vol.npy")):
                V, meta = figdata.load(fd)
                dx, dz, centre = meta["dx"], meta["dz"], meta["centre"]
                h5 = os.path.join(fd, "phase_vol.npy")          # config only needs an existing path
                af = os.path.join(an, f"atomfind_{label}")
                break
    if V is None:
        return None
    cfg = figdata.config_for(h5, GT_DIR, dz, dx, centre)
    pos, Z = align.load_gt(cfg)
    found = gtidx = None
    if af and os.path.exists(os.path.join(af, "found_atoms.npy")):
        found = np.load(os.path.join(af, "found_atoms.npy"))
        al = align.refine_with_atoms(align.register(V, dx, pos, Z, cfg), found, pos, Z, cfg)
        gtidx = al.site_to_index(pos[:, 0], pos[:, 1], pos[:, 2])
    return dict(V=V, dx=dx, dz=dz, pos=pos, Z=Z, found=found, gt=gtidx, cfg=cfg)


def insitu(ax, lg, half=10.0):
    """The depth-summed phase over the central field."""
    V, dx = lg["V"], lg["dx"]
    s = V.sum(0); ny, nx = s.shape
    ax.imshow(s, cmap="magma", extent=[-nx * dx / 2, nx * dx / 2, ny * dx / 2, -ny * dx / 2],
              vmin=np.percentile(s, 1), vmax=np.percentile(s, 99.7), interpolation="bilinear")
    ax.set_xlim(-half, half); ax.set_ylim(half, -half)
    clean(ax)


def section_row_A(base):
    """The physical row (A from the crop's top) of the B-site column the earlier figures cut through."""
    mep = _load_mod("mep", os.path.join(HERE, "make_mep_volumes_fig.py"))
    r = mep.pick_row(base["found"], 22, base["V"].shape[1], base["dx"]) if base["found"] is not None else None
    return (base["V"].shape[1] / 2 if r is None else r) * base["dx"]


def placeholder(ax, text):
    clean(ax); ax.set_facecolor("#f4f3ee")
    ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=9, color=INK2, transform=ax.transAxes, wrap=True)


# ------------------------------------------------------------------------------------------------ fig 4
def fig4_scan(a):
    R = rows()
    fig = plt.figure(figsize=(17.0, 9.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1], hspace=0.34, wspace=0.34)
    L, N = 220.0, 3584; px = L / N; c = N // 2
    for j, al in enumerate((70, 80)):
        old, new = R[f"ceosopt_a{al:03d}"], R[f"ceosopt_a{al:03d}_f20"]
        P = mmf.build_probe(al, aberrations(old), float(old["c1"]), L=L, N=N)
        win = window_for(al, P, L)[0]; half = win / 2; h = int(round(half / px))
        logI = np.log10(np.abs(P) ** 2 / (np.abs(P) ** 2).max() + 1e-7)
        ax = fig.add_subplot(gs[0, j])
        ax.imshow(logI[c - h:c + h, c - h:c + h], cmap="magma", vmin=-6, vmax=0,
                  extent=[-half, half, -half, half], interpolation="bilinear")
        for rr, col in ((old, "#bdbcb4"), (new, "#ffd24a")):
            w = float(rr["win"])
            ax.add_patch(plt.Rectangle((-w / 2, -w / 2), w, w, fill=False, color=col, lw=1.4))
        clean(ax)
        ax.set_title(f"{'ab'[j]}   {al} mrad: the probe in its {win:g} Å window", loc="left", fontsize=10.5)
        ax.set_xlabel(f"grey: first sweep, {old['win']} Å field, {old['step']} Å step\n"
                      f"yellow: {new['win']} Å field, {new['step']} Å step (1600 positions each)", fontsize=8.5,
                      color=INK2)
        # inset: the probe core with both grids of scan positions at true scale
        ins = ax.inset_axes([0.60, 0.60, 0.38, 0.38]); z = CORE_A / 1.5; hz = int(round(z / px))
        cy, cx = np.unravel_index((np.abs(P) ** 2).argmax(), P.shape)
        ins.imshow(logI[cy - hz:cy + hz, cx - hz:cx + hz], cmap="magma", vmin=-4, vmax=0,
                   extent=[-z, z, -z, z], interpolation="bilinear")
        for rr, col, ms in ((old, "#e8e7e0", 3.4), (new, "#ffd24a", 1.6)):
            st = float(rr["step"]); g = np.arange(-z, z + 1e-9, st)
            gx, gy = np.meshgrid(g, g)
            ins.plot(gx.ravel(), gy.ravel(), ls="none", marker="o", ms=ms, color=col, mew=0)
        ins.set_xlim(-z, z); ins.set_ylim(-z, z); clean(ins)
        for sp in ins.spines.values():
            sp.set_visible(True); sp.set_color("white")
    ax = fig.add_subplot(gs[0, 2:])
    for al, col in ((70, INK), (80, MUTED)):
        for ls, d in (((0, (4, 2)), glob.glob(os.path.join(FIRST_HEAVY[0], f"recon_af_ceosopt_a{al:03d}_lab_NL*"))),
                      ("-", glob.glob(os.path.join(B1[0], f"recon_af_ceosopt_a{al:03d}_f20_lab_NL*")))):
            tr = full_trace(d[0]) if d else None
            if tr is None:
                continue
            it, e = tr
            ax.plot(it, e, color=col, ls=ls, lw=1.9)
            ax.annotate(f"{al} mrad, {'first-sweep scan' if ls != '-' else '20 Å scan'}", (it[-1], e[-1]),
                        xytext=(5, 0), textcoords="offset points", fontsize=8.5, va="center", color=col)
    ax.set_yscale("log"); ax.set_xlim(0, 262); ax.set_ylim(0.008, 0.7)
    ax.set_yticks([0.01, 0.03, 0.1, 0.3]); ax.set_yticklabels(["0.01", "0.03", "0.1", "0.3"]); ax.minorticks_off()
    ax.set_xlabel("full-engine iteration"); ax.set_ylabel("misfit to the data, ε")
    ax.set_title("c   Misfit through the full solve (same slices, same window; only the scan differs)", loc="left",
                 fontsize=10.5)
    for j, (lab, title) in enumerate((("ceosopt_a070", "70 mrad, first-sweep scan"), ("ceosopt_a070_f20", "70 mrad, 20 Å scan"),
                                      ("ceosopt_a080", "80 mrad, first-sweep scan"), ("ceosopt_a080_f20", "80 mrad, 20 Å scan"))):
        ax = fig.add_subplot(gs[1, j])
        lg = leg(lab)
        if lg is None:
            placeholder(ax, "figure data not pulled yet"); ax.set_title(f"{'defg'[j]}   {title}", fontsize=10, loc="left")
            continue
        insitu(ax, lg)
        ax.set_title(f"{'defg'[j]}   {title}", loc="left", fontsize=10)
        ax.set_xlabel("depth-summed phase, central 20 Å", fontsize=8.5, color=INK2)
    return save(fig, "fig4_scan.png")


# ------------------------------------------------------------------------------------------------ fig 5
def fig5_recons(a):
    sfig = _load_mod("make_simple_figs", os.path.join(HERE, "make_simple_figs.py"))
    T = {t["label"]: t for t in results_table()}
    cols = [al for al in ALPHAS if any(r[1] == al for r in RESULTS)]
    ceos = {al: lab for inst, al, lab, _ in RESULTS if inst == "CEOS"}
    rnd = {al: lab for inst, al, lab, _ in RESULTS if inst == "round"}
    base = leg("round_a070")
    rowA = section_row_A(base)
    fig, axes = plt.subplots(3, len(cols), figsize=(3.05 * len(cols), 10.2),
                             gridspec_kw=dict(height_ratios=[0.8, 1, 1]))
    seen = set()
    for i, al in enumerate(cols):
        for k, lab in ((0, ceos.get(al)), (1, ceos.get(al)), (2, rnd.get(al))):
            ax = axes[k, i]
            lg = leg(lab) if lab else None
            if lg is None:
                placeholder(ax, "figure data\nnot pulled yet"); continue
            if k == 0:
                insitu(ax, lg)
                t = T[lab]
                ax.set_title(f"{al} mrad\nCEOS · Pb {100 * t['Pb']:.0f} / Ti {100 * t['Ti']:.0f} / O {100 * t['O']:.0f} %",
                             fontsize=9.5, pad=5)
                continue
            nx = lg["V"].shape[2]
            sfig.xz_panel(ax, lg, rowA / lg["dx"], half_A=0.25,
                          xlim=(nx * lg["dx"] / 2 - 8, nx * lg["dx"] / 2 + 8), show_y=(i == 0))
            ax.grid(False)
            seen |= lg.get("species_seen", set())
            t = T[lab]
            who = "CEOS" if k == 1 else "round control"
            flag = "" if t["confusion"] <= CONF_MAX else "\nspecies labels unreliable"
            ax.set_title(f"{who} · depth error {t['z_rms']:.2f} Å{flag}", fontsize=9, pad=4)
    fig.suptitle("The reconstructions: the depth-summed phase, and depth sections through the same row of B-site columns",
                 x=0.01, ha="left", fontsize=12, fontweight="semibold")
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    sfig.legend_axes(fig, y=0.0, species=seen or None)
    return save(fig, "fig5_recons.png")


# ------------------------------------------------------------------------------------------------ the explorer
# Krivanek C_nm -> the symbols a CEOS tableau uses, and what each is
SYMBOL = {"C10": ("C1", "defocus"), "C12": ("A1", "two-fold astigmatism"), "C21": ("B2", "axial coma"),
          "C23": ("A2", "three-fold astigmatism"), "C30": ("C3", "spherical aberration"), "C32": ("S3", "star"),
          "C34": ("A3", "four-fold astigmatism"), "C41": ("B4", "fourth-order axial coma"),
          "C43": ("D4", "three-lobe"), "C45": ("A4", "five-fold astigmatism"),
          "C50": ("C5", "fifth-order spherical"), "C56": ("A5", "six-fold astigmatism")}
ROLE = {"C10": ("fought", "set for the smallest probe"), "C30": ("fought", "set for the smallest probe"),
        "C21": ("fought", "set against B4, its only partner"), "C50": ("held", "held at 1 mm"),
        "C41": ("hardware", "no knob (parasitic)"), "C45": ("hardware", "no knob (parasitic)"),
        "C56": ("hardware", "no knob (intrinsic)")}                       # the rest: retuned to 0.1 waves
CEOS_ADJUSTABLE = {"C10", "C12", "C21", "C23", "C30", "C32", "C34", "C43"}   # CEOS: C1 A1 B2 A2 C3 S3 A3 D4 (and C5)


def _unit(term, v):
    n = int(term[1])
    return (f"{v:.2f} Å" if abs(v) < 10 else f"{v:.0f} Å") if n == 1 else \
           (f"{v / 10:.1f} nm" if n == 2 else (f"{v / 1e4:.2f} µm" if n in (3, 4) else f"{v / 1e7:.2f} mm"))


def _uri(fig, fmt="jpeg"):
    import base64, io
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=100, bbox_inches="tight", pad_inches=0,
                **({"pil_kwargs": {"quality": 86}} if fmt == "jpeg" else {}))
    plt.close(fig)
    return f"data:image/{fmt};base64," + base64.b64encode(buf.getvalue()).decode()


def setup_entry(key, r, al):
    import make_ronchigram_fig as mrf
    ab = aberrations(r); c1 = float(r["c1"])
    terms = dict(ab); terms["C10"] = -c1                               # abTEM defocus = -C10
    out = []
    for t in sorted((k for k in terms if k.startswith("C") and k in SYMBOL), key=lambda k: (int(k[1]), int(k[2]))):
        v = float(terms[t]); sym, name = SYMBOL[t]
        role, why = ROLE.get(t, ("retuned", "retuned at this aperture to 0.1 waves"))
        out.append(dict(k=t, sym=sym, name=name, C=v, phi=float(terms.get("phi" + t[1:], 0.0)),
                        shown=_unit(t, -c1 if t == "C10" else v), deg=(None if t[2] == "0" else
                        round(float(np.degrees(terms.get("phi" + t[1:], 0.0))), 1)),
                        waves=round(abs(aw.waves(t, v, al)), 3), role=role, why=why, ceos=t in CEOS_ADJUSTABLE))
    L, N = 220.0, 3584; px = L / N; c = N // 2
    P = mmf.build_probe(al, ab, c1, L=L, N=N)
    win, d90, d99 = window_for(al, P, L)
    half = win / 2; h = int(round(half / px))
    fig, ax = plt.subplots(figsize=(3.2, 3.2))
    ax.imshow((np.abs(P[c - h:c + h, c - h:c + h]) ** 2) ** 0.4, cmap="magma", extent=[-half, half, -half, half],
              interpolation="bilinear")
    cy, cx = np.unravel_index((np.abs(P) ** 2).argmax(), P.shape)
    ax.add_artist(plt.Circle(((cx - c) * px, (cy - c) * px), d90 / 2, fill=False, color="white", lw=0.9, ls=(0, (3, 2))))
    bar = 0.2 * win
    ax.plot([-half * 0.9, -half * 0.9 + bar], [-half * 0.85] * 2, color="white", lw=2.2, solid_capstyle="butt")
    ax.text(-half * 0.9, -half * 0.8, f"{bar:g} Å", color="white", fontsize=9, va="bottom")
    ax.set_position([0, 0, 1, 1]); clean(ax)
    probe = _uri(fig)
    P2 = mmf.build_probe(al, ab, c1, L=160.0, N=2560)
    Rg, hr = mrf.ronchigram(P2, 160.0, al)
    fig, ax = plt.subplots(figsize=(3.2, 3.2))
    ax.imshow(Rg, cmap="gray", extent=[-hr, hr, -hr, hr], interpolation="bilinear")
    ax.add_artist(plt.Circle((0, 0), al, fill=False, color="#ffd24a", lw=1.0, ls=":"))
    ax.set_position([0, 0, 1, 1]); clean(ax)
    ronchi = _uri(fig)
    return dict(key=key, alpha=al, run=r["label"], terms=out, d90=round(d90, 1), window=win,
                probe=probe, ronchi=ronchi)


def explorer(a):
    """The page's setup explorer as an HTML fragment: per setup, the coefficient table (CEOS-adjustable terms marked)
    and, computed live in the page from that table, the round part, the non-round part and their sum across the
    aperture -- one colour cycle per wave, the same map as fig 1 -- beside the Ronchigram and the probe."""
    import matplotlib.cm as cm
    R = rows()
    setups = []
    for al in ALPHAS:
        setups.append(setup_entry(f"ceos-{al}", R[f"ceosopt_a{al:03d}"], al))
        setups.append(setup_entry(f"round-{al}", R[ROUND_CONTROL[al]], al))
        print(f"  explorer {al} mrad", flush=True)
    lut = (np.array([cm.twilight(i / 255.0)[:3] for i in range(256)]) * 255).round().astype(int).tolist()
    data = json.dumps(dict(lam=LAM, lut=lut, setups=setups), separators=(",", ":"))
    tabs = "".join(f'<button type="button" role="tab" id="xp-a{al}" data-a="{al}" aria-selected="false" tabindex="-1">'
                   f'{al}</button>' for al in ALPHAS)
    frag = XP_TEMPLATE.replace("{{TABS}}", tabs).replace("{{DATA}}", data.replace("</", "<\\/"))
    p = os.path.join(REPO, "aberration_experiment", "results", datetime.date.today().strftime("%G-W%V"),
                     "ceos_explorer.html")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, "w").write(frag)
    print(f"wrote {p} ({os.path.getsize(p) / 1e6:.2f} MB)")
    return p


XP_TEMPLATE = """<div class="xp" id="xp">
  <div class="xp-bar">
    <div class="xp-tabs" role="tablist" aria-label="Aperture semi-angle, mrad">{{TABS}}<span class="xp-unit">mrad</span></div>
    <div class="xp-seg" role="group" aria-label="Which column">
      <button type="button" id="xp-col-ceos" data-col="ceos" aria-pressed="true">CEOS column</button>
      <button type="button" id="xp-col-round" data-col="round" aria-pressed="false">round control</button>
    </div>
  </div>
  <div class="xp-body">
    <div class="xp-left">
      <div class="xp-meta" id="xp-meta"></div>
      <div class="xp-tablewrap"><table class="xp-table">
        <thead><tr><th scope="col"><span class="xp-sr">include</span></th><th scope="col">term</th><th scope="col">value</th>
          <th scope="col">angle</th><th scope="col">waves at edge</th><th scope="col">who sets it</th></tr></thead>
        <tbody id="xp-rows"></tbody>
      </table></div>
      <p class="xp-note">Highlighted rows are the terms a CEOS corrector lets the operator adjust. Untick a term to
        take it out of the three wavefront maps; the Ronchigram and probe are always the full setup.</p>
    </div>
    <div class="xp-maps">
      <figure class="xp-f"><canvas id="xp-round" width="280" height="280" role="img"></canvas>
        <figcaption><b>round part</b><span id="xp-round-pv"></span></figcaption></figure>
      <figure class="xp-f"><canvas id="xp-nonround" width="280" height="280" role="img"></canvas>
        <figcaption><b>non-round part</b><span id="xp-nonround-pv"></span></figcaption></figure>
      <figure class="xp-f"><canvas id="xp-total" width="280" height="280" role="img"></canvas>
        <figcaption><b>both together</b><span id="xp-total-pv"></span></figcaption></figure>
      <figure class="xp-f"><img id="xp-ronchi" alt=""><figcaption><b>Ronchigram</b><span>full setup</span></figcaption></figure>
      <figure class="xp-f"><img id="xp-probe" alt=""><figcaption><b>probe</b><span id="xp-probe-d"></span></figcaption></figure>
      <p class="xp-key">Maps: one colour cycle per wave across the aperture, the same picture as fig. 1.
        Probe: intensity in its reconstruction window; dashed circle holds 90 %.</p>
    </div>
  </div>
</div>
<script type="application/json" id="xp-data">{{DATA}}</script>
<script>
(function () {
  var D = JSON.parse(document.getElementById("xp-data").textContent);
  var byKey = {}; D.setups.forEach(function (s) { byKey[s.key] = s; });
  var ROUND = { C10: 1, C30: 1, C50: 1 };
  var st = { a: 80, col: "ceos", off: {} };
  try { var sv = JSON.parse(localStorage.getItem("xp-state") || "null"); if (sv && byKey[sv.col + "-" + sv.a]) { st.a = sv.a; st.col = sv.col; } } catch (e) {}

  function waves(s, part) {
    var n = 280, a = s.alpha * 1e-3, out = new Float32Array(n * n), lo = Infinity, hi = -Infinity;
    for (var i = 0; i < n; i++) {
      var tx = -a + 2 * a * i / (n - 1);                 // rows: theta_x, as abTEM's first array axis
      for (var j = 0; j < n; j++) {
        var ty = -a + 2 * a * j / (n - 1), th = Math.hypot(tx, ty), v = NaN;
        if (th <= a) {
          var ph = Math.atan2(ty, tx); v = 0;
          for (var k = 0; k < s.terms.length; k++) {
            var t = s.terms[k], isR = !!ROUND[t.k];
            if (st.off[t.k] || (part === "round" && !isR) || (part === "nonround" && isR)) continue;
            var nn = +t.k[1], m = +t.k[2];
            v += t.C * Math.pow(th, nn + 1) / (nn + 1) * Math.cos(m * (ph - t.phi));
          }
          v /= D.lam;
          if (v < lo) lo = v; if (v > hi) hi = v;
        }
        out[i * n + j] = v;
      }
    }
    return { w: out, pv: hi - lo };
  }
  function paint(id, s, part) {
    var cv = document.getElementById(id), ctx = cv.getContext("2d"), n = cv.width;
    var r = waves(s, part), img = ctx.createImageData(n, n), d = img.data;
    for (var p = 0; p < n * n; p++) {
      var v = r.w[p], q = p * 4;
      if (v !== v) { d[q + 3] = 0; continue; }
      var c = D.lut[Math.floor((v - Math.floor(v)) * 255.999)];
      d[q] = c[0]; d[q + 1] = c[1]; d[q + 2] = c[2]; d[q + 3] = 255;
    }
    ctx.putImageData(img, 0, 0);
    var pv = isFinite(r.pv) ? r.pv : 0;
    document.getElementById(id + "-pv").textContent = pv < 0.05 ? "none" : (pv < 10 ? pv.toFixed(1) : pv.toFixed(0)) + " waves peak to valley";
    cv.setAttribute("aria-label", part + " part of the wavefront at " + s.alpha + " mrad, " + pv.toFixed(1) + " waves peak to valley");
  }
  function rowsFor(s) {
    var tb = document.getElementById("xp-rows"); tb.textContent = "";
    s.terms.forEach(function (t) {
      var tr = document.createElement("tr"); if (t.ceos) tr.className = "xp-ceos";
      var id = "xp-t-" + t.k, cb = '<input type="checkbox" id="' + id + '"' + (st.off[t.k] ? "" : " checked") + ' aria-label="include ' + t.sym + '">';
      var chip = { fought: "c-open", retuned: "c-ok", hardware: "c-stop", held: "c-none" }[t.role];
      tr.innerHTML = "<td>" + cb + "</td><td><label for='" + id + "'><b>" + t.sym + "</b> <span class='xp-k'>" + t.k +
        "</span><span class='xp-name'>" + t.name + "</span></label></td><td class='xp-num'>" + t.shown + "</td><td class='xp-num'>" +
        (t.deg === null ? "–" : t.deg.toFixed(1) + "°") + "</td><td class='xp-num'>" + (t.waves < 10 ? t.waves.toFixed(2) : t.waves.toFixed(1)) +
        "</td><td><span class='chip " + chip + "'>" + t.role + "</span><span class='xp-why'>" + t.why + "</span></td>";
      tb.appendChild(tr);
      tr.querySelector("input").addEventListener("change", function (e) { st.off[t.k] = !e.target.checked; maps(); });
    });
  }
  function maps() { var s = byKey[st.col + "-" + st.a]; paint("xp-round", s, "round"); paint("xp-nonround", s, "nonround"); paint("xp-total", s, "total"); }
  function show() {
    var s = byKey[st.col + "-" + st.a];
    document.querySelectorAll("#xp .xp-tabs button").forEach(function (b) {
      var on = +b.dataset.a === st.a; b.setAttribute("aria-selected", on); b.tabIndex = on ? 0 : -1; });
    document.querySelectorAll("#xp .xp-seg button").forEach(function (b) { b.setAttribute("aria-pressed", b.dataset.col === st.col); });
    document.getElementById("xp-meta").innerHTML = "<b>" + s.alpha + " mrad · " + (st.col === "ceos" ? "CEOS column" : "round control") +
      "</b><span>run <code>" + s.run + "</code> · probe d90 " + s.d90.toFixed(1) + " Å · window " + s.window + " Å</span>";
    rowsFor(s); maps();
    var ro = document.getElementById("xp-ronchi"), pr = document.getElementById("xp-probe");
    ro.src = s.ronchi; ro.alt = "Ronchigram at " + s.alpha + " mrad";
    pr.src = s.probe; pr.alt = "Probe intensity at " + s.alpha + " mrad, d90 " + s.d90 + " Å";
    document.getElementById("xp-probe-d").textContent = "d90 " + s.d90.toFixed(1) + " Å, window " + s.window + " Å";
    try { localStorage.setItem("xp-state", JSON.stringify({ a: st.a, col: st.col })); } catch (e) {}
  }
  var tabs = Array.prototype.slice.call(document.querySelectorAll("#xp .xp-tabs button"));
  tabs.forEach(function (b, i) {
    b.addEventListener("click", function () { st.a = +b.dataset.a; st.off = {}; show(); });
    b.addEventListener("keydown", function (e) {
      var k = e.key === "ArrowRight" ? 1 : (e.key === "ArrowLeft" ? -1 : 0); if (!k) return;
      e.preventDefault(); var nb = tabs[(i + k + tabs.length) % tabs.length]; nb.focus(); nb.click(); });
  });
  document.querySelectorAll("#xp .xp-seg button").forEach(function (b) {
    b.addEventListener("click", function () { st.col = b.dataset.col; st.off = {}; show(); }); });
  show();
})();
</script>
"""


FIGS = {1: fig1_aberration, 2: fig2_probes, 3: fig3_growth, 4: fig4_scan, 5: fig5_recons, 6: fig6_numbers,
        7: fig7_stability, 8: explorer}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figs", type=int, nargs="+", default=sorted(FIGS))
    a = ap.parse_args()
    mmf.style()
    for k in a.figs:
        FIGS[k](a)


if __name__ == "__main__":
    main()
