"""
@file make_mono_page.py
@brief The Coherence Wall page (300 kV monochromated microscope, phase 1): every figure, judged table and tab viewer
       from the run outputs, then the page.

The ARM page's builder, reused unchanged (analysis/make_arm_page.py): runs found by tag in the search folders, the
honest atom scores (right element within 0.5 A in depth against a guessing level), the verdict rule, the depth sections,
the tab viewer, and the page's CSS (taken from the ARM template at build time, so the two pages stay one system).
Runs are listed in aberration_experiment/page/mono_results.json; one that has not landed is shown as a marked
placeholder carrying its 'pending' text, so the page can be republished at any point and fills in as runs are pulled.

Per run it adds two numbers: the share of the beam in its strongest coherent mode (the exact eigen-spectrum of the
simulator's own 16 configurations, from aberration_experiment/mono300/mode_spectrum.py), and the vacuum loading (rms
phase per slice in the two vacuum bands over that in the slab: the partial-coherence artefact as one number).

    ~/hyperspy-bundle/bin/python analysis/make_mono_page.py          # figures + page -> page/built/mono_page_built.html
    ~/hyperspy-bundle/bin/python analysis/make_mono_page.py --list   # what has landed, what is pending

Publishing and filling in: aberration_experiment/page/MONO_PAGE.md.
"""
from __future__ import annotations

import argparse
import contextlib
import datetime
import html
import importlib.util
import io
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location("make_arm_page", os.path.join(HERE, "make_arm_page.py"))
arm = importlib.util.module_from_spec(spec); sys.modules["make_arm_page"] = arm; spec.loader.exec_module(arm)
plt = arm.plt
PAGE = os.path.join(REPO, "aberration_experiment", "page")
MANIFEST = os.path.join(PAGE, "mono_results.json")
TEMPLATE = os.path.join(PAGE, "mono_2026-10-12.html")
OUT = os.path.join(PAGE, "built", "mono_page_built.html")
arm.FIGS = os.path.join(REPO, "aberration_experiment", "figs", "mono_page")      # this page's figures (print copies)
MONO = os.path.join(REPO, "aberration_experiment", "mono300")
FS_PER_EV_MM = None                     # A rms of focal spread per eV (FWHM) x mm of Cc at 300 kV, set by focal_spread()
OK, WARN, STOP = "#00795c", "#a97400", "#c2490a"                                 # the page's verdict colours
VCOL = {"good": OK, "okay": WARN, "off": STOP}
PB, TI = arm.SP[82][1], arm.SP[22][1]


# ====================================================================================== physics, read from the code
def focal_spread(dE_eV, cc_mm):
    """rms focal spread [A] of an energy spread (FWHM, eV) through Cc (mm) at 300 kV -- simulate_4dstem's own formula."""
    e0, E = 510998.95, 300e3
    return cc_mm * 1e7 * (dE_eV / (2 * np.sqrt(2 * np.log(2))) / E) * (1 + E / e0) / (1 + E / (2 * e0))


FS_PER_EV_MM = focal_spread(1.0, 1.0)


def _mode_spectrum_module():
    """mode_spectrum.py's definitions without its print loop (the loop runs 1000-configuration states for every point)."""
    src = open(os.path.join(MONO, "mode_spectrum.py")).read().split('\nprint(f"lambda')[0]
    ns = {"__name__": "mode_spectrum_defs"}
    exec(compile(src, os.path.join(MONO, "mode_spectrum.py"), "exec"), ns)
    return ns


_MS = {}


def strongest_mode(fs, src):
    """Share of the beam in its strongest mode: the sim's own 16 configurations (seed 7), as in mode_spectrum.txt."""
    if fs == 0 and src == 0:
        return 1.0
    if "ns" not in _MS:
        _MS["ns"] = _mode_spectrum_module()
    ns = _MS["ns"]
    return float(ns["spectrum"](*ns["offsets_sim"](fs, src))[0])


def vacuum_loading(lg):
    """rms phase per slice in the two 4 A vacuum bands over that in the slab (coherent reference ~0.3)."""
    V, dz = lg["V"], lg["dz"]
    zc = (np.arange(V.shape[0]) + 0.5) * dz
    s = np.array([np.std(V[i]) for i in range(V.shape[0])])
    vac = (s[zc < arm.sfig.ZVAC].mean() + s[zc > arm.sfig.BOXZ - arm.sfig.ZVAC].mean()) / 2
    slab = s[(zc > arm.sfig.ZVAC) & (zc < arm.sfig.BOXZ - arm.sfig.ZVAC)].mean()
    return float(vac / slab)


# ====================================================================================== one record per run
def record(R, key):
    """One run: its scores, verdict, vacuum loading and probe record, and its state --
    in (own kernels) · provisional (atoms read with the fixed-probe control's kernels, a '_kern-' tag) · volume (a
    reconstruction but no atoms) · crashed (an analysis row but no reconstruction) · pending (nothing landed).
    found = scored, so every table and figure that needs atoms reads it directly."""
    info = R.find(key)
    spec = info["spec"]
    tag = info["tag"]
    leg = spec["legs"][spec["tag"].index(tag)] if tag and "legs" in spec else spec["leg"]
    provisional = bool(tag) and "_kern-" in tag
    earlier = bool(tag) and not info["primary"] and not provisional and bool(spec.get("score_fallbacks"))
    if provisional or earlier:
        info = dict(info, primary=True)           # its atoms are the point of the tag: score them, labelled so
    r = dict(key=key, info=info, spec=spec, leg=leg, tag=tag, fs=float(spec["fs"]), src=float(spec["src"]),
             w1=strongest_mode(float(spec["fs"]), float(spec["src"])), sc=None, v=None, vac=None, kick=None,
             found=False, vol=False, provisional=provisional, state="pending")
    if info["found"] and info["primary"] and info["af"]:
        r["sc"] = arm.atom_scores(info)
        r["found"] = r["sc"] is not None
        r["v"] = arm.verdict(r["sc"])
    if info["found"] and info["fd"]:
        lg = arm.load_leg(info)
        r["vol"] = lg is not None
        r["vac"] = vacuum_loading(lg) if lg is not None else None
    if info["found"]:
        r["kick"], _ = R.kick(leg)
        r["state"] = ("provisional" if provisional else "earlier" if earlier else "in") if r["found"] else \
                     ("volume" if r["vol"] else "crashed")
    return r


STATE_NOTE = {"volume": "reconstruction in, no atoms yet", "crashed": "crashed before a reconstruction",
              "provisional": "provisional: atoms read with the fixed-probe control's kernels",
              "earlier": "earlier run, solver before the 2026-10-16 fix; its rerun is running"}


def pct(sc, z, which="right"):
    return 100 * sc[which][z]


# ====================================================================================== figures
def fig_wall(recs, man):
    """Pb and Ti in the right place against each single effect, with the guessing level; probe-mode legs when landed."""
    fig, axs = plt.subplots(1, 2, figsize=(11.6, 4.5), sharey=True, gridspec_kw=dict(wspace=0.07))
    ref = recs["known_a080"]
    for ax, grp, xkey, xlab in ((axs[0], "focal", "fs", "focal spread (Å rms)"), (axs[1], "source", "src", "source size (Å FWHM)")):
        pts = [ref] + sorted([r for r in recs.values() if r["spec"]["group"] == grp and r["found"]], key=lambda r: r[xkey])
        x = np.array([r[xkey] for r in pts])
        for z, col, nm in ((82, PB, "Pb"), (22, TI, "Ti")):
            ax.plot(x, [pct(r["sc"], z) for r in pts], "-o", color=col, lw=1.8, ms=6, label=f"{nm}, right place", zorder=4)
            mix = [m for m in recs.values() if m["spec"]["group"] == "modes" and m["found"] and m[xkey] > 0
                   and (m["fs"] == 0 if xkey == "src" else m["src"] == 0)]
            for prov, lab, al in ((False, f"{nm}, probe modes", 1.0), (True, f"{nm}, probe modes (provisional)", .45)):
                mm = [m for m in mix if m["provisional"] == prov]
                if mm:
                    dx = 0.03 * (7.6 if xkey == "fs" else 0.64)      # drawn just right of their point, clear of its dots
                    ax.plot([m[xkey] + dx for m in mm], [pct(m["sc"], z) for m in mm], "D", mfc="none", mec=col, mew=1.6, ms=8,
                            alpha=al, label=lab, zorder=5)
        good = [r[xkey] for r in pts if r["v"] == "good"]
        bad = [r[xkey] for r in pts if r["v"] == "off"]
        if good and bad:
            lo, hi = max(good), min(b for b in bad if b > max(good))
            ax.axvspan(lo, hi, ymax=100 / 112, color=WARN, alpha=.13, lw=0, zorder=0)
            ax.text((lo + hi) / 2, 101, "the wall", ha="center", va="bottom", fontsize=9, color=WARN)
        for r in pts:                                   # the verdict of each point, a coloured square above the plot
            ax.plot(r[xkey], 108, "s", ms=7, color=VCOL[r["v"]], zorder=6)
        ax.set_xlabel(xlab); ax.set_ylim(0, 112); ax.set_yticks(range(0, 101, 20)); ax.grid(alpha=.25, lw=.6)
        ax.spines[["top", "right"]].set_visible(False)
    rs = man["real_source_A"]
    axs[1].axvspan(rs[0], rs[1], ymax=100 / 112, facecolor="none", edgecolor="#8a96a3", hatch="////", lw=0, zorder=0)
    axs[1].text(0.415, 62, "typical source size\non real instruments", fontsize=8.5, color=arm.INK2, va="center",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none"))
    axs[1].set_xlim(-0.02, rs[1] + 0.02)
    axs[0].set_xlim(-0.2, 7.95)
    top = axs[0].secondary_xaxis("top", functions=(lambda f: f / (FS_PER_EV_MM * 1.4), lambda e: e * FS_PER_EV_MM * 1.4))
    top.set_xlabel("energy spread through Cc 1.4 mm (eV FWHM)", fontsize=9.5)
    axs[0].set_ylabel("atoms in the right place (%)")
    for ax, t in zip(axs, ("Focal spread alone", "Source size alone")):
        ax.text(0, 1.2, t, transform=ax.transAxes, fontsize=11, ha="left", va="bottom", color=arm.INK)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=len(l), fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.06))
    return arm.save(fig, "fig_wall.png")


def fig_instruments(recs, man):
    """Each instrument class as the range of focal spread its energy spread and Cc give, on the measured wall."""
    focal = sorted([recs["known_a080"]] + [r for r in recs.values() if r["spec"]["group"] == "focal" and r["found"]],
                   key=lambda r: r["fs"])
    good = max(r["fs"] for r in focal if r["v"] == "good")
    off = min(r["fs"] for r in focal if r["v"] == "off" and r["fs"] > good)
    ins = man["instruments"]
    fig, ax = plt.subplots(figsize=(11.6, 0.62 * len(ins) + 1.5))
    lo, hi = 0.5, 40
    ax.axvspan(lo, good, color=OK, alpha=.10, lw=0); ax.axvspan(good, off, color=WARN, alpha=.16, lw=0)
    ax.axvspan(off, hi, color=STOP, alpha=.08, lw=0)
    for x, t, c in ((np.sqrt(lo * good), "good", OK), (np.sqrt(good * off), "the wall", WARN), (np.sqrt(off * hi), "off", STOP)):
        ax.text(x, len(ins) - 0.35, t, ha="center", va="bottom", fontsize=9.5, color=c, fontweight="bold")
    for i, s in enumerate(ins[::-1]):
        a = focal_spread(s["dE"][0], s["cc"][0]); b = focal_spread(s["dE"][1], s["cc"][1])
        if abs(b - a) < 1e-9:
            ax.plot(a, i, "o", color=arm.INK, ms=7, zorder=4)
        else:
            ax.plot([a, b], [i, i], "-", color=arm.INK, lw=6, solid_capstyle="butt", alpha=.85, zorder=4)
        ax.text(np.sqrt(a * b), i + 0.3, f"{s['dE'][0]:g}{'' if s['dE'][0] == s['dE'][1] else '–' + format(s['dE'][1], 'g')} eV × "
                f"{s['cc'][0]:g}{'' if s['cc'][0] == s['cc'][1] else '–' + format(s['cc'][1], 'g')} mm",
                ha="center", va="bottom", fontsize=8.5, color=arm.INK2)
    for r in focal:                                     # the tested points
        if r["fs"] > 0:
            ax.plot(r["fs"], -0.8, "v", color=VCOL[r["v"]], ms=8, clip_on=False)
    ax.text(lo * 1.03, -0.8, "tested:", va="center", fontsize=8.5, color=arm.INK2)
    ax.set_xscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(-1.2, len(ins) - 0.1)
    ax.set_yticks(range(len(ins))); ax.set_yticklabels([s["name"] for s in ins[::-1]], fontsize=9.5)
    ax.set_xticks([0.5, 1, 2, 3.5, 5, 10, 20, 40]); ax.set_xticklabels(["0.5", "1", "2", "3.5", "5", "10", "20", "40"])
    ax.set_xlabel("focal spread (Å rms) = energy spread × Cc × " + f"{FS_PER_EV_MM:.1f} Å per eV·mm at 300 kV")
    ax.spines[["top", "right", "left"]].set_visible(False); ax.tick_params(axis="y", length=0); ax.grid(axis="x", alpha=.25, lw=.6)
    return arm.save(fig, "fig_instruments.png")



WHERE_ROWS = [("focal spread only (§02)", ["known_a080", "fs0p75", "fs1p5", "fs2p5", "fs3p5", "fs5", "fs5_mix"]),
              ("+ source 0.1 Å", ["src0p1", "fs2p5_src0p1", "ceil_s01", "fs3p5_src0p1"]),
              ("+ source 0.2 Å", ["src0p2", "src0p2_mix", "fs2p5_src0p2", "ceil_s02", "fs3p5_src0p2"]),
              ("+ source 0.1 Å + 100 K", ["cool100_fixed", "cool100_mix"]),
              ("+ source 0.1 Å + 20 K", ["cool20_fixed", "cool20_mix"]),
              ("+ source 0.1 Å, probe searched (C)", ["searched"])]


def fig_where(recs, man):
    """Where it works: the combined runs at their focal spread, under the Spectra 300 monochromators' ranges (no off-scale
    instruments). Square colour = verdict, hollow = true probe fixed, filled = probe modes, grey = not scored yet."""
    focal = sorted([recs["known_a080"]] + [r for r in recs.values() if r["spec"]["group"] == "focal" and r["found"]],
                   key=lambda r: r["fs"])
    good = max(r["fs"] for r in focal if r["v"] == "good")
    off = min(r["fs"] for r in focal if r["v"] == "off" and r["fs"] > good)
    ins = {s["name"]: s for s in man["instruments"]}
    mono, ulti = ins["Spectra 300, X-FEG Mono"], ins["Spectra 300, X-FEG UltiMono"]
    bars = [("Spectra 300, X-FEG Mono, at the 90 meV quoted in use", mono["dE"][0], mono["dE"][0], mono["cc"]),
            ("Spectra 300, X-FEG Mono, at its < 0.2 eV specification", mono["dE"][1], mono["dE"][1], mono["cc"]),
            ("Spectra 300, X-FEG UltiMono", ulti["dE"][0], ulti["dE"][1], ulti["cc"])]
    lo, hi = -0.45, 6.0
    fig, (at, ab) = plt.subplots(2, 1, figsize=(11.6, 6.4), sharex=True, gridspec_kw=dict(height_ratios=[1.25, 3.4], hspace=.08))
    for ax in (at, ab):
        ax.axvspan(lo, good, color=OK, alpha=.10, lw=0); ax.axvspan(good, off, color=WARN, alpha=.16, lw=0)
        ax.axvspan(off, hi, color=STOP, alpha=.08, lw=0)
        ax.spines[["top", "right", "left"]].set_visible(False); ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", alpha=.25, lw=.6)
        ax.axvline(2.5, color=arm.INK2, lw=1, ls="--", alpha=.7, zorder=2)
    at.text(2.55, -0.5, "the setting used: 2.5 Å (90 meV, Cc ≈ 1.6 mm)", ha="left", va="bottom", fontsize=8.5, color=arm.INK2)
    for x, t, c in (((lo + good) / 2, "good", OK), ((good + off) / 2, "the wall", WARN), ((off + hi) / 2, "off", STOP)):
        at.text(x, len(bars) - 0.35, t, ha="center", va="bottom", fontsize=9.5, color=c, fontweight="bold")
    for i, (name, d0, d1, cc) in enumerate(bars[::-1]):
        a, b = focal_spread(d0, cc[0]), focal_spread(d1, cc[1])
        at.plot([a, min(b, hi)], [i, i], "-", color=arm.INK, lw=6, solid_capstyle="butt", alpha=.85, zorder=4)
        lab = f"{d0:g}{'' if d0 == d1 else '–' + format(d1, 'g')} eV × {cc[0]:g}–{cc[1]:g} mm: {a:.1f}–{b:.1f} Å"
        at.text(min((a + min(b, hi)) / 2, hi - 0.9), i + 0.22, lab + (" (runs off the plot)" if b > hi else ""),
                ha="center", va="bottom", fontsize=8.5, color=arm.INK2)
    at.set_ylim(-0.6, len(bars) + 0.1); at.set_yticks(range(len(bars)))
    at.set_yticklabels([b[0] for b in bars[::-1]], fontsize=9.5)
    n = len(WHERE_ROWS)
    for i, (lab, keys) in enumerate(WHERE_ROWS):
        y = n - 1 - i
        pts = [recs[k] for k in keys if k in recs]
        byx = {}
        for r in pts:
            byx.setdefault(r["fs"], []).append(r)
        for x, rs in byx.items():
            rs.sort(key=lambda r: (r["spec"]["group"] in ("modes", "ceiling") or r["spec"].get("kind") == "modes"
                                   or r["key"] == "searched", r["key"]))
            for j, r in enumerate(rs):
                xx = x + (j - (len(rs) - 1) / 2) * 0.26
                mode = r["spec"]["group"] in ("modes", "ceiling") or r["spec"].get("kind") == "modes" or r["key"] == "searched"
                if r["found"] and not r["provisional"]:
                    c = VCOL[r["v"]]
                    ab.plot(xx, y, "s", ms=11, mec=c, mew=2, mfc=c if mode else "white", zorder=5)
                    ab.text(xx, y + 0.24, f"{pct(r['sc'], 82):.0f}/{pct(r['sc'], 22):.0f}", ha="center", va="bottom",
                            fontsize=8, color=arm.INK)
                else:
                    ab.plot(xx, y, "s", ms=11, mec="#9aa3ad", mew=2, mfc="#9aa3ad" if mode else "white", zorder=5)
                    ab.text(xx, y + 0.24, "running", ha="center", va="bottom", fontsize=8, color="#7d8a99")
    ab.set_ylim(-0.6, n - 0.2); ab.set_yticks(range(n)); ab.set_yticklabels([r[0] for r in WHERE_ROWS[::-1]], fontsize=9.5)
    ab.set_xlim(lo, hi); ab.set_xticks([0, 1, 2, 2.5, 3.5, 5, 6]); ab.set_xticklabels(["0", "1", "2", "2.5", "3.5", "5", "6"])
    ab.set_xlabel("focal spread (Å rms) = energy spread × Cc × " + f"{FS_PER_EV_MM:.1f} Å per eV·mm at 300 kV")
    h = [plt.Line2D([], [], ls="", marker="s", ms=9, mec=arm.INK2, mew=2, mfc="white", label="true probe held fixed"),
         plt.Line2D([], [], ls="", marker="s", ms=9, mec=arm.INK2, mew=2, mfc=arm.INK2, label="probe modes"),
         plt.Line2D([], [], ls="", marker="s", ms=9, mec="#9aa3ad", mew=2, mfc="#9aa3ad", label="running, not scored yet")]
    ab.legend(handles=h, loc="lower right", fontsize=8.5, frameon=True, framealpha=.9, ncol=3,
              title="colour: verdict · numbers: Pb / Ti % at the right depth", title_fontsize=8.5)
    return arm.save(fig, "fig_where_it_works.png")

def fig_strip(name, keys, recs, cuts, title_of):
    """Depth sections side by side: one column per run, the Pb row and the B-site row of the coherent reference's cuts."""
    keys = [k for k in keys if recs[k]["vol"]]
    n = len(keys)
    fig, axs = plt.subplots(2, n, figsize=(2.05 * n + 0.4, 8.7), squeeze=False)
    for j, k in enumerate(keys):
        r = recs[k]; lg = arm.load_leg(r["info"]); sc = r["sc"]
        nx = lg["V"].shape[2]
        for i, (cut, rA) in enumerate(cuts):
            ax = axs[i, j]
            arm.xz_section(ax, lg, sc, rA / lg["dx"], xlim=(nx * lg["dx"] / 2 - 8, nx * lg["dx"] / 2 + 8), show_y=(j == 0))
            if i == 0:
                body = (f"\nPb {pct(sc, 82):.0f} vs {pct(sc, 82, 'guess'):.0f}\nTi {pct(sc, 22):.0f} vs {pct(sc, 22, 'guess'):.0f}"
                        f"\n{r['v']}{' (provisional)' if r['provisional'] else ' (earlier solver)' if r['state'] == 'earlier' else ''}") \
                    if sc else "\nno atoms yet\n\n"
                ax.set_title(f"{title_of(r)}\nstrongest mode {100 * r['w1']:.0f} %{body}", fontsize=8.5, loc="left", color=arm.INK)
            ax.set_xlabel(f"x (Å) · {cut}", fontsize=8)
    fig.tight_layout()
    return arm.save(fig, name)


def fig_modes(mix, recs):
    """Mode powers of each landed probe-mode leg: the exact state it started from against what it recovered."""
    m = [r for r in mix if r["kick"] and r["kick"].get("mode_power")]
    if not m:
        return None
    fig, axs = plt.subplots(1, len(m), figsize=(2.9 * len(m) + 0.6, 3.2), squeeze=False, sharey=True)
    for ax, r in zip(axs[0], m):
        k = r["kick"]; rec = np.array(k["mode_power"]); st = np.array(k.get("start_mode_power") or [])
        i = np.arange(1, len(rec) + 1)
        if len(st):
            ax.bar(i - 0.2, st, 0.4, color="#9aa6b3", label="start: a real gun's mixed state")
        ax.bar(i + 0.2, rec, 0.4, color=arm.INK, label="recovered")
        ax.axhline(r["w1"], color=WARN, lw=1.2, ls=(0, (4, 3)), label="truth: the data's strongest mode")
        ax.set_title(r["spec"]["short"], fontsize=10, loc="left"); ax.set_xticks(i)
        ax.set_xlabel("mode"); ax.spines[["top", "right"]].set_visible(False)
    axs[0][0].set_ylabel("share of the beam")
    h, l = axs[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.1))
    fig.tight_layout()
    return arm.save(fig, "fig_modes.png")


# ====================================================================================== HTML parts
def setting(r):
    parts = []
    if r["fs"]:
        parts.append(f"focal spread {r['fs']:g} Å")
    if r["src"]:
        parts.append(f"source {r['src']:g} Å")
    return " + ".join(parts) or "coherent"


def judged_table(rows, caption):
    head = ["setting", "strongest mode (%)", *arm.ATOM_HEAD, "vacuum / slab", "verdict"]
    trs = []
    for r, label, cls in rows:
        if not r["found"]:
            trs.append(f"<tr><td>{html.escape(label)}</td><td>{100 * r['w1']:.0f}</td><td colspan='6'>{arm.pending_box(r['info'])}</td></tr>")
            continue
        trs.append(f"<tr class='{cls}'><td>{html.escape(label)}</td><td>{100 * r['w1']:.0f}</td>{arm.atoms_cells(r['sc'])}"
                   f"<td>{r['vac']:.2f}</td>{arm.vcell(r['v'])}</tr>")
    return arm.table(head, trs, caption)


def series_tabs(tid, keys, recs, cuts):
    """One tab per run: its in-plane phase and the two depth sections, big (arm.compare_fig, one row)."""
    panes = []
    for k in keys:
        r = recs[k]
        if not r["vol"]:
            continue
        kind = r["spec"]["name"].split(", ")[-1] + ("" if r["state"] == "in" else f"  ·  {STATE_NOTE[r['state']]}")
        p = arm.compare_fig(f"row_{k}.png", [(f"{setting(r)}  ·  {kind}", r["info"])], CUTS[0])
        panes.append((r["spec"]["short"], arm.img(p, r["spec"]["name"], maxw=1500, jpeg=80)))
    return arm.tabs(tid, panes) if panes else "<p>Nothing landed yet.</p>"


def headline(recs):
    """The summary's opening sentence: the last good focal spread and source size alone, and together."""
    def wall(grp, key):
        pts = sorted([recs["known_a080"]] + [r for r in recs.values() if r["spec"]["group"] == grp and r["found"]], key=lambda r: r[key])
        g = max(r[key] for r in pts if r["v"] == "good")
        o = min(r[key] for r in pts if r["v"] == "off" and r[key] > g)
        mid = [r for r in pts if g < r[key] < o]
        return g, o, mid
    fg, fo, _ = wall("focal", "fs")
    sg, so, smid = wall("source", "src")
    joint = [r for r in recs.values() if r["spec"]["group"] == "joint" and r["found"]]
    jgood = sorted([r for r in joint if r["v"] == "good"], key=lambda r: (r["fs"], r["src"]))
    line = (f"With the true probe given, the depth survives up to <b>{fg:g} Å</b> of focal spread alone, a <b>{sg:g} Å</b> "
            f"source alone")
    line += (f", or <b>{jgood[-1]['fs']:g} Å + {jgood[-1]['src']:g} Å</b> together." if jgood else ".")
    return dict(WALL_LINE=line)

ORDER = {"off": 0, "okay": 1, "good": 2}


# ====================================================================================== the Spectra Mono point (A, B, C)
P2_KEYS = ["fs2p5_src0p1", "ceil_s01", "fs2p5_src0p2", "ceil_s02", "cool100_fixed", "cool100_mix", "cool20_fixed",
           "cool20_mix", "searched"]
LAM_A, ALPHA = 0.019687, 0.080


def p2_table(recs):
    """A, B and C beside their controls: the ceiling at the Spectra Mono point, cooled, and with the probe found by search."""
    head = ["run, all at 2.5 Å focal spread", "probe", "strongest mode (%)", *arm.ATOM_HEAD, "vacuum / slab", "verdict"]
    rows = [("A", "fs2p5_src0p1", "fixed (step 2A)"), ("A", "ceil_s01", "8 modes"), ("A", "fs2p5_src0p2", "fixed (step 2A)"),
            ("A", "ceil_s02", "8 modes"), ("B", "ph_rt", "fixed"), ("B", "ph_100", "fixed"), ("B", "cool100_fixed", "fixed"),
            ("B", "cool100_mix", "8 modes"), ("B", "cool20_fixed", "fixed"), ("B", "cool20_mix", "8 modes"),
            ("C", "searched", "searched, 8 modes")]
    trs = []
    for part, k, probe in rows:
        r = recs[k]
        label = {"ph_rt": "coherent beam, room temperature", "ph_100": "coherent beam, 100 K",
                 "searched": f"source {r['src']:g} Å, probe searched"}.get(
            k, f"source {r['src']:g} Å" + (f" + {r['spec']['temp']} K" if r["spec"].get("kind") in ("fixed", "modes") else ""))
        cls = "ref" if k in ("ph_rt", "ph_100") else ""
        lead = f"<td><b>{part}</b> {html.escape(label)}</td><td>{html.escape(probe)}</td><td>{100 * r['w1']:.0f}</td>"
        if r["found"]:
            trs.append(f"<tr class='{cls}'>{lead}{arm.atoms_cells(r['sc'])}<td>{r['vac']:.2f}</td>{arm.vcell(r['v'])}</tr>")
        elif r["vol"]:
            why = ("failed: the search left C1 32 Å off, the reconstruction drifted into the vacuum, and the finder's Pb "
                   "and Ti templates could not be made. A rerun with a third search stage is running." if k == "searched"
                   else "reconstruction in, no atoms yet")
            trs.append(f"<tr>{lead}<td colspan='4' class='wrap'><span class='chip c-stop'>no atoms</span> "
                       f"<span class='muted-s'>{html.escape(why)}</span></td><td>{r['vac']:.2f}</td><td>–</td></tr>")
        else:
            trs.append(f"<tr>{lead}<td colspan='6'>{arm.pending_box(r['info'])}</td></tr>")
    return arm.table(head, trs, "All at the Spectra Mono's setting: focal spread 2.5 Å (about 90 meV through C<sub>c</sub> 1.6 mm), "
                     "the ARM aberrations, 80 mrad, no shot noise. <b>A</b>: the source size, 0.1 or 0.2 Å. <b>B</b>: plus atomic "
                     "vibration at room temperature, 100 K (liquid nitrogen) and 20 K (liquid helium, where only the zero-point "
                     "motion is left). <b>C</b>: the probe found from the data (below). Probe: <b>fixed</b> = the true probe given "
                     "and held; <b>8 modes</b> = refined from the data as 8 modes. The first two B rows are references with a "
                     "perfectly coherent beam, vibration only. " + COLS_NOTE)


def _search_dir(R, name):
    for d in R.dirs:
        for c in (os.path.join(d, "results", name), os.path.join(d, name)):
            if os.path.exists(os.path.join(c, "start.json")):
                return c
    return None


def search_table(R, recs):
    """The outer search: start drawn from the corrector's precision, each stage's answer, against the truth (record only)."""
    spec = recs["searched"]["spec"]
    found = [(n, _search_dir(R, n)) for n in spec.get("search_dirs", [])]
    found = [(n, d) for n, d in found if d]
    if not found:
        return arm.pending_box(recs["searched"]["info"])
    name, d = found[0]
    st = json.load(open(os.path.join(d, "start.json")))
    tr, s0 = st["truth_record_only"], st["start"]
    stages = [json.load(open(os.path.join(d, f"stage{k}_best.json"))) for k in (1, 2, 3)
              if os.path.exists(os.path.join(d, f"stage{k}_best.json"))]
    vec = lambda c12, phi: np.array([c12 * np.cos(2 * phi), c12 * np.sin(2 * phi)])
    def errs(x):
        e1 = abs(x["c1_A"] - tr["c1_A"]); e3 = abs(x["c3_A"] - tr["c3_A"])
        ea = float(np.hypot(*(vec(x["c12_A"], x["phi12"]) - vec(tr["c12_A"], tr["phi12"]))))
        w = lambda e, n: e * ALPHA ** n / (n * LAM_A)           # waves at the aperture edge for a term of order n in alpha
        return [(e1, w(e1, 2)), (e3, w(e3, 4)), (ea, w(ea, 2))]
    cols = [("drawn start", s0)] + [(f"after stage {x['stage']}", x) for x in stages]
    head = ["term", "truth"] + [c for c, _ in cols]
    fmt = [lambda v: f"{v['c1_A']:.1f} Å", lambda v: f"{v['c3_A'] / 1e4:.3f} µm",
           lambda v: f"{v['c12_A']:.1f} Å @ {np.degrees(v['phi12']):.0f}°"]
    names = ["C1 (defocus)", "C3 (spherical)", "A1 (two-fold astigmatism)"]
    trs = []
    for i in range(3):
        cells = "".join(f"<td>{fmt[i](v)}<br><span class='muted-s'>off by {errs(v)[i][0]:.3g} {'µm' if i == 1 else 'Å'}"
                        .replace(f"{errs(v)[i][0]:.3g} µm", f"{errs(v)[i][0] / 1e4:.3g} µm")
                        + f" · {errs(v)[i][1]:.1f} waves</span></td>" for _, v in cols)
        trs.append(f"<tr><td>{names[i]}</td><td>{fmt[i](tr)}</td>{cells}</tr>")
    note = (" A third stage, C1 × C3 again at the found A1, is running." if len(stages) < 3 else "")
    return arm.table(head, trs, f"Each row is one aberration the search had to find. <b>Drawn start</b>: what an operator would "
                     f"begin with, the true value plus an error drawn at the corrector's own quoted 95 % precision (C1 39.5 Å, "
                     f"C3 0.80 µm, A1 59 Å; seed 1). <b>Stage 1</b> tries a 7 × 7 grid of C1 and C3 with A1 left at the start; "
                     f"<b>stage 2</b> a 5 × 5 grid of A1 with C1 and C3 at stage 1's answer. Each trial is a short reconstruction "
                     f"with that probe held fixed; the answer is the lowest data error, refined by a smooth fit through its "
                     f"neighbours. <b>Waves</b>: the phase error the "
                     f"remaining mistake leaves at the edge of the 80 mrad aperture, in wavelengths; a quarter wave is the usual "
                     f"tolerance. The truth is shown only for the record; the search never sees it.{note}")


COLS_NOTE = ("Columns: <b>strongest mode</b>, how coherent the beam is (the share of its intensity in its single most "
             "coherent component; 100 % = perfectly coherent). <b>Right place</b>, the score: atoms the finder puts as the right "
             "element within 0.5 Å of their true depth. <b>Guessing</b>, the same found atoms slid to random depths along their "
             "columns: what the score would be with no depth information. <b>Wrong element</b>, found atoms given the wrong "
             "element. <b>95 % depth bar</b>, the finder's own uncertainty on an atom's depth (median). <b>Vacuum / slab</b>, how "
             "much false structure the reconstruction puts in the empty vacuum layers, relative to the specimen: 0.31 for a "
             "perfect beam, near 1 when the vacuum looks as strong as the specimen.")


def status(recs, man):
    """What has been done: one row per step, what it changes, how many runs are in, and its state."""
    rows = [("Part 1, the scan: how coherent must the beam be? The true probe is given, no other physics.",),
            ("Step 1", "one knob at a time: six focal spreads (0.75–7.29 Å), four source sizes (0.05–0.4 Å)",
             "true probe, held fixed", "focal source"),
            ("Step 2A", "both knobs together: 2.5 and 3.5 Å, each with a 0.1 and a 0.2 Å source", "true probe, held fixed", "joint"),
            ("Step 2B", "probe modes at four settings: 5 and 7.29 Å focal spread, 0.2 and 0.4 Å source",
             "refined from the data as 7–8 modes, started from a real gun's mixed beam", "modes"),
            ("Part 2, the physics put back, one piece at a time, at the Spectra Mono's setting (2.5 Å focal spread).",),
            ("A", "the source size: 0.1 and 0.2 Å", "refined as 8 modes (corrected solver)", "ceiling"),
            ("B", "add atomic vibration: phonons at 100 K and 20 K (0.1 Å source)",
             "true probe held fixed, and refined as 8 modes", "cooled"),
            ("C", "add not knowing the probe: start from the corrector's quoted precision, find it by an outer search",
             "searched, then refined as 8 modes", "searched")]
    trs = []
    for row in rows:
        if len(row) == 1:
            trs.append(f"<tr class='grp'><td colspan='5'>{html.escape(row[0])}</td></tr>")
            continue
        step, what, probe, gs = row
        rs = [r for r in recs.values() if r["spec"]["group"] in gs.split()]
        n = sum(r["state"] == "in" for r in rs)
        chip = '<span class="chip c-ok">done</span>' if n == len(rs) else '<span class="chip c-warn">running</span>'
        trs.append(f"<tr><td><b>{step}</b></td><td class='wrap'>{html.escape(what)}</td><td class='wrap'>{html.escape(probe)}</td>"
                   f"<td>{n} of {len(rs)}</td><td>{chip}</td></tr>")
    trs.append("<tr><td><b>D</b></td><td class='wrap'>add shot noise: a dose series, 1e5–1e8 e/Å², with the searched start, at "
               "the point that survives</td><td class='wrap'>searched, then modes</td><td>–</td><td><span class='chip c-none'>next"
               "</span></td></tr>")
    return arm.table(["", "what changes", "probe", "runs in", ""], trs,
                     "Every run: 300 kV, 80 mrad, no shot noise; atomic vibration only in B. Running now: step 2B's four points, "
                     "with a solver bug fixed (§05), and C with a third search stage (§07).")


def joint_grid(recs):
    """Focal spread (rows) x source size (columns): step-1 edges and the step-2A joint points, judged."""
    by = {(r["fs"], r["src"]): r for r in recs.values() if r["spec"]["group"] in ("ref", "focal", "source", "joint")
          and r["key"] != "cold_feg"}
    fss, srcs = [0.0, 2.5, 3.5], [0.0, 0.1, 0.2]
    trs = []
    for f in fss:
        tds = []
        for s in srcs:
            r = by.get((f, s))
            if r is None:
                tds.append("<td>–</td>"); continue
            if not r["found"]:
                tds.append(f"<td class='wrap'><span class='chip c-warn'>pending</span><br><span class='muted-s'>strongest mode "
                           f"{100 * r['w1']:.0f} %</span></td>"); continue
            v = r["v"]
            tds.append(f"<td class='{arm.CHIP[v][1]} wrap'><span class='chip {arm.CHIP[v][0]}'>{v}</span><br>Pb {pct(r['sc'], 82):.0f} · "
                       f"Ti {pct(r['sc'], 22):.0f} %<br><span class='muted-s'>guessing {pct(r['sc'], 82, 'guess'):.0f} · "
                       f"{pct(r['sc'], 22, 'guess'):.0f} · strongest mode {100 * r['w1']:.0f} %</span></td>")
        trs.append(f"<tr><td><b>{'no spread' if f == 0 else f'{f:g} Å'}</b></td>{''.join(tds)}</tr>")
    head = ["focal spread ↓ · source →", *[("no source" if s == 0 else f"{s:g} Å") for s in srcs]]
    return arm.table(head, trs, "Rows: focal spread; columns: source size. The top row and the left column are step 1, one "
                     "effect alone; the four inner cells are step 2A, both at once. Each cell: the verdict, the share of Pb and Ti "
                     "atoms at the right depth, what guessing gives, and how coherent the beam is (strongest mode). True probe "
                     "given and held fixed throughout.")


def modes_table(recs):
    head = ["setting", "modes", "probe fixed (step 1): Pb / Ti (%)", "verdict", "probe modes: Pb / Ti (%)", "verdict",
            "strongest mode: truth / start / recovered (%)", "vacuum / slab: fixed → modes"]
    trs = []
    for r in [x for x in recs.values() if x["spec"]["group"] == "modes"]:
        c = recs[r["spec"]["control"]]
        cpart = f"<td>{pct(c['sc'], 82):.0f} / {pct(c['sc'], 22):.0f}</td>{arm.vcell(c['v'])}" if c["found"] else "<td>–</td><td>–</td>"
        k = r["kick"] or {}
        st, rec_ = k.get("start_mode_power") or [], k.get("mode_power") or []
        mp = (f"{100 * r['w1']:.0f} / {100 * st[0]:.0f} / {100 * rec_[0]:.0f}" if st and rec_ else f"{100 * r['w1']:.0f} / – / –")
        vac = f"{c['vac']:.2f} → {r['vac']:.2f}" if r["vac"] is not None and c["vac"] is not None else "–"
        if r["found"]:
            tag = (" <span class='chip c-warn'>provisional</span>" if r["provisional"] else
                   " <span class='chip c-none'>earlier solver</span>" if r["state"] == "earlier" else "")
            mpart = f"<td>{pct(r['sc'], 82):.0f} / {pct(r['sc'], 22):.0f}{tag}</td>{arm.vcell(r['v'])}"
        else:
            mpart = (f"<td colspan='2' class='wrap'><span class='chip c-warn'>{'no atoms yet' if r['state'] == 'volume' else 'pending'}</span> "
                     f"<span class='muted-s'>{html.escape(r['spec'].get('pending', ''))}</span></td>")
        trs.append(f"<tr><td>{html.escape(setting(r))}</td><td>{r['spec']['modes']}</td>{cpart}{mpart}<td>{mp}</td><td>{vac}</td></tr>")
    return arm.table(head, trs, "No shot noise, no atomic vibration. Fixed: the true probe given and held. Modes: the beam "
                     "refined from the data as 7–8 modes. Strongest mode, three numbers: <b>truth</b>, the beam that made the data "
                     "(the simulator's 16 configurations); <b>start</b>, what the reconstruction was given (the mixed beam a real gun "
                     "with this setting would have); <b>recovered</b>, what it ended with after 250 iterations. A recovered number "
                     "near the truth means the reconstruction kept the beam mixed rather than collapsing it to one probe. The finder's "
                     "templates for a lone Pb and a lone Ti atom are reconstructed the same way as each run.")


def instruments_table(man, recs):
    trs = []
    for s in man["instruments"]:
        a = focal_spread(s["dE"][0], s["cc"][0]); b = focal_spread(s["dE"][1], s["cc"][1])
        rng = f"{a:.1f}" if abs(a - b) < 0.05 else f"{a:.1f}–{b:.1f}"
        trs.append(f"<tr><td class='wrap'>{html.escape(s['name'])}</td><td>{rng}</td><td class='wrap muted'>{html.escape(s['src'])}</td></tr>")
    return arm.table(["instrument", "focal spread (Å rms)", "energy spread and Cc: where they come from"], trs,
                     f"Focal spread = Cc · (σ<sub>E</sub>/E) · (1 + E/E<sub>0</sub>)/(1 + E/2E<sub>0</sub>), σ<sub>E</sub> = FWHM/2.355: "
                     f"{FS_PER_EV_MM:.2f} Å per eV·mm at 300 kV. High-tension ripple and objective-current noise add in "
                     f"quadrature, so every range is a lower bound.")


def renders(recs):
    panes = []
    for r in recs.values():
        p = r["info"]["found"] and arm._RUNS[0].render(r["leg"])
        if p:
            panes.append((r["spec"]["short"], arm.img(p, r["spec"]["name"], maxw=1100, jpeg=78)))
    return arm.tabs("tabs-renders", panes) if panes else "<p>No renders pulled yet.</p>"


def repro():
    rows = [("Step 1: ten single-effect points, probe fixed (sims kept)", "bash aberration_experiment/mono300/phase1_fixed_submit.sh"),
            ("Step 2: joint points and probe-mode legs", "bash aberration_experiment/mono300/phase1_step2_submit.sh"),
            ("Step 2: the four mode legs that stopped at a NaN, resubmitted unchanged", "bash aberration_experiment/mono300/phase1_step2_retry.sh"),
            ("One point by hand (here 2.5 Å focal spread, probe fixed)",
             'RTAG=_fixed LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv GROUPING="32;16" PACK_H5=0 FOCAL_SPREAD=2.5 bash campaign/run_thin_atomfind.sh'),
            ("Its probe-mode leg from the mixed state",
             'RECON_ONLY=1 RTAG=_mix LABELS=armf_a080 TSV=campaign/ceos_sweep.tsv FOCAL_SPREAD=5 KICK_ALL=1 KICK_START=mixed '
             'KICK_MODES=7 KICK_NITER=250 KICK_GROUPING="4;2" RMEM=96G PACK_H5=0 bash campaign/run_thin_atomfind.sh'),
            ("Mode spectra of every point", "python aberration_experiment/mono300/mode_spectrum.py"),
            ("This page", "~/hyperspy-bundle/bin/python analysis/make_mono_page.py")]
    trs = [f"<tr><td class='wrap'>{html.escape(a)}</td><td class='wrap'><code>{html.escape(b)}</code></td></tr>" for a, b in rows]
    return arm.table(["what", "command (repository root; the submit scripts run on Blythe)"], trs,
                     "Partial coherence in the simulator (sim/simulate_4dstem.py): 16 configurations, each with its own defocus "
                     "offset (stratified Gaussian quantiles, rescaled to exactly the asked rms) and its own probe-position offset "
                     "(the source), seed 7; the diffraction patterns are averaged. The recorded positions stay nominal.")


# ====================================================================================== main
CUTS = [None]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="only print what has landed and what is pending")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    man = json.load(open(MANIFEST))
    R = arm.Runs(man); arm._RUNS[0] = R
    print(f"searching {len(R.dirs)} folders, {len(R.analyses)} analyses")
    for k in man["runs"]:
        f = R.find(k)
        print(f"  {'IN ' if f['found'] and f['primary'] else ('part' if f['found'] else 'pending'):7s}  {k:15s} {f['tag'] or man['runs'][k]['tag'][0]}")
    if a.list:
        return 0
    with contextlib.redirect_stdout(io.StringIO()):
        recs = {k: record(R, k) for k in man["runs"]}
    CUTS[0] = arm.cuts_from(arm.load_leg(recs["known_a080"]["info"]))
    focal = ["known_a080", "fs0p75", "fs1p5", "fs2p5", "fs3p5", "fs5", "fs7p29", "cold_feg"]
    source = ["known_a080", "src0p05", "src0p1", "src0p2", "src0p4", "cold_feg"]
    joint = ["fs2p5", "fs2p5_src0p1", "fs2p5_src0p2", "fs3p5", "fs3p5_src0p1", "fs3p5_src0p2"]
    mix = [k for k in man["runs"] if man["runs"][k]["group"] == "modes"]
    short = lambda r: (r["spec"]["short"] if r["key"] not in ("known_a080", "cold_feg") else
                       ("coherent" if r["key"] == "known_a080" else "cold FEG\n7.29 Å + 0.4 Å"))
    frag = dict(BUILT=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    frag.update(headline(recs))
    frag["FRAG_STATUS"] = status(recs, man)
    frag["FIG_WALL"] = arm.img(fig_wall(recs, man), "atoms in the right place against focal spread and source size")
    ref_row = lambda k: (recs[k], "coherent beam (reference)" if k == "known_a080" else "ARM200F cold FEG: 7.29 + 0.4 Å", "ref")
    frag["FRAG_FOCAL_TABLE"] = judged_table([ref_row("known_a080")] + [(recs[k], setting(recs[k]), "") for k in focal[1:-1]]
                                            + [ref_row("cold_feg")], "No shot noise, no atomic vibration; the reconstruction is given "
                                            "the true probe and holds it fixed. " + COLS_NOTE)
    frag["FIG_FOCAL_STRIP"] = arm.img(fig_strip("fig_strip_focal.png", focal, recs, CUTS[0], short), "depth sections, focal spread")
    frag["FRAG_FOCAL_TABS"] = series_tabs("tabs-focal", focal[1:-1], recs, CUTS[0])
    frag["FRAG_SOURCE_TABLE"] = judged_table([ref_row("known_a080")] + [(recs[k], setting(recs[k]), "") for k in source[1:-1]]
                                             + [ref_row("cold_feg")], "As in §02: no shot noise, no atomic "
                                             "vibration, the true probe given and held fixed. " + COLS_NOTE)
    frag["FIG_SOURCE_STRIP"] = arm.img(fig_strip("fig_strip_source.png", source, recs, CUTS[0], short), "depth sections, source size")
    frag["FRAG_SOURCE_TABS"] = series_tabs("tabs-source", source[1:-1], recs, CUTS[0])
    frag["FRAG_JOINT_GRID"] = joint_grid(recs)
    jl = [k for k in joint if recs[k]["spec"]["group"] == "joint" and recs[k]["found"]]
    frag["FRAG_JOINT_FIG"] = (f'<figure class="wide">{arm.img(fig_strip("fig_strip_joint.png", joint, recs, CUTS[0], short), "depth sections, both at once")}'
                              f'<figcaption><span class="fig-name">fig_strip_joint.png</span> Read as the strips above; each '
                              f'focal spread first alone, then with a 0.1 and a 0.2 Å source.</figcaption></figure>'
                              + series_tabs("tabs-joint", jl, recs, CUTS[0])) if jl else ""
    frag["FRAG_MODES_TABLE"] = modes_table(recs)
    ml = [r for r in (recs[k] for k in mix) if r["vol"]]
    fm = fig_modes([recs[k] for k in mix], recs)
    pairs = [k for k2 in mix if recs[k2]["vol"] for k in (recs[k2]["spec"]["control"], k2)]
    frag["FRAG_MODES_FIG"] = ((f'<figure class="wide">{arm.img(fm, "mode powers, start against recovered")}<figcaption>'
                               f'<span class="fig-name">fig_modes.png</span><div class="howto"><b>How to read it</b> Grey: the '
                               f'share of the beam in each mode at the start (a real gun with this setting). Black: what the '
                               f'reconstruction ended with. Dashed: the strongest mode of the state that made the data. Black bars '
                               f'piling onto mode 1 mean the reconstruction is leaving the mixed state.</div></figcaption></figure>'
                               if fm else "")
                              + (f'<figure class="wide">{arm.img(fig_strip("fig_strip_modes.png", pairs, recs, CUTS[0], short), "fixed probe against probe modes")}'
                                 f'<figcaption><span class="fig-name">fig_strip_modes.png</span> Each setting twice: the true probe '
                                 f'fixed (step 1), then the probe modes started from the mixed state.</figcaption></figure>'
                                 + series_tabs("tabs-modes", [r["key"] for r in ml], recs, CUTS[0]) if ml else ""))
    frag["FRAG_P2_TABLE"] = p2_table(recs)
    p2 = [k for k in P2_KEYS if recs[k]["vol"]]
    p2short = lambda r: {"fs2p5_src0p1": "2.5 + 0.1\nfixed", "fs2p5_src0p2": "2.5 + 0.2\nfixed"}.get(r["key"], r["spec"]["short"].replace(", ", "\n"))
    frag["FIG_P2_STRIP"] = arm.img(fig_strip("fig_strip_p2.png", p2, recs, CUTS[0], p2short), "depth sections, the Spectra Mono point")
    frag["FRAG_P2_TABS"] = series_tabs("tabs-p2", p2, recs, CUTS[0])
    frag["FRAG_SEARCH_TABLE"] = search_table(R, recs)
    frag["FIG_INSTRUMENTS"] = arm.img(fig_instruments(recs, man), "instruments on the focal-spread wall")
    frag["FIG_WHERE"] = arm.img(fig_where(recs, man), "the combined runs placed under the Spectra 300 monochromators")
    frag["FRAG_INSTR_TABLE"] = instruments_table(man, recs)
    frag["FRAG_RENDERS"] = renders(recs)
    frag["FRAG_REPRO"] = repro()
    css = re.search(r"<style>(.*?)</style>", open(arm.TEMPLATE).read(), re.S).group(1)
    frag["ARM_CSS"] = css
    for k in frag:
        if k.startswith("FRAG_"):
            frag[k] = re.sub(r"<caption>(.*?)</caption></table></div>", r'</table></div><p class="tnote">\1</p>', frag[k], flags=re.S)
    page = open(TEMPLATE).read()
    for k, v in frag.items():
        if "{{%s}}" % k not in page:
            sys.exit(f"template has no placeholder {{{{{k}}}}}")
        page = page.replace("{{%s}}" % k, v)
    left = re.findall(r"\{\{([A-Z_]+)\}\}", page)
    if left:
        sys.exit("unfilled placeholders: " + ", ".join(left))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    open(a.out, "w").write(page)
    mb = os.path.getsize(a.out) / 1e6
    print(f"wrote {a.out} ({mb:.1f} MB{' -- OVER the 16 MB page limit' if mb > 16 else ''})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
