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
    if provisional:
        info = dict(info, primary=True)           # its atoms are the point of the tag: score them, labelled provisional
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
        r["state"] = ("provisional" if provisional else "in") if r["found"] else ("volume" if r["vol"] else "crashed")
    return r


STATE_NOTE = {"volume": "reconstruction in, no atoms yet", "crashed": "crashed before a reconstruction",
              "provisional": "provisional: atoms read with the fixed-probe control's kernels"}


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
            ax.plot(x, [pct(r["sc"], z, "guess") for r in pts], "--", color=col, lw=1.1, alpha=.65, label=f"{nm}, guessing", zorder=3)
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
    axs[1].text(0.415, 62, "real sources\nat the specimen\n(CEOS 2006)", fontsize=8.5, color=arm.INK2, va="center",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none"))
    axs[1].set_xlim(-0.02, rs[1] + 0.02)
    axs[0].set_xlim(-0.2, 7.6)
    top = axs[0].secondary_xaxis("top", functions=(lambda f: f / (FS_PER_EV_MM * 1.4), lambda e: e * FS_PER_EV_MM * 1.4))
    top.set_xlabel("energy spread through Cc 1.4 mm (eV FWHM)", fontsize=9.5)
    axs[0].set_ylabel("bulk atoms in the right place (%)")
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
                        f"\n{r['v']}{' (provisional)' if r['provisional'] else ''}") if sc else "\nno atoms yet\n\n"
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
    def wall(grp, key):
        pts = sorted([recs["known_a080"]] + [r for r in recs.values() if r["spec"]["group"] == grp and r["found"]], key=lambda r: r[key])
        g = max(r[key] for r in pts if r["v"] == "good")
        o = min(r[key] for r in pts if r["v"] == "off" and r[key] > g)
        mid = [r for r in pts if g < r[key] < o]
        return g, o, mid
    fg, fo, _ = wall("focal", "fs")
    sg, so, smid = wall("source", "src")
    dE14, dE24 = fg / (FS_PER_EV_MM * 1.4), fg / (FS_PER_EV_MM * 2.4)
    cold = recs["cold_feg"]
    alone = {(r["fs"], r["src"]): r for r in recs.values() if r["spec"]["group"] in ("focal", "source") and r["found"]}
    joint = [r for r in recs.values() if r["spec"]["group"] == "joint" and r["found"]]
    worse = [r for r in joint if r["v"] == "off" and alone.get((r["fs"], 0.0), {}).get("v") in ("good", "okay")
             and alone.get((0.0, r["src"]), {}).get("v") in ("good", "okay")]
    jgood = sorted([r for r in joint if r["v"] == "good"], key=lambda r: (r["fs"], r["src"]))
    lifted = [r for r in recs.values() if r["spec"]["group"] == "modes" and r["found"] and not r["provisional"]
              and ORDER[r["v"]] > ORDER[recs[r["spec"]["control"]]["v"]]]
    extra = []
    for r in worse[:1]:
        a, b = alone[(r["fs"], 0.0)], alone[(0.0, r["src"])]
        extra.append((f"{r['fs']:g} Å + {r['src']:g} Å: off", f"both at once, though {r['fs']:g} Å alone is {a['v']} and {r['src']:g} Å "
                      f"alone {b['v']}: Pb {pct(r['sc'], 82):.0f} %, Ti {pct(r['sc'], 22):.0f} %. The losses compound."))
    for r in lifted:
        c = recs[r["spec"]["control"]]
        extra.append((f"{setting(r)}: {c['v']} → {r['v']}", f"with probe modes: Pb {pct(c['sc'], 82):.0f} → {pct(r['sc'], 82):.0f} %, "
                      f"vacuum / slab {c['vac']:.2f} → {r['vac']:.2f}; Ti {pct(r['sc'], 22):.0f} % (good needs 45)"))
    tiles = [(f"{fg:g} → {fo:g} Å", "focal-spread wall (rms): the last good point and the first off"),
             (f"{sg:g} → {so:g} Å", "source-size wall (FWHM)" + "".join(
                 f"; {m['src']:g} Å is {m['v']}: Pb {pct(m['sc'], 82):.0f} %, Ti {pct(m['sc'], 22):.0f} %" for m in smid)),
             (f"≤ {dE14:.2f} eV", f"energy spread for good at Cc 1.4 mm (≤ {dE24:.2f} eV at 2.4 mm)"),
             ("both at once", f"the cold FEG ({cold['fs']:g} Å + {cold['src']:g} Å) failed on each count: Pb {pct(cold['sc'], 82):.0f} % "
                              f"against {pct(cold['sc'], 82, 'guess'):.0f} % by guessing")]
    tiles = tiles[:3] + extra + tiles[3:]
    grid = "".join(f"<div><b>{html.escape(a)}</b><span>{html.escape(b)}</span></div>" for a, b in tiles)
    line = f"alone: focal spread ≤ <em>{fg:g} Å</em>, source ≤ <em>{sg:g} Å</em>"
    if jgood:
        line += f" · together: <em>{jgood[-1]['fs']:g} Å + {jgood[-1]['src']:g} Å</em>"
    return dict(WALL_LINE=line, FRAG_HEAD_GRID=f'<div class="spec-grid">{grid}</div>')


ORDER = {"off": 0, "okay": 1, "good": 2}


def status(recs, man):
    groups = [("Step 1: one effect at a time, true probe fixed", "focal source ref"),
              ("Step 2A: both at once, true probe fixed", "joint"),
              ("Step 2B: probe modes started from the mixed state", "modes")]
    out = []
    for name, gs in groups:
        rs = [r for r in recs.values() if r["spec"]["group"] in gs.split()]
        pend = [r for r in rs if r["state"] != "in"]
        ok = len(rs) - len(pend)
        chip = '<span class="chip c-ok">all in</span>' if not pend else (
            f'<span class="chip c-warn">{len(pend)} pending</span>' if ok else '<span class="chip c-warn">running</span>')
        det = "".join(f"<li><b>{html.escape(r['spec']['name'])}</b>"
                      f"{(' (' + STATE_NOTE[r['state']] + ')') if r['state'] in STATE_NOTE else ''}: "
                      f"{html.escape(r['spec'].get('pending', 'not landed yet'))}</li>" for r in pend)
        out.append(f'<div><span class="lg-t">{html.escape(name)}</span>{chip}<span class="lg-d">{ok} of {len(rs)} in'
                   f'{f"<ul>{det}</ul>" if det else ""}</span></div>')
    out.append(f'<div><span class="lg-t">Phase 2: temperature</span><span class="chip c-none">not designed</span>'
               f'<span class="lg-d">{html.escape(man["phase2"])}</span></div>')
    return f'<div class="ledger">{"".join(out)}</div>'


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
    return arm.table(head, trs, "Edges: step 1 (one effect at a time). Inside: step 2A, both at once. Bulk atoms in the right "
                     "place; verdict as in §02. The strongest-mode share is computed for every cell, landed or not.")


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
            tag = " <span class='chip c-warn'>provisional</span>" if r["provisional"] else ""
            mpart = f"<td>{pct(r['sc'], 82):.0f} / {pct(r['sc'], 22):.0f}{tag}</td>{arm.vcell(r['v'])}"
        else:
            mpart = (f"<td colspan='2' class='wrap'><span class='chip c-warn'>{'no atoms yet' if r['state'] == 'volume' else 'pending'}</span> "
                     f"<span class='muted-s'>{html.escape(r['spec'].get('pending', ''))}</span></td>")
        trs.append(f"<tr><td>{html.escape(setting(r))}</td><td>{r['spec']['modes']}</td>{cpart}{mpart}<td>{mp}</td><td>{vac}</td></tr>")
    return arm.table(head, trs, "Noiseless, no phonons, bulk atoms in the right place. Kernels from Pb and Ti grids reconstructed "
                     "the same way as their own leg, except where marked provisional (the fixed-probe control's kernels). "
                     "Strongest mode: truth, the data's own state (the simulator's 16 configurations); start, the state of a real "
                     "gun with this setting, which the leg was given; recovered, after 250 iterations.")


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
    ref_row = lambda k: (recs[k], "coherent (reference)" if k == "known_a080" else "the cold FEG: 7.29 Å + 0.4 Å source", "ref")
    frag["FRAG_FOCAL_TABLE"] = judged_table([ref_row("known_a080")] + [(recs[k], setting(recs[k]), "") for k in focal[1:-1]]
                                            + [ref_row("cold_feg")], "Noiseless, no phonons, true coherent probe fixed; kernels with "
                                            "the same physics. Vacuum / slab: rms phase per slice in the vacuum bands over that in the slab.")
    frag["FIG_FOCAL_STRIP"] = arm.img(fig_strip("fig_strip_focal.png", focal, recs, CUTS[0], short), "depth sections, focal spread")
    frag["FRAG_FOCAL_TABS"] = series_tabs("tabs-focal", focal[1:-1], recs, CUTS[0])
    frag["FRAG_SOURCE_TABLE"] = judged_table([ref_row("known_a080")] + [(recs[k], setting(recs[k]), "") for k in source[1:-1]]
                                             + [ref_row("cold_feg")], "As above.")
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
    frag["FIG_INSTRUMENTS"] = arm.img(fig_instruments(recs, man), "instruments on the focal-spread wall")
    frag["FRAG_INSTR_TABLE"] = instruments_table(man, recs)
    frag["FRAG_RENDERS"] = renders(recs)
    frag["FRAG_REPRO"] = repro()
    css = re.search(r"<style>(.*?)</style>", open(arm.TEMPLATE).read(), re.S).group(1)
    frag["ARM_CSS"] = css
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
