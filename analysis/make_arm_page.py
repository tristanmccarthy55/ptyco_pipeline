#!/usr/bin/env python
"""
@file make_arm_page.py
@brief The ARM200F results page: every figure, judged table and tab viewer from the run outputs, then the page.

One command rebuilds the whole page from whatever has been pulled. Runs are listed in
aberration_experiment/page/arm_results.json (the manifest): each is found by its TAG (an analysis tarball's summary.csv
label, figdata/<tag>, atomfind_<tag>) and its LEG (recon_af_<leg>_NL<n>: phase render, probe record, error trace) in
the manifest's search_dirs. A run that has not landed is drawn as a marked placeholder carrying its 'pending' text,
so the page can be published before every run is in and filled in later without touching this script.

Figures follow the campaign's conventions: one colour per species everywhere (Pb blue, Ti orange, O green), real atoms
as dots and found atoms as small rings, the vacuum edges dashed, no number typed in (all read from the outputs, the
sweep table or the measured tableau). Comparison tables are judged good / okay / off with a colour AND a word, by the
rule in verdict().

    ~/hyperspy-bundle/bin/python analysis/make_arm_page.py            # figures + page -> page/built/arm_page_built.html
    ~/hyperspy-bundle/bin/python analysis/make_arm_page.py --list     # what has landed, what is pending (no figures)

Publishing (Artifact tool): republish the built file to the SAME url recorded in aberration_experiment/page/ARM_PAGE.md.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import glob
import html
import io
import json
import os
import re
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
PAGE = os.path.join(REPO, "aberration_experiment", "page")
MANIFEST = os.path.join(PAGE, "arm_results.json")
TEMPLATE = os.path.join(PAGE, "arm_2026-10-09.html")
FIGS = os.path.join(REPO, "aberration_experiment", "figs", "arm_page")          # the PNGs (print copies)
CACHE = os.path.join(REPO, "aberration_experiment", "results", "arm_page_cache")  # region GT caches (gitignored)
OUT = os.path.join(PAGE, "built", "arm_page_built.html")                         # the page to publish (gitignored)
TSV = os.path.join(REPO, "campaign", "ceos_sweep.tsv")
TABLEAU = os.path.join(REPO, "campaign", "arm200f_tableau.tsv")
PLAN = [os.path.join(REPO, "aberration_experiment", "results", "2026-W40", "block2_plan.json"),
        os.path.join(REPO, "aberration_experiment", "results", "2026-W40", "plan_runs", "armf_a100_sizes.json")]
KICK = 1.05                    # run_thin_atomfind.sh KICK default: C1, C3, A1 x 1.05 (the kicked start)

import importlib.util


def _load_mod(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


cf = _load_mod("make_ceos_figs", os.path.join(HERE, "make_ceos_figs.py"))      # insitu, xz_small, clean, rows
sfig = _load_mod("make_simple_figs", os.path.join(HERE, "make_simple_figs.py"))  # species colours, box geometry
mep = _load_mod("make_mep_volumes_fig", os.path.join(HERE, "make_mep_volumes_fig.py"))
bp = _load_mod("build_page", os.path.join(PAGE, "build_page.py"))             # encode(): PNG/JPEG data URIs
import figdata                                                                   # noqa: E402
cf.mmf.style()
INK, INK2, MUTED = cf.INK, cf.INK2, cf.MUTED
SP = sfig.SP                                                                     # {82: (Pb, blue), 22: (Ti, ...), 8: (O, ...)}
LINES = {"known": INK, "k200": "#9a9a92", "k500": "#a8438a", "m6": "#5b4b8a"}    # never a species colour
LAM = cf.LAM


# =========================================================================================== finding the runs
class Runs:
    """Every analysis and result folder under the manifest's search_dirs, and lookups by tag / leg."""

    def __init__(self, man):
        self.man = man
        self.dirs = [d for d in (os.path.expanduser(x) for x in man["search_dirs"]) if os.path.isdir(d)]
        an = [a for d in self.dirs for a in glob.glob(os.path.join(d, "analysis", "analysis_*")) if os.path.isdir(a)]
        self.analyses = sorted(set(an), key=os.path.getmtime)                       # newest last
        self.summary = {}
        for a in self.analyses:                                                     # a newer analysis wins
            p = os.path.join(a, "summary.csv")
            if os.path.exists(p):
                for r in csv.DictReader(open(p)):
                    self.summary[r["label"]] = (r, a)

    def _newest(self, paths):
        paths = [p for p in paths if os.path.exists(p)]
        return max(paths, key=os.path.getmtime) if paths else None

    def find(self, key):
        """dict(found, tag, row, fd, af, primary) for a manifest run; found=False when no tag has landed."""
        spec = self.man["runs"][key]
        for i, tag in enumerate(spec.get("tag", [])):
            row, an = self.summary.get(tag, (None, None))
            fd = self._newest([os.path.join(a, "figdata", tag) for a in self.analyses])
            af = self._newest([os.path.join(a, f"atomfind_{tag}") for a in self.analyses
                               if os.path.exists(os.path.join(a, f"atomfind_{tag}", "found_atoms.npy"))])
            if row is not None or fd is not None:
                return dict(found=True, key=key, tag=tag, row=row, fd=fd, af=af, primary=(i == 0), spec=spec)
        return dict(found=False, key=key, tag=None, row=None, fd=None, af=None, primary=False, spec=spec)

    def render(self, leg):
        return self._newest([p for a in self.analyses for p in glob.glob(os.path.join(a, "phase", f"recon_af_{leg}_NL*.png"))])

    def kick(self, leg):
        j = self._newest([p for a in self.analyses for p in glob.glob(os.path.join(a, "kick", f"recon_af_{leg}_NL*.json"))])
        return (json.load(open(j)), j[:-5] + ".png") if j else (None, None)

    def trace(self, leg):
        t = self._newest([p for d in self.dirs for p in glob.glob(os.path.join(d, "results", f"recon_af_{leg}_NL*", "analysis",
                                                                               "**", "*_error_trace.csv"), recursive=True)])
        if not t:
            return None
        it, e = [], []
        for l in open(t):
            if l[:1].isdigit():
                a, b = l.strip().split(",")[:2]; it.append(int(a)); e.append(float(b))
        n = int(re.search(r"(\d{3,4})x\d{3,4}", os.path.basename(t)).group(1))
        return np.array(it), np.array(e) * n / 2e5                                  # eps = error x N / 2e5

    def search(self, name):
        for d in self.dirs:
            for p in (os.path.join(d, name), os.path.join(d, "results", name)):
                if os.path.exists(os.path.join(p, "search_best.json")):
                    return p
        return None


def num(r, k, scale=1.0):
    try:
        return float(r[k]) * scale if r and r.get(k) not in (None, "", "nan") else None
    except ValueError:
        return None


def verdict(r):
    """good / okay / off from the atom finder's numbers (stated on the page):
    good  Pb >= 95 %, Ti and O >= 60 %, depth error <= 0.6 A, precision >= 0.9  (the known-probe 80 mrad standard)
    okay  Pb >= 80 %, O >= 25 %, depth error <= 1.1 A  (lattice and Pb right; species, depth or false hits weaker)
    off   anything else.  None when the finder did not run."""
    pb, ti, o, z, pr = (num(r, "Pb_recall_bulk", 100), num(r, "Ti_recall_bulk", 100), num(r, "O_recall_bulk", 100),
                        num(r, "z_rms"), num(r, "precision"))
    if None in (pb, ti, o, z, pr):
        return None
    if pb >= 95 and ti >= 60 and o >= 60 and z <= 0.6 and pr >= 0.9:
        return "good"
    if pb >= 80 and o >= 25 and z <= 1.1:
        return "okay"
    return "off"


# =========================================================================================== shared geometry
def gt_dir(region):
    """Region ground truth (atomfind.make_gt_cache, the simulator's own builder), cached in the repo (gitignored)."""
    d = os.path.join(CACHE, f"gt_region{region}")
    if not os.path.exists(os.path.join(d, "gt_prepared.npz")):
        os.makedirs(d, exist_ok=True)
        subprocess.run([sys.executable, "-m", "atomfind.make_gt_cache", "--thin-cells", "5", "--z-vacuum", "4",
                        "--region-side", str(region), "--out", os.path.join(d, "gt_prepared.npz")], cwd=HERE, check=True,
                       stdout=subprocess.DEVNULL)
    return d


_LEGS = {}


def load_leg(info):
    """Volume (finder's frame), ground truth and found atoms, as make_ceos_figs.leg does, from a run's figure data."""
    if not info["found"] or not info["fd"]:
        return None
    if info["fd"] in _LEGS:
        return _LEGS[info["fd"]]
    from atomfind import align
    V, meta = figdata.load(info["fd"])
    dx, dz, centre = meta["dx"], meta["dz"], meta["centre"]
    region = int(round(2 * centre[0]))                       # the scan sits at the region's centre
    cfg = figdata.config_for(os.path.join(info["fd"], "phase_vol.npy"), gt_dir(region), dz, dx, centre)
    pos, Z = align.load_gt(cfg)
    found = gtidx = None
    if info["af"]:
        found = np.load(os.path.join(info["af"], "found_atoms.npy"))
        al = align.refine_with_atoms(align.register(V, dx, pos, Z, cfg), found, pos, Z, cfg)
        gtidx = al.site_to_index(pos[:, 0], pos[:, 1], pos[:, 2])
    lg = dict(V=V, dx=dx, dz=dz, pos=pos, Z=Z, found=found, gt=gtidx, cfg=cfg)
    _LEGS[info["fd"]] = lg
    return lg


def cuts_from(lg):
    """The two depth sections every figure cuts, in A from the crop's top, from the 80 mrad known-probe run: a Pb row
    (Pb every 3.9 A down the column) and a B-site row (Ti and O alternating every 1.95 A)."""
    ny = lg["V"].shape[1]
    out = []
    for name, sp in (("Pb row", 82), ("B-site row", 22)):
        r = mep.pick_row(lg["found"], sp, ny, lg["dx"]) if lg["found"] is not None else None
        out.append((name, (ny / 2 if r is None else r) * lg["dx"]))
    return out


def numbers_line(r):
    if r is None:
        return ""
    pb, ti, o, z = num(r, "Pb_recall_bulk", 100), num(r, "Ti_recall_bulk", 100), num(r, "O_recall_bulk", 100), num(r, "z_rms")
    if None in (pb, ti, o, z):
        return "atoms: finder not run"
    pr = num(r, "precision")
    return f"Pb {pb:.0f} / Ti {ti:.0f} / O {o:.0f} %, depth {z:.2f} Å" + (f", precision {pr:.2f}" if pr is not None else "")


def run_row(fig, gs, i, info, title, cuts, seen, ncols=3):
    """One run across a figure row: depth-summed phase (central 20 A, cuts marked) and the two depth sections."""
    lg = load_leg(info)
    if lg is None:                                            # one panel across the row, the reason written out
        import textwrap
        ax = fig.add_subplot(gs[i, :])
        txt = ("Pending\n\n" + textwrap.fill(info["spec"].get("pending", "not landed yet"), 70)) if not info["found"] else \
              "No volume was shipped for this run.\nIts full reconstruction render is in the appendix."
        cf.placeholder(ax, ""); ax.text(0.5, 0.5, txt, ha="center", va="center", fontsize=11, color=INK2, transform=ax.transAxes)
        ax.set_title(title, fontsize=10.5, loc="left")
        return
    ax = fig.add_subplot(gs[i, 0]); cf.insitu(ax, lg)
    top = lg["V"].shape[1] * lg["dx"] / 2
    for cut, rA in cuts:
        ax.axhline(rA - top, color="white", lw=0.8, ls=(0, (4, 3)), alpha=0.8)
    ax.plot([-9.2, -4.2], [8.9, 8.9], color="white", lw=3, solid_capstyle="butt")      # 5 A scale bar
    ax.text(-6.7, 8.4, "5 Å", color="white", ha="center", va="bottom", fontsize=9)
    nums = numbers_line(info["row"] if info["primary"] else None)
    ax.set_title(title + (f"\n{nums}" if nums else ""), fontsize=10.5, loc="left")
    nx = lg["V"].shape[2]
    for j, (cut, rA) in enumerate(cuts):
        ax = fig.add_subplot(gs[i, 1 + j])
        seen |= cf.xz_small(ax, lg, rA / lg["dx"], sfig, xlim=(nx * lg["dx"] / 2 - 8, nx * lg["dx"] / 2 + 8), show_y=(j == 0))
        ax.set_title(f"{cut}  ·  {lg['V'].shape[0]} slices of {lg['dz']:.2f} Å" if j == 0 else cut, fontsize=9.5, loc="left")
        ax.set_xlabel("x (Å)", fontsize=9)


def legend_handles(seen):
    h = [Line2D([], [], ls="none", marker=".", color=c, ms=11, label=f"{n} in the structure") for z, (n, c) in SP.items() if z in seen]
    h += [Line2D([], [], ls="none", marker="o", mfc="none", mec=INK2, mew=1.3, ms=8, label="atom the finder reported"),
          Line2D([], [], color="#7fd4ff", ls=(0, (4, 3)), lw=1.2, label="vacuum edge"),
          Line2D([], [], color=INK2, ls=(0, (4, 3)), lw=1.0, label="where the sections cut (left image)")]
    return h


def compare_fig(name, items, cuts):
    """Rows of runs, each: depth-summed phase + Pb-row and B-site-row sections. items: [(title, info), ...]."""
    n = len(items)
    H = 4.45 * n + 1.1                                        # inches: rows plus a strip for the legend
    fig = plt.figure(figsize=(10.6, H))
    gs = fig.add_gridspec(n, 3, width_ratios=[1.35, 0.78, 0.78], wspace=0.16, hspace=0.42, top=1 - 0.55 / H, bottom=1.25 / H)
    seen = set()
    for i, (title, info) in enumerate(items):
        run_row(fig, gs, i, info, title, cuts, seen)
    fig.legend(handles=legend_handles(seen or {82, 22, 8}), loc="lower center", ncol=3, fontsize=9.5, bbox_to_anchor=(0.5, 0.0))
    return save(fig, name)


def save(fig, name):
    os.makedirs(FIGS, exist_ok=True)
    p = os.path.join(FIGS, name)
    fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"  wrote {os.path.relpath(p, REPO)}", flush=True)
    return p


def img(path, alt, maxw=1700, jpeg=None):
    """An embedded image. Figures go as PNG unless heavy (build_page.encode); tab images and renders (jpeg=quality)
    always as JPEG at a moderate width, which keeps the page well under the 16 MB limit with every run in."""
    if jpeg is None:
        uri, _ = bp.encode(path, maxw)
    else:
        import base64
        from PIL import Image
        im = Image.open(path).convert("RGB")
        if im.width > maxw:
            im = im.resize((maxw, round(im.height * maxw / im.width)), Image.LANCZOS)
        buf = io.BytesIO(); im.save(buf, "JPEG", quality=jpeg, optimize=True, progressive=True)
        uri = "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode()
    return f'<img src="{uri}" alt="{html.escape(alt)}" loading="lazy">'


def pending_box(info, what=""):
    t = info["spec"].get("pending", "not landed yet")
    return (f'<div class="pend"><span class="chip c-warn">pending</span> <b>{html.escape(info["spec"]["name"])}</b>'
            f'{(" — " + html.escape(what)) if what else ""}<br><span>{html.escape(t)}</span></div>')


def tabs(tid, panes, default=0, unit=""):
    """A tab viewer fragment (the CEOS page's .xp classes): panes = [(button label, inner html)]."""
    default = min(default, len(panes) - 1)
    btn = "".join(f'<button type="button" role="tab" aria-selected="{str(i == default).lower()}" tabindex="{0 if i == default else -1}">'
                  f'{html.escape(lab)}</button>' for i, (lab, _) in enumerate(panes))
    body = "".join(f'<div class="tab-pane" data-i="{i}"{"" if i == default else " hidden"}>{inner}</div>' for i, (_, inner) in enumerate(panes))
    return f"""<div class="xp rx" id="{tid}">
  <div class="xp-bar"><div class="xp-tabs" role="tablist">{btn}{f'<span class="xp-unit">{unit}</span>' if unit else ''}</div></div>
  <div class="rx-img">{body}</div>
</div>
<script>
(function () {{
  var root = document.getElementById("{tid}");
  var tabs = Array.prototype.slice.call(root.querySelectorAll(".xp-tabs button"));
  var panes = Array.prototype.slice.call(root.querySelectorAll(".tab-pane"));
  function show(i) {{
    tabs.forEach(function (b, k) {{ b.setAttribute("aria-selected", k === i); b.tabIndex = k === i ? 0 : -1; }});
    panes.forEach(function (p, k) {{ p.hidden = k !== i; }});
  }}
  tabs.forEach(function (b, i) {{
    b.addEventListener("click", function () {{ show(i); }});
    b.addEventListener("keydown", function (e) {{
      var k = e.key === "ArrowRight" ? 1 : (e.key === "ArrowLeft" ? -1 : 0); if (!k) return;
      e.preventDefault(); var j = (i + k + tabs.length) % tabs.length; tabs[j].focus(); show(j); }});
  }});
}})();
</script>"""


CHIP = {"good": ("c-ok", "v-good"), "okay": ("c-warn", "v-ok"), "off": ("c-stop", "v-bad")}


def vcell(v):
    if v is None:
        return '<td><span class="chip c-none">not measured</span></td>'
    return f'<td class="{CHIP[v][1]}"><span class="chip {CHIP[v][0]}">{v}</span></td>'


def atoms_cells(r):
    pb, ti, o, z, pr = (num(r, "Pb_recall_bulk", 100), num(r, "Ti_recall_bulk", 100), num(r, "O_recall_bulk", 100),
                        num(r, "z_rms"), num(r, "precision"))
    if None in (pb, ti, o, z):
        return "<td>–</td><td>–</td><td>–</td>"
    return f"<td>{pb:.0f} / {ti:.0f} / {o:.0f}</td><td>{z:.2f}</td><td>{pr:.2f}</td>"


def table(head, rows_html, caption="", cls=""):
    th = "".join(f"<th>{h}</th>" for h in head)
    return (f'<div class="tablewrap"><table class="{cls}"><thead><tr>{th}</tr></thead><tbody>{"".join(rows_html)}</tbody>'
            f'{f"<caption>{caption}</caption>" if caption else ""}</table></div>')


# =========================================================================================== physics helpers
def tsv_rows():
    return cf.rows()


def loss_from_note(note):
    m = re.search(r"window ([\d.]+) A loses ([\d.]+) %", note or "")
    return (float(m.group(1)), float(m.group(2))) if m else (None, None)


def phase_error_waves(alpha_mrad, d_defocus, dC30, dC12, phi12):
    """rms and peak-to-valley (waves) of a probe-start error over the aperture, piston removed. d_defocus is in abTEM's
    defocus (= -C10), so dC10 = -d_defocus."""
    a = alpha_mrad / 1e3
    t = np.linspace(-a, a, 601); X, Y = np.meshgrid(t, t); R = np.hypot(X, Y); m = R <= a
    r, p = R[m], np.arctan2(Y, X)[m]
    w = (0.5 * (-d_defocus) * r ** 2 + 0.25 * dC30 * r ** 4 + 0.5 * dC12 * r ** 2 * np.cos(2 * (p - phi12))) / LAM
    w = w - w.mean()
    return float(w.std()), float(np.ptp(w))


def tableau_ci(term):
    """run3's 95 % interval of a measured term, in A (the corrector's own precision)."""
    unit = {"nm": 10.0, "um": 1e4, "mm": 1e7}
    for r in csv.DictReader((l for l in open(TABLEAU) if not l.startswith("#")), delimiter="\t"):
        if r["source"] == "run3" and r["term"] == term:
            return float(r["uncertainty"]) * unit[r["unit"]]
    return None


# =========================================================================================== the page's parts
def part_probes(R, T):
    """Probe size against aperture with each run's window, and the share of the probe outside it."""
    plan = {}
    for p in PLAN:
        for r in json.load(open(p)):
            plan[r["label"]] = r
    al, d90, d99, win, loss = [], [], [], [], []
    for lab in sorted(plan, key=lambda k: plan[k]["alpha"]):
        r = T.get(lab)
        if r is None:
            continue
        w, l = loss_from_note(r["note"])
        al.append(float(r["alpha"])); d90.append(plan[lab]["d90"]); d99.append(plan[lab]["d99"])
        win.append(float(r["side"]) / float(r["bin"]) if r["side"] not in ("-", "") else w); loss.append(l)
    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    ax.plot(al, d90, "o-", color=INK, label="probe d90 (90 % of the intensity inside)")
    ax.plot(al, d99, "s--", color=MUTED, label="probe d99")
    ax.step(al, win, where="mid", color="#a8438a", lw=2.2, label="reconstruction window used")
    for a_, w_, l_ in zip(al, win, loss):
        if l_ is not None:                                   # large losses below the line, clear of the d99 curve
            big = l_ >= 1
            ax.annotate(f"{l_:.1f} %\noutside" if big else f"{l_:.2f} %", (a_, w_), xytext=(0, -7 if big else 7),
                        textcoords="offset points", ha="center", va="top" if big else "bottom", fontsize=8.5, color="#a8438a",
                        bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.6))
    fr = T.get("armf_a090_w140")
    if fr:
        w, l = loss_from_note(fr["note"])
        ax.plot([90], [w], "D", color="#a8438a", mfc="white", ms=8, mew=1.8)
        ax.annotate(f"final run: {w:g} Å window,\n{l:.1f} % outside", (90, w), xytext=(-70, 30), textcoords="offset points",
                    ha="right", va="bottom", fontsize=8.5, color="#a8438a", arrowprops=dict(arrowstyle="-", color="#a8438a", lw=0.8))
    ax.set_yscale("log"); ax.set_xlabel("probe semi-angle α (mrad)"); ax.set_ylabel("size (Å)")
    ax.set_xticks(al); ax.legend(loc="upper left", fontsize=9)
    return save(fig, "fig_probe_window.png")


def part_known(R, T):
    keys = [k for k in R.man["runs"] if k.startswith("known_")]
    infos = {k: R.find(k) for k in keys}
    rows = {k: infos[k]["row"] for k in keys}
    # chart: atoms and depth against aperture
    al = [float(R.man["runs"][k]["short"]) for k in keys]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.9))
    for sp, key in ((82, "Pb_recall_bulk"), (22, "Ti_recall_bulk"), (8, "O_recall_bulk")):
        v = [num(rows[k], key, 100) for k in keys]
        a1.plot([a for a, x in zip(al, v) if x is not None], [x for x in v if x is not None], "o-", color=SP[sp][1], label=SP[sp][0])
    a1.set_ylim(0, 108); a1.set_xlabel("α (mrad)"); a1.set_ylabel("atoms found, bulk (%)"); a1.legend(loc="lower center", ncol=3)
    z = [num(rows[k], "z_rms") for k in keys]
    a2.plot([a for a, x in zip(al, z) if x is not None], [x for x in z if x is not None], "o-", color=INK)
    a2.set_xlabel("α (mrad)"); a2.set_ylabel("depth error of found atoms, rms (Å)"); a2.set_ylim(0, None)
    for a in (a1, a2):
        a.set_xticks(al)
        a.axvspan(85, 105, color="#c2490a", alpha=0.06, lw=0)
    for a in (a1, a2):
        a.text(95, 0.98, "window\ntoo small", ha="center", va="top", fontsize=8.5, color="#c2490a", transform=a.get_xaxis_transform())
    fchart = save(fig, "fig_known_numbers.png")
    # judged table
    trs = []
    for k in keys:
        r, info = rows[k], infos[k]
        t = T.get(info["tag"] or "", {})
        w, l = loss_from_note(t.get("note"))
        tr = R.trace(R.man["runs"][k]["leg"])
        eps = f"{tr[1][-1]:.4f}" if tr else "–"
        cls = ' class="ref"' if k == R.man["reference"] else ""
        trs.append(f"<tr{cls}><td>{R.man['runs'][k]['short']}</td><td>{float(t['side']) / float(t['bin']) if t else 0:g}</td>"
                   f"<td>{'–' if l is None else f'{l:.2f}'}</td><td>{num(r, 'NL') and int(num(r, 'NL')) or '–'}</td>"
                   f"{atoms_cells(r)}<td>{eps}</td>{vcell(verdict(r))}</tr>")
    tab = table(["α (mrad)", "window (Å)", "probe outside (%)", "slices", "Pb / Ti / O (%)", "depth (Å)", "precision", "ε", "verdict"],
                trs, "ε = final residual × pattern width / 2·10⁵ (comparable across runs of different size). The tinted row is "
                     "the reference every later section is measured against.")
    # tabs: one figure per aperture
    ref = load_leg(infos[R.man["reference"]])
    cuts = cuts_from(ref)
    panes = []
    for k in keys:
        info = infos[k]
        if not info["found"]:
            panes.append((R.man["runs"][k]["short"], pending_box(info))); continue
        p = compare_fig(f"tab_known_{k}.png", [(R.man["runs"][k]["name"], info)], cuts)
        panes.append((R.man["runs"][k]["short"], img(p, R.man["runs"][k]["name"], maxw=1300, jpeg=86)))
    return dict(FIG_KNOWN=img(fchart, "atoms found and depth error against aperture, probe known"), FRAG_KNOWN_TABLE=tab,
                FRAG_KNOWN_TABS=tabs("tabs-known", panes, default=keys.index(R.man["reference"]), unit="mrad")), cuts


def part_window(R, T, cuts):
    keys = ["win_105", "win_070", "win_052"]
    infos = {k: R.find(k) for k in keys}
    p = compare_fig("fig_window.png", [(R.man["runs"][k]["name"], infos[k]) for k in keys], cuts)
    trs = []
    for k in keys + ["win_a090", "win_a100"]:
        info = R.find(k); t = T.get(info["tag"] or "", {})
        w, l = loss_from_note(t.get("note"))
        trs.append(f"<tr><td>{html.escape(R.man['runs'][k]['short'])}</td><td>{'–' if l is None else f'{l:.1f}'}</td>"
                   f"{atoms_cells(info['row'])}{vcell(verdict(info['row']))}</tr>")
    tab = table(["run", "probe outside the window (%)", "Pb / Ti / O (%)", "depth (Å)", "precision", "verdict"], trs,
                "The share outside the window is plan_probe.window_loss for that row's probe, about the beam axis.")
    panes = []
    for k in keys + ["win_a090", "win_a100"]:
        r = R.render(R.man["runs"][k]["leg"])
        panes.append((R.man["runs"][k]["short"], img(r, R.man["runs"][k]["name"], maxw=1250, jpeg=82) if r else pending_box(R.find(k), "no render")))
    return dict(FIG_WINDOW=img(p, "80 mrad probe in three windows: depth-summed phase and depth sections"),
                FRAG_WINDOW_TABLE=tab, FRAG_WINDOW_TABS=tabs("tabs-window", panes, default=0))


def part_kicks(R, T, cuts):
    # table: the start's error in waves, how far the probe came back, the atoms
    trs = []
    for k in ["kick_a040", "kick_a060", "kick_a080", "kick_a100", "long_m3", "long_m6"]:
        info = R.find(k); spec = R.man["runs"][k]
        lab = "armf_a" + re.search(r"(\d{3})", spec.get("leg", "")).group(1) if re.search(r"a(\d{3})", spec.get("leg", "")) else None
        t = T.get(lab, {})
        ab = cf.aberrations(t) if t else {}
        rms = pv = None
        if t:
            rms, pv = phase_error_waves(float(t["alpha"]), (KICK - 1) * float(t["c1"]), (KICK - 1) * ab.get("C30", 0.0),
                                        (KICK - 1) * ab.get("C12", 0.0), ab.get("phi12", 0.0))
        kj, _ = R.kick(spec["leg"])
        ov = f"{kj['ov_start']:.2f} → {kj['ov_final']:.2f}" if kj and "ov_final" in kj else "–"
        m1 = f"{kj['mode_power'][0]:.2f}" if kj and kj.get("mode_power") else "–"
        trs.append(f"<tr><td>{html.escape(spec['name'])}</td><td>{'–' if rms is None else f'{rms:.2f} / {pv:.1f}'}</td>"
                   f"<td>{ov}</td><td>{m1}</td>{atoms_cells(info['row'])}{vcell(verdict(info['row']))}</tr>")
    ref = R.find(R.man["reference"])
    trs.insert(0, f'<tr class="ref"><td>80 mrad, probe known (reference)</td><td>0</td><td>–</td><td>–</td>{atoms_cells(ref["row"])}'
                  f'{vcell(verdict(ref["row"]))}</tr>')
    tab = table(["run", "start error, waves (rms / p-v)", "probe overlap, start → end", "mode 1 power", "Pb / Ti / O (%)",
                 "depth (Å)", "precision", "verdict"], trs,
                "Start error: the phase the kicked start gets wrong across the aperture (C1, C3, A1 × 1.05), piston removed. "
                "Overlap: |⟨start or recovered mode 1, true probe⟩|, best over shifts; it falls to ~0.2 for any smooth ~1-wave "
                "error, so read it with the start error. 40 mrad: 4 slices, too few for the atom finder.")
    # convergence: the residual against iteration
    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    for leg, lab, col, ls in (("armf_a080_lab", "probe known", LINES["known"], "--"), ("armf_a080_kick", "kick, 3 modes, 200 it", LINES["k200"], "-"),
                              ("armf_a080_lab_long_m3", "kick, 3 modes, 500 it", LINES["k500"], "-"),
                              ("armf_a080_lab_long_m6", "kick, 6 modes, 250 it", LINES["m6"], "-")):
        tr = R.trace(leg)
        if tr:
            ax.plot(tr[0], tr[1], ls=ls, color=col, label=f"{lab} (end {tr[1][-1]:.4f})")
    ax.set_xlabel("iteration (full-resolution engine)"); ax.set_ylabel("ε = residual × N / 2·10⁵")
    ax.legend(fontsize=9); ax.set_ylim(0.008, None); ax.set_yscale("log")
    fconv = save(fig, "fig_kick_convergence.png")
    pc = compare_fig("fig_kicks.png", [("80 mrad, probe known", R.find("known_a080")), ("5 % kick, 3 modes, 200 iterations", R.find("kick_a080")),
                                        ("5 % kick, 3 modes, 500 iterations", R.find("long_m3")), ("5 % kick, 6 modes, 250 iterations", R.find("long_m6"))], cuts)
    # probe-recovery tabs
    panes = []
    for k in ["kick_a040", "kick_a060", "kick_a080", "kick_a100", "long_m3", "long_m6"]:
        kj, png = R.kick(R.man["runs"][k]["leg"])
        panes.append((R.man["runs"][k]["short"], img(png, R.man["runs"][k]["name"], maxw=1100, jpeg=84) if png and os.path.exists(png) else pending_box(R.find(k))))
    return dict(FRAG_KICK_TABLE=tab, FIG_KICK_CONV=img(fconv, "residual against iteration for the known probe and three kicked runs"),
                FIG_KICKS=img(pc, "80 mrad: probe known against kicked runs"), FRAG_KICK_TABS=tabs("tabs-kick", panes, default=2, unit="mrad"))


def part_physics(R, T, cuts):
    items = [("coherent, probe known (reference)", R.find("known_a080")), ("full physics, 1e7, true probe fixed", R.find("fixed_1e7")),
             ("full physics, 1e8, true probe fixed", R.find("fixed_1e8")), ("full physics, 1e7, 5 % kick, 200 it", R.find("hail_1e7")),
             ("full physics, 1e7, 5 % kick, 500 it, no search", R.find("nosearch_1e7"))]
    p = compare_fig("fig_physics.png", items, cuts)
    trs = []
    for title, info in items:
        kj, _ = R.kick(info["spec"].get("leg", ""))
        ov = f"{kj['ov_start']:.2f} → {kj['ov_final']:.2f}" if kj and "ov_final" in kj else "–"
        m1 = " / ".join(f"{x:.2f}" for x in kj["mode_power"]) if kj and kj.get("mode_power") else "–"
        note = "" if info["primary"] or not info["found"] else " <span class='chip c-none'>volume only</span>"
        trs.append(f"<tr><td>{html.escape(title)}{note}</td><td>{ov}</td><td>{m1}</td>{atoms_cells(info['row'] if info['primary'] else None)}"
                   f"{vcell(verdict(info['row']) if info['primary'] else None)}</tr>")
    tab = table(["run", "probe overlap", "mode powers", "Pb / Ti / O (%)", "depth (Å)", "precision", "verdict"], trs,
                "Fixed-probe runs use their own matched kernels (same physics, same operator). The 200-iteration kick is read "
                "with block 2's coherent known-probe kernels (its own grids failed); the 500-iteration one with the fixed-probe "
                "1e7 kernels once analysis 1310612 is pulled.")
    return dict(FIG_PHYSICS=img(p, "full physics against the coherent reference"), FRAG_PHYSICS_TABLE=tab)


def part_search(R, T):
    t = T["armf_a080"]; ab = cf.aberrations(t)
    c1k, c3k = KICK * float(t["c1"]), KICK * ab["C30"]
    ci1, ci3 = tableau_ci("C1"), tableau_ci("C3")
    fig, axs = plt.subplots(1, 2, figsize=(11.6, 4.7))
    trs = []
    for ax, key in zip(axs, ("search_1e7", "search_1e8")):
        d = R.search(R.man["runs"][key]["search"])
        dose = R.man["runs"][key]["short"]
        if d is None:
            cf.placeholder(ax, "pending"); trs.append(f"<tr><td>{dose}</td>" + "<td>–</td>" * 6 + "</tr>"); continue
        best = json.load(open(os.path.join(d, "search_best.json")))
        rows_ = list(csv.DictReader(open(os.path.join(d, "search_trials.csv"))))
        u1 = sorted({float(r["f_c1"]) for r in rows_}); u3 = sorted({float(r["f_c3"]) for r in rows_})
        E = np.full((len(u1), len(u3)), np.nan)
        for r in rows_:
            try:
                E[u1.index(float(r["f_c1"])), u3.index(float(r["f_c3"]))] = float(r["final_error"])
            except ValueError:
                pass
        d1, d3 = np.diff(u1).mean() / 2, np.diff(u3).mean() / 2
        im = ax.imshow(E, origin="lower", cmap="viridis", aspect="auto", extent=[u3[0] - d3, u3[-1] + d3, u1[0] - d1, u1[-1] + d1])
        fig.colorbar(im, ax=ax, label="final residual of the fixed-probe trial", shrink=0.9)
        sim = best.get("sim_record_only") or {}
        if sim:
            ax.add_patch(Rectangle((sim["f_c3"] - ci3 / abs(c3k), sim["f_c1"] - ci1 / abs(c1k)), 2 * ci3 / abs(c3k), 2 * ci1 / abs(c1k),
                                   fill=False, ec="white", lw=1.4, ls=(0, (4, 3))))
            ax.plot(sim["f_c3"], sim["f_c1"], "x", color="#ff5a5a", ms=10, mew=2.2)
        ax.plot(best["f_c3"], best["f_c1"], "+", color="white", ms=15, mew=2.2)
        ax.set_xlabel("trial C3 / kicked C3"); ax.set_ylabel("trial C1 / kicked C1"); ax.grid(False)
        ax.set_title(f"{dose} e/Å²: {best['n_finished']} of {best['n_trials']} trials", fontsize=11, loc="left")
        e1, e3 = best["c1_A"] - sim.get("c1_A", np.nan), best["c3_A"] - sim.get("c3_A", np.nan)
        rms, pv = phase_error_waves(float(t["alpha"]), best["c1_A"] - sim.get("c1_A", 0), best["c3_A"] - sim.get("c3_A", 0),
                                    (KICK - 1) * ab.get("C12", 0.0), ab.get("phi12", 0.0))
        trs.append(f"<tr><td>{dose}</td><td>{best['n_finished']} / {best['n_trials']}</td><td>{best['c1_A']:.1f} ({e1:+.1f})</td>"
                   f"<td>{best['c3_A'] / 1e4:.3f} ({e3 / 1e4:+.3f})</td><td>{rms:.2f} / {pv:.1f}</td><td>{html.escape(best['status'])}</td>"
                   f"<td>{html.escape(best['method'])}</td></tr>")
    axs[1].legend(handles=[Line2D([], [], ls="none", marker="+", color=INK, ms=11, mew=2, label="chosen start"),
                           Line2D([], [], ls="none", marker="x", color="#ff5a5a", ms=9, mew=2, label="the simulation's truth"),
                           Line2D([], [], color=INK2, ls=(0, (4, 3)), lw=1.4, label="the corrector's own 95 % interval\n(run3: C1 and C3)")],
                  loc="upper left", bbox_to_anchor=(1.28, 1.0), fontsize=9)
    p = save(fig, "fig_search.png")
    kick_rms, kick_pv = phase_error_waves(float(t["alpha"]), (KICK - 1) * float(t["c1"]), (KICK - 1) * ab["C30"],
                                          (KICK - 1) * ab.get("C12", 0.0), ab.get("phi12", 0.0))
    trs.insert(0, f'<tr class="ref"><td>the kicked start</td><td>–</td><td>{c1k:.1f} ({c1k - float(t["c1"]):+.1f})</td>'
                  f'<td>{c3k / 1e4:.3f} ({(c3k - ab["C30"]) / 1e4:+.3f})</td><td>{kick_rms:.2f} / {kick_pv:.1f}</td><td>–</td><td>–</td></tr>')
    trs.append(f'<tr><td>corrector, run3 95 % interval</td><td>–</td><td>± {ci1:.0f}</td><td>± {ci3 / 1e4:.2f}</td><td>–</td><td>–</td><td>–</td></tr>')
    tab = table(["dose (e/Å²)", "trials", "C1 found, Å (error)", "C3 found, µm (error)", "start error, waves (rms / p-v)", "status", "method"],
                trs, f"True values: C1 {float(t['c1']):.1f} Å (abTEM defocus), C3 {ab['C30'] / 1e4:.2f} µm. The start keeps A1 at its "
                     "kicked value (× 1.05); the search moves C1 and C3 only.")
    return dict(FIG_SEARCH=img(p, "outer search error surfaces at 1e7 and 1e8"), FRAG_SEARCH_TABLE=tab)


def part_refined(R, T, cuts):
    items = [("full physics, 1e7, true probe fixed", R.find("fixed_1e7")), ("1e7, kicked start, no search, 500 it", R.find("nosearch_1e7")),
             ("1e7, searched start, refined, 500 it", R.find("refined_1e7")), ("1e8, searched start, refined, 500 it", R.find("refined_1e8")),
             ("1e8, TRUE start, refined, 500 it", R.find("true4_1e8"))]
    p = compare_fig("fig_refined.png", items, cuts)
    trs = []
    for title, info in items + [("1e7, searched start, its own kernels", R.find("refined_1e7_own"))]:
        if not info["found"]:
            trs.append(f"<tr><td>{html.escape(title)}</td><td colspan='6'>{pending_box(info)}</td></tr>"); continue
        kj, _ = R.kick(info["spec"].get("leg", ""))
        ov = f"{kj['ov_start']:.2f} → {kj['ov_final']:.2f}" if kj and "ov_final" in kj else "–"
        m1 = " / ".join(f"{x:.2f}" for x in kj["mode_power"]) if kj and kj.get("mode_power") else "–"
        r = info["row"] if info["primary"] else None
        trs.append(f"<tr><td>{html.escape(title)}</td><td>{ov}</td><td>{m1}</td>{atoms_cells(r)}{vcell(verdict(r))}</tr>")
    tab = table(["run", "probe overlap", "mode powers", "Pb / Ti / O (%)", "depth (Å)", "precision", "verdict"], trs,
                "All rows read with the same fixed-probe 1e7 kernels where that analysis exists (like for like); "
                "'its own kernels' uses the refined Pb/Ti grid legs.")
    panes = []
    for k in ["hail_1e7", "nosearch_1e7", "refined_1e7", "refined_1e8", "true4_1e8"]:
        kj, png = R.kick(R.man["runs"][k]["leg"])
        panes.append((R.man["runs"][k]["short"], img(png, R.man["runs"][k]["name"], maxw=1100, jpeg=84) if png and os.path.exists(png) else pending_box(R.find(k))))
    return dict(FIG_REFINED=img(p, "searched and true starts against the fixed probe and the plain kick"), FRAG_REFINED_TABLE=tab,
                FRAG_REFINED_TABS=tabs("tabs-refined", panes, default=2))


def part_final(R, T, cuts):
    items = [("90 mrad, 140 Å window, true probe fixed", R.find("final90_fixed")), ("90 mrad, 140 Å window, true start, 4 modes", R.find("final90"))]
    p = compare_fig("fig_final90.png", items, cuts)
    trs = []
    for title, info in items:
        if not info["found"]:
            trs.append(f"<tr><td>{html.escape(title)}</td><td colspan='5'>{pending_box(info)}</td></tr>"); continue
        trs.append(f"<tr><td>{html.escape(title)}</td>{atoms_cells(info['row'])}{vcell(verdict(info['row']))}<td></td></tr>")
    tab = table(["run", "Pb / Ti / O (%)", "depth (Å)", "precision", "verdict", ""], trs)
    return dict(FIG_FINAL=img(p, "90 mrad in a 140 Å window with the full physics"), FRAG_FINAL_TABLE=tab)


def part_status(R):
    groups = [("Probe known, 40–100 mrad", [k for k in R.man["runs"] if k.startswith("known_")]),
              ("Window test", ["win_070", "win_052"]), ("Coherent kicks", ["kick_a040", "kick_a060", "kick_a080", "kick_a100", "long_m3", "long_m6"]),
              ("Full physics, fixed probe and plain kick", ["fixed_1e7", "fixed_1e8", "hail_1e7", "nosearch_1e7"]),
              ("Outer search and refined probe", ["search_1e7", "search_1e8", "refined_1e7", "refined_1e7_own", "refined_1e8", "true4_1e8"]),
              ("Final run, 90 mrad in a 140 Å window", ["final90", "final90_fixed"])]
    out = []
    for name, keys in groups:
        ok, pend = 0, []
        for k in keys:
            spec = R.man["runs"][k]
            landed = R.search(spec["search"]) is not None if "search" in spec else (R.find(k)["found"] and R.find(k)["primary"])
            if landed:
                ok += 1
            else:
                pend.append(spec)
        chip = '<span class="chip c-ok">all in</span>' if not pend else (f'<span class="chip c-warn">{len(pend)} pending</span>' if ok else '<span class="chip c-stop">pending</span>')
        det = "".join(f"<li><b>{html.escape(s['name'])}</b>: {html.escape(s.get('pending', 'not landed yet'))}</li>" for s in pend)
        out.append(f'<div><span class="lg-t">{html.escape(name)}</span>{chip}<span class="lg-d">{ok} of {len(keys)} landed'
                   f'{f"<ul>{det}</ul>" if det else ""}</span></div>')
    return dict(FRAG_STATUS=f'<div class="ledger">{"".join(out)}</div>')


def part_renders(R):
    """Appendix: every landed run's full reconstruction render (render_phase.py), grouped in one viewer."""
    panes = []
    for k, spec in R.man["runs"].items():
        if "leg" not in spec or k.startswith("win_") or k == "refined_1e7_own":
            continue
        r = R.render(spec["leg"])
        if r:
            panes.append((spec["short"] if not k.startswith("known_") else f"{spec['short']} mrad", img(r, spec["name"], maxw=1100, jpeg=78)))
    return dict(FRAG_RENDERS=tabs("tabs-renders", panes, default=0) if panes else "<p>No renders pulled yet.</p>")


# =========================================================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="only print what has landed and what is pending")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    man = json.load(open(MANIFEST))
    R = Runs(man)
    print(f"searching {len(R.dirs)} folders, {len(R.analyses)} analyses")
    for k, spec in man["runs"].items():
        if "search" in spec:
            print(f"  {'IN ' if R.search(spec['search']) else 'pending'}  {k:20s} {spec['search']}")
        else:
            f = R.find(k)
            st = "IN " if f["found"] and f["primary"] else ("part" if f["found"] else "pending")
            print(f"  {st:7s}  {k:20s} {f['tag'] or spec['tag'][0]}")
    if a.list:
        return 0
    T = tsv_rows()
    frag = dict(BUILT=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    frag["FIG_PROBES"] = img(part_probes(R, T), "probe size and reconstruction window against aperture")
    known, cuts = part_known(R, T); frag.update(known)
    frag.update(part_window(R, T, cuts))
    frag.update(part_kicks(R, T, cuts))
    frag.update(part_physics(R, T, cuts))
    frag.update(part_search(R, T))
    frag.update(part_refined(R, T, cuts))
    frag.update(part_final(R, T, cuts))
    frag.update(part_status(R))
    frag.update(part_renders(R))
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
