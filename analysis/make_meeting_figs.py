#!/usr/bin/env python
"""
@file make_meeting_figs.py
@brief The meeting figure set for the thin-slab aberration experiment: one clean figure per result,
       every number read from the run data (no values typed into this file).

Figures written to aberration_experiment/figs/<ISO-week>/meeting/:
  fig1_probe.png        what opening the aperture does to the probe (probe diameter d90 vs alpha)
  fig2_depth.png        depth recovery vs alpha with a known probe (the round-alpha result)
  fig3_volume.png       a depth section of the reconstruction with the blindly found atoms on it
  fig4_c1.png           the C1 search: data-fit error vs trial defocus, and the defocus/depth trade
  fig5_probeupdate.png  what the free probe update does with an offset start (it does not move C1)
  fig6_dose.png         shot noise: depth error, recall and kernel quality vs dose

Sources (all produced by the pipeline, nothing hand-entered):
  campaign/round_sweep.tsv                              probe plan per alpha (d90 column)
  <af_final>/out/atomfind_a<A>/report.json              step 0, per alpha
  <af_final>/recon_af_a<A>_lab_NL<NL>/**/*_recons.h5    the phase volumes
  results/<week>/c1_objective.csv + _summary.json       step 1
  results/<week>/probe_refit.csv                        stage 2.5
  results/<week>/relaxation_ladder.csv                  step 0 + step 2 rows
  <relax>/dose_analysis/psf/psf_<el>_a<A>_dose<D>_vol.npy   matched kernels per dose

    ~/hyperspy-bundle/bin/python analysis/make_meeting_figs.py            # all six
    ~/hyperspy-bundle/bin/python analysis/make_meeting_figs.py --only 4 5
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
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)

LAM = 0.0196877                      # 300 keV electron wavelength [A]
BOXZ, ZVAC = 27.525, 4.0             # thin-slab recon box and its z-vacuum band [A]
# dataviz reference palette (light surface): categorical slots 1-3, then ink
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRIDC = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
SPECIES = {"Pb": S1, "Ti": S2, "O": S3}
ALPHA_C = {70: S1, 90: S2, 50: MUTED, 100: INK2}
DEFAULT_AF = "~/Desktop/thin_ab_af_final"
DEFAULT_RELAX = "~/Desktop/relax_0917"


def week_dir():
    d = os.path.join(REPO, "aberration_experiment", "figs", datetime.date.today().strftime("%G-W%V"), "meeting")
    os.makedirs(d, exist_ok=True)
    return d


def style():
    plt.rcParams.update({
        "figure.facecolor": "white", "savefig.facecolor": "white", "savefig.dpi": 200,
        "savefig.bbox": "tight", "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 10.5,
        "axes.titlesize": 11.5, "axes.labelsize": 10.5, "axes.titleweight": "semibold",
        "axes.labelcolor": INK2, "axes.edgecolor": "#c3c2b7", "axes.linewidth": 0.9,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
        "grid.color": GRIDC, "grid.linewidth": 0.7, "xtick.color": MUTED, "ytick.color": MUTED,
        "xtick.labelsize": 9.5, "ytick.labelsize": 9.5, "text.color": INK,
        "legend.frameon": False, "legend.fontsize": 9.5, "lines.linewidth": 2.0,
        "lines.markersize": 6,
    })


def tag(ax, letter, title):
    ax.set_title(f"{letter}   {title}", loc="left", pad=8)


def save(fig, name):
    p = os.path.join(week_dir(), name)
    fig.savefig(p)
    plt.close(fig)
    print(f"wrote {p}")
    return p


# ----------------------------------------------------------------- data readers
def read_round_sweep():
    """(alpha, C3_um, C1_A, d90_A, bin, commented) per planned probe, including commented rows."""
    rows = []
    for line in open(os.path.join(REPO, "campaign", "round_sweep.tsv")):
        if re.match(r"#?a\d{3}\t", line):
            f = line.lstrip("#").rstrip("\n").split("\t")
            rows.append(dict(alpha=int(f[1]), c3_um=float(f[3]) / 1e4, c1=float(f[4]), bin=int(f[6]),
                             d90=float(f[10]), note=f[12], skipped=line.startswith("#")))
    return sorted(rows, key=lambda r: r["alpha"])


def read_reports(af_root):
    """{alpha: finder-v3 report} for the step-0 runs."""
    out = {}
    for d in sorted(glob.glob(os.path.expanduser(os.path.join(af_root, "out", "atomfind_a*")))):
        m = re.search(r"atomfind_a(\d+)$", d)
        rp = os.path.join(d, "report.json")
        if m and os.path.isfile(rp):
            out[int(m.group(1))] = json.load(open(rp))["finder"]["v3"]
    return out


def read_csv_rows(path):
    with open(path) as f:
        return [r for r in csv.DictReader(line for line in f if not line.startswith("#"))]


def results_dir():
    """The newest results/<week> holding the step-1 CSVs."""
    cands = sorted(glob.glob(os.path.join(REPO, "aberration_experiment", "results", "2026-W*")))
    for d in reversed(cands):
        if os.path.isfile(os.path.join(d, "c1_objective.csv")):
            return d
    raise SystemExit("no results/<week>/c1_objective.csv found")


def kernel_quality(npy):
    """(peak phase, peak/background) of a matched kernel, exactly as extract_psf.py reports it:
    background = std over the layer of the peak, outside 0.75 x the crop half-width. extract_psf saves
    exp(i*phase), so the phase has to be taken first -- checked against the pipeline's own printed
    numbers (a70 265/159/93, a90 784/381/76 at 1e7/1e6/1e5)."""
    V = np.load(npy); K = np.angle(V) if np.iscomplexobj(V) else V
    l, r, c = np.unravel_index(np.argmax(K), K.shape)
    W = K.shape[1] // 2
    rad = np.hypot(*np.meshgrid(np.arange(K.shape[2]) - c, np.arange(K.shape[1]) - r))
    bg = np.std(K[l][rad > 0.75 * W])
    return float(K[l, r, c]), float(K[l, r, c] / bg)


def _load_mod(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def build_probe(alpha, aberrations, c1, L=140.0, N=2048):
    """Probe on a converged box for any aberration set (round or not). abTEM applies the
    convention, so nothing here hand-derives the six-fold term."""
    import abtem
    return np.asarray(abtem.Probe(energy=300e3, semiangle_cutoff=alpha, extent=L, gpts=N,
                                  defocus=c1, aberrations=aberrations).build().compute().array)


def aperture_phase(P, alpha, L=140.0, keep=1.12):
    """The wave aberration across the aperture, read back from the probe abTEM built: the phase of
    its Fourier transform inside the illuminated disc. Returns (wrapped phase, half-width in mrad)."""
    # ifftshift first: the probe sits at the array centre, and transforming it there puts a half-pixel
    # phase ramp on every pixel, which displays as a checkerboard rather than the wavefront.
    F = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(P)))
    A = np.abs(F)
    ph = np.where(A > 0.35 * A.max(), np.angle(F), np.nan)
    px = LAM / L * 1e3                                   # mrad per pixel
    h = int(round(keep * alpha / px)); c = P.shape[0] // 2
    return ph[c - h:c + h, c - h:c + h], h * px


def probe_crop(P, L=140.0, half_A=13.0):
    I = np.abs(P) ** 2
    cy, cx = np.unravel_index(I.argmax(), I.shape); h = int(half_A / (L / P.shape[0]))
    return I[cy - h:cy + h, cx - h:cx + h], half_A


# ----------------------------------------------------------------- fig 7: what the operator sees
def fig7_ronchigram(a):
    """Ronchigram and aperture phase as the aperture opens -- the round sweep, cleaned up: the two
    rows that carry the story, without the ray-landing profile or the twin-axis probe-size panel."""
    mrf = _load_mod("make_ronchigram_fig", os.path.join(HERE, "make_ronchigram_fig.py"))
    plan = {r["alpha"]: r for r in read_round_sweep()}
    show = [50, 70, 90, 100]
    rng = np.random.default_rng(3)
    fig, axes = plt.subplots(2, len(show), figsize=(12.4, 6.6))
    for i, alpha in enumerate(show):
        r = plan[alpha]; c3 = r["c3_um"] * 1e4
        P = build_probe(alpha, {"C30": c3, "C50": 1e7}, r["c1"])
        R, hr = mrf.ronchigram(P, 140.0, alpha, rng)
        ax = axes[0, i]
        ax.imshow(R, cmap="gray", extent=[-hr, hr, -hr, hr], interpolation="bilinear")
        ax.add_artist(plt.Circle((0, 0), alpha, fill=False, color="#ffd24a", lw=1.2, ls=":"))
        ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
        ax.set_title(f"{alpha} mrad\nC3 {r['c3_um']:+.0f} µm, focus {r['c1']:+.0f} Å",
                     pad=6, color=INK, fontsize=10.5)
        if i == 0:
            ax.set_ylabel("Ronchigram\n(what the operator sees)", fontsize=10, color=INK2)
        ph, hp = aperture_phase(P, alpha)
        ax = axes[1, i]
        ax.imshow(ph, cmap="twilight_shifted", vmin=-np.pi, vmax=np.pi,
                  extent=[-hp, hp, -hp, hp], interpolation="bilinear")
        ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
        pv = mrf.PLAN.nondefocus_pv(alpha, c3, 1e7)
        ax.set_xlabel(f"wavefront error {pv:.1f} rad", fontsize=9,
                      color=(S3 if pv <= np.pi / 2 else INK2), labelpad=4)
        if i == 0:
            ax.set_ylabel("aberration across\nthe aperture", fontsize=10, color=INK2)
    fig.suptitle("a   Opening a corrected aperture: the fringes are the aberration the corrector no longer holds",
                 x=0.012, ha="left", fontsize=11.5, fontweight="semibold", color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return save(fig, "fig7_ronchigram.png")


# ----------------------------------------------------------------- fig 1: the probe
def fig1_probe(args):
    plan = read_round_sweep()
    spec = importlib.util.spec_from_file_location("plan_probe", os.path.join(REPO, "campaign", "plan_probe.py"))
    pp = importlib.util.module_from_spec(spec); spec.loader.exec_module(pp)
    import abtem

    show_a = (50, 70, 90, 100)
    fig = plt.figure(figsize=(13.0, 4.0))
    gs = fig.add_gridspec(1, 1 + len(show_a), width_ratios=[2.0] + [1] * len(show_a), wspace=0.10)
    ax = fig.add_subplot(gs[0])

    run = [r for r in plan if not r["skipped"]]
    ok = [r for r in run if r["alpha"] <= 100]; bad = [r for r in run if r["alpha"] > 100]
    ax.plot([r["alpha"] for r in run], [r["d90"] for r in run], "-", color=MUTED, lw=1.6, zorder=1)
    ax.plot([r["alpha"] for r in ok], [r["d90"] for r in ok], "o", color=S1, ms=9, zorder=3,
            label="usable")
    ax.plot([r["alpha"] for r in bad], [r["d90"] for r in bad], "o", color="white", mec="#d03b3b",
            mew=2.2, ms=9, zorder=3, label="too big for the scan field")
    ax.axhline(4.0, color=INK2, ls="--", lw=1.2, zorder=0)
    ax.annotate("4 Å target: one unit cell", (118, 4.0), xytext=(0, -16), textcoords="offset points",
                color=INK2, fontsize=9, ha="right")
    for r in run:
        dy = -19 if r["d90"] < 5 else 10
        ax.annotate(f"{r['d90']:.1f} Å", (r["alpha"], r["d90"]), xytext=(0, dy),
                    textcoords="offset points", ha="center", fontsize=9, color=INK)
    ax.set_yscale("log"); ax.set_yticks([2, 4, 10, 25, 60])
    ax.set_yticklabels(["2", "4", "10", "25", "60"]); ax.minorticks_off()
    ax.set_ylim(2.2, 95); ax.set_xlim(43, 127)
    ax.set_xticks([r["alpha"] for r in run])
    ax.set_xlabel("aperture semi-angle α (mrad)"); ax.set_ylabel("probe diameter (Å)")
    ax.legend(loc="upper left")
    tag(ax, "a", "Opening the aperture eventually wrecks the probe")

    ext, N, HALF = 140.0, 1024, 12.0
    for i, alpha in enumerate(show_a):
        r = [q for q in run if q["alpha"] == alpha][0]
        P = pp.build(abtem, alpha, r["c3_um"] * 1e4, r["c1"], 1e7, ext, N)
        I = np.abs(P) ** 2
        cy, cx = np.unravel_index(I.argmax(), I.shape); h = int(HALF / (ext / N))
        crop = I[cy - h:cy + h, cx - h:cx + h]
        a2 = fig.add_subplot(gs[1 + i])
        a2.imshow(crop ** 0.28, cmap="magma", extent=[-HALF, HALF, -HALF, HALF], interpolation="bilinear")
        a2.add_artist(plt.Circle((0, 0), r["d90"] / 2, fill=False, color="white", lw=1.5, alpha=0.95))
        a2.set_xticks([]); a2.set_yticks([]); a2.grid(False)
        a2.set_title(f"{alpha} mrad", fontsize=10.5, pad=5, color=INK)
        a2.text(0.5, 0.04, f"{r['d90']:.1f} Å across", transform=a2.transAxes, ha="center",
                color="white", fontsize=9.5)
        if i == 0:
            a2.plot([-HALF + 1.2, -HALF + 6.2], [HALF - 2.0] * 2, "-", color="white", lw=3)
            a2.text(-HALF + 3.7, HALF - 3.4, "5 Å", color="white", ha="center", va="top", fontsize=9)
            fig.text(a2.get_position().x0, 0.995,
                     "b   The same probes, one scale (circle = the diameter on the left)",
                     fontsize=11.5, fontweight="semibold", color=INK, va="top")
    return save(fig, "fig1_probe.png")


# ----------------------------------------------------------------- fig 2: depth vs alpha
def fig2_depth(args):
    rep = read_reports(args.af)
    al = sorted(rep)
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.9))

    ax = axes[0]
    # atomfind's own health check: above 5% species confusion the labels are not trustworthy
    def confusion(a):
        cf = rep[a].get("confusion") or {}
        off = sum(cf.get(f"{x}->{y}", 0) for x in (82, 22, 8) for y in (82, 22, 8) if x != y)
        dia = sum(cf.get(f"{x}->{x}", 0) for x in (82, 22, 8))
        return off / max(off + dia, 1)
    bad = [a for a in al if confusion(a) > 0.05]
    for sp, col in SPECIES.items():
        y = [100 * rep[a][sp]["recall_bulk"] for a in al]
        ax.plot(al, y, "-", color=col, label=sp)
        good_m = [a not in bad for a in al]
        ax.plot([a for a, g in zip(al, good_m) if g], [v for v, g in zip(y, good_m) if g], "o", color=col)
        ax.plot([a for a, g in zip(al, good_m) if not g], [v for v, g in zip(y, good_m) if not g],
                "o", mfc="white", mec=col, mew=1.8)
    if bad:
        ax.annotate("hollow: species labels\nnot trustworthy here", (bad[0], 30), xytext=(6, 0),
                    textcoords="offset points", fontsize=9, color=INK2, va="center")
    ax.set_ylim(0, 105); ax.set_xticks(al); ax.set_xlabel("aperture semi-angle α (mrad)")
    ax.set_ylabel("atoms found (%)"); ax.legend(loc="lower right", ncol=3)
    tag(ax, "a", "Atoms found, by species")

    ax = axes[1]
    ax.plot(al, [rep[a]["z_rms_A"] for a in al], "-o", color=S1, label="measured error")
    aa = np.linspace(min(al) - 3, max(al) + 3, 100)
    ax.plot(aa, LAM / (aa / 1000) ** 2, "--", color=MUTED, lw=1.4, label="classical depth limit λ/α²")
    ax.set_yscale("log"); ax.set_yticks([0.3, 1, 3, 8]); ax.set_yticklabels(["0.3", "1", "3", "8"])
    ax.minorticks_off(); ax.set_xticks(al)
    ax.set_xlabel("aperture semi-angle α (mrad)"); ax.set_ylabel("depth error (Å)")
    ax.legend(loc="upper right")
    tag(ax, "b", "Depth error falls as the aperture opens")

    ax = axes[2]
    ax.plot(al, [rep[a]["xy_rms_A"] for a in al], "-o", color=S1, label="measured error")
    ax.axhline(0.24, color=MUTED, ls="--", lw=1.4)
    ax.annotate("real polar displacements in this\nmaterial: 0.24 Å typical", (al[0], 0.24),
                xytext=(2, 8), textcoords="offset points", fontsize=9, color=INK2)
    ax.set_ylim(0, 0.32); ax.set_xticks(al); ax.set_xlabel("aperture semi-angle α (mrad)")
    ax.set_ylabel("in-plane error (Å)")
    tag(ax, "c", "In-plane error stays far below the signal")
    fig.tight_layout()
    return save(fig, "fig2_depth.png")


# ----------------------------------------------------------------- fig 3: the volume itself
def fig3_volume(args):
    """A depth section of the reconstruction with the blindly found atoms on it, alpha 70 vs 90."""
    from atomfind import config, align                      # noqa: E402
    from dataclasses import replace
    from matplotlib.lines import Line2D
    os.environ.setdefault("ATOMFIND_DATA", os.path.expanduser("~/Desktop/thin_ab_af/gtdata"))
    show = [(70, 14), (90, 23)]
    fig, axes = plt.subplots(1, len(show), figsize=(11.8, 4.8))
    for ax, (alpha, nl) in zip(np.atleast_1d(axes), show):
        h5 = sorted(glob.glob(os.path.expanduser(os.path.join(
            args.af, f"recon_af_a{alpha}_lab_NL{nl}", "**", "*_recons.h5")), recursive=True))
        af = os.path.expanduser(os.path.join(args.af, "out", f"atomfind_a{alpha}"))
        cfg = replace(config.preset("thin"), recon_vol=h5[-1], dz=BOXZ / nl)
        V, dx = align.load_phase(cfg); V = align.crop_to_fov(V, dx, cfg)
        pos, Z = align.load_gt(cfg)
        found = np.load(os.path.join(af, "found_atoms.npy"))
        a1 = align.refine_with_atoms(align.register(V, dx, pos, Z, cfg), found, pos, Z, cfg)
        gr, gc, gl = a1.site_to_index(pos[:, 0], pos[:, 1], pos[:, 2])
        fr, fc, fl = found["row"], found["col"], found["layer"]
        row0 = int(np.median(fr)); halfpx = int(1.0 / dx)          # a 2 A thick slice of the volume
        sec = V[:, max(row0 - halfpx, 0):row0 + halfpx, :].mean(1)
        zrec = (np.arange(V.shape[0]) + 0.5) * cfg.dz
        xmax = (V.shape[2] - 1) * dx
        ax.imshow(sec, cmap="magma", vmin=0, vmax=np.percentile(sec, 99.6), aspect="auto",
                  interpolation="bilinear", extent=[0, xmax, zrec[-1], zrec[0]])
        inb = lambda r, c: (np.abs(r - row0) < halfpx) & (c * dx >= 0) & (c * dx <= xmax)
        m = inb(gr, gc)
        ax.plot(gc[m] * dx, gl[m] * cfg.dz, ".", color="#bdbdbd", ms=3.0, alpha=0.75, zorder=2)
        mf = inb(fr, fc)
        for sp, zz in (("Pb", 82), ("Ti", 22), ("O", 8)):
            s2 = mf & (found["species"] == zz)
            ax.plot(fc[s2] * dx, fl[s2] * cfg.dz, "o", mfc="none", mec=SPECIES[sp], mew=1.9, ms=11, zorder=3)
        for z0, z1 in ((zrec[0] - cfg.dz, ZVAC), (BOXZ - ZVAC, zrec[-1] + cfg.dz)):
            ax.axhspan(z0, z1, color="white", alpha=0.22, lw=0, zorder=1)
        ax.set_xlim(0, xmax); ax.set_ylim(zrec[-1] + cfg.dz / 2, zrec[0] - cfg.dz / 2)
        ax.set_xlabel("position across the sample (Å)")
        ax.set_ylabel("depth into the sample (Å)")
        ax.grid(False)
        rep = json.load(open(os.path.join(af, "report.json")))["finder"]["v3"]
        tag(ax, "a" if alpha == show[0][0] else "b",
            f"α = {alpha} mrad — depth error {rep['z_rms_A']:.2f} Å")
    h = [Line2D([], [], ls="", marker="o", mfc="none", mec=c, mew=1.9, ms=9) for c in SPECIES.values()]
    h += [Line2D([], [], ls="", marker=".", color="#bdbdbd", ms=9)]
    np.atleast_1d(axes)[0].legend(h, [f"found: {k}" for k in SPECIES] + ["true positions"],
                                  loc="upper left", bbox_to_anchor=(0.0, -0.16), ncol=4,
                                  handletextpad=0.3, columnspacing=1.4)
    fig.text(0.53, -0.01, "pale bands: the empty space either side of the sample, where surface artefacts are parked",
             fontsize=9, color=INK2)
    fig.tight_layout()
    return save(fig, "fig3_volume.png")


# ----------------------------------------------------------------- fig 4: the C1 search
def fig4_c1(args):
    rows = read_csv_rows(os.path.join(results_dir(), "c1_objective.csv"))
    summ = json.load(open(os.path.join(results_dir(), "c1_objective_summary.json")))
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.1))

    for ax, alpha in zip(axes[:2], (70, 90)):
        for niter, col, lab in ((50, S1, "50 iterations"), (200, S2, "200 iterations")):
            pts = sorted([(float(r["dc1"]), float(r["final_error"])) for r in rows
                          if int(r["alpha"]) == alpha and int(r["niter"]) == niter])
            if not pts:
                continue
            d = sorted({p[0] for p in pts})
            mean = [np.mean([e for x, e in pts if x == v]) for v in d]
            ax.plot(d, mean, "-o", color=col, ms=5, label=lab)
        E = summ.get(f"a{alpha}_lab_n50")
        if E:
            ax.axvspan(E["lo"] - E["c1_true"], E["hi"] - E["c1_true"], color=S1, alpha=0.12, lw=0)
            ax.errorbar([E["c1_fit"] - E["c1_true"]], [min(mean)], xerr=[E["halfwidth"]],
                        fmt="D", ms=8, color=INK, mfc="white", mew=1.6, capsize=3, zorder=5,
                        label=f"best fit {E['c1_fit'] - E['c1_true']:+.0f} ± {E['halfwidth']:.0f} Å")
        ax.axvline(0, color=INK2, ls="--", lw=1.2)
        ax.annotate("true focus", (0, 0.42), xycoords=("data", "axes fraction"), xytext=(-5, 0),
                    textcoords="offset points", fontsize=9, color=INK2, va="center", ha="right",
                    rotation=90)
        ax.set_yscale("log")
        ax.set_yticks([5, 10, 20, 50, 100, 200]); ax.set_yticklabels(["5", "10", "20", "50", "100", "200"])
        ax.minorticks_off()
        ax.set_xlabel("focus error of the trial probe (Å)")
        ax.set_ylabel("mismatch to the measured data")
        ax.legend(loc="upper center")
        tag(ax, "a" if alpha == 70 else "b", f"α = {alpha} mrad")

    ax = axes[2]
    # only where the reconstruction itself survives: past +2 A at a70 the object collapses into the
    # vacuum layers (fig 4a's error jump), so its depth centroid stops meaning anything
    valid = {70: (-4, 2), 90: (-4, 4)}
    for alpha, col in ((70, S1), (90, S2)):
        lo, hi = valid[alpha]
        pts = {float(r["dc1"]): float(r["z_centroid"]) for r in rows
               if int(r["alpha"]) == alpha and int(r["niter"]) == 50 and r["z_centroid"]
               and lo <= float(r["dc1"]) <= hi and int(r["rep"]) == 1}
        if 0.0 in pts:
            d = sorted(pts)
            ax.plot(d, [pts[v] - pts[0.0] for v in d], "-o", color=col, ms=7, label=f"α = {alpha} mrad")
    xx = np.array([-4.4, 4.4]); ax.plot(xx, xx, "--", color=MUTED, lw=1.4, label="one-for-one")
    ax.set_xlabel("focus error of the trial probe (Å)")
    ax.set_ylabel("shift of the sample in depth (Å)")
    ax.legend(loc="upper left")
    tag(ax, "c", "Why the answer has a ± 2 Å floor")
    fig.tight_layout()
    return save(fig, "fig4_c1.png")


# ----------------------------------------------------------------- fig 5: the probe update
def fig5_probeupdate(args):
    rows = read_csv_rows(os.path.join(results_dir(), "probe_refit.csv"))
    # a crashed run has no probe record, so its own row carries no truth: take the true C1 per alpha
    # from the runs that did finish, and the start from the dir name (which always has it)
    truth = {int(r["alpha"]): float(r["c1_true"]) for r in rows if r["c1_true"]}
    def start_off(r):
        return float(r["c1_start"]) - truth[int(r["alpha"])]
    num = lambda v: float(v) if v not in (None, "") else np.nan

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.3))
    ax = axes[0]
    for alpha, col in ((70, S1), (90, S2)):
        up = [r for r in rows if int(r["alpha"]) == alpha and r["variant"] != "fixed"]
        good = [(start_off(r), num(r["dc1_final"])) for r in up
                if r["status"] == "ok" and num(r["ov_best_c1"]) > 0.5]
        junk = [start_off(r) for r in up if r["status"] != "ok" or not (num(r["ov_best_c1"]) > 0.5)]
        if good:
            ax.plot([g[0] for g in good], [g[1] for g in good], "o", color=col, ms=9,
                    label=f"α = {alpha} mrad")
        for x in junk:
            ax.plot([x], [0], "x", color="#d03b3b", ms=11, mew=2.4)
    lim = 36
    ax.plot([-lim, lim], [-lim, lim], "--", color=MUTED, lw=1.4)
    ax.annotate("probe did not move", (lim * 0.60, lim * 0.60), rotation=39, fontsize=9,
                color=INK2, ha="center", va="bottom")
    ax.axhline(0, color=S3, lw=1.8)
    ax.annotate("what success would look like", (-lim * 0.95, 1.6), fontsize=9, color=S3)
    ax.plot([], [], "x", color="#d03b3b", ms=10, mew=2.2, label="probe destroyed, or run crashed")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("focus error handed to the solver (Å)")
    ax.set_ylabel("focus error left in the refined probe (Å)")
    ax.legend(loc="lower right")
    tag(ax, "a", "The free probe update does not find the focus")

    ax = axes[1]
    xi, lab, st, fi, cols, div = [], [], [], [], [], None
    for alpha, col in ((70, S1), (90, S2)):
        up = [r for r in rows if int(r["alpha"]) == alpha and r["variant"] != "fixed" and r["status"] == "ok"]
        for off in sorted({round(start_off(r)) for r in up}):
            grp = [r for r in up if round(start_off(r)) == off]      # the two release schedules, averaged
            xi.append(len(xi)); lab.append(f"{off:+d}")
            st.append(np.mean([num(r["ov_start"]) for r in grp]))
            fi.append(np.mean([num(r["ov_final"]) for r in grp])); cols.append(col)
        if div is None:
            div = len(xi) - 0.5
    xi = np.array(xi)
    ax.bar(xi - 0.2, st, 0.38, color="#c9d8ee", label="probe handed in")
    ax.bar(xi + 0.2, fi, 0.38, color=cols, label="probe after refinement")
    ax.axvline(div, color=GRIDC, lw=1.2)
    ax.text(div / 2, 1.02, "α = 70 mrad", ha="center", fontsize=9.5, color=S1)
    ax.text((div + len(xi) - 0.5) / 2, 1.02, "α = 90 mrad", ha="center", fontsize=9.5, color=S2)
    ax.set_xticks(xi); ax.set_xticklabels(lab)
    ax.set_ylim(0, 1.16); ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("focus error handed to the solver (Å)")
    ax.set_ylabel("match to the true probe")
    ax.legend(loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.17))
    tag(ax, "b", "Refinement barely improves the probe")
    fig.tight_layout()
    return save(fig, "fig5_probeupdate.png")


# ----------------------------------------------------------------- fig 6: dose
def fig6_dose(args):
    lad = read_csv_rows(os.path.join(results_dir(), "relaxation_ladder.csv"))
    base = {int(r["alpha"]): r for r in lad if r["step"] == "0"}
    dose_rows = [r for r in lad if r["step"] == "2"]
    dose_of = lambda r: float(re.search(r"([\d.]+e[\d+-]+)", r["label"]).group(1))
    done = sorted({dose_of(r) for r in dose_rows}, reverse=True)
    psf = os.path.expanduser(os.path.join(args.relax, "dose_analysis", "psf"))
    # every dose that was RUN, including the one with no result, from the kernel dir
    ran = sorted({float(m.group(1)) for f in glob.glob(os.path.join(psf, "psf_Pb_a*_dose*_vol.npy"))
                  for m in [re.search(r"_dose([\d.]+e[\d+-]+)_", os.path.basename(f))] if m}
                 | set(done) | {1e4}, reverse=True)
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.1))

    def get(alpha, dose, key):
        for r in dose_rows:
            if int(r["alpha"]) == alpha and dose_of(r) == dose:
                return float(r[key]) if r[key] else np.nan
        return np.nan

    def dose_axis(ax):
        ax.set_xscale("log"); ax.invert_xaxis()
        ax.set_xticks(ran); ax.set_xticklabels([f"$10^{{{int(np.log10(d))}}}$" for d in ran])
        ax.minorticks_off(); ax.set_xlabel("dose (electrons per Å²)")

    ax = axes[0]
    for alpha, col in ((70, S1), (90, S2)):
        ax.plot(done, [get(alpha, d, "z_rms") for d in done], "-o", color=col, label=f"α = {alpha} mrad")
        b = float(base[alpha]["z_rms"]); ax.axhline(b, color=col, ls=":", lw=1.4)
        ax.annotate(f"noiseless {b:.2f} Å", (done[-1], b), xytext=(0, 4), textcoords="offset points",
                    fontsize=8.5, color=col, ha="left")
    ax.axvspan(1e4 / 1.6, 1e4 * 1.6, color="#d03b3b", alpha=0.07, lw=0)
    ax.set_ylim(0, 1.25); dose_axis(ax); ax.set_ylabel("depth error (Å)")
    ax.legend(loc="upper left")
    tag(ax, "a", "Depth error grows slowly as dose falls")

    ax = axes[1]
    for alpha, ls in ((70, "-"), (90, "--")):
        for sp, col in SPECIES.items():
            ax.plot(done, [100 * get(alpha, d, f"{sp}_recall_bulk") for d in done], ls, marker="o",
                    ms=5, color=col)
    ax.axvspan(1e4 / 1.6, 1e4 * 1.6, color="#d03b3b", alpha=0.07, lw=0)
    ax.set_ylim(0, 105); dose_axis(ax); ax.set_ylabel("atoms found (%)")
    h = [Line2D([], [], color=c, marker="o", ms=5) for c in SPECIES.values()]
    h += [Line2D([], [], color=MUTED, ls="-"), Line2D([], [], color=MUTED, ls="--")]
    ax.legend(h, list(SPECIES) + ["α = 70", "α = 90"], loc="lower left", ncol=2)
    tag(ax, "b", "Oxygen is the first to go")

    ax = axes[2]
    for alpha, col in ((70, S1), (90, S2)):
        xs, ys = [], []
        for d in ran:
            f = os.path.join(psf, f"psf_Pb_a{alpha}_dose{d:.0e}".replace("e+0", "e") + "_vol.npy")
            if os.path.isfile(f):
                xs.append(d); ys.append(kernel_quality(f)[1])
        ax.plot(xs, ys, "-o", color=col, label=f"α = {alpha} mrad")
        ax.plot([min(ran)], [min(ys) * 0.45], "x", color="#d03b3b", ms=13, mew=2.6)
    ax.axvspan(1e4 / 1.6, 1e4 * 1.6, color="#d03b3b", alpha=0.07, lw=0)
    ax.annotate("reference could not be\nmeasured at all", (min(ran), 40), xytext=(-6, 26),
                textcoords="offset points", fontsize=9, color="#d03b3b", ha="right")
    ax.set_yscale("log"); ax.set_yticks([30, 100, 300, 800]); ax.set_yticklabels(["30", "100", "300", "800"])
    ax.minorticks_off(); dose_axis(ax); ax.set_ylabel("reference-signal quality")
    ax.legend(loc="lower left")
    tag(ax, "c", "The calibration reference fails first")
    fig.tight_layout()
    return save(fig, "fig6_dose.png")


FIGS = {1: fig1_probe, 2: fig2_depth, 3: fig3_volume, 4: fig4_c1, 5: fig5_probeupdate,
        6: fig6_dose, 7: fig7_ronchigram}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--af", default=DEFAULT_AF, help="canonical step-0 set (recons + out/atomfind_a*)")
    ap.add_argument("--relax", default=DEFAULT_RELAX, help="relaxation tarball root (dose_analysis inside)")
    ap.add_argument("--only", nargs="+", type=int, default=sorted(FIGS))
    a = ap.parse_args()
    style()
    for n in a.only:
        FIGS[n](a)


if __name__ == "__main__":
    main()
