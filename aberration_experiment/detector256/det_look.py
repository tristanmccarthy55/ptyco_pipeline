"""@file det_look.py
@brief detector256: every run's honest numbers, depth sections, phase images and the limits figure, saved in results/.

Same machinery as aberration_experiment/mono300/phase1_look.py (analysis/make_arm_page.py, unchanged): right place (right
element within 0.5 A in depth) against guessing (atom_scores), the verdict rule (good: Pb >= 80 and Ti >= 45 %; okay:
Pb >= 50 %), the vacuum loading (rms phase per slice in the vacuum bands over that in the slab; 0.31 for a clean run) and
the depth sections. Misfit eps = final Fourier error x N / 2e5 with N the width of the DATA each run fitted (its own
sim_meta.mat); eps compares like with like only on the same data grid. Measurements per unknown = positions x data pixels
over (scan 20 A + modelled window)^2 / dx^2 x 18 slices x 2 real numbers.

Each run's kernels are checked in its extraction logs: one whose brightest voxel is NOT its own atom ('argmax is off the
crop centre': a neighbour inside the crop) or averaged from fewer than 16 atoms is flagged, and its scores are not about
the reconstruction (the first R1 / R2 analysis, 2026-10-23, before PSF_HALF_A=2.22).

    ~/hyperspy-bundle/bin/python aberration_experiment/detector256/det_look.py [more pull folders]
Pull folders may be laid out flat (the tarballs unpacked in place: analysis_* / recon_af_* / sim_out_af_* side by side) or
as <d>/analysis + <d>/results; a view is built in a temporary folder, nothing is moved. Defaults: ~/Desktop/mono_phase1
(the reference), ~/Desktop/detector256_{r1r2, r1r2_kh, r3, limits}. For a run in several analyses the newest wins.
Writes aberration_experiment/detector256/results/: det_scores.json (every number + where it came from), det_scores.csv,
det_sections_cameras.png, det_sections_limits.png, det_limits.png, phase/<run>_lab.jpg (a 75 % copy).
"""
import csv
import glob
import importlib.util
import json
import os
import re
import shutil
import sys
import tempfile

import numpy as np
from scipy.io import loadmat

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(HERE, "results")
spec = importlib.util.spec_from_file_location("arm", os.path.join(REPO, "analysis", "make_arm_page.py"))
arm = importlib.util.module_from_spec(spec); sys.modules["arm"] = arm; spec.loader.exec_module(arm)
spec = importlib.util.spec_from_file_location("p1", os.path.join(REPO, "aberration_experiment", "mono300", "phase1_look.py"))
p1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(p1)
plt = arm.plt
arm.FIGS = OUT
PULLS = ["~/Desktop/mono_phase1", "~/Desktop/detector256_r1r2", "~/Desktop/detector256_r1r2_kh", "~/Desktop/detector256_r3",
         "~/Desktop/detector256_limits"] + sys.argv[1:]
C = "_fs2p5_src0p1"
SCAN_A, NL = 20.0, 18
# key, name, label, tag after the label, data dir (the lab leg's data), model grid px, iterations per engine, what changes vs R2
RUNS = [
    ("ref", "reference: 1420 px, ±133 mrad", "armf_a080", f"{C}_fixed", f"armf_a080_lab{C}", 1420, 200, "full camera"),
    ("R1", "R1: 1024 px, ±96 mrad", "armf_a080", f"{C}_detc96", f"armf_a080_lab{C}_detc96", 1024, 200, "full pixels"),
    ("C512", "C512: 512 px, ±96 mrad", "armf_a080", f"{C}_detc96b2_up1", f"armf_a080_lab{C}_detc96b2", 1024, 200, "512 px"),
    ("R2", "R2: 256 px, ±96 mrad", "armf_a080", f"{C}_detc96b4_up2", f"armf_a080_lab{C}_detc96b4", 1024, 200, "-"),
    ("C128", "C128: 128 px, ±96 mrad", "armf_a080", f"{C}_detc96b8_up3", f"armf_a080_lab{C}_detc96b8", 1024, 200, "128 px"),
    ("R3", "R3: 256 px, not modelled", "armf_a080", f"{C}_detc96b4", f"armf_a080_lab{C}_detc96b4", 256, 200, "no upsampling"),
    ("R4", "R4: 354 px, ±133 mrad", "armf_a080", f"{C}_detc133b4_up2", f"armf_a080_lab{C}_detc133b4", 1416, 200, "±133 mrad"),
    ("S04", "S04: 256 px, 0.4 Å scan", "armf_a080_s04", f"{C}_detc96b4_up2", f"armf_a080_s04_lab{C}_detc96b4", 1024, 200, "0.4 Å step"),
    ("R2_n400", "R2, 400 + 400 it", "armf_a080", f"{C}_detc96b4_up2_n400", f"armf_a080_lab{C}_detc96b4", 1024, 400, "2x iterations"),
    ("R2_n1000", "R2, 1000 + 1000 it", "armf_a080", f"{C}_detc96b4_up2_n1000", f"armf_a080_lab{C}_detc96b4", 1024, 1000, "5x iterations"),
]
SHEETS = (("cameras", ["ref", "R1", "C512", "R2", "C128", "R3"]), ("limits", ["R2", "R2_n400", "R2_n1000", "R4", "S04"]))
BLUE, GREY, INK, INK2 = "#2a78d6", "#8a8984", "#0b0b0b", "#52514e"
LABEL_OFFSET = {"ref": (-24, 4), "R1": (7, -12)}          # the two full cameras share an x: label them apart
KEYS = ("ref 1420 px ±133 mrad · R1 1024 px ±96 · C512 512 px · R2 256 px · C128 128 px (all ±96, 0.5 Å scan) · "
        "R4 354 px ±133 · S04 256 px on a 0.4 Å scan; all but ref and R1 summed pixels modelled by upsampling")


def build_view(folders):
    """A temporary <view>/analysis + <view>/results of symlinks over the pull folders (flat or structured)."""
    v = tempfile.mkdtemp(prefix="det_view_")
    for sub in ("analysis", "results"):
        os.makedirs(os.path.join(v, sub))
    seen = set()
    for f in (os.path.expanduser(x) for x in folders):
        if not os.path.isdir(f):
            continue
        cands = [(p, "analysis") for p in glob.glob(os.path.join(f, "analysis_*")) + glob.glob(os.path.join(f, "analysis", "analysis_*"))]
        cands += [(p, "results") for pat in ("recon_af_*", "sim_out_af_*")
                  for p in glob.glob(os.path.join(f, pat)) + glob.glob(os.path.join(f, "results", pat))]
        for p, sub in cands:
            if os.path.isdir(p) and (sub, os.path.basename(p)) not in seen:
                seen.add((sub, os.path.basename(p)))
                os.symlink(p, os.path.join(v, sub, os.path.basename(p)))
    return v


def kernel_check(info):
    """'' or why the run's Pb / Ti kernels are suspect, from the extraction logs of the analysis the run was found in."""
    if not info.get("af"):
        return ""
    bad = []
    for el in ("Pb", "Ti"):
        lg = os.path.join(os.path.dirname(info["af"]), "logs", f"extract_{el}_{info['tag']}.log")
        if os.path.exists(lg):
            t = open(lg).read()
            if "off the crop centre" in t:
                bad.append(f"{el} kernel off its atom")
            m = re.search(r"\[avg\] (\d+) grid atoms", t)
            if m and int(m.group(1)) < 16:
                bad.append(f"{el} from {m.group(1)} atoms")
    return "; ".join(bad)


def trace_of(R, leg):
    """(final Fourier error, iterations per engine from the trace header, trace path) or (None, None, None)."""
    t = R._newest([p for d in R.dirs for p in glob.glob(os.path.join(d, "results", f"recon_af_{leg}_NL*", "analysis", "**",
                                                                     "*_error_trace.csv"), recursive=True)])
    if not t:
        return None, None, None
    head = open(t).readline()
    it = re.search(r"NITER=(\d+)", head)
    e = [float(l.split(",")[1]) for l in open(t) if l[:1].isdigit()]
    return e[-1], int(it.group(1)) if it else None, t


def meta_of(R, data_dir):
    m = R._newest([os.path.join(d, "results", f"sim_out_af_{data_dir}", "01", "sim_meta.mat") for d in R.dirs])
    return (loadmat(m, simplify_cells=True)["meta"], m) if m else (None, None)


def score_all(R):
    out = []
    for key, name, lab, tag, data, model_px, niter, change in RUNS:
        info = R.find(key)
        leg = f"{lab}_lab{tag}"
        r = dict(key=key, run=name, change_vs_R2=change, label=lab, recon_leg=f"recon_af_{leg}", data_dir=f"sim_out_af_{data}",
                 model_px=model_px, iterations_per_engine=niter, found=info["found"])
        if info["found"]:
            sc = arm.atom_scores(info)
            err, it, tpath = trace_of(R, leg)
            meta, mpath = meta_of(R, data)
            r.update(analysis=os.path.basename(os.path.dirname(info["af"])) if info.get("af") else
                     (os.path.basename(os.path.dirname(info["fd"])) if info.get("fd") else None),
                     kernels=kernel_check(info), vacuum_ratio=p1.vacuum_ratio(info.get("fd")),
                     error_trace=os.path.basename(tpath) if tpath else None, final_error=err)
            if it and it != niter:
                r["kernels"] = (r["kernels"] + "; " if r["kernels"] else "") + f"trace says NITER={it}, expected {niter}"
            if meta is not None:
                n, dx, step = int(meta["Ndpx"]), float(meta["dx_object_A"]), float(meta["scan_step_A"])
                npos = int(round(SCAN_A / step)) ** 2
                unknowns = ((SCAN_A + model_px * dx) / dx) ** 2 * NL * 2
                r.update(data_px=n, theta_max_mrad=round(n / 2 * float(meta["d_alpha_rad"]) * 1e3, 2), object_px_A=dx,
                         scan_step_A=step, positions=npos, measurements_per_unknown=npos * n * n / unknowns,
                         eps=err * n / 2e5 if err is not None else None)
            if sc:
                r.update(right={k: sc["right"][z] for k, z in (("Pb", 82), ("Ti", 22), ("O", 8))},
                         guess={k: sc["guess"][z] for k, z in (("Pb", 82), ("Ti", 22), ("O", 8))},
                         wrong=sc["wrong"], n_found=sc["n"], verdict=arm.verdict(sc))
            else:
                r["verdict"] = "failed (no atoms)" if r.get("analysis") else None
        out.append(r)
    return out


def print_table(out):
    f3 = lambda d: "/".join(f"{100 * d[s]:3.0f}" for s in ("Pb", "Ti", "O"))
    print(f"{'run':30s}{'right Pb/Ti/O %':>16s} {'guess Pb/Ti/O %':>16s} {'wrong':>5s} {'vac':>5s} {'eps (N)':>14s} {'m/u':>6s}  verdict")
    for r in out:
        if not r["found"]:
            print(f"{r['run']:30s}{'not pulled':>16s}"); continue
        vr = f"{r['vacuum_ratio']:5.2f}" if r.get("vacuum_ratio") is not None else "    -"
        ep = f"{r['eps']:.4f} ({r['data_px']})" if r.get("eps") is not None else "-"
        mu = f"{r['measurements_per_unknown']:6.1f}" if r.get("measurements_per_unknown") else "     -"
        sc = (f"{f3(r['right']):>16s} {f3(r['guess']):>16s} {100 * r['wrong']:4.0f}%") if "right" in r else f"{'no atoms':>16s} {'':>16s} {'':>5s}"
        print(f"{r['run']:30s}{sc} {vr} {ep:>14s} {mu}  {r.get('verdict')}"
              + (f"   KERNELS SUSPECT: {r['kernels']}" if r.get("kernels") else ""))


def save_tables(out):
    json.dump(dict(written_by="aberration_experiment/detector256/det_look.py", rules=dict(
        right="right element within 0.6 A in plane and 0.5 A in depth (make_arm_page.atom_scores)",
        guess="the same found atoms slid along their columns: the score of depths that are guesses",
        verdict="good: Pb >= 80 and Ti >= 45 %; okay: Pb >= 50 %", eps="final Fourier error x data width / 2e5",
        measurements_per_unknown="positions x data px^2 / (((20 A + model px x dx) / dx)^2 x 18 slices x 2)"),
        runs=out), open(os.path.join(OUT, "det_scores.json"), "w"), indent=1, ensure_ascii=False)
    cols = ["key", "run", "change_vs_R2", "verdict", "Pb_right", "Ti_right", "O_right", "Pb_guess", "Ti_guess", "O_guess", "wrong",
            "vacuum_ratio", "eps", "data_px", "theta_max_mrad", "model_px", "object_px_A", "scan_step_A", "positions",
            "iterations_per_engine", "measurements_per_unknown", "kernels", "analysis", "recon_leg", "data_dir"]
    with open(os.path.join(OUT, "det_scores.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in out:
            row = {k: r.get(k) for k in cols}
            for s in ("Pb", "Ti", "O"):
                row[f"{s}_right"] = round(100 * r["right"][s], 1) if "right" in r else None
                row[f"{s}_guess"] = round(100 * r["guess"][s], 1) if "guess" in r else None
            for k in ("wrong",):
                row[k] = round(100 * r[k], 1) if r.get(k) is not None else None
            for k, nd in (("vacuum_ratio", 3), ("eps", 5), ("object_px_A", 5), ("measurements_per_unknown", 2)):
                row[k] = round(r[k], nd) if r.get(k) is not None else None
            w.writerow(row)


def copy_phase(R, out):
    os.makedirs(os.path.join(OUT, "phase"), exist_ok=True)
    for r in out:
        if r.get("analysis"):
            src = R.render(r["recon_leg"].replace("recon_af_", "", 1))
            if src:                       # a reduced copy for the repo (the full-size one stays in the pulled analysis)
                from PIL import Image
                im = Image.open(src).convert("RGB")
                im = im.resize((int(im.width * 0.75), int(im.height * 0.75)), Image.LANCZOS)
                # JPEG, not a palette PNG: 256 colours contoured the faint x-z shading (vacuum loading), ~200 KB a run
                im.save(os.path.join(OUT, "phase", f"{r['key']}_lab.jpg"), "JPEG", quality=92)


def sections(R, out, cuts):
    for series, keys in SHEETS:
        keys = [k for k in keys if arm.load_leg(R.find(k)) is not None]
        if not keys:
            continue
        fig, axs = plt.subplots(2, len(keys), figsize=(2.3 * len(keys), 8.6), squeeze=False)
        for j, k in enumerate(keys):
            info = R.find(k); lg = arm.load_leg(info); sc = arm.atom_scores(info)
            r = next(x for x in out if x["key"] == k)
            for i, (cut, rA) in enumerate(cuts):
                ax = axs[i, j]
                nx = lg["V"].shape[2]
                arm.xz_section(ax, lg, sc, rA / lg["dx"], xlim=(nx * lg["dx"] / 2 - 8, nx * lg["dx"] / 2 + 8), show_y=(j == 0))
                if i == 0:
                    ttl = r["run"].replace(": ", ":\n")
                    if "right" in r:
                        ttl += f"\nPb {100 * r['right']['Pb']:.0f} vs {100 * r['guess']['Pb']:.0f}" \
                               f"\nTi {100 * r['right']['Ti']:.0f} vs {100 * r['guess']['Ti']:.0f}\n{r['verdict']}"
                    else:
                        ttl += "\n(no atoms)"
                    ax.set_title(ttl, fontsize=8.5, loc="left")
                ax.set_xlabel(f"x (Å) · {cut}", fontsize=8)
        fig.tight_layout()
        arm.save(fig, f"det_sections_{series}.png")


def limits_figure(out):
    """Top: right place (filled) and guessing (hollow) vs measurements per unknown, runs at 200 + 200 iterations, labelled.
    Bottom: the same for R2 against iterations (redundancy fixed at R2's), R1 as a dashed reference."""
    base = [r for r in out if "right" in r and r["iterations_per_engine"] == 200 and r.get("measurements_per_unknown")]
    iters = [r for r in out if "right" in r and r["key"].startswith("R2")]
    r1 = next((r for r in out if r["key"] == "R1" and "right" in r), None)
    fig, axs = plt.subplots(2, 2, figsize=(10, 7.4), gridspec_kw=dict(height_ratios=[1.35, 1]))
    for c, (el, bar) in enumerate((("Pb", 80), ("Ti", 45))):
        for row, (rows, xkey, xlab) in enumerate(((base, "measurements_per_unknown", "measurements per unknown (log)"),
                                                  (iters, "iterations_per_engine", "iterations per engine (R2, 256 px)"))):
            ax = axs[row, c]
            for r in rows:
                x = r[xkey]; y, g = 100 * r["right"][el], 100 * r["guess"][el]
                ax.plot([x, x], [g, y], color=GREY, lw=1.2, zorder=1)
                ax.plot(x, g, "o", ms=8, mfc="white", mec=GREY, mew=1.6, zorder=2)
                ax.plot(x, y, "o", ms=9, color=BLUE, mec="white", mew=1.5, zorder=3)
                ax.annotate(r["key"], (x, y), xytext=LABEL_OFFSET.get(r["key"], (6, 4)) if row == 0 else (6, 4),
                            textcoords="offset points", fontsize=8.5, color=INK2)
            ax.axhline(bar, color=INK2, lw=0.8, ls=(0, (4, 3)))
            ax.text(1.0, bar, f" good: {el} ≥ {bar} %", transform=ax.get_yaxis_transform(), fontsize=8, color=INK2,
                    va="center", ha="left", clip_on=False)
            if row == 1 and r1:
                ax.axhline(100 * r1["right"][el], color=BLUE, lw=0.8, ls=(0, (1, 2)))
                ax.text(0.98, 100 * r1["right"][el] + 2, "R1 (1024 px, 200 it)", transform=ax.get_yaxis_transform(), fontsize=8,
                        color=INK2, ha="right")
            ax.set_xscale("log")
            ticks = [0.5, 1, 2, 5, 10, 20, 50] if row == 0 else [200, 400, 1000]
            ax.set_xticks(ticks); ax.set_xticklabels([f"{t:g}" for t in ticks]); ax.minorticks_off()
            ax.set_xlim((0.3, 70) if row == 0 else (150, 1400))
            ax.set_ylim(0, 100); ax.set_xlabel(xlab, fontsize=9, color=INK2); ax.set_ylabel(f"{el} right place (%)", fontsize=9, color=INK2)
            ax.grid(True, color="#e6e5e1", lw=0.6); ax.set_axisbelow(True)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            ax.tick_params(colors=INK2, labelsize=8.5)
    axs[0, 0].plot([], [], "o", ms=9, color=BLUE, mec="white", label="right element, right depth (±0.5 Å)")
    axs[0, 0].plot([], [], "o", ms=8, mfc="white", mec=GREY, mew=1.6, label="guessing (depths slid along the column)")
    axs[0, 0].legend(fontsize=8, frameon=False, loc="lower right")
    fig.suptitle("detector256: what a smaller camera costs (80 mrad, ARM tableau, fs 2.5 Å + src 0.1 Å, true probe)",
                 fontsize=10.5, color=INK, x=0.01, ha="left")
    fig.text(0.01, 0.005, KEYS, fontsize=7.5, color=INK2, ha="left", va="bottom", wrap=True)
    fig.tight_layout(rect=(0, 0.03, 0.93, 0.97))
    arm.save(fig, "det_limits.png")


def main():
    os.makedirs(OUT, exist_ok=True)
    view = build_view(PULLS)
    man = json.load(open(arm.MANIFEST))
    man["search_dirs"] = man["search_dirs"] + [view]
    for key, name, lab, tag, *_ in RUNS:
        man["runs"][key] = dict(tag=[f"{lab}{tag}"], leg=f"{lab}_lab{tag}", short=key, name=name)
    R = arm.Runs(man); arm._RUNS[0] = R
    out = score_all(R)
    print_table(out)
    save_tables(out)
    copy_phase(R, out)
    sections(R, out, arm.cuts_from(arm.load_leg(R.find("known_a080"))))
    limits_figure(out)
    shutil.rmtree(view, ignore_errors=True)
    print(f"saved in {os.path.relpath(OUT, REPO)}/: det_scores.json, det_scores.csv, det_sections_*.png, det_limits.png, phase/")


if __name__ == "__main__":
    main()
