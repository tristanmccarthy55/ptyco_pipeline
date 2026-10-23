"""@file det_look.py
@brief detector256: the honest numbers and depth sections for R1-R4 beside their reference, from the pulled runs.

Same machinery as aberration_experiment/mono300/phase1_look.py (analysis/make_arm_page.py, unchanged): right place (right
element within 0.5 A in depth) against guessing (atom_scores), the verdict rule (good: Pb >= 80 and Ti >= 45 %; okay:
Pb >= 50 %), the vacuum loading (rms phase per slice in the vacuum bands over that in the slab; 0.31 for a clean run) and
the depth sections. The misfit eps = final Fourier error x N / 2e5 with N the width of the DATA each run fitted (its own
sim_meta.mat: 1420 reference, 1024 R1, 256 R2/R3, 354 R4), not the width in the error-trace file name, which for an
upsampled run is the model's. eps compares like with like only on the same data grid (R2 with R3, R1 with R4's crop...).

    ~/hyperspy-bundle/bin/python aberration_experiment/detector256/det_look.py
Runs are looked for in ~/Desktop/mono_phase1 (the reference) and ~/Desktop/detector256_{r1r2,r3,r4} (submit.sh's pull
folders). Writes aberration_experiment/detector256/det_scores.json and aberration_experiment/figs/detector256/det_sections.png.
"""
import glob
import importlib.util
import json
import os
import re
import sys

import numpy as np
from scipy.io import loadmat

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
spec = importlib.util.spec_from_file_location("arm", os.path.join(REPO, "analysis", "make_arm_page.py"))
arm = importlib.util.module_from_spec(spec); sys.modules["arm"] = arm; spec.loader.exec_module(arm)
spec = importlib.util.spec_from_file_location("p1", os.path.join(REPO, "aberration_experiment", "mono300", "phase1_look.py"))
p1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(p1)
plt = arm.plt
arm.FIGS = os.path.join(REPO, "aberration_experiment", "figs", "detector256")
SEARCH = ["~/Desktop/mono_phase1", "~/Desktop/detector256_r1r2", "~/Desktop/detector256_r3", "~/Desktop/detector256_r4"]
CORE = "fs2p5_src0p1"
# (key, display name, recon tag after the label, camera, model grid)
RUNS = [("ref", "reference: +-133 mrad, 1420 px", f"{CORE}_fixed", "1420 px, +-133 mrad", "1420 px"),
        ("R1", "R1: +-96 mrad, full pixels", f"{CORE}_detc96", "1024 px, +-96 mrad", "1024 px"),
        ("R2", "R2: 256 px, upsampling 2", f"{CORE}_detc96b4_up2", "256 px, +-96 mrad", "1024 px (4x4 sums)"),
        ("R3", "R3: 256 px, no upsampling", f"{CORE}_detc96b4", "256 px, +-96 mrad", "256 px (26 A window)"),
        ("R4", "R4: 354 px, upsampling 2", f"{CORE}_detc133b4_up2", "354 px, +-133 mrad", "1416 px (4x4 sums)")]


def eps_of(R, leg):
    """(final eps, data width N) from the leg's error trace and the sim_meta.mat of the data it fitted."""
    t = R._newest([p for d in R.dirs for p in glob.glob(os.path.join(d, "results", f"recon_af_{leg}_NL*", "analysis", "**",
                                                                     "*_error_trace.csv"), recursive=True)])
    sim = re.sub(r"_up\d+$", "", leg.replace(f"_{CORE}_fixed", f"_{CORE}"))       # the data dir (a copy, or the kept sim)
    m = R._newest([p for d in R.dirs for p in (os.path.join(d, "results", f"sim_out_af_{sim}", "01", "sim_meta.mat"),)])
    if not t or not m:
        return None, None
    e = [float(l.split(",")[1]) for l in open(t) if l[:1].isdigit()]
    n = int(loadmat(m, simplify_cells=True)["meta"]["Ndpx"])
    return e[-1] * n / 2e5, n


def main():
    man = json.load(open(arm.MANIFEST))
    man["search_dirs"] = man["search_dirs"] + SEARCH
    for key, name, tag, *_ in RUNS:
        man["runs"][key] = dict(tag=[f"armf_a080_{tag}"], leg=f"armf_a080_lab_{tag}", short=key, name=name)
    R = arm.Runs(man); arm._RUNS[0] = R
    cuts = arm.cuts_from(arm.load_leg(R.find("known_a080")))
    out = []
    print(f"{'run':34s}{'right Pb/Ti/O %':>17s}  {'guess Pb/Ti/O %':>17s}  {'wrong':>5s}  {'vac/slab':>8s}  {'eps (N)':>13s}  verdict")
    f3 = lambda d: "/".join(f"{100 * d[s]:3.0f}" for s in ("Pb", "Ti", "O"))
    for key, name, tag, cam, grid in RUNS:
        info = R.find(key)
        sc = arm.atom_scores(info) if info["found"] else None
        eps, n = eps_of(R, f"armf_a080_lab_{tag}") if info["found"] else (None, None)
        r = dict(key=key, run=name, camera=cam, model_grid=grid, found=info["found"], eps=eps, data_px=n,
                 vacuum_ratio=p1.vacuum_ratio(info.get("fd")))
        if sc:
            r.update(right={k: sc["right"][z] for k, z in (("Pb", 82), ("Ti", 22), ("O", 8))},
                     guess={k: sc["guess"][z] for k, z in (("Pb", 82), ("Ti", 22), ("O", 8))},
                     wrong=sc["wrong"], bar=sc["bar"], n=sc["n"], verdict=arm.verdict(sc))
        out.append(r)
        vr = f"{r['vacuum_ratio']:8.2f}" if r["vacuum_ratio"] is not None else "       -"
        ep = f"{eps:.4f} ({n})" if eps is not None else "-"
        if sc:
            print(f"{name:34s}{f3(r['right']):>17s}  {f3(r['guess']):>17s}  {100 * r['wrong']:4.0f}%  {vr}  {ep:>13s}  {r['verdict']}")
        else:
            print(f"{name:34s}{'not pulled' if not info['found'] else 'no atoms (a leg or kernel missing)':>37s}         {vr}  {ep:>13s}")
    json.dump(out, open(os.path.join(HERE, "det_scores.json"), "w"), indent=1)
    keys = [k for k, *_ in RUNS if arm.load_leg(R.find(k)) is not None]
    fig, axs = plt.subplots(2, len(keys), figsize=(2.3 * len(keys), 8.6), squeeze=False)
    for j, k in enumerate(keys):
        info = R.find(k); lg = arm.load_leg(info); sc = arm.atom_scores(info)
        r = next(x for x in out if x["key"] == k)
        for i, (cut, rA) in enumerate(cuts):
            ax = axs[i, j]
            nx = lg["V"].shape[2]
            arm.xz_section(ax, lg, sc, rA / lg["dx"], xlim=(nx * lg["dx"] / 2 - 8, nx * lg["dx"] / 2 + 8), show_y=(j == 0))
            if i == 0:
                ttl = r["run"].replace(": ", ":\n") + f"\nmodel {r['model_grid']}"
                if "right" in r:
                    ttl += f"\nPb {100 * r['right']['Pb']:.0f} vs {100 * r['guess']['Pb']:.0f}" \
                           f"\nTi {100 * r['right']['Ti']:.0f} vs {100 * r['guess']['Ti']:.0f}\n{r['verdict']}"
                else:
                    ttl += "\n(no atoms)"
                ax.set_title(ttl, fontsize=8.5, loc="left")
            ax.set_xlabel(f"x (Å) · {cut}", fontsize=8)
    fig.tight_layout()
    arm.save(fig, "det_sections.png")


if __name__ == "__main__":
    main()
