"""@file phase1_look.py
@brief Phase 1 step 1 (fixed-probe controls): honest atom scores per point and depth-section sheets, from the pulled runs.

Reuses the ARM page's machinery unchanged (analysis/make_arm_page.py): runs found by tag in the search folders, the
"right place vs guessing" scores (atom_scores), the verdict rule (good: Pb >= 80 and Ti >= 45 %; okay: Pb >= 50 %), the
depth sections (xz_section, the two cuts of the coherent reference). Each point's coherence is also given as the share
of the beam in its first mode (the sim's 16-configuration state, read from mode_spectrum.txt), so the focal-spread and
source-size series can be read on one axis.

    ~/hyperspy-bundle/bin/python aberration_experiment/mono300/phase1_look.py     # pulled to ~/Desktop/mono_phase1
Writes aberration_experiment/mono300/phase1_scores.json and aberration_experiment/figs/mono300/phase1_*.png.
"""
import importlib.util
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
spec = importlib.util.spec_from_file_location("arm", os.path.join(REPO, "analysis", "make_arm_page.py"))
arm = importlib.util.module_from_spec(spec); sys.modules["arm"] = arm; spec.loader.exec_module(arm)
plt = arm.plt
arm.FIGS = os.path.join(REPO, "aberration_experiment", "figs", "mono300")

POINTS = [("fs", "7.29", "fs 7.29"), ("fs", "5", "fs 5"), ("fs", "3.5", "fs 3.5"), ("fs", "2.5", "fs 2.5"), ("fs", "1.5", "fs 1.5"),
          ("fs", "0.75", "fs 0.75"), ("src", "0.05", "src 0.05"), ("src", "0.1", "src 0.1"), ("src", "0.2", "src 0.2"),
          ("src", "0.4", "src 0.4")]


def mode1_share():
    """{point name in mode_spectrum.txt: first-mode share of the sim's 16-configuration state}."""
    out = {}
    for line in open(os.path.join(HERE, "mode_spectrum.txt")):
        m = re.match(r"^(\S+ \S+|_coh \S+)\s+sim16\s+([0-9.]+)", line)
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


def main():
    man = json.load(open(arm.MANIFEST))
    man["search_dirs"] = man["search_dirs"] + ["~/Desktop/mono_phase1"]
    for kind, v, _ in POINTS:
        t = f"{kind}{v.replace('.', 'p')}"
        man["runs"][t] = {"tag": [f"armf_a080_{t}_fixed"], "leg": f"armf_a080_lab_{t}_fixed", "short": t, "name": t}
    R = arm.Runs(man); arm._RUNS[0] = R
    w1 = mode1_share()
    ref = R.find("known_a080")
    cuts = arm.cuts_from(arm.load_leg(ref))
    rows = [("coherent", 0.0, 0.0, 1.0, "known_a080")]
    rows += [(name, float(v) if k == "fs" else 0.0, float(v) if k == "src" else 0.0, w1.get(name), f"{k}{v.replace('.', 'p')}")
             for k, v, name in POINTS]
    rows += [("_coh (fs 7.29 + src 0.4)", 7.29, 0.4, w1.get("_coh 7.29+0.4"), "split_coh")]
    out = []
    print(f"{'point':26s}{'mode-1':>7s}  {'right Pb/Ti/O %':>17s}  {'guess Pb/Ti/O %':>17s}  {'wrong':>5s}  {'bar A':>5s}  "
          f"{'eps':>7s}  verdict")
    for name, fs, src, w, key in rows:
        info = R.find(key)
        sc = arm.atom_scores(info) if info["found"] else None
        tr = R.trace(info["spec"]["leg"])
        eps = float(tr[1][-1]) if tr is not None else None
        r = dict(point=name, focal_spread_A=fs, source_fwhm_A=src, mode1_share=w, key=key, found=info["found"], eps=eps)
        if sc:
            r.update(right={k: sc["right"][z] for k, z in (("Pb", 82), ("Ti", 22), ("O", 8))},
                     guess={k: sc["guess"][z] for k, z in (("Pb", 82), ("Ti", 22), ("O", 8))},
                     wrong=sc["wrong"], bar=sc["bar"], n=sc["n"], verdict=arm.verdict(sc))
            f = lambda d: "/".join(f"{100 * d[s]:3.0f}" for s in ("Pb", "Ti", "O"))
            print(f"{name:26s}{w if w is not None else float('nan'):7.3f}  {f(r['right']):>17s}  {f(r['guess']):>17s}  "
                  f"{100 * r['wrong']:4.0f}%  {r['bar']:5.2f}  {eps if eps is not None else float('nan'):7.4f}  {r['verdict']}")
        else:
            print(f"{name:26s}  not scored (found={info['found']})")
        out.append(r)
    json.dump(out, open(os.path.join(HERE, "phase1_scores.json"), "w"), indent=1)
    # depth-section sheets: one column per run, the Pb row and the B-site row of the coherent reference's cuts
    for series, keys in (("focal", ["known_a080", "fs0p75", "fs1p5", "fs2p5", "fs3p5", "fs5", "fs7p29", "split_coh"]),
                         ("source", ["known_a080", "src0p05", "src0p1", "src0p2", "src0p4", "split_coh"])):
        n = len(keys)
        fig, axs = plt.subplots(2, n, figsize=(2.05 * n, 8.6), squeeze=False)
        for j, k in enumerate(keys):
            info = R.find(k); lg = arm.load_leg(info); sc = arm.atom_scores(info)
            r = next(x for x in out if x["key"] == k)
            for i, (cut, rA) in enumerate(cuts):
                ax = axs[i, j]
                nx = lg["V"].shape[2]
                arm.xz_section(ax, lg, sc, rA / lg["dx"], xlim=(nx * lg["dx"] / 2 - 8, nx * lg["dx"] / 2 + 8), show_y=(j == 0))
                if i == 0:
                    ttl = r["point"].replace(" (fs 7.29 + src 0.4)", "\nfs 7.29 + src 0.4")
                    ttl += f"\nmode 1 {r['mode1_share']:.2f}" if r["mode1_share"] is not None else ""
                    if "right" in r:
                        ttl += f"\nPb {100 * r['right']['Pb']:.0f} vs {100 * r['guess']['Pb']:.0f}" \
                               f"\nTi {100 * r['right']['Ti']:.0f} vs {100 * r['guess']['Ti']:.0f}\n{r['verdict']}"
                    ax.set_title(ttl, fontsize=8.5, loc="left")
                ax.set_xlabel(f"x (Å) · {cut}", fontsize=8)
        fig.tight_layout()
        arm.save(fig, f"phase1_{series}.png")


if __name__ == "__main__":
    main()
