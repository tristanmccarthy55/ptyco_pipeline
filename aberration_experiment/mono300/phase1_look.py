"""@file phase1_look.py
@brief Phase 1 (focal spread x source size): honest atom scores per run and depth-section sheets, from the pulled runs.

Reuses the ARM page's machinery unchanged (analysis/make_arm_page.py): runs found by tag in the search folders, the
"right place vs guessing" scores (atom_scores), the verdict rule (good: Pb >= 80 and Ti >= 45 %; okay: Pb >= 50 %), the
depth sections (xz_section, the two cuts of the coherent reference). Each run's coherence is also given as the share of
the beam in its first mode (the sim's own 16-configuration state, read from mode_spectrum.txt), the vacuum loading as
the rms phase per slice in the two vacuum bands over that in the slab, and for a probe-mode leg the start and recovered
mode powers (kick/*.json).
  step 1  ten fixed-probe controls, one effect at a time (_fixed)
  step 2  joint points (_fs<X>_src<Y>_fixed) and probe-mode legs started from the mixed state (_mix, old engine)
  _mixF   the four single wall points' mode legs rerun with the engine fix (all-modes illumination; rerunF.sh)
  A/B/C   the ceiling, cooling and searched start at the Spectra Mono settings (C: first 2-stage search, then 3-stage _s3)
Runs not pulled (or without atoms) are listed as such and drawn without rings.

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
BOXZ, ZVAC = 27.525, 4.0

STEP1 = [("fs", "7.29"), ("fs", "5"), ("fs", "3.5"), ("fs", "2.5"), ("fs", "1.5"), ("fs", "0.75"),
         ("src", "0.05"), ("src", "0.1"), ("src", "0.2"), ("src", "0.4")]
JOINT = [("2.5", "0.1"), ("2.5", "0.2"), ("3.5", "0.1"), ("3.5", "0.2")]
MODES = [("fs", "7.29"), ("fs", "5"), ("src", "0.2"), ("src", "0.4")]
# the ceiling and the searched start at the realistic Spectra Mono settings (ceiling_submit.sh, search_submit.sh; 2026-10-16):
# (section, display name, key, recon tag core, mode_spectrum name, fs, src) -- tags armf_a080_<core>, legs armf_a080_lab_<core>
CEIL = [("A ceiling", "fs 2.5 + src 0.1, modes", "fs2p5_src0p1_mixF", "fs2p5_src0p1_mixF", "fs 2.5+src 0.1", 2.5, 0.1),
        ("A ceiling", "fs 2.5 + src 0.2, modes", "fs2p5_src0p2_mixF", "fs2p5_src0p2_mixF", "fs 2.5+src 0.2", 2.5, 0.2),
        ("B cooled", "100 K, fixed", "T100_fixed", "ph16_T100_fs2p5_src0p1_mixF_fixed", "fs 2.5+src 0.1", 2.5, 0.1),
        ("B cooled", "100 K, modes", "T100_mixF", "ph16_T100_fs2p5_src0p1_mixF", "fs 2.5+src 0.1", 2.5, 0.1),
        ("B cooled", "20 K, fixed", "T20_fixed", "ph16_T20_fs2p5_src0p1_mixF_fixed", "fs 2.5+src 0.1", 2.5, 0.1),
        ("B cooled", "20 K, modes", "T20_mixF", "ph16_T20_fs2p5_src0p1_mixF", "fs 2.5+src 0.1", 2.5, 0.1),
        ("C searched", "fs 2.5 + src 0.1, searched start, modes", "fs2p5_src0p1_ps", "fs2p5_src0p1_ps", "fs 2.5+src 0.1", 2.5, 0.1),
        ("C searched", "3-stage search, modes", "fs2p5_src0p1_s3_ps", "fs2p5_src0p1_s3_ps", "fs 2.5+src 0.1", 2.5, 0.1)]
p_ = lambda v: v.replace(".", "p")


def mode1_share():
    """{point name in mode_spectrum.txt: first-mode share of the sim's 16-configuration state}."""
    out = {}
    for line in open(os.path.join(HERE, "mode_spectrum.txt")):
        m = re.match(r"^(.+?)\s+sim16\s+([0-9.]+)", line)
        if m:
            out[m.group(1).strip()] = float(m.group(2))
    return out


def vacuum_ratio(fd):
    """Mean rms phase per slice in the two 4 A vacuum bands over that in the slab, from a run's figure-data volume."""
    if not fd or not os.path.exists(os.path.join(fd, "phase_vol.npy")):
        return None
    V = np.load(os.path.join(fd, "phase_vol.npy")); dz = json.load(open(os.path.join(fd, "phase_vol.json")))["dz"]
    zc = (np.arange(V.shape[0]) + 0.5) * dz
    s = np.array([np.std(V[i]) for i in range(V.shape[0])])
    return float((s[zc < ZVAC].mean() + s[zc > BOXZ - ZVAC].mean()) / 2 / s[(zc > ZVAC) & (zc < BOXZ - ZVAC - 0.025)].mean())


def runs():
    """[(section, display name, manifest key, run spec, mode_spectrum name, fs, src)]"""
    out = [("step 1", "coherent", "known_a080", None, None, 0.0, 0.0)]
    for k, v in STEP1:
        t = f"{k}{p_(v)}"
        out.append(("step 1", f"{k} {v}", t, {"tag": [f"armf_a080_{t}_fixed"], "leg": f"armf_a080_lab_{t}_fixed"},
                    f"{k} {v}", float(v) if k == "fs" else 0.0, float(v) if k == "src" else 0.0))
    out.append(("step 1", "_coh (fs 7.29 + src 0.4)", "split_coh", None, "_coh 7.29+0.4", 7.29, 0.4))
    for f, s in JOINT:
        t = f"fs{p_(f)}_src{p_(s)}"
        out.append(("step 2 joint", f"fs {f} + src {s}", t, {"tag": [f"armf_a080_{t}_fixed"], "leg": f"armf_a080_lab_{t}_fixed"},
                    f"fs {f}+src {s}", float(f), float(s)))
    for k, v in MODES:
        t = f"{k}{p_(v)}"
        out.append(("step 2 modes", f"{k} {v}, modes", f"{t}_mix", {"tag": [f"armf_a080_{t}_mix"], "leg": f"armf_a080_lab_{t}_mix"},
                    f"{k} {v}", float(v) if k == "fs" else 0.0, float(v) if k == "src" else 0.0))
    for k, v in MODES:
        t = f"{k}{p_(v)}"
        out.append(("step 2 modes, engine fixed", f"{k} {v}, modes", f"{t}_mixF",
                    {"tag": [f"armf_a080_{t}_mixF"], "leg": f"armf_a080_lab_{t}_mixF"},
                    f"{k} {v}", float(v) if k == "fs" else 0.0, float(v) if k == "src" else 0.0))
    out.append(("B cooled", "phonons only, RT (ARM page)", "split_ph", None, None, 0.0, 0.0))
    out.append(("B cooled", "phonons only, 100 K (ARM page)", "split_ph100", None, None, 0.0, 0.0))
    for sec, name, key, core, wkey, fs, src in CEIL:
        out.append((sec, name, key, {"tag": [f"armf_a080_{core}"], "leg": f"armf_a080_lab_{core}"}, wkey, fs, src))
    return out


def main():
    man = json.load(open(arm.MANIFEST))
    man["search_dirs"] = man["search_dirs"] + ["~/Desktop/mono_phase1"]
    table = runs()
    for _, name, key, spec_, *_ in table:
        if spec_ is not None:
            man["runs"][key] = dict(spec_, short=key, name=name)
    R = arm.Runs(man); arm._RUNS[0] = R
    w1 = mode1_share()
    cuts = arm.cuts_from(arm.load_leg(R.find("known_a080")))
    out, sec = [], None
    f3 = lambda d: "/".join(f"{100 * d[s]:3.0f}" for s in ("Pb", "Ti", "O"))
    for section, name, key, _, wkey, fs, src in table:
        if section != sec:
            sec = section
            print(f"\n== {section}\n{'point':26s}{'mode-1':>7s}  {'right Pb/Ti/O %':>17s}  {'guess Pb/Ti/O %':>17s}  "
                  f"{'wrong':>5s}  {'vac/slab':>8s}  {'eps':>7s}  verdict   mode powers (start -> recovered)")
        info = R.find(key)
        sc = arm.atom_scores(info) if info["found"] else None
        tr = R.trace(info["spec"]["leg"]) if info["found"] else None
        kj, _ = R.kick(info["spec"]["leg"]) if info["found"] else (None, None)
        r = dict(section=section, point=name, key=key, focal_spread_A=fs, source_fwhm_A=src,
                 mode1_share=1.0 if key == "known_a080" else w1.get(wkey), found=info["found"],
                 eps=float(tr[1][-1]) if tr is not None else None, vacuum_ratio=vacuum_ratio(info.get("fd")))
        if kj and kj.get("mode_power"):
            r.update(start_mode_power=kj.get("start_mode_power"), mode_power=kj["mode_power"])
        if sc:
            r.update(right={k: sc["right"][z] for k, z in (("Pb", 82), ("Ti", 22), ("O", 8))},
                     guess={k: sc["guess"][z] for k, z in (("Pb", 82), ("Ti", 22), ("O", 8))},
                     wrong=sc["wrong"], bar=sc["bar"], n=sc["n"], verdict=arm.verdict(sc))
        out.append(r)
        w = f"{r['mode1_share']:7.3f}" if r["mode1_share"] is not None else "      -"
        vr = f"{r['vacuum_ratio']:8.2f}" if r["vacuum_ratio"] is not None else "       -"
        ep = f"{r['eps']:7.4f}" if r["eps"] is not None else "      -"
        mp = ""
        if r.get("mode_power"):
            mp = (" ".join(f"{x:.2f}" for x in (r.get("start_mode_power") or [])[:4]) + " -> "
                  + " ".join(f"{x:.2f}" for x in r["mode_power"][:4]))
        if sc:
            print(f"{name:26s}{w}  {f3(r['right']):>17s}  {f3(r['guess']):>17s}  {100 * r['wrong']:4.0f}%  {vr}  {ep}  "
                  f"{r['verdict']:8s}  {mp}")
        else:
            why = "not pulled" if not info["found"] else "no atoms (a leg or kernel missing)"
            print(f"{name:26s}{w}  {why:>37s}         {vr}  {ep}  {'':8s}  {mp}")
    json.dump(out, open(os.path.join(HERE, "phase1_scores.json"), "w"), indent=1)
    sheets = (("focal", ["known_a080", "fs0p75", "fs1p5", "fs2p5", "fs3p5", "fs5", "fs7p29", "split_coh"]),
              ("source", ["known_a080", "src0p05", "src0p1", "src0p2", "src0p4", "split_coh"]),
              ("joint", ["known_a080", "fs2p5", "fs2p5_src0p1", "fs2p5_src0p2", "fs3p5", "fs3p5_src0p1", "fs3p5_src0p2",
                         "src0p1", "src0p2"]),
              ("modes", ["fs7p29", "fs7p29_mix", "fs5", "fs5_mix", "src0p2", "src0p2_mix", "src0p4", "src0p4_mix"]),
              ("modesF", ["fs7p29", "fs7p29_mix", "fs7p29_mixF", "fs5", "fs5_mixF", "src0p2", "src0p2_mixF", "src0p4",
                          "src0p4_mixF"]),
              ("ceiling", ["fs2p5_src0p1", "fs2p5_src0p1_mixF", "fs2p5_src0p2", "fs2p5_src0p2_mixF", "T100_fixed", "T100_mixF",
                           "T20_fixed", "T20_mixF", "fs2p5_src0p1_ps", "fs2p5_src0p1_s3_ps"]))
    for series, keys in sheets:
        keys = [k for k in keys if arm.load_leg(R.find(k)) is not None]
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
                    ttl = r["point"].replace(" (fs 7.29 + src 0.4)", "\nfs 7.29 + src 0.4").replace(", modes", "\nprobe modes")
                    ttl += " (old engine)" if series == "modesF" and k.endswith("_mix") else ""
                    ttl += "\nfixed probe" if (series in ("modes", "modesF") and not k.endswith(("_mix", "_mixF"))) or \
                                               (series == "ceiling" and k in ("fs2p5_src0p1", "fs2p5_src0p2")) else ""
                    ttl += f"\nmode 1 {r['mode1_share']:.2f}" if r["mode1_share"] is not None else ""
                    if "right" in r:
                        ttl += f"\nPb {100 * r['right']['Pb']:.0f} vs {100 * r['guess']['Pb']:.0f}" \
                               f"\nTi {100 * r['right']['Ti']:.0f} vs {100 * r['guess']['Ti']:.0f}\n{r['verdict']}"
                    else:
                        ttl += "\n(no atoms)"
                    ax.set_title(ttl, fontsize=8.5, loc="left")
                ax.set_xlabel(f"x (Å) · {cut}", fontsize=8)
        fig.tight_layout()
        arm.save(fig, f"phase1_{series}.png")


if __name__ == "__main__":
    main()
