#!/usr/bin/env python
"""
@file relaxation_ladder.py
@brief Append one row per (step, label, alpha) to the relaxation-ladder table, so every step of
       aberration_experiment/NEXT_PHASE.md records the same numbers in one place.

Per atomfind run (a dir holding report.json): precision; whole-slab and bulk recall per species; xy-RMS;
z-RMS; species confusion; false positives in total and by species; and, for a fitted probe, the recovered
C1 from analysis/c1_objective.py's summary JSON. The atomfind numbers come through
collate_atomfind_depth.load_run, the same reader the round-alpha figure uses, so the two cannot drift.
False positives by species are the confusion matrix's "<species>->none" entries (found as that species,
matched to no ground-truth atom); their total is n_found - n_matched.

A row with the same (step, label, alpha) is replaced, never duplicated, so re-running after a fix is safe.

    python analysis/relaxation_ladder.py --step 0 --label "known probe, noiseless" \\
        --atomfind ~/Desktop/thin_ab_af_final/out/atomfind_a70 ~/Desktop/thin_ab_af_final/out/atomfind_a90
    python analysis/relaxation_ladder.py --step 1 --label "C1 fitted, C3/C5 fixed" \\
        --atomfind <out>/atomfind_a70 --c1-summary aberration_experiment/results/<week>/c1_objective_summary.json
"""
from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from collate_atomfind_depth import SPECIES, load_run, week_dir   # noqa: E402

COLS = (["step", "label", "alpha", "date", "precision"]
        + [f"{k}_recall" for k in SPECIES] + [f"{k}_recall_bulk" for k in SPECIES]
        + ["xy_rms", "z_rms"] + [f"{k}_zrms" for k in SPECIES]
        + ["confusion", "n_found", "n_matched", "fp_total", "fp_Pb", "fp_Ti", "fp_O",
           "c1_true", "c1_fit", "c1_halfwidth", "c1_within", "c1_source", "note", "report"])
FP_KEY = {"Pb": "82->none", "Ti": "22->none", "O": "8->none"}


def run_alpha(run_dir, row):
    """Alpha from the atomfind OUT-DIR name when it carries one (atomfind_a90_dose1e5). collate's
    find_alpha scans the whole recon path, which for a combined tarball (atomfind_results_a70-90_...)
    matches the tarball's first alpha instead -- it labelled every a90 dose run a70 (2026-09-18)."""
    m = re.search(r"atomfind_a0*(\d{2,3})", os.path.basename(os.path.abspath(run_dir)))
    if not m:
        return row["alpha"]
    a = int(m.group(1))
    if row["alpha"] != a:
        print(f"  NOTE {os.path.basename(run_dir)}: alpha {a} from the out-dir name, not {row['alpha']} from the recon path")
    return a


def ladder_row(run_dir, step, label, note, c1):
    row = load_run(run_dir)
    row["alpha"] = run_alpha(run_dir, row)
    rp = run_dir if run_dir.endswith(".json") else os.path.join(run_dir, "report.json")
    v3 = json.load(open(rp)).get("finder", {}).get("v3", {})
    cf = v3.get("confusion") or {}
    row.update(step=step, label=label, date=datetime.date.today().isoformat(), note=note or "",
               report=os.path.abspath(rp), n_matched=v3.get("n_matched"),
               fp_total=(v3["n_found"] - v3["n_matched"]) if "n_matched" in v3 else None,
               **{f"fp_{k}": cf.get(key) for k, key in FP_KEY.items()})
    if c1:
        row.update(c1_true=c1.get("c1_true"), c1_fit=c1.get("c1_fit"), c1_halfwidth=c1.get("halfwidth"),
                   c1_within=c1.get("within"), c1_source=c1.get("_key"))
    return row


def pick_c1(summary, alpha, mode, niter):
    """The summary entry for this alpha: a<alpha>_<mode>_n<NITER>, the given NITER or else the smallest."""
    cands = {k: v for k, v in summary.items()
             if v.get("alpha") == alpha and v.get("mode") == mode and "niter" in v
             and (niter is None or v["niter"] == niter)}
    if not cands:
        return None
    k = min(cands, key=lambda k: cands[k]["niter"])
    return dict(cands[k], _key=k)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--step", required=True, help="ladder step (0 = baseline, 1 = C1 fit, ...)")
    ap.add_argument("--label", required=True, help="what this step relaxes, e.g. 'C1 fitted, C3/C5 fixed'")
    ap.add_argument("--atomfind", nargs="+", required=True, help="atomfind out dirs (report.json inside)")
    ap.add_argument("--c1-summary", default=None, help="c1_objective_summary.json (fitted-probe steps)")
    ap.add_argument("--c1-mode", default="lab"); ap.add_argument("--c1-niter", type=int, default=None)
    ap.add_argument("--note", default=None)
    ap.add_argument("--csv", default=None, help="default aberration_experiment/results/<week>/relaxation_ladder.csv")
    a = ap.parse_args()

    summary = json.load(open(a.c1_summary)) if a.c1_summary else {}
    new = []
    for d in a.atomfind:
        d = os.path.expanduser(d)
        alpha = run_alpha(d, load_run(d))
        c1 = pick_c1(summary, alpha, a.c1_mode, a.c1_niter) if summary else None
        if summary and c1 is None:
            print(f"  WARNING a{alpha}: no C1 estimate in {a.c1_summary} for mode {a.c1_mode}")
        new.append(ladder_row(d, str(a.step), a.label, a.note, c1))

    path = a.csv or os.path.join(week_dir("results"), "relaxation_ladder.csv")
    old = list(csv.DictReader(open(path))) if os.path.isfile(path) else []
    keys = {(r["step"], r["label"], str(r["alpha"])) for r in new}
    kept = [r for r in old if (r["step"], r["label"], str(r["alpha"])) not in keys]
    rows = sorted(kept + new, key=lambda r: (float(r["step"]), int(r["alpha"]), r["label"]))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, extrasaction="ignore"); w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c) for c in COLS})
    print(f"wrote {path}: {len(new)} row(s) for step {a.step} ({len(old) - len(kept)} replaced), {len(rows)} total")
    pct = lambda x: "-" if x in (None, "") else f"{100 * float(x):.0f}"
    for r in new:
        fp = ", ".join(f"{k} {r[f'fp_{k}']}" for k in SPECIES if r.get(f"fp_{k}"))
        c1s = (f" | C1 {r['c1_fit']:g} +- {r['c1_halfwidth']:g} (true {r['c1_true']:g})" if r.get("c1_fit") is not None else "")
        print(f"  a{r['alpha']}: precision {float(r['precision']):.2f} | bulk recall Pb/Ti/O "
              f"{pct(r['Pb_recall_bulk'])}/{pct(r['Ti_recall_bulk'])}/{pct(r['O_recall_bulk'])}% | "
              f"xy {float(r['xy_rms']):.2f} A z {float(r['z_rms']):.2f} A | confusion {100 * float(r['confusion']):.1f}% | "
              f"FP {r['fp_total']} of {r['n_found']} ({fp}){c1s}")


if __name__ == "__main__":
    main()
