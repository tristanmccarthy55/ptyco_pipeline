#!/usr/bin/env python
"""@file collate_atomfind_depth.py
@brief Depth-localisation-vs-alpha figure for the aberration atomfind sweep (the headline result).

run_atomfind.py writes one report.json per run; run it once per convergence angle (each into its
own out-dir), then this reads them all and plots the blind finder (v3) per-species recall and
depth (z) RMS against alpha -- showing depth recovery IMPROVING as the aperture opens past the
Cs-corrector spec (A-site by ~50 mrad, the BO2/O column only by ~100), and degrading again once
the probe blows up. Also drops a small CSV summary.

  python analysis/collate_atomfind_depth.py results/2026-W37/atomfind_a50 results/.../atomfind_a100 ...
  python analysis/collate_atomfind_depth.py --glob 'results/2026-W37/atomfind_a*'

Each argument is a run dir holding report.json (or a report.json path). Alpha is read from
report['vol'] (…a<NN>_lab_NL…) or the dir name (a<NN>). Figure defaults to
aberration_experiment/figs/<ISO-week>/atomfind_depth_vs_alpha.png and the CSV to
aberration_experiment/results/<ISO-week>/atomfind_depth_summary.csv (committed, dated).
"""
import argparse, csv, datetime, glob, json, os, re

LAMBDA_A = 0.0196877          # 300 keV electron wavelength [Å]
SPECIES = ("Pb", "Ti", "O")   # A-site / B-site cation / oxygen
SP_COLOR = {"Pb": "#1f77b4", "Ti": "#d62728", "O": "#2ca02c"}


def week_dir(sub):
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    d = os.path.join(repo, "aberration_experiment", sub, datetime.date.today().strftime("%G-W%V"))
    os.makedirs(d, exist_ok=True)
    return d


def find_alpha(report, path):
    """Alpha (mrad) from the recon-volume name or the run-dir path (a050 / a50 / a100)."""
    for s in (str(report.get("vol", "")), path):
        m = re.search(r"a0*(\d{2,3})", s)
        if m:
            return int(m.group(1))
    return None


def load_run(d):
    p = d if d.endswith(".json") else os.path.join(d, "report.json")
    if not os.path.isfile(p):
        raise SystemExit(f"no report.json at {p}")
    with open(p) as f:
        r = json.load(f)
    frep = r.get("finder", {}).get("v3", {})
    row = dict(alpha=find_alpha(r, p), dz=r.get("dz"),
               precision=frep.get("precision"),
               xy_rms=frep.get("xy_rms_A"), z_rms=frep.get("z_rms_A"),
               n_found=frep.get("n_found"))
    for k in SPECIES:
        sp = frep.get(k, {})
        row[f"{k}_recall"] = sp.get("recall")
        row[f"{k}_recall_bulk"] = sp.get("recall_bulk")
        row[f"{k}_zrms"] = sp.get("z_rms_A")
    if row["alpha"] is None:
        raise SystemExit(f"could not parse alpha from {p} (vol={r.get('vol')!r})")
    return row


def write_csv(rows, path):
    cols = ["alpha", "dz", "precision", "xy_rms", "z_rms", "n_found"]
    cols += [f"{k}_{m}" for k in SPECIES for m in ("recall", "recall_bulk", "zrms")]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c) for c in cols})


def _xy(rows, xk, yk):
    """(alpha, value) pairs where value is present -- so a missing/failed alpha just drops out."""
    xs, ys = [], []
    for r in rows:
        if r.get(yk) is not None:
            xs.append(r[xk]); ys.append(r[yk])
    return xs, ys


def make_fig(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (axR, axZ) = plt.subplots(1, 2, figsize=(12.5, 5.2), constrained_layout=True)
    alphas = [r["alpha"] for r in rows]

    # --- recall vs alpha (bulk = interior planes, the honest metric; total dashed) ----
    for k in SPECIES:
        x, y = _xy(rows, "alpha", f"{k}_recall_bulk")
        if x:
            axR.plot(x, [v * 100 for v in y], "-o", color=SP_COLOR[k], label=f"{k} (bulk)")
        xt, yt = _xy(rows, "alpha", f"{k}_recall")
        if xt:
            axR.plot(xt, [v * 100 for v in yt], "--", color=SP_COLOR[k], alpha=0.45, lw=1)
    axR.set_xlabel("convergence semi-angle α (mrad)")
    axR.set_ylabel("depth recall (%)")
    axR.set_ylim(-3, 103)
    axR.set_title("atom recovery vs α  (solid = bulk, dashed = all)")
    axR.grid(alpha=0.3)
    axR.legend(loc="lower right", fontsize=9)

    # --- z-RMS vs alpha, with the diffraction depth-resolution reference δz = λ/α² ----
    for k in SPECIES:
        x, y = _xy(rows, "alpha", f"{k}_zrms")
        if x:
            axZ.plot(x, y, "-o", color=SP_COLOR[k], label=f"{k} z-RMS")
    x, y = _xy(rows, "alpha", "z_rms")
    if x:
        axZ.plot(x, y, "-s", color="0.25", lw=2, label="all z-RMS")
    if alphas:
        aa = sorted(alphas)
        dz_res = [LAMBDA_A / (a * 1e-3) ** 2 for a in aa]
        axZ.plot(aa, dz_res, ":", color="0.5", label="δz = λ/α² (depth res.)")
    axZ.set_xlabel("convergence semi-angle α (mrad)")
    axZ.set_ylabel("depth error / resolution (Å)")
    axZ.set_title("depth localisation vs α")
    axZ.grid(alpha=0.3)
    axZ.legend(loc="upper right", fontsize=9)

    fig.suptitle("Multislice ptychography depth recovery vs aperture (Cs-corrected@30 opened up)",
                 fontsize=13)
    fig.savefig(path, dpi=130)
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="*", help="run dirs (each with report.json) or report.json paths")
    ap.add_argument("--glob", default=None, help="glob for run dirs, e.g. 'results/2026-W37/atomfind_a*'")
    ap.add_argument("--fig", default=None, help="figure path (default figs/<week>/atomfind_depth_vs_alpha.png)")
    ap.add_argument("--csv", default=None, help="csv path (default results/<week>/atomfind_depth_summary.csv)")
    a = ap.parse_args()

    paths = list(a.runs) + (sorted(glob.glob(a.glob)) if a.glob else [])
    if not paths:
        raise SystemExit("no runs given (positional dirs or --glob)")
    rows = sorted((load_run(p) for p in paths), key=lambda r: r["alpha"])
    print("alphas:", [r["alpha"] for r in rows])
    for r in rows:
        print(f"  a{r['alpha']:<3} dz={r['dz']}  Pb={r['Pb_recall_bulk']} Ti={r['Ti_recall_bulk']} "
              f"O={r['O_recall_bulk']}  z-RMS={r['z_rms']}")

    csv_path = a.csv or os.path.join(week_dir("results"), "atomfind_depth_summary.csv")
    write_csv(rows, csv_path)
    print(f"wrote {csv_path}")
    make_fig(rows, a.fig or os.path.join(week_dir("figs"), "atomfind_depth_vs_alpha.png"))


if __name__ == "__main__":
    main()
