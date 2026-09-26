#!/usr/bin/env python
"""
@file analyse_sweep.py
@brief One command for a sweep's analysis: triage -> matched kernels -> atomfind -> one summary table.

Written for runs too big to download (a BIN-1-class recon h5 is ~1 GB): run it where the reconstructions are --
on Blythe, as a CPU job (campaign/run_analysis.sh) -- and bring home only the small outputs. It drives the
pipeline's own CLIs (extract_psf.py, run_atomfind.py) and triage/ladder functions, so nothing is re-implemented.

Every geometric constant is READ, never typed, because a region leg's box, pixel and scan centre differ from the
old 70 A box's (the Rules: "geometry-coupled constants derived, not typed"):
  NL          the reconstructed object's layer count (h5)
  dx          the object pixel, reconstruction/p/dx_spec (h5)  -- 0.0492 A in the 70 A box, 0.074 A at 80 mrad
  dz          sim_meta beam_thickness_A / NL;  zdrop = round(z_vacuum / dz)  (PSF_KERNELS.md)
  scan centre (region_side_A / 2, same) for a --region-side leg (sim_meta), else the preset's (40, 20)
  alpha       sim_meta convergence_mrad

    python analysis/analyse_sweep.py --root <dir holding recon_af_*> --gt <dir holding gt_prepared.npz> \\
        --labels round_a040 ceosopt_a080 --out <out dir>

Writes <out>/psf/ (kernels + extractor logs), <out>/atomfind_<label>/ (report.json, figures), <out>/logs/,
<out>/summary.csv (one row per label) and prints the table. A label whose recon fails triage (saturated phase, no
h5) is reported and skipped, never analysed. The region GT is built with
`python -m atomfind.make_gt_cache --thin-cells 5 --z-vacuum 4 --region-side 210 --out <gt>/gt_prepared.npz`.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import triage_recon  # noqa: E402  (object_health, WRAP_FRAC_MAX, STD_MAX)

MODES = ("lab", "Pb", "Ti")
COLS = ["label", "alpha", "NL", "dx_A", "dz_A", "zdrop", "scan_centre", "wrapped", "phase_std", "status",
        "precision", "Pb_recall_bulk", "Ti_recall_bulk", "O_recall_bulk", "xy_rms", "z_rms", "confusion",
        "n_found", "fp_total"]


def newest(paths):
    paths = [p for p in paths if os.path.exists(p)]
    return max(paths, key=os.path.getmtime) if paths else None


def leg_dirs(root, label):
    """{mode: recon dir}: recon_af_<label>_<mode>_NL<n>, the newest when several NL dirs exist."""
    out = {}
    for m in MODES:
        d = newest(glob.glob(os.path.join(root, f"recon_af_{label}_{m}_NL*")))
        if d:
            out[m] = d
    return out


def recon_h5(d):
    return newest(glob.glob(os.path.join(d, "analysis", "**", "*_recons.h5"), recursive=True))


def geometry(recon_dir, h5, z_vacuum):
    """Everything the kernels and the finder need, from the leg's own h5 and sim_meta.mat."""
    import h5py
    from scipy.io import loadmat
    with h5py.File(h5, "r") as f:
        nl = int(f["reconstruction/object"].shape[0])
        dx = float(np.ravel(f["reconstruction/p/dx_spec"][...])[0]) * 1e10
    # on Blythe the recon's 01/ links to its sim; a packed tarball carries sim_out_af_<leg>/01/ instead
    leg = os.path.basename(os.path.normpath(recon_dir)).replace("recon_af_", "", 1).rsplit("_NL", 1)[0]
    meta_path = newest([os.path.join(recon_dir, "01", "sim_meta.mat"),
                        os.path.join(os.path.dirname(os.path.normpath(recon_dir)), f"sim_out_af_{leg}", "01", "sim_meta.mat")])
    if meta_path is None:
        raise FileNotFoundError(f"no sim_meta.mat for {recon_dir}")
    meta = loadmat(meta_path, simplify_cells=True)["meta"]
    dz = float(meta["beam_thickness_A"]) / nl
    region = float(meta.get("region_side_A", 0.0) or 0.0)
    return dict(NL=nl, dx_A=dx, dz_A=dz, zdrop=int(round(z_vacuum / dz)), alpha=float(meta["convergence_mrad"]),
                scan_centre=(region / 2, region / 2) if region else None)


def run(cmd, log):
    with open(log, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    return r.returncode


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, help="directory holding the recon_af_<label>_<mode>_NL<n> dirs")
    ap.add_argument("--labels", nargs="+", required=True, help="sweep labels (tsv column 1)")
    ap.add_argument("--gt", required=True, help="directory holding the matching gt_prepared.npz")
    ap.add_argument("--out", required=True)
    ap.add_argument("--z-vacuum", type=float, default=4.0, help="the sims' Z_VACUUM (campaign: 4 A)")
    ap.add_argument("--python", default=sys.executable, help="interpreter for the pipeline CLIs")
    a = ap.parse_args()
    for sub in ("psf", "logs"):
        os.makedirs(os.path.join(a.out, sub), exist_ok=True)
    extract = os.path.join(HERE, "atomfind", "extract_psf.py")
    finder = os.path.join(HERE, "atomfind", "run_atomfind.py")
    import relaxation_ladder

    rows = []
    for label in a.labels:
        row = dict(label=label, status="")
        dirs = leg_dirs(a.root, label)
        h5s = {m: recon_h5(d) for m, d in dirs.items()}
        missing = [m for m in MODES if not h5s.get(m)]
        if missing:
            row["status"] = "NO H5: " + " ".join(missing); rows.append(row); print(f"{label}: {row['status']}"); continue
        health = {m: triage_recon.object_health(h5s[m]) for m in MODES}
        bad = [m for m in MODES if health[m][0] > triage_recon.WRAP_FRAC_MAX or health[m][1] > triage_recon.STD_MAX]
        row["wrapped"], row["phase_std"] = f"{health['lab'][0]:.2e}", f"{health['lab'][1]:.3f}"
        if bad:
            row["status"] = "SATURATED: " + " ".join(bad); rows.append(row); print(f"{label}: {row['status']}"); continue
        g = geometry(dirs["lab"], h5s["lab"], a.z_vacuum)
        row.update(alpha=g["alpha"], NL=g["NL"], dx_A=f"{g['dx_A']:.5f}", dz_A=f"{g['dz_A']:.4f}", zdrop=g["zdrop"],
                   scan_centre="preset" if g["scan_centre"] is None else f"{g['scan_centre'][0]:g},{g['scan_centre'][1]:g}")
        psf = {}
        for el in ("Pb", "Ti"):
            name = f"{el}_{label}"
            rc = run([a.python, extract, dirs[el], name, "--zdrop", str(g["zdrop"]), "--dx", f"{g['dx_A']:.6f}",
                      "--out", os.path.join(a.out, "psf")], os.path.join(a.out, "logs", f"extract_{name}.log"))
            psf[el] = os.path.join(a.out, "psf", f"psf_{name}_vol.npy")
            if rc or not os.path.exists(psf[el]):
                row["status"] = (row["status"] + " " if row["status"] else "KERNEL FAILED:") + f" {el} (logs/extract_{name}.log)"
        if row["status"]:
            rows.append(row); print(f"{label}: {row['status']}"); continue
        outd = os.path.join(a.out, f"atomfind_{label}")
        cmd = [a.python, finder, "--preset", "thin", "--recon", h5s["lab"], "--dz", f"{g['dz_A']:.6f}",
               "--data-dir", a.gt, "--single-atom-vol", psf["Pb"], "--ti-kernel-vol", psf["Ti"], "--out", outd,
               "--set", f"dx={g['dx_A']:.6f}"]
        if g["scan_centre"] is not None:
            cmd += ["--set", f"scan_center_xy={g['scan_centre'][0]:g},{g['scan_centre'][1]:g}"]
        rc = run(cmd, os.path.join(a.out, "logs", f"atomfind_{label}.log"))
        if rc or not os.path.exists(os.path.join(outd, "report.json")):
            row["status"] = f"ATOMFIND FAILED (logs/atomfind_{label}.log)"; rows.append(row); print(f"{label}: {row['status']}"); continue
        r = relaxation_ladder.ladder_row(outd, "-", label, "", None, int(round(g["alpha"])))
        row.update(status="ok", **{k: r.get(k) for k in COLS if k in r and k not in ("label", "alpha")})
        rows.append(row)
        pct = lambda k: "-" if r.get(k) is None else f"{100 * float(r[k]):.0f}"
        print(f"{label}: ok | Pb/Ti/O bulk {pct('Pb_recall_bulk')}/{pct('Ti_recall_bulk')}/{pct('O_recall_bulk')} % "
              f"| z-RMS {float(r['z_rms']):.3f} A | precision {float(r['precision']):.3f}", flush=True)

    path = os.path.join(a.out, "summary.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    print(f"wrote {path} ({sum(r.get('status') == 'ok' for r in rows)} of {len(rows)} labels ok)")
    return 0 if all(r.get("status") == "ok" for r in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
