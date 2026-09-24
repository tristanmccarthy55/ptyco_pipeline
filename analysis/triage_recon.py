#!/usr/bin/env python
"""
@file triage_recon.py
@brief Triage a batch of reconstructions BEFORE any of them is analysed.

The old triage was "is there an h5, and does the log say NaN?". That is not enough. On 2026-09-22
the a70 Ti leg of the 70 A batch passed both checks -- COMPLETED, an h5 on disk, no NaN in the log --
and its object is saturated junk: every layer, vacuum included, sits at a phase standard deviation of
1.5 rad with 3 % of pixels wrapped against pi. Feeding that to extract_psf produces "0 grid atoms
found"; feeding it to the finder would have produced numbers.

So this also looks at the object. Measured across known-good, merely-degraded and saturated legs, the
separation is total rather than marginal:

    leg                                   |phase|>3 rad     phase std
    thin a70 / a90 lab (published)            0.0000      0.12 / 0.03
    six-fold 0.45 waves (mis-oriented probe,
      badly degraded but not saturated)       0.0000           0.13
    70 A a70 Ti (saturated)                   0.0314           1.53

A real reconstruction wraps almost nothing: the round a70 control has 41 pixels over 3 rad (5e-6), the
fixed six-fold legs 2e-5 and 5e-5, all on the edge of the illuminated field, where it is poorly
constrained. Saturation is two to three orders above that, so WRAP_FRAC_MAX = 1e-4 separates them; the
phase standard deviation is quoted beside it as the size of the failure.

    ~/hyperspy-bundle/bin/python analysis/triage_recon.py --root ~/Desktop/thin18_0922 --box 77.9347
    ~/hyperspy-bundle/bin/python analysis/triage_recon.py --root ~/Desktop/relax_0922
"""
from __future__ import annotations

import argparse
import glob
import os
import re

import numpy as np

WRAP_RAD = 3.0        # |phase| above this is as good as wrapped against pi
WRAP_FRAC_MAX = 1e-4  # genuine legs measured at <= 5e-5 (edge pixels), saturated at 3e-2
STD_MAX = 0.5         # rad: an order of magnitude above the worst genuine leg measured


def object_health(h5_path):
    """(wrapped fraction, deramped phase std, layer count) of a reconstruction's object."""
    import h5py
    with h5py.File(h5_path, "r") as f:
        obj = np.asarray(f["reconstruction/object"]).squeeze()
    ph = np.angle(obj)
    dr = ph - ph.mean(axis=tuple(range(1, ph.ndim)), keepdims=True)
    return float((np.abs(ph) > WRAP_RAD).mean()), float(dr.std()), int(ph.shape[0])


def triage(d):
    row = {"dir": os.path.basename(d.rstrip("/"))}
    h5 = sorted(glob.glob(os.path.join(d, "analysis", "**", "*_recons.h5"), recursive=True))
    logs = sorted(glob.glob(os.path.join(d, "slurm_*.out")))
    row["h5"] = len(h5)
    if logs:
        txt = open(logs[-1], errors="replace").read()
        row["nan_in_log"] = "contains NaNs" in txt
        m = re.search(r"Number of probe positions:\s*(\d+)", txt)
        row["positions"] = int(m.group(1)) if m else None
    if not h5:
        row["verdict"] = "NO H5 (MATLAB crashed)"
        return row
    row["wrapped"], row["phase_std"], row["layers"] = object_health(h5[-1])
    if row["wrapped"] > WRAP_FRAC_MAX or row["phase_std"] > STD_MAX:
        row["verdict"] = "SATURATED (unusable)"
    elif row.get("nan_in_log"):
        row["verdict"] = "NaN in log"
    else:
        row["verdict"] = "ok"
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", nargs="+", required=True)
    ap.add_argument("--glob", default="recon_*", help="recon dir pattern (default recon_*)")
    a = ap.parse_args()

    dirs = []
    for r in a.root:
        dirs += sorted(glob.glob(os.path.join(os.path.expanduser(r), "**", a.glob), recursive=True))
    dirs = [d for d in dirs if os.path.isdir(d) and glob.glob(os.path.join(d, "01")) or os.path.isdir(d)]
    if not dirs:
        raise SystemExit(f"no {a.glob} directories under {a.root}")

    print(f"{'recon':<52} {'h5':>3} {'pos':>5} {'layers':>6} {'wrapped':>8} {'std':>7}  verdict")
    bad = 0
    for d in dirs:
        r = triage(d)
        if r["verdict"] != "ok":
            bad += 1
        f = lambda k, fmt: "-" if r.get(k) is None else format(r[k], fmt)
        print(f"{r['dir']:<52} {r['h5']:>3} {f('positions','d'):>5} {f('layers','d'):>6} "
              f"{f('wrapped','.4f'):>8} {f('phase_std','.3f'):>7}  {r['verdict']}")
    print(f"\n{len(dirs) - bad} of {len(dirs)} usable")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
