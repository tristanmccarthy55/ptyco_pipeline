#!/usr/bin/env python
"""
@file probe_refit_check.py
@brief Stage 2.5 check: did the probe update, started from the corrector's C3/C5 and an offset C1,
       converge to the true probe? Per trial, from the recon h5 and the in-job probe record.

For every recon_<camp>_a<A>_<mode>_df<C1>[_rN][_ps<P>[x<P2>][dfo[o<OS2>][c<cap>]]]_n<NITER>_NL<NL> dir
under the given roots:
  - probe_alive     final probe finite, with its power not collapsed (ratio to the start probe)
  - ov_start/final  overlap |<P, P_true>| of the start and final probe with the TRUE probe, invariant to a
                    real-space shift and a global phase (a shift is a gauge the object can absorb)
  - c1_final        the C1 whose probe (C3/C5 held at the sim's values) best matches the final probe, and that
                    best overlap: how defocus-like the refined probe is. Calibrated 2026-09-17 on known-probe
                    h5s: returns -60.0 A (a70) and -160.0 A (a90) exactly from a reference 7 A away, overlap
                    0.9997 / 0.9994; the engine stores the probe as abTEM builds it (not conjugated).
  - final_error     last value of the error trace; slab_z0/z1 the slab edges (c1_objective.h5_slab_profile)
Fixed-probe controls (no _ps) at the same start show what the update bought.

    python analysis/probe_refit_check.py --root ~/Desktop/c1fit_a70 ~/Desktop/c1fit_a90
"""
from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import json
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
# The defocus-only tag is dfo[o<OSTART2>][c<cap>] (campaign/run_c1_search.sh): "oinf" = object frozen in the
# full engine, "c0.05" = focus step capped at 0.05 A/iteration. Both appeared on 2026-09-22.
NAME_RE = re.compile(r"recon_(?P<camp>[A-Za-z0-9]+)_a(?P<alpha>\d+)_(?P<mode>lab|Pb|Ti)_df(?P<c1>-?\d+(?:\.\d+)?)"
                     r"(?:_r(?P<rep>\d+))?(?:_ps(?P<ps>\d+)(?:x(?P<ps2>\d+))?"
                     r"(?P<dfo>dfo)?(?:o(?P<os2>[A-Za-z0-9]+))?(?:c(?P<cap>[0-9.]+))?)?_n(?P<niter>\d+)_NL(?P<nl>\d+)$")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def shift_invariant_overlap(a, Fb):
    """max over real-space shifts of |<a, b>| / (||a|| ||b||), with b given by its unnormalised fft2 Fb."""
    x = np.fft.ifft2(np.fft.fft2(a) * np.conj(Fb))
    return float(np.abs(x).max() * np.sqrt(a.size) / (np.linalg.norm(a) * np.linalg.norm(Fb)))


class ProbeModel:
    """P(C1) on the recon grid, C3/C5 fixed: abTEM probe at a reference C1 times exp(i phi(k) dC1) in Fourier space,
    phi from a 1 A finite difference (0.8 rad at the 70 mrad edge, no wrapping)."""

    def __init__(self, info, n_b):
        sim = _load("simulate_4dstem", os.path.join(REPO, "sim", "simulate_4dstem.py"))
        sim.CONVERGENCE_MRAD = float(info["alpha_mrad"]); sim.BIN_FACTOR = int(info["bin_factor"])
        sim.ABERRATIONS = {k: float(v) for k, v in info["aberrations_Cnm_A_phi_rad"].items()}
        self.c1_ref = float(info["c1_sim_A"])
        build = lambda c1: (setattr(sim, "DEFOCUS_A", c1),
                            sim.build_initial_probe(n_b, float(info["box_A"]), aberrated=True))[1].astype(np.complex128)
        self.F0 = np.fft.fft2(build(self.c1_ref)); F1 = np.fft.fft2(build(self.c1_ref + 1.0))
        inside = np.abs(self.F0) > 0.5 * np.abs(self.F0).max()
        self.phi = np.where(inside, np.angle(F1 * np.conj(self.F0)), 0.0)

    def F(self, c1):
        return self.F0 * np.exp(1j * self.phi * (c1 - self.c1_ref))

    def fit_c1(self, P, half_range=80.0, step=0.5):
        d = np.arange(-half_range, half_range + step / 2, step)
        ov = np.array([shift_invariant_overlap(P, self.F(self.c1_ref + x)) for x in d])
        i = int(np.argmax(ov)); best = d[i]
        if 0 < i < len(d) - 1:                           # parabolic refinement
            y0, y1, y2 = ov[i - 1], ov[i], ov[i + 1]; den = y0 - 2 * y1 + y2
            if den < 0:
                best = d[i] + 0.5 * step * (y0 - y2) / den
        return self.c1_ref + float(best), float(ov[i]), bool(i in (0, len(d) - 1))


def variant_label(m):
    """How the probe was treated, from the directory name. The 2026-09-22 batch added two
    variants of the defocus-only update: the object frozen in the full engine (OBJECT_START2=inf,
    tag "oinf") and a cap on the focus step (PROBE_DFO_MAX_STEP, tag "c<cap>")."""
    if not m["ps"]:
        return "fixed"
    kind = "C1-only" if m["dfo"] else "full-pixel"
    where = f"{m['ps']}" + (f"+{m['ps2']}" if m["ps2"] else " (presolve only)")
    extra = []
    if m["os2"]:
        extra.append("object frozen" if m["os2"].lower().startswith("inf") else f"object from {m['os2']}")
    if m["cap"]:
        extra.append(f"step <= {m['cap']} A")
    return f"{kind} {where}" + (", " + ", ".join(extra) if extra else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", nargs="+", required=True)
    ap.add_argument("--camp", nargs="+", default=["c1fit", "c1"])
    ap.add_argument("--out-csv", default=None)
    a = ap.parse_args()
    obj = _load("c1_objective", os.path.join(HERE, "c1_objective.py"))

    rows, models = [], {}
    for root in a.root:
        for camp in a.camp:
            for d in sorted(glob.glob(os.path.join(os.path.expanduser(root), "**", f"recon_{camp}_a*"), recursive=True)):
                m = NAME_RE.match(os.path.basename(d.rstrip("/")))
                if not m or not os.path.isdir(d):
                    continue
                h5 = sorted(glob.glob(os.path.join(d, "analysis", "**", "*_recons.h5"), recursive=True))
                pj = os.path.join(d, "01", "probe_initial.json")
                variant = variant_label(m)
                row = dict(alpha=int(m["alpha"]), mode=m["mode"], variant=variant, c1_start=float(m["c1"]),
                           rep=int(m["rep"] or 1), niter=int(m["niter"]), dir=d)
                log = sorted(glob.glob(os.path.join(d, "slurm_*.out")))
                if log:
                    txt = open(log[-1], errors="replace").read()
                    row["nan_in_log"] = "contains NaNs" in txt
                if not h5 or not os.path.isfile(pj):
                    row["status"] = "NO H5 (crashed?)" if not h5 else "no probe_initial.json"
                    rows.append(row); continue
                info = json.load(open(pj)); truth = float(info["c1_sim_A"])
                import h5py
                with h5py.File(h5[-1], "r") as f:
                    Pf = np.asarray(f["reconstruction/probes"]).squeeze().astype(np.complex128)
                    P0 = np.asarray(f["reconstruction/p/probe_initial"]).squeeze().astype(np.complex128)
                if Pf.ndim == 3:                         # several probe modes: the first is the dominant one
                    Pf = Pf[..., 0] if Pf.shape[-1] < Pf.shape[0] else Pf[0]
                key = (row["alpha"], Pf.shape[0], truth)
                if key not in models:
                    models[key] = ProbeModel(info, Pf.shape[0])
                M = models[key]
                alive = bool(np.isfinite(Pf).all() and np.linalg.norm(Pf) > 1e-6 * np.linalg.norm(P0))
                row.update(c1_true=truth, dc1_start=row["c1_start"] - truth, status="ok" if alive else "PROBE DEAD",
                           power_ratio=float((np.abs(Pf) ** 2).sum() / (np.abs(P0) ** 2).sum()) if np.isfinite(Pf).all() else np.nan)
                trace = sorted(glob.glob(os.path.join(d, "analysis", "**", "*_error_trace.csv"), recursive=True))
                if trace:
                    tr = obj.read_csv_rows(trace[-1]); row["final_error"] = float(tr[-1]["fourier_error"])
                    for line in open(trace[-1]):                     # the engine's own cumulative defocus shift
                        if line.startswith("# probe_defocus_shift_A="):
                            row["engine_dc1"] = float(line.split("=")[1])
                if alive:
                    row["ov_start"] = shift_invariant_overlap(P0, M.F(truth))
                    row["ov_final"] = shift_invariant_overlap(Pf, M.F(truth))
                    c1f, ovb, edge = M.fit_c1(Pf)
                    row.update(c1_final=c1f, dc1_final=c1f - truth, ov_best_c1=ovb, c1_fit_at_edge=edge)
                    row["slab_centroid"], row["slab_z0"], row["slab_z1"] = obj.h5_slab_profile(h5[-1])
                rows.append(row)

    if not rows:
        raise SystemExit(f"no recon_{{{','.join(a.camp)}}}_* dirs under {a.root}")
    rows.sort(key=lambda r: (r["alpha"], r["variant"], r["c1_start"], r["rep"]))
    f = lambda v, fmt: "-" if v is None or (isinstance(v, float) and not np.isfinite(v)) else format(v, fmt)
    print(f"{'alpha':>5} {'variant':<34} {'dC1 start':>9} {'status':<10} {'error':>8} {'ov start':>8} {'ov final':>8} "
          f"{'dC1 final':>9} {'ov@C1':>6} {'slab (A)':>12}")
    for r in rows:
        line = (f"{r['alpha']:>5} {r['variant']:<26} {f(r.get('dc1_start'), '+.0f'):>9} {r['status']:<10} "
                f"{f(r.get('final_error'), '.3f'):>8} {f(r.get('ov_start'), '.3f'):>8} {f(r.get('ov_final'), '.3f'):>8} "
                f"{f(r.get('dc1_final'), '+.1f'):>9} {f(r.get('ov_best_c1'), '.3f'):>6} "
                f"{f(r.get('slab_z0'), '.1f'):>5}-{f(r.get('slab_z1'), '.1f'):<5}")
        if r.get("engine_dc1") is not None:
            line += f"  engine says dC1 {r['engine_dc1']:+.1f}"
        if r.get("c1_fit_at_edge"):
            line += "  C1 FIT AT GRID EDGE"
        if r.get("nan_in_log"):
            line += "  NaN in log"
        print(line)
    out = a.out_csv or os.path.join(obj.week_dir("results"), "probe_refit.csv")
    cols = ["alpha", "mode", "variant", "niter", "rep", "c1_true", "c1_start", "dc1_start", "status", "nan_in_log",
            "power_ratio", "final_error", "ov_start", "ov_final", "c1_final", "dc1_final", "engine_dc1", "ov_best_c1", "c1_fit_at_edge",
            "slab_centroid", "slab_z0", "slab_z1", "dir"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
