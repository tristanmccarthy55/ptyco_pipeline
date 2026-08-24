#!/usr/bin/env python
"""@file noise_sweep.py
@brief Noise-robustness harness: inject Gaussian noise into the phase volume, re-run the finder.

Regenerates the noise table (RESULTS.md sec.3 / the report's Table D1). Everything except the
noise is held fixed -- same volume, same kernels, same alignment, same config -- so the only
variable is sigma. Noise is added to the mean-subtracted PHASE volume, which is where the
reconstruction's own noise lives, and the draw is seeded so the table is reproducible.

For scale: the volume's intrinsic vacuum noise floor is sigma ~ 0.023 and an oxygen peak is
~0.06, so the sweep brackets the point where oxygen falls below the noise. The collapse past
sigma ~ 0.02 is that physical limit, not a thresholding artefact -- see RESULTS.md sec.3.

    python atomfind/noise_sweep.py                          # the standard ladder
    python atomfind/noise_sweep.py --sigmas 0 0.02 --seed 7
    python atomfind/noise_sweep.py --json noise.json
"""
from __future__ import annotations
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from atomfind import config, align, psf as psfmod, find, validate, uncertainty

SIGMAS = (0.0, 0.01, 0.02, 0.04, 0.08)


def run_one(V0, dx, pos, Z, olab, cfg, sigma, seed=0, verbose=True):
    """One rung of the ladder: V0 + N(0, sigma), then the full blind finder."""
    V = V0 if sigma == 0 else V0 + np.random.default_rng(seed).normal(0.0, sigma, V0.shape)
    al = align.register(V, dx, pos, Z, cfg)
    kernels = psfmod.species_kernels(cfg, dx)
    found, _ = find.find_atoms_v3(V, cfg, dx, kernels)
    al = align.refine_with_atoms(al, found, pos, Z, cfg)
    rep, m = validate.finder_report(found, pos, Z, al, cfg, olabel=olab)
    cov = {}
    try:                                   # coverage is the canary the results section names
        qtab = uncertainty.calibrate(found, m, cfg, alphas=cfg.uq_alphas,
                                     min_n=cfg.uq_min_stratum)
        _rows, cov = uncertainty.coverage_table(found, m, cfg, qtab, 0.32)
    except Exception as e:
        if verbose:
            print(f"    [uq] {type(e).__name__}: {e}")
    if verbose:
        print(f"  sigma={sigma:<5g} found {rep['n_found']:>5}  prec {rep['precision']:.2f}  "
              f"Pb {rep['Pb']['recall']:.0%}  Ti {rep['Ti']['recall']:.0%}  "
              f"O {rep['O']['recall']:.0%}  conf {validate.confusion_rate(rep):.1%}")
        for w in validate.health_warnings(rep):
            print(f"    ! {w}")
    return rep, cov


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", default="NL70_coherent")
    ap.add_argument("--sigmas", nargs="*", type=float, default=list(SIGMAS))
    ap.add_argument("--seed", type=int, default=0, help="noise draw seed (table is seeded)")
    ap.add_argument("--json", default=None)
    ap.add_argument("--data-dir", default=None, help="see run_atomfind.py --data-dir")
    a = ap.parse_args()
    if a.data_dir:
        config.set_data_dir(a.data_dir)

    cfg = config.preset(a.preset).resolve()
    print(f"[noise_sweep] {cfg.name}  {cfg.recon_vol}  seed={a.seed}")
    V0, dx = align.load_phase(cfg)
    pos, Z = align.load_gt(cfg)
    olab = validate.classify_oxygen(pos, Z, cfg)
    floor = float(np.std(V0[V0 < np.percentile(V0, 50)]))
    print(f"[noise_sweep] intrinsic sub-median phase spread ~ {floor:.3f} "
          f"(the injected sigmas are relative to this)")

    rows, blob = [], {}
    for s in a.sigmas:
        rep, cov = run_one(V0, dx, pos, Z, olab, cfg, s, seed=a.seed)
        ov = rep.get("O_axial_overlap", {}).get("recall", float("nan"))
        rows.append((s, rep["n_found"], rep["precision"], rep["Pb"]["recall"],
                     rep["Ti"]["recall"], rep["O"]["recall"], ov,
                     validate.confusion_rate(rep)))
        blob[f"{s:g}"] = {"report": rep, "coverage": cov}

    print("\n" + "=" * 74)
    print(f"{'sigma':>7} {'found':>6} {'prec':>6} {'Pb':>6} {'Ti':>6} {'O':>6} "
          f"{'O ovlp':>7} {'conf':>7}")
    for s, n, p, pb, ti, o, ov, cf in rows:
        print(f"{s:>7g} {n:>6} {p:>6.2f} {pb:>5.0%} {ti:>5.0%} {o:>5.0%} {ov:>6.0%} {cf:>6.1%}")
    print("\nPrecision is NOT a safety indicator here: it stays ~0.95 while recall collapses,\n"
          "because the surviving bright Pb are still placed correctly. Read the confusion rate\n"
          "and the sigma-coverage instead.")

    if a.json:
        with open(a.json, "w") as f:
            json.dump(blob, f, indent=2, default=float)
        print(f"\nwrote {a.json}")


if __name__ == "__main__":
    main()
