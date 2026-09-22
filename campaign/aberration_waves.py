#!/usr/bin/env python
"""
@file aberration_waves.py
@brief Convert between an aberration coefficient and "waves at the aperture edge", for any order.

The campaign quotes non-round aberrations in waves at the edge of the aperture, because that is the
quantity a corrector's tableau is judged on and the only way terms of different order can be put on
one axis. With the Krivanek convention the wave aberration is

    chi(theta, phi) = (2 pi / lambda) * sum_nm  C_nm * theta^(n+1) / (n+1) * cos(m (phi - phi_nm))

so a term of order n reaches W = C_nm * theta^(n+1) / (n+1) angstroms of path difference at the edge,
and one WAVE there costs

    C_nm(1 wave) = (n + 1) * lambda / theta^(n+1).

That falls fast with order and with aperture: at 70 mrad one wave is 8 A of two-fold astigmatism but
1.0e6 A of six-fold, and at 90 mrad the six-fold figure is 2.2e5 A. It is also why the tolerance
should NOT be assumed to transfer between orders at equal waves: for the same edge phase, a term of
order n displaces the edge ray by (n+1) lambda / theta, which is 0.56 A at n = 1 and 1.69 A at n = 5.

    python campaign/aberration_waves.py --alpha 70 90                  # the scale per order
    python campaign/aberration_waves.py --alpha 70 --waves 0.1 --terms C12 C23 C34 C43 C56
    python campaign/aberration_waves.py --alpha 70 --tsv campaign/nonround_sweep.tsv
"""
from __future__ import annotations

import argparse
import csv
import json

LAMBDA_A = 0.0196877          # 300 keV electron wavelength [A]

#: abTEM 1.0.5 accepts all of these, with a phi<nm> angle for every m > 0
ORDER = {"C10": 1, "C12": 1, "C21": 2, "C23": 2, "C30": 3, "C32": 3, "C34": 3,
         "C41": 4, "C43": 4, "C45": 4, "C50": 5, "C52": 5, "C54": 5, "C56": 5}
NAME = {"C12": "two-fold astigmatism", "C21": "coma", "C23": "three-fold astigmatism",
        "C30": "spherical (Cs)", "C32": "star", "C34": "four-fold astigmatism",
        "C41": "fourth-order coma", "C43": "three-lobe", "C45": "five-fold astigmatism",
        "C50": "fifth-order spherical", "C52": "fifth-order star", "C54": "rosette",
        "C56": "six-fold astigmatism", "C10": "defocus"}


def per_wave(term: str, alpha_mrad: float) -> float:
    """C_nm [A] that puts exactly one wave at the aperture edge."""
    n = ORDER[term]
    return (n + 1) * LAMBDA_A / (alpha_mrad * 1e-3) ** (n + 1)


def waves(term: str, cnm: float, alpha_mrad: float) -> float:
    return cnm / per_wave(term, alpha_mrad)


def edge_ray_A(term: str, alpha_mrad: float) -> float:
    """Edge-ray displacement [A] for ONE wave of this term: dW/dtheta = (n+1) lambda / theta."""
    return (ORDER[term] + 1) * LAMBDA_A / (alpha_mrad * 1e-3)


def tableau_waves(ab: dict, alpha_mrad: float) -> dict:
    """{term: waves at the edge} for every non-round term of an abTEM aberration dict."""
    return {k: waves(k, float(v), alpha_mrad)
            for k, v in ab.items() if k in ORDER and not k.startswith("phi") and k[-1] != "0"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alpha", nargs="+", type=float, default=[70.0])
    ap.add_argument("--terms", nargs="+", default=["C12", "C23", "C34", "C43", "C45", "C56"])
    ap.add_argument("--waves", type=float, default=None,
                    help="print the C_nm that gives this many waves, per term and aperture")
    ap.add_argument("--tsv", default=None, help="a sweep tsv: report the waves each row actually carries")
    a = ap.parse_args()

    if a.tsv:
        with open(a.tsv) as f:
            for r in csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"):
                if r["aber_json"] in ("-", "", None):
                    continue
                w = tableau_waves(json.loads(r["aber_json"]), float(r["alpha"]))
                terms = ", ".join(f"{k} {v:.3f} waves" for k, v in sorted(w.items()))
                print(f"{r['label']:<22} alpha {r['alpha']:>3} mrad   {terms or '(round only)'}")
        return

    for alpha in a.alpha:
        print(f"\nalpha = {alpha:g} mrad")
        print(f"  {'term':<6} {'aberration':<24} {'order':>5} {'C_nm for 1 wave':>16} {'edge ray / wave':>16}"
              + (f" {'C_nm at ' + format(a.waves, 'g') + ' waves':>22}" if a.waves else ""))
        for t in a.terms:
            line = (f"  {t:<6} {NAME.get(t, ''):<24} {ORDER[t]:>5} {per_wave(t, alpha):>16.4g} "
                    f"{edge_ray_A(t, alpha):>15.2f} A")
            if a.waves:
                line += f" {a.waves * per_wave(t, alpha):>22.4g}"
            print(line)


if __name__ == "__main__":
    raise SystemExit(main())
