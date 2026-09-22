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


#: The instrument the campaign assumes -- sim/simulate_4dstem.py's ABERRATIONS ("a real JEOL ARM, hexapole
#: corrector, corrected to 3rd order, measured to 5th") -- plus order-of-magnitude post-tuning residuals
#: from published hexapole-corrector tableaus for the terms the sim leaves out. "tuned" = the operator
#: nulls it at the design aperture and it drifts; "hardware" = it is what the corrector leaves behind.
ASSUMED_TABLEAU = {
    "C12": (5.0,  "two-fold astigmatism A1",  "tuned",    "sim: 0.5 nm residual"),
    "C21": (3e2,  "coma B2",                  "tuned",    "typical ~30 nm"),
    "C23": (3e2,  "three-fold astigmatism A2","tuned",    "typical ~30 nm"),
    "C32": (5e3,  "star S3",                  "tuned",    "typical ~0.5 um"),
    "C34": (1e4,  "four-fold astigmatism A3", "tuned",    "typical ~1 um"),
    "C41": (2e5,  "fourth-order coma B4",     "hardware", "typical ~20 um"),
    "C43": (1e6,  "three-lobe D4",            "hardware", "typical ~100 um; the known 4th-order hexapole residual"),
    "C45": (5e5,  "five-fold astigmatism A4", "hardware", "typical ~50 um"),
    "C52": (3e6,  "fifth-order star S5",      "hardware", "typical ~0.3 mm"),
    "C54": (3e6,  "rosette R5",               "hardware", "typical ~0.3 mm"),
    "C56": (1e7,  "six-fold astigmatism A5",  "hardware", "sim: 1 mm, the hexapole signature"),
}


def growth_figure(out, alphas=(30, 40, 50, 60, 70, 80, 90, 100), limit=0.1):
    """Every non-round term of the assumed instrument, in waves at the edge, against aperture."""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    a = np.linspace(min(alphas), max(alphas), 200)
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    cols = {"tuned": "#898781", "hardware": "#0072b2"}
    ends = []
    for t, (v, name, kind, _) in ASSUMED_TABLEAU.items():
        w = np.array([waves(t, v, x) for x in a])
        lw = 2.4 if t == "C56" else 1.4
        col = "#c2490a" if t == "C56" else cols[kind]
        ax.plot(a, w, color=col, lw=lw, alpha=0.95 if kind == "hardware" else 0.7)
        short = {"C12": "two-fold", "C21": "coma", "C23": "three-fold", "C32": "star", "C34": "four-fold",
                 "C41": "4th coma", "C43": "three-lobe", "C45": "five-fold", "C52": "5th star",
                 "C54": "rosette", "C56": "six-fold"}[t]
        ends.append([np.log10(w[-1]), f"{t} {short}", col])
    # labels at the right edge, pushed apart so none overprints its neighbour (log-spaced)
    ends.sort(key=lambda e: e[0]); gap = 0.115
    for i in range(1, len(ends)):
        if ends[i][0] - ends[i - 1][0] < gap:
            ends[i][0] = ends[i - 1][0] + gap
    for y, lab, col in ends:
        ax.annotate(lab, (a[-1], 10 ** y), xytext=(5, 0), textcoords="offset points",
                    fontsize=8, color=col, va="center")
    ax.axhline(limit, color="#00795c", lw=1.8)
    ax.annotate(f"{limit:g} waves: the six-fold limit measured at 70 mrad", (min(alphas) + 1, limit),
                xytext=(0, 5), textcoords="offset points", fontsize=9, color="#00795c")
    ax.axvline(30, color="#898781", lw=1, ls=(0, (3, 3)))
    ax.annotate("design aperture", (30.6, 0.012), fontsize=8.5, color="#52514e", ha="left", va="bottom")
    ax.set_yscale("log"); ax.set_ylim(0.01, 500); ax.set_xlim(min(alphas), max(alphas) + 18)
    ax.set_yticks([0.01, 0.1, 1, 10, 100]); ax.set_yticklabels(["0.01", "0.1", "1", "10", "100"])
    ax.set_xlabel("aperture semi-angle (mrad)")
    ax.set_ylabel("waves at the aperture edge")
    ax.set_title("How every non-round residual grows as the corrector is opened past its design aperture",
                 loc="left", fontsize=11, pad=10)
    ax.text(0.01, -0.16, "Blue: what the hardware leaves behind. Grey: what the operator tunes out at the design "
            "aperture, and which drifts. Red: six-fold, the term measured in this campaign. Values are the "
            "campaign's own assumed instrument plus published order-of-magnitude residuals, not a measured tableau.",
            transform=ax.transAxes, fontsize=8.5, color="#52514e", va="top", wrap=True)
    ax.grid(True, which="both", color="#e1e0d9", lw=0.6)
    fig.tight_layout()
    fig.savefig(out, dpi=190, bbox_inches="tight")
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fig", default=None, help="write the growth-with-aperture figure to this path")
    ap.add_argument("--alpha", nargs="+", type=float, default=[70.0])
    ap.add_argument("--terms", nargs="+", default=["C12", "C23", "C34", "C43", "C45", "C56"])
    ap.add_argument("--waves", type=float, default=None,
                    help="print the C_nm that gives this many waves, per term and aperture")
    ap.add_argument("--tsv", default=None, help="a sweep tsv: report the waves each row actually carries")
    a = ap.parse_args()

    if a.fig:
        growth_figure(a.fig)
        return
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
