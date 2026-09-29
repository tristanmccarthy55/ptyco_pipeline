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
#: corrector, corrected to 3rd order, measured to 5th") -- plus order-of-magnitude residuals for the terms the
#: sim leaves out. "tuned" = the operator nulls it at the design aperture and it drifts; "hardware" = what the
#: corrector leaves behind. The split follows CEOS for its CESCOR hexapole probe corrector: adjustable C1, A1,
#: B2, A2, C3, S3, A3, D4, C5; limiting A5 (intrinsic) and A4, B4 (parasitic) -- ceos-gmbh.de, CESCOR page.
#: SOURCED VALUE: A5 = 1.0 mm (CEOS, residual aberrations of hexapole-type correctors). The others are
#: UNSOURCED orders of magnitude from an earlier session; do not quote them as a measured tableau.
#: (2026-09-24: D4 was filed as "hardware" at 100 um -- it is adjustable, and at 100 um it dominated the probe.)
ASSUMED_TABLEAU = {
    "C12": (5.0,  "two-fold astigmatism A1",  "tuned",    "sim: 0.5 nm residual"),
    "C21": (3e2,  "coma B2",                  "tuned",    "typical ~30 nm"),
    "C23": (3e2,  "three-fold astigmatism A2","tuned",    "typical ~30 nm"),
    "C32": (5e3,  "star S3",                  "tuned",    "typical ~0.5 um"),
    "C34": (1e4,  "four-fold astigmatism A3", "tuned",    "typical ~1 um"),
    "C41": (2e5,  "fourth-order coma B4",     "hardware", "typical ~20 um"),
    "C43": (1e6,  "three-lobe D4",            "tuned",    "unsourced ~100 um; ADJUSTABLE on a CESCOR (CEOS)"),
    "C45": (5e5,  "five-fold astigmatism A4", "hardware", "typical ~50 um"),
    "C52": (3e6,  "fifth-order star S5",      "hardware", "typical ~0.3 mm"),
    "C54": (3e6,  "rosette R5",               "hardware", "typical ~0.3 mm"),
    "C56": (1e7,  "six-fold astigmatism A5",  "hardware", "1.0 mm: CEOS, hexapole-type corrector"),
}


#: THE CEOS-APPROX INSTRUMENT (2026-09-24, the user's choice): an older, widespread hexapole-corrected column --
#: CEOS CESCOR-like, the pre-DELTA generation -- tuned at a 30 mrad design aperture and opened past it with the
#: probe known ("tanking" the aberrations). Split per CEOS: operator-tuned A1 B2 A2 S3 A3 D4; hardware A4, B4
#: (parasitic) and A5 (intrinsic). A5 = 1.0 mm is CEOS's figure (and consistent: it reaches the pi/4 tuning limit at
#: ~34 mrad at 300 kV, ~35 at 200 kV -- the range these correctors were used to; DCOR/ASCOR, built to beat it, hold
#: |A5| < 0.2 mm). No other term has a sourced value. Tuned terms sit at CEOS_DESIGN_WAVES at the aperture they were
#: tuned at (A1 0.44 nm, B2/A2 22 nm, S3/A3 1 um, D4 40 um at 30 mrad). Parasitic A4, B4: CEOS names A5 as THE
#: limiting residual of an aligned corrector, so A4/B4 must reach pi/4 no sooner than A5 does (~34 mrad), i.e.
#: <= ~27 um; CEOS_PARASITIC_A = 10 um is the central value (2026-09-24; 40 um, the first choice, was the worst an
#: in-spec column could be -- 0.1 waves at 30 mrad -- and broke that ordering). Hardware terms are FIXED as the
#: aperture opens (each grows as alpha^(n+1)).
#: S5/R5 (C52/C54) are left out: CEOS names A5, A4, B4 as the limiting residuals. Orientations are one seeded
#: uniform draw per term (the angle of an m-fold term is only defined modulo 2pi/m); probe sizes vary < 1 A across
#: draws. The simulation energy stays 300 keV (the campaign's), though the target era's column is 200 keV.
#: ASSUMPTION -- probably OPTIMISTIC (flagged 2026-09-29): with alpha given, ceos_tableau() sets every operator-tuned
#: term (A1 A2 S3 A3 D4) to CEOS_DESIGN_WAVES = 0.1 waves at the edge of THAT aperture, i.e. the operator re-runs the
#: corrector tuning after opening the aperture and reaches the same accuracy as at the 30 mrad design aperture. CEOS
#: lists these terms as adjustable, but how well they can be measured and nulled past the design aperture is NOT known,
#: and holding 0.1 waves means the coefficients must shrink as the aperture opens (D4 9.6 um at 40 mrad -> 0.30 um at 80,
#: ~30x finer). The other bound is alpha=None ("as built": tuned to 0.1 waves at 30 mrad and left, each term then grows
#: as alpha^(n+1): at 80 mrad A1 0.7, A2 1.9, S3/A3 5, D4 13 waves) -- probably too pessimistic, since an operator
#: would retune something. The real instrument lies between. Every ceosopt_* row uses the optimistic retuned case;
#: ceosbuilt_a080 is the pessimistic bound at 80 mrad (2026-09-29).
CEOS_DESIGN_MRAD, CEOS_DESIGN_WAVES, CEOS_A5_A, CEOS_SEED = 30.0, 0.1, 1.0e7, 0
CEOS_PARASITIC_A = 1.0e5        # A4, B4 = 10 um (central; see above)
CEOS_TUNED = ["C12", "C21", "C23", "C32", "C34", "C43"]
CEOS_HARDWARE = ["C41", "C45", "C56"]


def ceos_tableau(c3: float, c5: float = 1.0e7, alpha: float | None = None, b2: float | None = None) -> dict:
    """abTEM aberration dict: the round part (C30, C50) plus the CEOS-approx non-round residuals.

    alpha=None: the tuned terms as left at the 30 mrad design aperture -- the operator opens the aperture and
    retunes only C1/C3 ("as built"). alpha given: the operator retunes every tunable term AT that aperture, to the
    same accuracy (CEOS_DESIGN_WAVES at alpha), and b2 [A, signed] sets coma B2 against fourth-order coma B4, the
    one hardware term a knob shares symmetry with (negative = opposed). A5 and A4 have no such partner: no knob
    touches them, so they set the probe's floor. Angles are drawn in a fixed order, so every variant shares them."""
    import numpy as np
    rng = np.random.default_rng(CEOS_SEED)
    ab = {"C30": float(c3), "C50": float(c5)}
    for t in CEOS_TUNED + CEOS_HARDWARE:
        if t == "C56":
            v = CEOS_A5_A
        elif t in CEOS_HARDWARE:
            v = CEOS_PARASITIC_A
        elif alpha is None:
            v = CEOS_DESIGN_WAVES * per_wave(t, CEOS_DESIGN_MRAD)
        else:
            v = CEOS_DESIGN_WAVES * per_wave(t, alpha)
        ab[t] = round(v, 4)
        ab["phi" + t[1:]] = round(float(rng.uniform(0, 2 * np.pi)) / int(t[2]), 6)
    if b2 is not None:
        ab["C21"] = round(abs(float(b2)), 4)
        ab["phi21"] = round((ab["phi41"] + (np.pi if b2 < 0 else 0.0)) % (2 * np.pi), 6)
    return ab


#: HAIDER (CEOS) NOTATION -> abTEM (2026-09-29). The CEOS software reports each coefficient as a magnitude and the angle
#: of the complex coefficient X in
#:   chi(w) = (2 pi / lambda) Re{ 1/2 C1 w w' + 1/2 A1 w'^2 + B2 w^2 w' + 1/3 A2 w'^3 + 1/4 C3 (w w')^2 + S3 w^3 w'
#:            + 1/4 A3 w'^4 + B4 w^3 w'^2 + D4 w^4 w' + 1/5 A4 w'^5 + 1/6 C5 (w w')^3 + S5 w^4 w'^2 + R5 w^5 w' + 1/6 A5 w'^6 }
#: with w = theta e^{i phi} and w' its conjugate (Uhlemann & Haider 1998; the same form in arXiv 2510.01493 eq. A1 and
#: arXiv 2603.23958 Table 1). Two consequences for abTEM's C_nm theta^(n+1)/(n+1) cos(m (phi - phi_nm)):
#:   * the mixed terms B2, S3, B4, D4, S5, R5 carry NO 1/(n+1): C21 = 3 B2, C32 = 4 S3, C41 = 5 B4, C43 = 5 D4,
#:     C52 = 6 S5, C54 = 6 R5. The A and C terms map one to one.
#:   * the angle is m x the azimuth: an A term (w'^m) points at phi_nm = +angle/m, a mixed term (w^p w'^q, p > q)
#:     at phi_nm = -angle/m. The sign matters only BETWEEN kinds: A1 against S3, A2 against D4 (B2 against B4 share it).
#: That the software's angle is the complex coefficient's argument (not the pattern's azimuth) is read off the data:
#: every term's angles fill -180..180 deg (A3 -159/-176/-169, A5 150), where azimuths would fill only 360/m.
#: A global rotation or mirror of every angle only rotates or mirrors the probe. check: haider_chi() evaluates the
#: formula above directly, so haider_to_abtem() is tested against it, and abTEM's probe against both (tests below).
HAIDER = {"C1": (1, 1, 1 / 2), "A1": (0, 2, 1 / 2), "B2": (2, 1, 1.0), "A2": (0, 3, 1 / 3), "C3": (2, 2, 1 / 4),
          "S3": (3, 1, 1.0), "A3": (0, 4, 1 / 4), "B4": (3, 2, 1.0), "D4": (4, 1, 1.0), "A4": (0, 5, 1 / 5),
          "C5": (3, 3, 1 / 6), "S5": (4, 2, 1.0), "R5": (5, 1, 1.0), "A5": (0, 6, 1 / 6)}     # name: (p, q, prefactor)
UNIT_A = {"nm": 10.0, "um": 1e4, "mm": 1e7}


def haider_to_abtem(name: str, value_A: float, angle_deg: float | None = None) -> dict:
    """One Haider-notation coefficient (magnitude [A], CEOS angle [deg]) as abTEM's {C_nm: .., phi_nm: ..}."""
    import math
    p, q, pref = HAIDER[name]
    n, m = p + q - 1, abs(p - q)
    key = f"C{n}{m}"
    if m == 0:
        return {key: float(value_A)}
    ang = math.radians(angle_deg or 0.0)
    phi = (ang if q > p else -ang) / m
    return {key: float(value_A) * (n + 1) * pref, "phi" + key[1:]: phi % (2 * math.pi / m)}


def haider_chi(tab: dict, tx, ty):
    """chi / (2 pi / lambda) [A] of a Haider-notation tableau {name: (magnitude A, angle deg)}, straight from the complex
    formula above -- the reference haider_to_abtem() is checked against."""
    import numpy as np
    w = np.asarray(tx) + 1j * np.asarray(ty)
    out = np.zeros(w.shape)
    for name, (v, ang) in tab.items():
        p, q, pref = HAIDER[name]
        out += np.real(pref * v * np.exp(1j * np.radians(ang or 0.0)) * w ** p * np.conj(w) ** q)
    return out


ARM_TSV = "arm200f_tableau.tsv"
#: THE ARM200F-CLASS INSTRUMENT (2026-09-29, the user's lab microscope: JEOL JEM-ARM200F, 2009 CEOS CESCOR-type probe
#: corrector). Measured terms from campaign/arm200f_tableau.tsv, held FIXED: run3 (the last of three tableaus, 23 Sep
#: 2026) for A1 A2 B2 A3 S3 A4 -- one real, internally consistent snapshot rather than an average of wandering vectors --
#: and the reference Table 3 ("manual") for B4 D4 A5, the only source for them. C5 is factory set and NOT a knob: +4 mm
#: (the user: a well-tuned instrument, a little generous, inside the reference's 6 +- 4 mm). C1 and C3 are the only
#: terms the planner moves. Coefficients are kept in length units at 300 kV (the ARM runs at 200 kV): the same length is
#: 1.27x more waves at 300 kV, so this is slightly pessimistic; scale=0.785 keeps the waves at the tuning aperture.
ARM_RUN, ARM_MANUAL_TERMS, ARM_C5_A = "run3", ("B4", "D4", "A5"), 4.0e7


def arm_measured(run: str = ARM_RUN) -> dict:
    """{name: (magnitude A, angle deg or None, uncertainty A, kind, source)} of the ARM200F-class tableau, from the tsv."""
    import csv
    import os
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), ARM_TSV)
    rows = [r for r in csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t")]
    out = {}
    for r in rows:
        name = r["term"]
        want = "manual" if name in ARM_MANUAL_TERMS else run
        if r["source"] != want or name in ("C1", "C3", "C5"):
            continue
        f = UNIT_A[r["unit"]]
        ang = float(r["angle_deg"]) if r["angle_deg"] else None
        out[name] = (float(r["value"]) * f, ang, float(r["uncertainty"]) * f, r["uncertainty_kind"], r["source"])
    return out


def arm_tableau(c3: float, c5: float = ARM_C5_A, override: dict | None = None, scale: float = 1.0,
                run: str = ARM_RUN) -> dict:
    """abTEM aberration dict of the ARM200F-class tableau: round C30 (fought) and C50 (fixed), plus every measured
    non-round term converted from Haider notation. override {name: |X| [A], Haider notation}: replace a magnitude,
    keeping the measured angle -- for the open decisions (A1 nulled on the day: {"A1": 0}; B4/D4 at other values).
    scale multiplies every measured term and C5 (0.785 = the same waves at 300 kV as at 200 kV)."""
    ab = {"C30": float(c3), "C50": float(c5) * scale}
    for name, (v, ang, *_rest) in arm_measured(run).items():
        if override and name in override:
            v = float(override[name])
        for k, x in haider_to_abtem(name, v * scale, ang).items():
            ab[k] = round(x, 6) if k.startswith("phi") else round(x, 4)
    return ab


def growth_figure(out, alphas=(30, 40, 50, 60, 70, 80, 90, 100)):
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
    # No tolerance line. The "below 0.1 waves" figure one used to carry was withdrawn on 2026-09-22
    # (mis-oriented probe) and drawing any horizontal limit here invites it being read as measured.
    ax.axvline(30, color="#898781", lw=1, ls=(0, (3, 3)))
    ax.annotate("design aperture", (30.6, 0.012), fontsize=8.5, color="#52514e", ha="left", va="bottom")
    ax.set_yscale("log"); ax.set_ylim(0.01, 500); ax.set_xlim(min(alphas), max(alphas) + 18)
    ax.set_yticks([0.01, 0.1, 1, 10, 100]); ax.set_yticklabels(["0.01", "0.1", "1", "10", "100"])
    ax.set_xlabel("aperture semi-angle (mrad)")
    ax.set_ylabel("waves at the aperture edge")
    ax.set_title("How every non-round residual grows as the corrector is opened past its design aperture",
                 loc="left", fontsize=11, pad=10)
    ax.text(0.01, -0.16, "Blue: what the hardware leaves behind. Grey: what the operator tunes out at the design "
            "aperture, and which drifts (split per CEOS for a hexapole probe corrector). Red: six-fold, 1 mm, the one "
            "sourced value (CEOS); the rest are unsourced orders of magnitude, not a measured tableau. "
            "This figure is analytic: it depends on no reconstruction.",
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
    ap.add_argument("--ceos-json", type=float, default=None, metavar="C3",
                    help="print the CEOS-approx tableau as a sweep-row aber_json, for this round C3 [A]")
    a = ap.parse_args()

    if a.ceos_json is not None:
        print(json.dumps(ceos_tableau(a.ceos_json), separators=(",", ":")))
        return
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
