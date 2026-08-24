#!/usr/bin/env python
"""@file analyze_elnes.py
@brief M4-M6: turn OptaDOS q-resolved ELNES into the along-beam-polarisation answer.

Three jobs:
  M4  dichroism   Delta(E) = S(q||c) - S(q_perp): the intrinsic along-beam signal.
  M5  calibration Delta metric vs Ti displacement (|P|) across the scan_* series.
  M6  detectability: average the intrinsic dichroism over a real EELS collection aperture
      at 300 keV (the sim geometry) -> how much SURVIVES vs the magic angle, and the SNR /
      dose / aperture needed to see it. This is the honest, pipeline-matched conclusion.

The geometry model (M6) is validated HERE on synthetic spectra (`--selftest`): the surviving
anisotropy must cross zero at the textbook magic angle ~3.97*theta_E, so the physics is
checked before any real spectrum exists. Real OptaDOS spectra drop in as
<seed>.qc_core_edge.dat / <seed>.qperp_core_edge.dat, pulled back from Blythe into
eels/runs/ (per-run subdirectories are searched too).

NB the .odi files set a LAB-frame q, so for a cell whose polar axis is not along the beam
the two files' crystal-frame meanings are exchanged -- see q_is_swapped().

    ~/hyperspy-bundle/bin/python analyze_elnes.py --selftest
    ~/hyperspy-bundle/bin/python analyze_elnes.py --seed tet_Pz_Oap --edge O_K
    ~/hyperspy-bundle/bin/python analyze_elnes.py --compare tet_Pz_Oap tet_Px_Oap
"""
from __future__ import annotations

import argparse
import glob
import os

import numpy as np
from scipy.integrate import quad

import config as C
from config import OPTICS

# NumPy 2.0 renamed np.trapz -> np.trapezoid (old name removed). Support both.
_trapz = getattr(np, "trapezoid", None) or np.trapz

# edge onsets (eV) for theta_E; refined against the M2 benchmark spectra.
EDGE_ONSET_eV = {"O_K": 532.0, "Ti_L23": 456.0, "Pb_M": 2484.0}


# ---------------------------------------------------------------- electron-optics geometry
def characteristic_angle_rad(edge_eV: float, E0_keV: float = 300.0) -> float:
    """Characteristic inelastic angle theta_E = dE / (gamma m0 v^2) (relativistic)."""
    mc2 = 510_998.9                              # electron rest energy, eV
    E0 = E0_keV * 1e3
    gamma = 1.0 + E0 / mc2
    beta2 = 1.0 - 1.0 / gamma**2                 # (v/c)^2
    return edge_eV / (gamma * mc2 * beta2)       # rad


def parallel_weight(beta_rad: float, theta_E: float) -> float:
    """Fraction of the collected dipole signal that probes S_parallel when the BEAM is along
    the polar axis c. q makes angle phi with c where cos^2 phi = theta_E^2/(theta^2+theta_E^2);
    the dipole angular distribution weights each theta by theta/(theta^2+theta_E^2). Returns
    W_par/(W_par+W_perp) in [0,1]; = 1 at beta->0 (pure q||c), -> 1/3 (isotropic) at the magic
    angle, < 1/3 beyond it (anisotropy inverts)."""
    def w_par(t):  return theta_E**2 / (t**2 + theta_E**2)**2 * t
    def w_perp(t): return t**2 / (t**2 + theta_E**2)**2 * t
    Wp = quad(w_par, 0, beta_rad)[0]
    Wperp = quad(w_perp, 0, beta_rad)[0]
    return Wp / (Wp + Wperp)


def surviving_anisotropy(beta_rad: float, theta_E: float) -> float:
    """Signed fraction of the intrinsic dichroism (S_par - S_perp) that survives collection at
    semi-angle beta, relative to the beam||c isotropic reference. f_par - 1/3, renormalised so
    beta->0 gives 1.0 (full dichroism) and the magic angle gives 0."""
    return (parallel_weight(beta_rad, theta_E) - 1.0 / 3.0) / (1.0 - 1.0 / 3.0)


def magic_angle_rad(theta_E: float) -> float:
    """Collection semi-angle where the anisotropic term vanishes (surviving_anisotropy=0)."""
    from scipy.optimize import brentq
    return brentq(lambda b: surviving_anisotropy(b, theta_E), 1e-6, 200 * theta_E)


# ---------------------------------------------------------------- spectra I/O + dichroism
def load_spectrum(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Two-column (energy, intensity) text; skips # / ! comments. For generic/synthetic files;
    for a real OptaDOS core-loss .dat use load_optados_core (multi-section)."""
    data = np.loadtxt(path, comments=("#", "!"))
    return data[:, 0], data[:, 1]


def load_optados_core(path: str, section: str = ":exc",
                      broadened: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Parse an OptaDOS `<seed>_core_edge.dat`. It holds ONE block per atom, each preceded by a
    header like `# O 1 K1` or `# O 1 K1 O:exc`; columns are (energy, raw, broadened). The
    PHYSICAL core-hole ELNES is the excited atom's block — the header containing `section`
    (default ':exc'). Non-excited blocks use the perturbed conduction states with un-cored
    atoms and are NOT the ELNES, so we must select, never sum. Returns (energy, intensity)
    for the matched block (broadened col by default). Verified on the M2(b) TiO2 O-K benchmark."""
    col = 2 if broadened else 1
    e, y, grabbing = [], [], False
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("#"):
                grabbing = section in s          # (re)start capture iff this header matches
                continue
            if grabbing and s:
                p = s.split()
                if len(p) > col:
                    try:
                        e.append(float(p[0])); y.append(float(p[col]))
                    except ValueError:
                        pass
    if not e:
        raise ValueError(f"no core-loss block matching {section!r} in {path}")
    return np.asarray(e), np.asarray(y)


def dichroism(e: np.ndarray, s_par: np.ndarray, s_perp: np.ndarray) -> np.ndarray:
    return s_par - s_perp


def dichroism_metric(e: np.ndarray, delta: np.ndarray, window: tuple | None = None) -> float:
    """Unnormalised scalar size of the dichroism: integral |Delta| dE over an (optional)
    energy window. For the FRACTIONAL metrics quoted in RESULTS.md (max|D|/max S and
    int|D|/int S) use dichroism_fractions(), which divides by the edge itself."""
    m = np.ones_like(e, bool) if window is None else (e >= window[0]) & (e <= window[1])
    return float(_trapz(np.abs(delta[m]), e[m]))


# ---------------------------------------------------------------- detectability
def required_counts(anisotropy_frac: float, snr: float = 3.0) -> float:
    """Poisson counts per spectral channel to detect a fractional anisotropy a at target SNR
    (difference of two spectra ~ a*I over noise ~ sqrt(2I))  ->  N >= 2 (snr/a)^2."""
    if anisotropy_frac <= 0:
        return float("inf")
    return 2.0 * (snr / anisotropy_frac) ** 2


def geometry_report(edge: str) -> None:
    """M6: surviving dichroism vs collection aperture at 300 keV, incl. the sim's 100 mrad."""
    tE = characteristic_angle_rad(EDGE_ONSET_eV[edge], OPTICS.energy_keV)
    ma = magic_angle_rad(tE)
    print(f"\n== M6 geometry ({edge}, {OPTICS.energy_keV:.0f} keV) ==")
    print(f"  theta_E = {tE*1e3:.2f} mrad | magic angle = {ma*1e3:.2f} mrad "
          f"(~{ma/tE:.2f} theta_E)")
    print(f"  {'beta(mrad)':>10} {'survive':>9} {'note':>24}")
    for b in OPTICS.collection_sweep_mrad:
        f = surviving_anisotropy(b * 1e-3, tE)
        note = "full-ish" if f > 0.7 else ("near magic (washed)" if abs(f) < 0.1
                                           else ("inverted" if f < 0 else "reduced"))
        print(f"  {b:>10.0f} {f:>9.2f} {note:>24}")
    fconv = surviving_anisotropy(OPTICS.convergence_mrad * 1e-3, tE)
    print(f"  -> your ptychography {OPTICS.convergence_mrad:.0f} mrad probe: surviving "
          f"anisotropy {fconv:+.2f} (dichroism is averaged away/inverted -> use a small "
          f"EELS collection aperture, ideally << {ma*1e3:.1f} mrad).")


# ---------------------------------------------------------------- self-test (validates M6)
def selftest() -> None:
    """Validate the geometry model on synthetic edges: (1) magic angle ~ 3.97 theta_E,
    (2) surviving anisotropy is 1 at beta->0 and negative beyond the magic angle,
    (3) required-counts blows up as anisotropy -> 0. No HPC output needed."""
    print("== analyze_elnes self-test (geometry model on synthetic spectra) ==")
    e = np.linspace(525, 575, 500)                       # synthetic O K window
    g = lambda c, w, a: a * np.exp(-0.5 * ((e - c) / w) ** 2)
    s_perp = g(533, 1.2, 1.0) + g(540, 3.0, 0.8)
    s_par = g(533, 1.2, 1.0) + g(538, 3.0, 1.1)          # anisotropic pre-peak shift
    d = dichroism(e, s_par, s_perp)
    metric = dichroism_metric(e, d, (530, 545))
    assert metric > 0, "synthetic dichroism should be nonzero"

    tE = characteristic_angle_rad(EDGE_ONSET_eV["O_K"], 300.0)
    ma = magic_angle_rad(tE)
    print(f"  theta_E(O K,300kV) = {tE*1e3:.3f} mrad ; magic angle = {ma/tE:.3f} theta_E "
          f"({ma*1e3:.2f} mrad)")
    assert 3.5 < ma / tE < 4.4, f"magic angle {ma/tE:.2f} theta_E off the ~3.97 expectation"
    assert abs(surviving_anisotropy(1e-5, tE) - 1.0) < 1e-3, "beta->0 must give full dichroism"
    assert surviving_anisotropy(20 * tE, tE) < 0, "should invert well past the magic angle"

    a_small = 0.02
    print(f"  required counts/channel for {a_small:.0%} anisotropy at SNR 3 = "
          f"{required_counts(a_small):.0f}")
    assert required_counts(0.2) < required_counts(0.02), "smaller anisotropy needs more counts"
    geometry_report("O_K")
    print("== self-test PASSED ==")


# ---------------------------------------------------------------- real-data driver
def _runs_dir() -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs")


def find_core_dat(seed: str, qtag: str) -> str | None:
    """@brief Locate an OptaDOS core-loss .dat for (seed, qtag) under eels/runs/.

    run_coreloss.slurm writes `<seed>.qc_core_edge.dat` / `<seed>.qperp_core_edge.dat` INSIDE a
    per-run subdirectory (runs/tetPz_Oap/...), so both levels are searched. Note the underscore:
    an earlier pattern of `<seed>.qc.*core*.dat` could never match `<seed>.qc_core_edge.dat`.
    Also accepts the hand-extracted `.exc.txt` form.
    """
    pats = [f"{seed}.{qtag}*core*.dat", f"{seed}.{qtag}*core*.exc.txt"]
    hits = sorted(h for p in pats
                  for h in glob.glob(os.path.join(_runs_dir(), p))
                  + glob.glob(os.path.join(_runs_dir(), "*", p)))
    return hits[0] if hits else None


def _load_edge(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Excited-atom (:exc) broadened ELNES from either an OptaDOS multi-block `_core_edge.dat`
    or an already-extracted 3-column `.exc.txt`. Never sums the per-atom blocks."""
    try:
        return load_optados_core(path)
    except (ValueError, IndexError):
        d = np.loadtxt(path)
        return d[:, 0], d[:, 2]


def cell_of(seed: str):
    """`tet_Px_Oap` -> the config.CELLS entry for `tet_Px` (excited-site suffix stripped).
    Returns None for the scan_* ladder, which is not in CELLS and is uniformly c||z."""
    for suffix in ("_Oap", "_Oeq", "_Ti", "_Pb", "_O"):
        if seed.endswith(suffix):
            seed = seed[: -len(suffix)]
            break
    return C.CELLS.get(seed)


def q_is_swapped(seed: str) -> bool:
    """@brief Does this cell's crystal frame invert the meaning of the two .odi files?

    The .odi files set a LAB-frame q: `coreloss_qc.odi` is `core_qdir 0 0 1` (along the beam,
    = C.Q_BEAM) and `coreloss_qperp.odi` is `1 0 0` (= C.Q_PERP). For a cell whose polar axis
    is c||z -- tet_Pz, cubic, the scan_* ladder -- those coincide with q||c and q-perp-c. For
    tet_Px, whose c lies along x, THE TWO ROLES ARE EXCHANGED: its `.qc` file is q-perp-c and
    its `.qperp` file is q||c.

    This matters for the M4 rotational-invariance cross-check, whose correct pairing is
        tet_Pz .qc    (q||c)      ==  tet_Px .qperp (q||c)
        tet_Pz .qperp (q-perp-c)  ==  tet_Px .qc    (q-perp-c)
    Pairing by filename instead compares a sigma* spectrum against a pi* one and reads as a
    catastrophic failure of a calculation that is in fact correct.

    `real` (polar_axis "off") is deliberately NOT swapped: its polar axis is tilted, so the
    physically meaningful split is the lab one, along-beam versus in-plane.
    """
    spec = cell_of(seed)
    return bool(spec and spec.polar_axis == "x")


def crystal_frame_spectra(seed: str, swap: bool | None = None):
    """(energy, S_par_to_c, S_perp_to_c) for a seed, with the lab->crystal q mapping applied.
    `swap=None` derives it from config.CELLS via q_is_swapped()."""
    p_qc, p_qp = find_core_dat(seed, "qc"), find_core_dat(seed, "qperp")
    if not p_qc or not p_qp:
        return None
    e, s_qc = _load_edge(p_qc)
    e2, s_qp = _load_edge(p_qp)
    if len(e2) != len(e) or not np.allclose(e2, e):
        s_qp = np.interp(e, e2, s_qp)                 # same SCF, but do not assume one grid
    if q_is_swapped(seed) if swap is None else swap:
        s_qc, s_qp = s_qp, s_qc
    return e, s_qc, s_qp


def dichroism_fractions(e, s_par, s_perp, window=None) -> dict:
    """@brief The two fractional metrics quoted for M4: max|Delta|/max S and int|Delta|/int S.

    S is the MEAN of the two orientations, so the fraction does not depend on which spectrum
    is nominated as parallel. `dichroism_metric` returns the UNnormalised int|Delta| dE.
    """
    m = np.ones_like(e, bool) if window is None else (e >= window[0]) & (e <= window[1])
    d = np.abs(np.asarray(s_par) - np.asarray(s_perp))[m]
    s = (0.5 * (np.asarray(s_par) + np.asarray(s_perp)))[m]
    return {"peak": float(d.max() / s.max()),
            "integral": float(_trapz(d, e[m]) / _trapz(s, e[m]))}


def analyze_seed(seed: str, edge: str, swap: bool | None = None) -> None:
    got = crystal_frame_spectra(seed, swap=swap)
    if got is None:
        print(f"[no OptaDOS output yet for {seed}] pull <seed>.qc_core_edge.dat / "
              f"<seed>.qperp_core_edge.dat from Blythe into eels/runs/. "
              f"Running --selftest instead so the pipeline is validated.")
        return selftest()
    e, s_par, s_perp = got
    win = (EDGE_ONSET_eV[edge] - 2, EDGE_ONSET_eV[edge] + 30)
    fr = dichroism_fractions(e, s_par, s_perp, win)
    sw = q_is_swapped(seed) if swap is None else swap
    print(f"{seed} ({edge}){'  [q axes swapped: c||x]' if sw else ''}")
    print(f"  max|D|/max S = {fr['peak']:.1%}   int|D|/int S = {fr['integral']:.1%}")
    geometry_report(edge)


def compare_seeds(seed_a: str, seed_b: str, edge: str) -> None:
    """@brief M4 rotational-invariance cross-check: two cells must agree in the CRYSTAL frame.

    tet_Pz and tet_Px are the same crystal in two orientations, so once each seed's lab q is
    mapped back onto its own c axis (see q_is_swapped) their q||c spectra must coincide, and
    likewise their q-perp-c spectra. The residual is quoted on the same scale as the M4
    dichroism, so it is directly comparable with the M3 cubic null test (0.02%).
    """
    ga, gb = crystal_frame_spectra(seed_a), crystal_frame_spectra(seed_b)
    missing = [s for s, g in ((seed_a, ga), (seed_b, gb)) if g is None]
    if missing:
        print(f"[invariance check] no OptaDOS output yet for: {', '.join(missing)}")
        return
    e, a_par, a_perp = ga
    eb, b_par, b_perp = gb
    if len(eb) != len(e) or not np.allclose(eb, e):
        b_par, b_perp = np.interp(e, eb, b_par), np.interp(e, eb, b_perp)
    win = (EDGE_ONSET_eV[edge] - 2, EDGE_ONSET_eV[edge] + 30)
    print(f"== M4 rotational invariance: {seed_a} vs {seed_b} ({edge}) ==")
    print(f"   pairing is crystal-frame; swapped: {seed_a}={q_is_swapped(seed_a)}, "
          f"{seed_b}={q_is_swapped(seed_b)}")
    worst = 0.0
    for lab, x, y in (("q||c    ", a_par, b_par), ("q-perp-c", a_perp, b_perp)):
        fr = dichroism_fractions(e, x, y, win)
        worst = max(worst, fr["peak"])
        print(f"   {lab}  residual max|D|/max S = {fr['peak']:.2%}   "
              f"int|D|/int S = {fr['integral']:.2%}")
    verdict = "PASS" if worst < 0.02 else ("MARGINAL" if worst < 0.05 else "FAIL")
    print(f"   -> {verdict} (M3 cubic null test sets the numerical floor at ~0.02%; a residual "
          f"comparable to the 78% dichroism means the q axes were paired by filename)")


def main() -> None:
    ap = argparse.ArgumentParser(description="ELNES dichroism + detectability (M4-M6).")
    ap.add_argument("--selftest", action="store_true", help="validate the geometry model now")
    ap.add_argument("--seed", help="OptaDOS seed, e.g. tet_Pz_Oap")
    ap.add_argument("--edge", default="O_K", choices=list(EDGE_ONSET_eV))
    ap.add_argument("--compare", nargs=2, metavar=("SEED_A", "SEED_B"),
                    help="M4 rotational-invariance cross-check, e.g. "
                         "--compare tet_Pz_Oap tet_Px_Oap")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--swap-q", dest="swap", action="store_true", default=None,
                   help="force the lab->crystal q swap (default: derived from config.CELLS)")
    g.add_argument("--no-swap-q", dest="swap", action="store_false",
                   help="force NO q swap")
    args = ap.parse_args()
    if args.compare:
        compare_seeds(args.compare[0], args.compare[1], args.edge)
    elif args.selftest or not args.seed:
        selftest()
    else:
        analyze_seed(args.seed, args.edge, swap=args.swap)


if __name__ == "__main__":
    main()
