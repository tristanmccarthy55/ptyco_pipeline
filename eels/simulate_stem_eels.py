#!/usr/bin/env python
"""@file simulate_stem_eels.py
@brief Full STEM-mode EELS forward simulator: HAADF + a background-included core-loss spectrum.

Combines the two codes (see LITERATURE.md / HANDOVER.md):
  - CASTEP/OptaDOS  -> the intrinsic ELNES edge SHAPE + q-anisotropy (Channel A), INJECTED.
  - abtem multislice -> the STEM dynamics CASTEP cannot do: convergent probe (alpha), channelling
    through a thick specimen, frozen-phonon TDS, the scan, and SIMULTANEOUS HAADF (Channel B).

Full spectrum in the CHANNELLING (local) approximation -- the standard atomic-resolution-EELS
framework, and crucially gpaw-FREE (uses only ELASTIC multislice):

    I(scan, E) = sum_atoms [ probe intensity at atom ]  x  sigma_species(E; beta)  +  background(E)

  * [probe intensity at atom]  = abtem elastic multislice (channelling, thickness, phonons) -> weights
  * sigma_species(E; beta)     = OptaDOS q-resolved ELNES, aperture-averaged over the collection
                                 angle beta (the analyze_elnes.py magic-angle model), placed at the
                                 edge onset on the eV-loss axis
  * background(E)              = power-law A*E^-r now (config); abtem plasmon low-loss later
The rigorous inelastic transition-potential route (abtem + gpaw) is an optional Phase-2 refinement
selected by config.STEMEELS.inelastic_model = "transition_potential".

    ~/hyperspy-bundle/bin/python simulate_stem_eels.py --selftest      # spectrum core, real M4 data
    <blythe abtem env>/python  simulate_stem_eels.py --run            # + abtem channelling/HAADF
"""
from __future__ import annotations

import argparse
import os

import numpy as np

import config as C
import analyze_elnes as A


# ---------------------------------------------------------------- spectroscopy core (gpaw-free)
def _load_elnes(path: str) -> tuple[np.ndarray, np.ndarray]:
    """(energy, broadened intensity) for an injected ELNES. Handles both the raw OptaDOS
    `_core_edge.dat` (multi-block, `:exc` header -> select that block) and an already-extracted
    `.exc.txt` (plain 3-column: energy, raw, broadened)."""
    try:
        return A.load_optados_core(path)                 # raw multi-section
    except (ValueError, IndexError):
        d = np.loadtxt(path)
        return d[:, 0], d[:, 2]                           # extracted .exc.txt -> col3 = broadened


def aperture_averaged_elnes(qc_path: str, qperp_path: str, onset_eV: float, beta_mrad: float,
                            energy_keV: float, eloss: np.ndarray) -> np.ndarray:
    """Injected edge shape sigma(E): the OptaDOS q||c and q_perp ELNES, combined by the fraction
    of the collection aperture that probes q||c (analyze_elnes.parallel_weight -> the magic-angle
    model), then shifted so its onset sits at `onset_eV` and interpolated onto the eV-loss axis.
    beam is assumed along the polar axis c (the tet_Pz geometry)."""
    e_rel, s_par = _load_elnes(qc_path)                  # relative energy, broadened col
    _, s_perp = _load_elnes(qperp_path)
    tE = A.characteristic_angle_rad(onset_eV, energy_keV)
    f_par = A.parallel_weight(beta_mrad * 1e-3, tE)      # 1 at beta->0, ->1/3 at the magic angle
    s_meas = f_par * s_par + (1.0 - f_par) * s_perp
    onset_rel = e_rel[np.argmax(s_meas > 0.02 * s_meas.max())]
    e_abs = e_rel - onset_rel + onset_eV                 # place the onset at the tabulated edge energy
    return np.interp(eloss, e_abs, s_meas, left=0.0, right=0.0)


def powerlaw_background(eloss: np.ndarray, onset_eV: float, edge_jump: float,
                        r: float, frac: float) -> np.ndarray:
    """Standard EELS power-law background A*E^-r, scaled so its height at the first edge onset is
    `frac` x the edge jump (how experimental spectra are modelled/subtracted)."""
    e = np.where(eloss > 0, eloss, np.nan)
    bg = e ** (-r)
    A_ = frac * edge_jump / (onset_eV ** (-r))
    return np.nan_to_num(A_ * bg)


def assemble_spectrum(eloss: np.ndarray, weights: dict, cfg: C.STEMEELS) -> dict:
    """Full spectrum = sum over edges of (channelling weight x aperture-averaged ELNES) + background.
    `weights` maps edge -> scalar channelling weight (sum of probe intensity over that species'
    atoms); default 1.0. Returns {'eloss','edges','background','total'}."""
    edges = {}
    for edge in cfg.edges:
        src = cfg.elnes_source.get(edge)
        if src is None:
            raise ValueError(f"no elnes_source for edge {edge!r} (set cfg.elnes_source[{edge!r}])")
        shape = aperture_averaged_elnes(src["qc"], src["qperp"], cfg.edge_onset_eV[edge],
                                        cfg.collection_mrad, cfg.energy_keV, eloss)
        edges[edge] = weights.get(edge, 1.0) * shape
    total_edges = np.sum(list(edges.values()), axis=0)
    bg = np.zeros_like(eloss)
    if cfg.background_model == "powerlaw":
        first = cfg.edge_onset_eV[cfg.edges[0]]
        bg = powerlaw_background(eloss, first, total_edges.max(), cfg.background_r, cfg.background_frac)
    total = total_edges + bg
    if cfg.dose_e_per_A2 is not None:                    # optional Poisson noise
        rng = np.random.default_rng(0)
        scale = cfg.dose_e_per_A2 * cfg.scan_step_A ** 2 / max(total.max(), 1e-30)
        total = rng.poisson(np.clip(total * scale, 0, None)) / scale
    return {"eloss": eloss, "edges": edges, "background": bg, "total": total}


def eloss_axis(cfg: C.STEMEELS) -> np.ndarray:
    return np.arange(cfg.eloss_min_eV, cfg.eloss_max_eV, cfg.eloss_dispersion_eV)


# ---------------------------------------------------------------- STEM dynamics (REUSE sim/, Blythe)
def _import_sim4d():
    """Import sim/simulate_4dstem.py so we REUSE its scattering (no duplication)."""
    import sys
    sim_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "sim"))
    if sim_dir not in sys.path:
        sys.path.insert(0, sim_dir)
    import simulate_4dstem as s4                        # noqa: E402
    return s4


def eels_collection_mrad(s4) -> float:
    """The EELS 'hole': the spectrometer takes the inner HALF of the pixelated-detector radius;
    the outer annulus [beta, max] is the simultaneous HAADF. beta = DETECTOR_MAX_ANGLE_MRAD / 2."""
    return s4.DETECTOR_MAX_ANGLE_MRAD / 2.0


def probe_channelling(cfg: C.STEMEELS, atoms=None):
    """REUSE sim/simulate_4dstem.py for ALL the scattering — probe (its overfocused 100 mrad
    probe), Lobato potential, frozen-phonon TDS, multislice, scan — and add the EELS geometry:
    a detector 'hole' of radius beta = DETECTOR_MAX_ANGLE_MRAD/2 sends the central cone to the
    spectrometer; the [beta, max] annulus is the SIMULTANEOUS HAADF. Returns
    (haadf_map, eels_weights_by_species, beta_mrad), where the weights are the channelling-summed
    core-loss coupling per species (each multiplies that edge's injected CASTEP ELNES in
    assemble_spectrum). Runs on the Blythe abtem env; the core-loss coupling uses abtem's
    transition_potential_scan (same multislice engine) -> needs gpaw. NOT tested locally."""
    import abtem
    from abtem.inelastic.core_loss import SubshellTransitions
    s4 = _import_sim4d()
    s4.DEVICE = cfg.device
    if atoms is None:
        atoms, _ = s4.load_and_prepare_atoms()          # the labyrinth, beam = z (as ptychography)
    beta = eels_collection_mrad(s4)                     # EELS hole = half the detector radius
    pot = s4.build_potential(atoms, announce=True)      # reuse sim's potential
    probe = s4.build_probe(pot)                         # reuse sim's overfocused 100 mrad probe
    scan, pos_xy, ny, _ = s4.make_scan(test_mode=False)

    # simultaneous HAADF = the outer annulus [beta, detector_max] of the SAME elastic scan
    haadf_det = abtem.AnnularDetector(inner=beta, outer=s4.DETECTOR_MAX_ANGLE_MRAD)
    haadf = probe.scan(pot, scan=scan, detectors=haadf_det).compute()

    # core-loss coupling per edge species, collected within the beta hole (abtem + gpaw).
    # Frozen phonons: wrap FrozenPhonons(atoms) as sim does; here one config for the interface.
    weights = {}
    for edge in cfg.edges:
        Z, n, l = _edge_nl(edge)                        # (Z, n, l) of the ionised subshell
        tp = SubshellTransitions(Z=Z, n=n, l=l).get_transition_potentials(
            extent=pot.extent, gpts=pot.gpts, energy=probe.energy)
        eels_det = abtem.FlexibleAnnularDetector()      # integrate within beta at reduce time
        m = probe.transition_potential_scan(potential=pot, transition_potentials=tp,
                                            scan=scan, detectors=eels_det).compute()
        weights[edge] = m                               # energy-integrated core-loss map (per scan)
    return haadf, weights, beta


def _edge_nl(edge):
    return {"O_K": (8, 1, 0), "Ti_L23": (22, 2, 1), "Pb_M": (82, 3, 2)}[edge]


# ---------------------------------------------------------------- self-test (spectroscopy core)
def selftest() -> None:
    """Build a full O-K spectrum (edge + background) from the REAL M4 OptaDOS data if present on
    the Desktop, else a synthetic ELNES -- and show the collection-angle (beta) dependence. No abtem."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    print("== simulate_stem_eels self-test (spectrum core, gpaw-free) ==")

    qc = os.path.expanduser("~/Desktop/tet_Pz_Oap.qc_core_edge.exc.txt")
    qp = os.path.expanduser("~/Desktop/tet_Pz_Oap.qperp_core_edge.exc.txt")
    tmp = None
    if not (os.path.exists(qc) and os.path.exists(qp)):
        import tempfile                                  # synthesise a 2-peak sigma/pi ELNES pair
        tmp = tempfile.mkdtemp(); e = np.linspace(-5, 40, 4500)
        g = lambda c, w, a: a * np.exp(-0.5 * ((e - c) / w) ** 2)
        spar = g(3, 1.5, 1.0) + g(11, 3, 0.8); sperp = g(0.8, 1.2, 1.0) + g(9, 3, 0.5)
        qc = os.path.join(tmp, "qc.dat"); qp = os.path.join(tmp, "qp.dat")
        for p, s in ((qc, spar), (qp, sperp)):
            with open(p, "w") as f:
                f.write(" #  O 1 K1 O:exc\n")
                for ei, si in zip(e, s):
                    f.write(f" {ei:.4f} 0 {si:.6e}\n")
        print("  (real M4 files not on Desktop -> using a synthetic sigma/pi ELNES)")
    else:
        print("  using REAL M4 data:", os.path.basename(qc))

    from dataclasses import replace
    cfg = replace(C.STEMEELS_CFG, edges=("O_K",), elnes_source={"O_K": {"qc": qc, "qperp": qp}},
                  eloss_min_eV=520.0, eloss_max_eV=580.0)
    E = eloss_axis(cfg)
    fig, ax = plt.subplots(figsize=(8, 5))
    for beta in (2, 10, 25, 100):                        # collection-angle sweep
        out = assemble_spectrum(E, weights={"O_K": 1.0}, cfg=replace(cfg, collection_mrad=beta))
        ax.plot(E, out["total"], label=f"β={beta} mrad")
    base = assemble_spectrum(E, weights={"O_K": 1.0}, cfg=cfg)
    ax.plot(E, base["background"], "k:", lw=1, label="power-law background")
    ax.set_xlabel("energy loss (eV)"); ax.set_ylabel("intensity"); ax.legend()
    ax.set_title("STEM-EELS O-K: full spectrum (edge + background) vs collection angle β")
    out_png = os.path.join(os.path.dirname(__file__), "..", "..", "stem_eels_selftest.png")
    out_png = os.path.abspath(os.path.join(os.path.dirname(__file__), "stem_eels_selftest.png"))
    plt.tight_layout(); plt.savefig(out_png, dpi=120)
    # sanity assertions
    assert base["total"].max() > base["background"].max(), "edge should rise above background"
    assert np.all(base["background"][E < cfg.edge_onset_eV["O_K"]] >= 0), "background non-negative"
    print(f"  built full O-K spectrum (edge+background) for β=2..100 mrad -> {out_png}")
    print("== self-test PASSED ==")
    if tmp:
        import shutil; shutil.rmtree(tmp)


def main() -> None:
    ap = argparse.ArgumentParser(description="Full STEM-EELS forward simulator (HAADF + spectrum).")
    ap.add_argument("--selftest", action="store_true", help="spectrum core on real M4 data (no abtem)")
    ap.add_argument("--run", action="store_true", help="full sim incl. abtem channelling (Blythe)")
    ap.add_argument("--cell", default="tet_Pz")
    args = ap.parse_args()
    if args.run:
        haadf, weights = probe_channelling(args.cell, C.STEMEELS_CFG)   # abtem, Blythe
        out = assemble_spectrum(eloss_axis(C.STEMEELS_CFG), weights, C.STEMEELS_CFG)
        print("HAADF + spectrum computed:", out["total"].shape)
    else:
        selftest()


if __name__ == "__main__":
    main()
