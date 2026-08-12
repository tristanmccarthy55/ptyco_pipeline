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

# abtem 1.0.9's core_loss.py calls np.trapz, which NumPy 2.0 REMOVED (renamed to np.trapezoid).
# Alias it back so abtem's transition-potential build works on numpy>=2 without downgrading numpy
# (which would risk the elastic sim / cupy / gpaw stack that already works).
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid

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


# ---------------------------------------------------------------- membrane sample (up/down P_z)
def build_membrane(n_lat: int = 8, n_thick: int = 4, vacuum_A: float = 6.0,
                   domains: str = "updown"):
    """A thin PbTiO3 membrane of the SIMPLE unit cell (tet_Pz: polar c || beam z) as a PERIODIC
    slab: n_lat x n_lat cells laterally (periodic in x,y), n_thick cells thick, free-standing
    (vacuum along the beam z). Polarisation DOMAINS along z: `updown` = left half P_z UP (+z) /
    right half DOWN (-z), split down the middle (periodic -> a wall at the middle and at the box
    edge); `up`/`down` = uniform. Beam = z (cell already oriented -> NO rotation, so the labyrinth
    load_and_prepare_atoms path is untouched). The data-fusion testbed: ptychography reads the
    domain (direction) from the depth-resolved DP; EELS reads |P| (magnitude) from the O-K ELNES.
    Scan its CENTRE (the run picks a small window across the middle wall)."""
    import build_cells as B
    from ase import Atoms
    up = B.make_cell(C.CELLS["tet_Pz"])                  # P4mm, polar along z, 5 atoms
    a, _, c = up.cell.lengths()
    down = up.copy()                                     # mirror z -> polar flips sign
    sp = down.get_scaled_positions(); sp[:, 2] = (-sp[:, 2]) % 1.0
    down.set_scaled_positions(sp)

    memb = Atoms(cell=[n_lat * a, n_lat * a, n_thick * c], pbc=True)   # periodic slab
    for ix in range(n_lat):
        base = up if (domains == "up" or (domains == "updown" and ix < n_lat // 2)) \
            else (down if domains in ("down", "updown") else up)
        for iy in range(n_lat):
            for iz in range(n_thick):
                cell_ij = base.copy(); cell_ij.translate([ix * a, iy * a, iz * c])
                memb += cell_ij
    memb.center(vacuum=vacuum_A, axis=2)                 # free-standing along the beam only
    print(f"[membrane] {len(memb)} atoms | {n_lat}x{n_lat} cells x {n_thick} thick | "
          f"box {np.round(memb.cell.lengths(),1)} A (periodic x,y; vacuum z) | domains={domains}")
    return memb, float(memb.cell.lengths()[0])


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


def run_membrane(out_dir: str, cfg: C.STEMEELS = None, n_lat: int = 8, n_thick: int = 6,
                 domains: str = "updown", scan_gpts: tuple = (5, 5), device: str = "gpu") -> str:
    """Membrane STEM-EELS test -> ITS OWN FOLDER: the pixelated 4D-STEM DP (for ptychography +
    other uses) AND a full O-K EELS spectrum-image (injected CASTEP ELNES + power-law background).
    Reuses sim/simulate_4dstem.py for the scattering (potential/probe); the EELS map uses abtem
    transition_potential_scan (needs gpaw). Small quick scan across the up/down domain wall.
    NOTE: first Blythe run of the abtem/gpaw path -- abtem-1.0.9 detector API may need a tweak."""
    import os
    import abtem
    from abtem.inelastic.core_loss import SubshellTransitions
    from abtem.multislice import transition_potential_multislice_and_detect   # abtem 1.0.9 API
    cfg = cfg or C.STEMEELS_CFG
    # default injected ELNES = the M4 tet_Pz O-K (on Blythe: runs/exc/)
    here = os.path.dirname(__file__)
    src = cfg.elnes_source.get("O_K") or {
        "qc": os.path.join(here, "runs/exc/tet_Pz_Oap.qc_core_edge.exc.txt"),
        "qperp": os.path.join(here, "runs/exc/tet_Pz_Oap.qperp_core_edge.exc.txt")}
    s4 = _import_sim4d(); s4.DEVICE = device
    os.makedirs(out_dir, exist_ok=True)

    memb, box = build_membrane(n_lat, n_thick, domains=domains)
    pot = s4.build_potential(memb, announce=True)        # reuse sim's Lobato potential
    probe = s4.build_probe(pot)                          # reuse sim's overfocused 100 mrad probe
    beta = eels_collection_mrad(s4)                      # EELS hole = det/2 = 100 mrad
    ctr = box / 2.0                                      # scan a small window ACROSS the mid wall
    scan = abtem.GridScan(start=(ctr - 6, ctr - 2), end=(ctr + 6, ctr + 2), gpts=scan_gpts)
    pos = np.asarray(scan.get_positions()).reshape(-1, 2)

    print(f"[membrane] pixelated 4D-STEM ({scan_gpts} scan) ...")
    dp = probe.scan(pot, scan=scan,
                    detectors=abtem.PixelatedDetector(max_angle=s4.DETECTOR_MAX_ANGLE_MRAD)).compute()
    dp_arr = np.asarray(dp.array).astype(np.float32)     # (Nscan..., Kx, Ky) pixelated DP

    print("[membrane] O-K core-loss (transition potentials + gpaw) ...")
    tp = SubshellTransitions(Z=8, n=1, l=0).get_transition_potentials(
        extent=pot.extent, gpts=pot.gpts, energy=probe.energy)
    # lazy=False -> eager numpy waves: abtem 1.0.9's loss-accumulation does an in-place += of a
    # dask array into a numpy array, which newer dask rejects; eager arrays avoid that path.
    waves = probe.build(scan, lazy=False)
    eels_out = transition_potential_multislice_and_detect(
        waves, pot, tp, detectors=[abtem.FlexibleAnnularDetector()])
    eels = eels_out[0] if isinstance(eels_out, (list, tuple)) else eels_out
    eels = eels.compute() if hasattr(eels, "compute") else eels
    weight = np.asarray(eels.integrate_radial(0.0, beta).array).reshape(-1)   # per-scan-pixel coupling

    # inject the CASTEP ELNES energy shape (aperture-averaged over β) + power-law background
    E = eloss_axis(cfg)
    shape = aperture_averaged_elnes(src["qc"], src["qperp"], cfg.edge_onset_eV["O_K"],
                                    beta, cfg.energy_keV, E)
    shape = shape / max(shape.max(), 1e-30)
    spectra = np.outer(weight, shape)                    # (Nscan, Neloss): edge scaled by channelling
    bg = powerlaw_background(E, cfg.edge_onset_eV["O_K"], spectra.max(), cfg.background_r, cfg.background_frac)
    spectra = (spectra + bg).astype(np.float32)

    np.save(os.path.join(out_dir, "dp.npy"), dp_arr)     # the pixelated DP you asked for
    np.savez(os.path.join(out_dir, "eels.npz"), eloss=E, spectra=spectra, weight=weight,
             positions=pos, background=bg.astype(np.float32), beta_mrad=beta,
             box_A=box, n_lat=n_lat, n_thick=n_thick, domains=domains)
    print(f"[membrane] DONE -> {out_dir}\n  dp.npy {dp_arr.shape} (pixelated DP)"
          f"\n  eels.npz: spectra {spectra.shape} (scan × e-loss), + weight/positions/background")
    return out_dir


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
    ap.add_argument("--membrane", action="store_true",
                    help="membrane test (Blythe/gpaw): pixelated DP + O-K EELS spectrum-image -> --out-dir")
    ap.add_argument("--out-dir", default=None, help="output folder (default stem_eels_out_<tag>/)")
    ap.add_argument("--n-lat", type=int, default=8, help="membrane lateral cells (x=y)")
    ap.add_argument("--n-thick", type=int, default=6, help="membrane thickness in cells (along beam)")
    ap.add_argument("--domains", default="updown", choices=["updown", "up", "down"])
    ap.add_argument("--scan", default="5x5", help="scan gpts as NxN (small = quick, e.g. 3x3)")
    ap.add_argument("--device", default="gpu", choices=["gpu", "cpu"])
    args = ap.parse_args()
    if args.membrane:
        gx, gy = (int(v) for v in args.scan.lower().split("x"))
        tag = f"membrane_{args.domains}_{args.n_lat}x{args.n_lat}x{args.n_thick}_{gx}x{gy}"
        out = args.out_dir or os.path.join(os.path.dirname(__file__), "runs", f"stem_eels_{tag}")
        run_membrane(out, n_lat=args.n_lat, n_thick=args.n_thick, domains=args.domains,
                     scan_gpts=(gx, gy), device=args.device)
    else:
        selftest()


if __name__ == "__main__":
    main()
