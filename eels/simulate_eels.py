#!/usr/bin/env python
"""@file simulate_eels.py
@brief M6b: dynamical multislice STEM-EELS forward model for the PbTiO3 cells (abtem).

Complements the CASTEP+OptaDOS route. Two physically distinct channels carry along-beam
polarisation into EELS:
  Channel A (spectroscopic): anisotropic unoccupied DOS  -> CASTEP+OptaDOS (analyze_elnes.py).
  Channel B (dynamical): convergent-probe channelling + thickness + atom depth -> THIS file.
abtem's transition potentials are isolated-atom (all-electron radial wavefunctions), so this
forward model carries Channel B but NOT Channel A's fine structure. Comparing the two routes
therefore reveals WHICH channel dominates the measurable along-beam signal (the user's ask).

For a handful of unit cells with accuracy as the goal, EXACT abtem multislice is the reference
(BiP-PRISM / scatterem, Pelz 2026 arXiv:2607.00756, is a speed/memory approximation TO this
exact multislice -- it only earns its keep at scale, so it adds no accuracy on unit cells).

Uses the same abtem engine as sim/simulate_4dstem.py. The core-loss backend needs `gpaw`
(all-electron atomic wavefunctions); the elastic machinery + --selftest run without it.

    ~/hyperspy-bundle/bin/python simulate_eels.py --selftest              # elastic path (no gpaw)
    ~/hyperspy-bundle/bin/python simulate_eels.py --edge O_K --thickness-series   # needs gpaw
"""
from __future__ import annotations

import argparse
import numpy as np
import abtem

import config as C
import build_cells as B

# EELS edge -> (Z, n, l) of the ionised subshell (K=1s, L2,3=2p, M4,5=3d).
EDGE_NL = {"O_K": (8, 1, 0), "Ti_L23": (22, 2, 1), "Pb_M": (82, 3, 2)}

# Forward-model geometry. NB: the pipeline's 100 mrad convergence needs a fine grid/GPU; the
# default here is a modest alpha for feasibility on CPU. Set --alpha 100 for the real run.
DEFAULTS = dict(sampling=0.10, slice_thickness=1.0, alpha_mrad=20.0, lateral=(2, 2))


def build_slab(cellname: str, nz: int, lateral=(2, 2)):
    """Repeat a unit cell along the beam (z) to a target thickness; pad + centre for abtem."""
    slab = B.make_cell(C.CELLS[cellname]).repeat((lateral[0], lateral[1], nz))
    slab.center()
    return slab


def _probe_and_potential(slab, alpha_mrad, sampling, slice_thickness):
    pot = abtem.Potential(slab, sampling=sampling, slice_thickness=slice_thickness,
                          projection="infinite")
    probe = abtem.Probe(energy=C.OPTICS.energy_keV * 1e3, semiangle_cutoff=alpha_mrad)
    probe.grid.match(pot)
    return probe, pot


def _on_column_scan(slab):
    """Single probe on the Pb/Ti/O column at the cell centre (cheapest diagnostic scan)."""
    return abtem.CustomScan(np.array([[slab.cell[0, 0] / 2, slab.cell[1, 1] / 2]]))


def eels_signal(cellname: str, edge: str, nz: int, *, alpha_mrad=DEFAULTS["alpha_mrad"],
                sampling=DEFAULTS["sampling"], slice_thickness=DEFAULTS["slice_thickness"],
                lateral=DEFAULTS["lateral"], scan=None) -> float:
    """Energy-filtered core-loss signal for one cell/edge/thickness (Channel B). Needs gpaw."""
    Z, n, l = EDGE_NL[edge]
    slab = build_slab(cellname, nz, lateral)
    probe, pot = _probe_and_potential(slab, alpha_mrad, sampling, slice_thickness)
    from abtem.inelastic.core_loss import SubshellTransitions          # gpaw import happens in build()
    tp = SubshellTransitions(Z=Z, n=n, l=l).get_transition_potentials(
        extent=pot.extent, gpts=pot.gpts, energy=probe.energy)
    det = abtem.FlexibleAnnularDetector()
    scan = scan if scan is not None else _on_column_scan(slab)
    m = probe.transition_potential_scan(potential=pot, transition_potentials=tp,
                                        scan=scan, detectors=det)
    m = m.compute() if hasattr(m, "compute") else m
    return float(np.asarray(m.array).sum())


def thickness_series(edge: str, cells=("tet_Pz", "tet_Px", "cubic"),
                     nz_list=(2, 4, 6, 8), **kw) -> None:
    """Channel-B along-beam contrast vs thickness: tet_Pz (P along beam) vs tet_Px (P perp)."""
    _require_gpaw()
    print(f"== M6b dynamical EELS ({edge}, {C.OPTICS.energy_keV:.0f} keV, alpha "
          f"{kw.get('alpha_mrad', DEFAULTS['alpha_mrad']):.0f} mrad) ==")
    print(f"  {'thick(A)':>9} " + " ".join(f"{c:>12}" for c in cells) + "   Pz-vs-Px")
    for nz in nz_list:
        vals = {c: eels_signal(c, edge, nz, **kw) for c in cells}
        t = build_slab(cells[0], nz, kw.get("lateral", DEFAULTS["lateral"])).cell.lengths()[2]
        pz, px = vals["tet_Pz"], vals["tet_Px"]
        contrast = 200 * (pz - px) / (pz + px) if (pz + px) else 0.0
        print(f"  {t:>9.1f} " + " ".join(f"{vals[c]:>12.5g}" for c in cells) +
              f"   {contrast:>+7.2f}%")
    print("  -> nonzero Pz-vs-Px = a DYNAMICAL along-beam signature (Channel B). Compare its "
          "size to the CASTEP dichroism (Channel A) from analyze_elnes.py.")


def _require_gpaw() -> None:
    try:
        import gpaw  # noqa: F401
    except ImportError:
        raise SystemExit(
            "gpaw not installed -> abtem's core-loss atomic wavefunctions are unavailable.\n"
            "  It is the accurate backend (all-electron), so it fits the accuracy goal.\n"
            "  Install into the bundle:  ~/hyperspy-bundle/bin/pip install gpaw && "
            "~/hyperspy-bundle/bin/gpaw install-data ~/hyperspy-bundle/share/gpaw-setups\n"
            "  (macOS may need libxc + a C toolchain; if it fights back, defer to scatterem's\n"
            "   hydrogenic backend on release.) Meanwhile: simulate_eels.py --selftest.")


def selftest() -> None:
    """Validate the elastic multislice path (channelling engine) WITHOUT gpaw: build cells ->
    potential -> probe -> scan -> finite pixelated output; and a HAADF proxy for the along-beam
    contrast (a Channel-B teaser that needs no core-loss backend)."""
    print("== simulate_eels self-test (elastic multislice path, no gpaw) ==")
    abtem.config.set({"device": "cpu"})

    def haadf(name, nz=6):
        slab = build_slab(name, nz)
        probe, pot = _probe_and_potential(slab, 20.0, 0.10, 1.0)
        det = abtem.PixelatedDetector(max_angle="valid")
        m = probe.scan(pot, scan=_on_column_scan(slab), detectors=det).compute()
        arr = np.asarray(m.array)
        assert np.isfinite(arr).all() and arr.sum() > 0, "elastic scan produced no signal"
        # HAADF proxy: high-angle fraction of the diffraction intensity
        ny, nx = arr.shape[-2:]
        yy, xx = np.mgrid[0:ny, 0:nx] - np.array([[ny / 2], [nx / 2]])[:, None]
        r = np.hypot(*np.fft.fftshift(np.stack([yy, xx]), axes=(1, 2)))
        return float(arr.reshape(-1, ny, nx)[0][r > 0.35 * min(ny, nx) / 2].sum())

    pz, px = haadf("tet_Pz"), haadf("tet_Px")
    print(f"  elastic HAADF-proxy on-column: tet_Pz={pz:.4g} tet_Px={px:.4g} "
          f"Δ={200*(pz-px)/(pz+px):+.2f}%")
    print("  elastic multislice path OK (channelling engine works). Core-loss EELS needs gpaw.")
    try:
        import gpaw  # noqa: F401
        print("  gpaw: PRESENT -> run `--edge O_K --thickness-series` for the real M6b.")
    except ImportError:
        print("  gpaw: MISSING -> install it (see --thickness-series message) to run M6b.")
    print("== self-test PASSED ==")


def main() -> None:
    ap = argparse.ArgumentParser(description="Dynamical multislice STEM-EELS forward model (M6b).")
    ap.add_argument("--selftest", action="store_true", help="elastic-path validation (no gpaw)")
    ap.add_argument("--edge", default="O_K", choices=list(EDGE_NL))
    ap.add_argument("--thickness-series", action="store_true", help="Pz vs Px vs thickness (gpaw)")
    ap.add_argument("--alpha", type=float, default=DEFAULTS["alpha_mrad"],
                    help="probe convergence semi-angle (mrad); use 100 to match the pipeline")
    ap.add_argument("--device", default="cpu", choices=["cpu", "gpu"])
    args = ap.parse_args()
    abtem.config.set({"device": args.device})
    if args.selftest or not args.thickness_series:
        selftest()
    else:
        thickness_series(args.edge, alpha_mrad=args.alpha)


if __name__ == "__main__":
    main()
