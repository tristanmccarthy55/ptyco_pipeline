#!/usr/bin/env python
"""@file config.py
@brief Central configuration for the PbTiO3 core-loss EELS (ELNES) polarisation study.

One place for every knob: the reference crystallography, the structure presets, the DFT
convergence/core-hole choices, and the momentum-transfer (q) directions used for the ELNES
anisotropy test. The physics question: does the ferroelectric polarisation COMPONENT ALONG
THE BEAM produce a distinguishable ELNES signal (a q||c vs q-perp dichroism) that projected
STEM imaging is blind to? See ../../.claude plan or eels/README.md for the milestone map.

Beam/orientation convention is inherited from the 4D-STEM sim (sim/simulate_4dstem.py,
ROTATE_DEG_Y=-90): the beam is Cartesian +z, so a polar axis c||z means "P along beam" and
c-perp-z (c||x) means "P in-plane / perpendicular to beam". Same viewing geometry as the
labyrinth PTO/STO example, so unit-cell results map onto the real acquisition.

Needs the hyperspy-bundle Python (ase + pymatgen + spglib):
    ~/hyperspy-bundle/bin/python build_cells.py --preset all
"""
from __future__ import annotations

from dataclasses import dataclass, field
import os

# ---------------------------------------------------------------- paths
# Repo root is two levels up (eels/ is peer to sim/ and analysis/). The labyrinth POSCAR is
# reused as the source of REALISTIC vortex displacements (see real_cell in build_cells.py).
_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
_VASP_CANDIDATES = [
    os.path.join(_REPO, "sim", "PTO6_STO6_18_18_labyrinthPoscar.vasp"),
    os.path.join(_REPO, "..", "PTO6_STO6_18_18_labyrinthPoscar.vasp"),
]


def find_vasp() -> str:
    for p in _VASP_CANDIDATES:
        if os.path.exists(p):
            return os.path.abspath(p)
    return _VASP_CANDIDATES[0]


# ---------------------------------------------------------------- reference crystallography
# Room-temperature tetragonal PbTiO3, space group P4mm (#99), origin on Pb (1a). z-fractional
# basis is a standard experimental refinement; the DERIVED off-centering |delta_Ti| ~ 0.3 A
# (Ti minus its O6 centroid) and c/a ~ 1.06 are the quantities build_cells.py VALIDATES, so a
# wrong z-set is caught rather than trusted. Cubic Pm-3m (#221) is the paraelectric reference.
@dataclass(frozen=True)
class RefCrystal:
    # lattice
    a_cubic: float = 3.970          # A, paraelectric Pm-3m
    a_tet: float = 3.904            # A
    c_tet: float = 4.152            # A  -> c/a ~ 1.0635
    # tetragonal P4mm fractional z (x,y fixed by symmetry). Pb at origin (z=0).
    z_Ti: float = 0.5377            # Ti  1b (1/2,1/2,z)
    z_O_ap: float = 0.1174          # apical O  1b (1/2,1/2,z)   -- the O along the polar axis
    z_O_eq: float = 0.6174          # equatorial O 2c (1/2,0,z) & (0,1/2,z)
    # validation targets (literature); build_cells.py asserts against these with tolerance
    target_c_over_a: float = 1.0635
    target_delta_Ti_A: float = 0.30    # |Ti - O6 centroid|, the polarisation proxy
    target_delta_Pb_A: float = 0.47    # |Pb - O12 / A-site cage centroid| (looser)
    tol_lattice_frac: float = 0.02     # 2% on a, c, c/a
    tol_delta_A: float = 0.08          # A, on the derived off-centering magnitudes


REF = RefCrystal()

# ---------------------------------------------------------------- M2(b) benchmark crystals
# Known, well-characterised ELNES references used to lock the core-hole method (OTFG string +
# final-state treatment) BEFORE trusting the PTO dichroism. SrTiO3 is the closest analogue
# (perovskite titanate, cubic, single O site); rutile TiO2 is the classic O-K/Ti-L benchmark.
BENCH = {
    "srtio3": dict(spacegroup=221, sg_symbol="Pm-3m", cellpar=(3.905, 3.905, 3.905, 90, 90, 90),
                   symbols=("Sr", "Ti", "O"), basis=((0, 0, 0), (0.5, 0.5, 0.5), (0.5, 0.5, 0.0))),
    "tio2_rutile": dict(spacegroup=136, sg_symbol="P4_2/mnm", cellpar=(4.593, 4.593, 2.959, 90, 90, 90),
                        symbols=("Ti", "O"), basis=((0, 0, 0), (0.305, 0.305, 0.0))),
}

# Atomic numbers (for the labyrinth extraction + species bookkeeping)
Z_PB, Z_TI, Z_O, Z_SR = 82, 22, 8, 38

# ---------------------------------------------------------------- edges (core-loss targets)
# Ranked by how defensible the single-particle DFT/OptaDOS result is. 'core' names the shell
# that carries the hole; 'note' is the honesty caveat surfaced in the writeup.
EDGES = {
    "O_K":    dict(element="O",  core="1s", note="primary: reliable + strongly polar-sensitive via Ti-O hybridisation"),
    "Ti_L23": dict(element="Ti", core="2p", note="secondary: use for e_g/t2g ANISOTROPY, not absolute lineshape (multiplet/SO missing)"),
    "Pb_M":   dict(element="Pb", core="3d", note="exploratory: 6s lone-pair covalency; deep core, needs care"),
}
# For O K in tetragonal PTO there are TWO inequivalent sites; the measured edge is their
# multiplicity-weighted sum (1 apical + 2 equatorial per formula unit).
O_SITE_MULTIPLICITY = {"O_ap": 1, "O_eq": 2}

# ---------------------------------------------------------------- momentum-transfer directions
# The ELNES anisotropy knob (OptaDOS core_qdir). Beam = +z. q||z == "P along beam"; q perp z
# (x) == "P in-plane". Fractional/Cartesian handling is verified against OptaDOS at M3/M4.
Q_BEAM = (0.0, 0.0, 1.0)      # q parallel to beam
Q_PERP = (1.0, 0.0, 0.0)      # q perpendicular to beam (in-plane)


# ---------------------------------------------------------------- structure presets
@dataclass
class CellSpec:
    """One structure to build + validate + emit as CASTEP .cell."""
    name: str
    kind: str                       # "cubic" | "tetragonal" | "scan" | "real"
    polar_axis: str = "z"           # "z" (P along beam) | "x" (P in-plane) | "off" (tilted)
    disp_scale: float = 1.0         # 0..1 fraction of the full ferroelectric displacement
    tilt_deg: float | None = None   # for kind=="real": polar axis angle from beam z
    supercell: tuple = (2, 2, 2)    # for the core-hole runs (M3+); (1,1,1) also emitted for M2a
    desc: str = ""


# The staircase's structures (see plan M1). All share the sim beam=z convention.
CELLS: dict[str, CellSpec] = {
    "cubic":   CellSpec("cubic", "cubic", polar_axis="z", disp_scale=0.0,
                        desc="Pm-3m paraelectric reference -> M3 null test (q||z must equal q-perp)"),
    "tet_Pz":  CellSpec("tet_Pz", "tetragonal", polar_axis="z", disp_scale=1.0,
                        desc="P4mm, polar c || beam z -> P ALONG BEAM (headline case, M4)"),
    "tet_Px":  CellSpec("tet_Px", "tetragonal", polar_axis="x", disp_scale=1.0,
                        desc="P4mm, polar c || x -> P PERP to beam (M4 rotational-invariance x-check)"),
    "real":    CellSpec("real", "real", polar_axis="off", disp_scale=1.0, tilt_deg=45.0,
                        desc="realistic tilted P from the labyrinth vortex (magnitude+angle) -> M6"),
    "srtio3":  CellSpec("srtio3", "srtio3",
                        desc="M2(b) benchmark: cubic perovskite titanate, known O-K/Ti-L ELNES"),
    "tio2":    CellSpec("tio2", "tio2_rutile",
                        desc="M2(b) benchmark: rutile TiO2, classic O-K/Ti-L ELNES reference"),
}
# Displacement-scaled series for the M5 calibration curve (dichroism vs |P|), all c||z.
SCAN_SCALES = (0.0, 0.25, 0.5, 0.75, 1.0)


# ---------------------------------------------------------------- DFT / core-hole choices
# SCHEMATIC values; the cutoff/k-points are LOCKED by the M2(a) convergence run and the
# core-hole treatment by the M2(b) benchmark against TiO2/SrTiO3. Do not treat as final.
@dataclass(frozen=True)
class DFT:
    xc: str = "PBE"
    cutoff_sweep_eV: tuple = (600, 700, 800, 900)     # M2a convergence
    cutoff_eV: int = 800                              # working value until M2a locks it
    kgrid_1cell_sweep: tuple = ((6, 6, 6), (8, 8, 8), (10, 10, 10))
    kgrid_1cell: tuple = (8, 8, 8)
    kgrid_2x2x2: tuple = (4, 4, 4)                    # scaled for the supercell
    nextra_bands: int = 100                           # empty bands to span the edge
    elec_energy_tol: float = 1e-8
    core_hole_charge: float = 1.0                     # full core hole + neutralising background
    corehole_label: str = "exc"                       # excited species suffix -> "O:exc" etc.


DFT_CFG = DFT()

# ---------------------------------------------------------------- experimental geometry (M6)
# Inherited from sim/simulate_4dstem.py so the detectability answer matches the real pipeline.
@dataclass(frozen=True)
class Optics:
    energy_keV: float = 300.0
    convergence_mrad: float = 100.0     # ptychography probe: HUGE for EELS -> averages anisotropy
    # collection semi-angle for a HYPOTHETICAL EELS acquisition; swept in analyze_elnes.py to
    # find where the dichroism survives vs the magic angle (~few mrad at 300 keV for O K).
    collection_sweep_mrad: tuple = (1, 2, 5, 10, 25, 50, 100)


OPTICS = Optics()

# ---------------------------------------------------------------- output
OUT_DIR = os.path.join(os.path.dirname(__file__), "structures")   # .cell/.param land here


def preset(name: str) -> CellSpec:
    if name not in CELLS:
        raise KeyError(f"unknown cell preset {name!r}; have {list(CELLS)} (+ 'all', 'scan')")
    return CELLS[name]
