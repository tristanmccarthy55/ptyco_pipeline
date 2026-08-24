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

# MEASURED final-state core-level shift between the two inequivalent oxygen sites, BE(apical) -
# BE(equatorial), in eV. Both core-hole runs use the SAME tet_Pz 2x2x2 supercell and therefore
# share a ground state, so the difference of their Delta-SCF total energies IS the shift:
#   E(hole on apical)     = -66672.67284 eV
#   E(hole on equatorial) = -66672.84345 eV   ->   +0.1706 eV
# OptaDOS writes each spectrum on an E_F-referenced axis and does NOT include the core binding
# energy, so combining the two sites into the full O K edge requires displacing the equatorial
# spectrum by -O_SITE_SHIFT_eV onto the apical frame. Each site's OWN dichroism is a within-run
# difference and is unaffected. (Measured 2026-08-24; previously assumed zero, which left the
# weighted edge uncertain over 22-37%.)
O_SITE_SHIFT_eV = 0.1706

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
# CASTEP 24.1 default on-the-fly (OTFG) pseudopotential strings, read from the M0 smoke-run
# pseudopotential report (runs/smoke/smoke.castep). Used to build COMPLETE core-hole SPECIES_POT
# blocks. Fill Pb/Sr from a run that includes them before doing PbTiO3/SrTiO3 core holes.
OTFG = {
    "O":  "2|1.1|17|20|23|20:21(qc=8)",
    "Ti": "3|1.8|9|10|11|30U:40:31:32(qc=5.5)",
    "Pb": "3|2.2|7|8|9|50U:60:51U:61:52",       # 5s5p5d6s6p semicore (22 e = ionic charge)
    # "Sr": "...",   # from an SrTiO3 run, if an SrTiO3 core hole is ever needed
}
# Core-hole occupancy suffix per edge: reduce the ionised shell by one electron. Appended to the
# excited species' OTFG string (Mizoguchi/NaGe convention: e.g. `...(qc=6){1s1}` = 1s K-hole).
COREHOLE_OCC = {"O_K": "{1s1}", "Ti_L23": "{2p5}", "Pb_M": "{3d9}", "Sr_L23": "{2p5}"}
# Map an excited-site tag (from corehole_sites) to its edge, hence its occupancy suffix.
SITE_EDGE = {"O": "O_K", "Oap": "O_K", "Oeq": "O_K", "Ti": "Ti_L23", "Pb": "Pb_M"}


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

# ---------------------------------------------------------------- STEM-EELS acquisition config
# The parameter home for the full forward simulator (simulate_stem_eels.py): a scanned STEM
# probe through a thick, phonon-averaged specimen producing SIMULTANEOUS HAADF + a full core-loss
# EELS spectrum (edges with the injected CASTEP/OptaDOS ELNES + a background). Optics conventions
# match sim/simulate_4dstem.py; phonon B values are the room-temperature isotropic Debye-Waller
# factors used there. NB the probe CONVERGENCE (alpha) and the EELS COLLECTION (beta) are separate
# angles: alpha forms the probe, beta is the spectrometer entrance aperture that sets how much of
# the anisotropy survives (magic angle ~4*theta_E; see analyze_elnes.py).
@dataclass(frozen=True)
class STEMEELS:
    # ---- beam / probe ----
    energy_keV: float = 300.0
    convergence_mrad: float = 20.0          # probe semi-angle (alpha). EELS: ~20-30 (NOT the 100 ptycho probe)
    defocus_A: float = 0.0
    # ---- EELS spectrometer collection (round aperture) ----
    # DETECTOR-HOLE geometry (simultaneous ptychography + EELS on ONE pixelated detector, per the
    # sim/simulate_4dstem.py setup): the spectrometer takes the inner HALF of the detector radius,
    # so beta = DETECTOR_MAX_ANGLE_MRAD / 2 = 100 mrad. simulate_stem_eels.eels_collection_mrad()
    # derives it from sim; this default mirrors it. (beta drives the OptaDOS q-anisotropy averaging,
    # M6: 100 mrad is far past the magic angle -> the along-beam dichroism is largely averaged out.)
    collection_mrad: float = 100.0          # = DETECTOR_MAX_ANGLE_MRAD / 2 (the EELS hole)
    # ---- simultaneous HAADF = the OUTER annulus [beta, detector_max] of the same scan ----
    haadf_inner_mrad: float = 100.0         # = beta (the hole edge)
    haadf_outer_mrad: float = 200.0         # = DETECTOR_MAX_ANGLE_MRAD
    # ---- specimen / multislice ----
    thickness_nm: float = 20.0              # specimen thickness along the beam (repeats the cell in z)
    slice_thickness_A: float = 2.0
    sampling_A: float = 0.05                # real-space potential sampling (band-limits the detector)
    # ---- scan ----
    scan_step_A: float = 0.2
    scan_window_A: float | None = None      # None -> full in-plane cell; else a square window (A)
    scan_center_A: tuple | None = None      # (x,y) window centre; None -> cell centre
    # ---- frozen phonons (thermal diffuse scattering) ----
    n_phonons: int = 8                      # 0 = coherent; 8-16 for realistic TDS
    phonon_seed: int = 1
    phonon_B: dict = field(default_factory=lambda: {"Pb": 0.90, "Sr": 0.55, "Ti": 0.45, "O": 0.80})  # A^2, RT
    # ---- energy-loss axis (the spectrum) ----
    eloss_min_eV: float = 500.0
    eloss_max_eV: float = 600.0
    eloss_dispersion_eV: float = 0.1        # eV per channel
    # ---- edges: onset (eV) + the OptaDOS ELNES source (a .exc.txt or processed shape) ----
    # The injected shape is the q-resolved OptaDOS ELNES, aperture-averaged over `collection_mrad`.
    edges: tuple = ("O_K",)                 # keys into EDGES; each contributes an edge to the spectrum
    edge_onset_eV: dict = field(default_factory=lambda: {"O_K": 532.0, "Ti_L23": 456.0, "Pb_M": 2484.0})
    elnes_source: dict = field(default_factory=dict)  # {"O_K": "runs/exc/tet_Pz_Oap.{q}.exc.txt"} injected shapes
    # ---- background ----
    background_model: str = "powerlaw"      # "powerlaw" (A*E^-r) now; "plasmon" (abtem low-loss) later; "none"
    background_r: float = 3.0               # power-law exponent
    background_frac: float = 0.5            # background height at the first edge onset, as a fraction of the edge jump
    # ---- dose / noise (optional Poisson) ----
    dose_e_per_A2: float | None = None      # None -> noiseless; else Poisson at this dose
    # ---- rigour / device ----
    inelastic_model: str = "channelling"    # "channelling" (elastic multislice x sigma, NO gpaw) |
    #                                         "transition_potential" (abtem rigorous, NEEDS gpaw)
    device: str = "gpu"                     # Blythe abtem env has cupy+CUDA


STEMEELS_CFG = STEMEELS()

# ---------------------------------------------------------------- output
OUT_DIR = os.path.join(os.path.dirname(__file__), "structures")   # .cell/.param land here


def preset(name: str) -> CellSpec:
    if name not in CELLS:
        raise KeyError(f"unknown cell preset {name!r}; have {list(CELLS)} (+ 'all', 'scan')")
    return CELLS[name]
