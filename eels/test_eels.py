#!/usr/bin/env python
"""@file test_eels.py
@brief Unit tests for everything that does NOT need CASTEP/OptaDOS -- run now to keep the
       pipeline honest while the HPC admin clears. Plain asserts (also pytest-discoverable):

    ~/hyperspy-bundle/bin/python test_eels.py         # runs all, exits nonzero on any failure

Covers: structure crystallography (M1), CASTEP .cell writer validity (ASE round-trip),
core-hole invariants, the M2(b) benchmark crystals, the M6 geometry model, the OptaDOS
spectrum parser, and a submit-readiness preflight (all template files present).
"""
from __future__ import annotations

import os
import tempfile
import warnings

import numpy as np
import ase.io

# ASE instantiates a Castep calculator to validate keywords on read; with no CASTEP binary it
# warns loudly but parses fine. Silence it so the test output is readable (CASTEP-validity of
# the .cell is still asserted via the round-trip below).
warnings.filterwarnings("ignore", module="ase")

import config as C
import build_cells as B
import analyze_elnes as A

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- M1: structures
def test_all_structures_validate():
    """Every built cell passes its crystallographic gate (0 total failures)."""
    built = B.cells_to_build("all")
    ladder = {}
    nfail = sum(len(B.validate(n, s, at, scan_record=ladder if n.startswith("scan_") else None))
                for n, (s, at) in built.items())
    assert nfail == 0, f"{nfail} structure checks failed"


def test_polar_axis_orientation():
    """tet_Pz puts the Ti displacement ALONG the beam; tet_Px puts it PERPENDICULAR."""
    dz = B.offcentering(B.make_cell(C.CELLS["tet_Pz"]), C.Z_TI, 6)[0]
    dx = B.offcentering(B.make_cell(C.CELLS["tet_Px"]), C.Z_TI, 6)[0]
    assert abs(dz[2]) > 0.9 * np.linalg.norm(dz), "tet_Pz should be along beam z"
    assert np.hypot(dx[0], dx[1]) > 0.9 * np.linalg.norm(dx), "tet_Px should be in-plane"


def test_scan_ladder_linear():
    """The M5 calibration ladder is linear + monotonic in the along-beam displacement."""
    d = [np.abs(B.offcentering(B.make_cell(
        C.CellSpec(f"s{s}", "tetragonal", "z", disp_scale=s)), C.Z_TI, 6)[0, 2])
        for s in C.SCAN_SCALES]
    assert np.all(np.diff(d) > 0), "ladder must be monotonic"
    assert np.corrcoef(C.SCAN_SCALES, d)[0, 1] > 0.9999, "ladder must be linear"


# ---------------------------------------------------------------- .cell writer validity
def test_cell_writer_ase_roundtrip():
    """The written CASTEP .cell parses back through ASE with matching geometry (proxy for
    CASTEP-validity)."""
    at = B.make_cell(C.CELLS["tet_Pz"])
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "x.cell")
        B.write_castep_cell(at, p, (4, 4, 4))
        back = ase.io.read(p, format="castep-cell")
    assert back.get_chemical_formula() == at.get_chemical_formula()
    assert np.allclose(sorted(back.cell.lengths()), sorted(at.cell.lengths()), atol=1e-6)
    assert np.allclose(np.sort(back.get_scaled_positions(), axis=0),
                       np.sort(at.get_scaled_positions(), axis=0), atol=1e-6)


def test_corehole_single_excited_atom():
    """Each core-hole .cell labels EXACTLY one atom X:exc, of the requested species."""
    spec = C.CELLS["tet_Pz"]
    sup = B.make_cell(spec).repeat(spec.supercell)
    for tag, site in B.corehole_sites(spec).items():
        ei = B.excited_index(sup, site)
        assert ei is not None, f"no {site} found"
        assert sup.get_array("site_label")[ei] == site
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, f"{tag}.cell")
            B.write_castep_cell(sup, p, (4, 4, 4), exc_index=ei)
            back = ase.io.read(p, format="castep-cell")
        cs = list(back.arrays["castep_custom_species"])
        assert sum(1 for x in cs if x.endswith(":exc")) == 1, f"{tag}: not exactly one :exc"


def test_corehole_species_pot():
    """Core-hole cells get a complete SPECIES_POT with the right core-hole occupancy suffix;
    a cell containing a species with no known OTFG string (Sr) falls back to None."""
    spec = C.CELLS["tio2"]
    sup = B.make_cell(spec).repeat(spec.supercell)
    ok = B.species_pot_block(sup, B.excited_index(sup, "O"), "O_K")
    assert ok and "O:exc" in ok and "{1s1}" in ok and "Ti " in ok      # O K = 1s hole
    ti = B.species_pot_block(sup, B.excited_index(sup, "Ti"), "Ti_L23")
    assert ti and "Ti:exc" in ti and "{2p5}" in ti                     # Ti L = 2p hole
    ps = C.CELLS["tet_Pz"]; psup = B.make_cell(ps).repeat(ps.supercell)  # Pb now known -> complete
    pb = B.species_pot_block(psup, B.excited_index(psup, "Pb"), "Pb_M")
    assert pb and "Pb:exc" in pb and "{3d9}" in pb and "Ti " in pb and "O " in pb
    sr = C.CELLS["srtio3"]; srsup = B.make_cell(sr).repeat(sr.supercell)  # Sr OTFG TBD -> None
    assert B.species_pot_block(srsup, B.excited_index(srsup, "O"), "O_K") is None


# ---------------------------------------------------------------- M2(b) benchmarks
def test_benchmarks_correct():
    """SrTiO3 (Pm-3m) and rutile TiO2 (P4_2/mnm) build with the right symmetry + lattice, and
    Ti is centrosymmetric (a good null reference for the method)."""
    for key in ("srtio3", "tio2"):
        spec = C.CELLS[key]
        at = B.make_cell(spec)
        assert not B.validate(key, spec, at), f"{key} failed its benchmark gate"
        assert np.linalg.norm(B.offcentering(at, C.Z_TI, 6), axis=1).mean() < REF_TOL, \
            f"{key} Ti should be centrosymmetric"


REF_TOL = C.REF.tol_delta_A


# ---------------------------------------------------------------- M6 geometry model
def test_geometry_magic_angle():
    """theta_E and the magic angle come out at their textbook values, and the surviving
    anisotropy is 1 at beta->0 and inverts past the magic angle."""
    tE = A.characteristic_angle_rad(A.EDGE_ONSET_eV["O_K"], 300.0)
    assert 0.9e-3 < tE < 1.3e-3, f"theta_E {tE} out of range for O K @300kV"
    ma = A.magic_angle_rad(tE)
    assert 3.5 < ma / tE < 4.4, f"magic angle {ma/tE:.2f} theta_E off ~3.97"
    assert abs(A.surviving_anisotropy(1e-6, tE) - 1.0) < 1e-3
    assert A.surviving_anisotropy(20 * tE, tE) < 0


def test_required_counts_monotonic():
    assert A.required_counts(0.2) < A.required_counts(0.02) < A.required_counts(0.002)
    assert A.required_counts(0.0) == float("inf")


# ---------------------------------------------------------------- OptaDOS parser
def test_optados_parser():
    """load_spectrum reads a two-column OptaDOS-style file, skipping # and ! comments."""
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "seed_core_edge.dat")
        with open(p, "w") as f:
            f.write("# OptaDOS core-loss\n! energy  intensity\n"
                    "530.0 0.0\n532.0 1.5\n534.0 0.7\n")
        e, s = A.load_spectrum(p)
    assert np.allclose(e, [530, 532, 534]) and np.allclose(s, [0.0, 1.5, 0.7])
    assert A.dichroism_metric(e, np.array([0.0, 1.0, -1.0])) > 0


def test_optados_core_section_parser():
    """load_optados_core selects the excited-atom (:exc) block only, returns the broadened col
    (matches the real OptaDOS <seed>_core_edge.dat format from the M2(b) TiO2 benchmark)."""
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "x_core_edge.dat")
        with open(p, "w") as f:
            f.write("# Core loss function\n")
            f.write(" #  O 1 K1\n 0.0 9.9 9.9\n 1.0 9.9 9.9\n")        # normal block (ignore)
            f.write(" #  O 1 K1 O:exc\n 0.0 1.0 2.0\n 1.0 3.0 4.0\n")  # excited block (physical)
            f.write(" #  O 2 K1\n 0.0 9.9 9.9\n")                      # another normal (ignore)
        e, y = A.load_optados_core(p)                         # default :exc, broadened (col3)
        assert list(e) == [0.0, 1.0] and list(y) == [2.0, 4.0]
        _, yr = A.load_optados_core(p, broadened=False)       # raw (col2)
        assert list(yr) == [1.0, 3.0]


# ---------------------------------------------------------------- submit-readiness preflight
def test_templates_present():
    """All HPC input templates exist so submission is turnkey once CASTEP lands."""
    need = ["templates/groundstate.param", "templates/geomopt.param", "templates/coreloss.param",
            "templates/species_pot_corehole.cellblock", "templates/coreloss_qc.odi",
            "templates/coreloss_qperp.odi", "submit_castep.slurm"]
    missing = [f for f in need if not os.path.exists(os.path.join(HERE, f))]
    assert not missing, f"missing templates: {missing}"


def test_q_swap_by_polar_axis():
    """The .odi files carry a LAB-frame q, so a cell with c||x has the two files' crystal-frame
    meanings exchanged. Getting this wrong turns the M4 invariance check into a false failure."""
    assert A.q_is_swapped("tet_Px_Oap"), "tet_Px has c||x -> its .qc file is q-perp-c"
    assert not A.q_is_swapped("tet_Pz_Oap"), "tet_Pz has c||z -> .qc really is q||c"
    assert not A.q_is_swapped("tet_Pz_Oeq"), "excited site must not change the q mapping"
    assert not A.q_is_swapped("cubic_Ti"), "cubic is isotropic; no swap"
    assert not A.q_is_swapped("scan_0.50_Oap"), "the scan ladder is uniformly c||z"
    assert not A.q_is_swapped("real_Oap"), "tilted P: the meaningful split stays the lab one"
    assert A.cell_of("tet_Px_Oeq").name == "tet_Px", "site suffix must be stripped"


def test_core_dat_glob_matches_slurm_naming():
    """run_coreloss.slurm writes `<seed>.qc_core_edge.dat` into a per-run subdirectory. The
    finder must match that underscore and descend one level (an earlier `<seed>.qc.*core*.dat`
    pattern could match neither)."""
    import shutil
    runs = os.path.join(HERE, "runs")
    d = os.path.join(runs, "_pytest_tmp")
    made_runs = not os.path.isdir(runs)
    os.makedirs(d, exist_ok=True)
    try:
        for tag in ("qc", "qperp"):
            with open(os.path.join(d, f"zz_probe.{tag}_core_edge.dat"), "w") as fh:
                fh.write("# O 1 K1 O:exc\n530.0 1.0 1.0\n531.0 2.0 2.0\n")
        assert A.find_core_dat("zz_probe", "qc"), "must find <seed>.qc_core_edge.dat in runs/*/"
        assert A.find_core_dat("zz_probe", "qperp")
        assert A.find_core_dat("zz_absent", "qc") is None
    finally:
        shutil.rmtree(d, ignore_errors=True)
        if made_runs:
            try: os.rmdir(runs)
            except OSError: pass


def test_dichroism_fractions_normalised():
    """dichroism_fractions returns FRACTIONS (unlike dichroism_metric, which is unnormalised),
    is symmetric in its two arguments, and reads zero for identical spectra."""
    e = np.linspace(525.0, 575.0, 400)
    base = np.exp(-0.5 * ((e - 535.0) / 3.0) ** 2)
    same = A.dichroism_fractions(e, base, base)
    assert same["peak"] == 0.0 and same["integral"] == 0.0, "null test must read exactly zero"
    other = base + 0.5 * np.exp(-0.5 * ((e - 545.0) / 3.0) ** 2)
    f1 = A.dichroism_fractions(e, base, other)
    f2 = A.dichroism_fractions(e, other, base)
    assert abs(f1["peak"] - f2["peak"]) < 1e-12, "must not depend on argument order"
    assert 0.0 < f1["peak"] < 2.0 and 0.0 < f1["integral"] < 2.0, "must read as a fraction"


# ---------------------------------------------------------------- runner
def _run():
    tests = sorted(k for k, v in globals().items() if k.startswith("test_") and callable(v))
    npass = 0
    for name in tests:
        try:
            globals()[name]()
            print(f"  PASS  {name}")
            npass += 1
        except Exception as e:
            print(f"  FAIL  {name}: {type(e).__name__}: {e}")
    print(f"== {npass}/{len(tests)} passed ==")
    return npass == len(tests)


if __name__ == "__main__":
    raise SystemExit(0 if _run() else 1)
