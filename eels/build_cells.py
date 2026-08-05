#!/usr/bin/env python
"""@file build_cells.py
@brief M1: build the PbTiO3 test cells (sim beam=z orientation) and VALIDATE them before DFT.

Emits CASTEP .cell files for the ELNES polarisation staircase and, crucially, checks each
structure is crystallographically correct (space group, c/a, and the Ti off-centering that
IS the polarisation proxy) so a bad structure is caught here, not 200 core-hole SCF hours
later. Pass criteria (from config.REF) are asserted and printed as a report.

Cells (config.CELLS + the scan series):
  cubic   Pm-3m paraelectric      -> M3 null test (q||z must equal q-perp by symmetry)
  tet_Pz  P4mm, polar c || beam z  -> P ALONG BEAM (headline, M4)
  tet_Px  P4mm, polar c || x       -> P PERP to beam (M4 rotational-invariance cross-check)
  scan_*  c||z, displacement scaled -> M5 calibration curve (dichroism vs |P|)

Orientation matches sim/simulate_4dstem.py (beam=+z). "P along beam" == polar axis c||z.

    ~/hyperspy-bundle/bin/python build_cells.py            # build + validate all, write .cell
    ~/hyperspy-bundle/bin/python build_cells.py --preset tet_Pz
    ~/hyperspy-bundle/bin/python build_cells.py --check-only   # validate, don't write files
"""
from __future__ import annotations

import argparse
import dataclasses
import os

import numpy as np
from scipy.spatial import cKDTree
from ase import Atoms

import config as C
from config import REF, DFT_CFG, Z_TI, Z_PB, Z_O

# Site labels carried through supercelling so the core-hole writer can pick a specific
# crystallographic O (apical vs equatorial), Ti, or Pb to excite. Order matches build below.
_SITE_LABELS = np.array(["Pb", "Ti", "O_ap", "O_eq", "O_eq"], dtype="<U5")


# ---------------------------------------------------------------- structure construction
def build_perovskite(a: float, c: float, s: float) -> Atoms:
    """5-atom ABO3 perovskite, Pb at origin, polar axis along z.

    s in [0,1] scales the ferroelectric internal displacements from the centrosymmetric
    (s=0) toward the full experimental P4mm (s=1). At a==c and s==0 this is cubic Pm-3m.
    Apical O sits at (1/2,1/2,z); the two equatorial O at (1/2,0,z) & (0,1/2,z).
    """
    z_ti = 0.5 + s * (REF.z_Ti - 0.5)
    z_oap = 0.0 + s * (REF.z_O_ap - 0.0)
    z_oeq = 0.5 + s * (REF.z_O_eq - 0.5)
    symbols = ["Pb", "Ti", "O", "O", "O"]
    scaled = [
        (0.0, 0.0, 0.0),        # Pb  (A-site, origin)
        (0.5, 0.5, z_ti),       # Ti  (B-site)
        (0.5, 0.5, z_oap),      # O apical  (along polar axis)
        (0.5, 0.0, z_oeq),      # O equatorial
        (0.0, 0.5, z_oeq),      # O equatorial
    ]
    at = Atoms(symbols=symbols, scaled_positions=scaled,
               cell=[[a, 0, 0], [0, a, 0], [0, 0, c]], pbc=True)
    at.set_array("site_label", _SITE_LABELS.copy())
    return at


def orient(atoms: Atoms, polar_axis: str, tilt_deg: float | None = None) -> Atoms:
    """Rotate so the polar axis (built along z) points where we want vs the beam (+z)."""
    at = atoms.copy()
    if polar_axis == "z":
        return at                                   # P along beam
    if polar_axis == "x":
        at.rotate(90, "y", rotate_cell=True)        # P in-plane, perpendicular to beam
        return at
    if polar_axis == "off":
        at.rotate(float(tilt_deg), "y", rotate_cell=True)   # P tilted from beam by tilt_deg
        return at
    raise ValueError(f"bad polar_axis {polar_axis!r}")


def build_benchmark(key: str) -> Atoms:
    """M2(b) reference crystal from its space group + Wyckoff basis (ase.spacegroup.crystal).
    Sets site_label = element symbol so the core-hole writer can pick a Ti or O to excite."""
    from ase.spacegroup import crystal
    b = C.BENCH[key]
    at = crystal(symbols=list(b["symbols"]), basis=list(b["basis"]),
                 spacegroup=b["spacegroup"], cellpar=list(b["cellpar"]))
    at.set_array("site_label", np.array(at.get_chemical_symbols(), dtype="<U5"))
    return at


def make_cell(spec: C.CellSpec) -> Atoms:
    """Realise a CellSpec into an oriented ASE Atoms (unit cell, before supercelling)."""
    if spec.kind in C.BENCH:
        return build_benchmark(spec.kind)                # benchmarks: no polar orientation
    if spec.kind == "cubic":
        base = build_perovskite(REF.a_cubic, REF.a_cubic, 0.0)
    else:
        # tetragonal lattice held fixed; disp_scale isolates the polar displacement so the
        # M5 scan separates crystal-field tetragonality (s=0) from ferroelectricity (s>0).
        base = build_perovskite(REF.a_tet, REF.c_tet, spec.disp_scale)
    return orient(base, spec.polar_axis, spec.tilt_deg)


# ---------------------------------------------------------------- validation (the M1 gate)
def offcentering(atoms: Atoms, center_Z: int, nn: int) -> np.ndarray:
    """|center - centroid(nn nearest O)| per center atom, with full PBC. The polarisation
    proxy from analysis/figures/pol_vortex.py (delta = r_center - O-cage centroid)."""
    Z = atoms.get_atomic_numbers()
    P = atoms.get_positions()
    cell = atoms.cell.array
    ctr = P[Z == center_Z]
    ox = P[Z == Z_O]
    shifts = np.array([i * cell[0] + j * cell[1] + k * cell[2]
                       for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    ox_t = (ox[None] + shifts[:, None]).reshape(-1, 3)
    _, idx = cKDTree(ox_t).query(ctr, k=nn)
    centroid = ox_t[idx].mean(1)
    return ctr - centroid                              # vector(s), Cartesian (beam = z)


def spacegroup(atoms: Atoms, symprec: float = 1e-3) -> str:
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    struct = AseAtomsAdaptor.get_structure(atoms)
    return SpacegroupAnalyzer(struct, symprec=symprec).get_space_group_symbol()


def validate(name: str, spec: C.CellSpec, atoms: Atoms,
             scan_record: dict | None = None) -> list[str]:
    """Return a list of FAILED checks (empty = pass). Prints a per-cell report line.

    Strict |dTi|~0.30 A is enforced only at FULL displacement (disp_scale==1); the s=0 cell
    must be centrosymmetric (|dTi|~0); intermediate scan cells are recorded for a post-loop
    linearity check rather than compared to the full target (they are meant to be scaled)."""
    fails: list[str] = []
    a, b, c = atoms.cell.lengths()
    sg = spacegroup(atoms)

    dti = offcentering(atoms, Z_TI, 6)                 # (n,3) Ti off-centering vector(s)
    dti_mag = float(np.linalg.norm(dti, axis=1).mean())
    has_pb = bool((atoms.get_atomic_numbers() == Z_PB).any())
    dpb_mag = float(np.linalg.norm(offcentering(atoms, Z_PB, 12), axis=1).mean()) if has_pb else float("nan")

    # projection of the Ti polar displacement onto the beam (+z) vs in-plane
    along = float(abs(dti[0, 2]))
    perp = float(np.hypot(dti[0, 0], dti[0, 1]))

    # ---- checks
    if spec.kind in C.BENCH:                            # M2(b) reference crystals
        want = C.BENCH[spec.kind]
        if sg != want["sg_symbol"] and not sg.startswith(want["sg_symbol"].split("/")[0]):
            fails.append(f"spacegroup {sg} != {want['sg_symbol']}")
        for got, exp, lab in zip((a, b, c), want["cellpar"][:3], "abc"):
            if abs(got - exp) > REF.tol_lattice_frac * exp:
                fails.append(f"{lab}={got:.3f} off ref {exp}")
    elif spec.kind == "cubic":
        if not sg.startswith("Pm"):
            fails.append(f"spacegroup {sg} != Pm-3m")
        if dti_mag > REF.tol_delta_A:
            fails.append(f"cubic |dTi|={dti_mag:.3f} should be ~0")
    elif spec.kind in ("tetragonal", "real"):
        c_over_a = max(a, b, c) / min(a, b, c)
        if abs(c_over_a - REF.target_c_over_a) > REF.tol_lattice_frac * REF.target_c_over_a:
            fails.append(f"c/a={c_over_a:.4f} off target {REF.target_c_over_a}")
        if spec.disp_scale == 0:
            if dti_mag > REF.tol_delta_A:
                fails.append(f"s=0 |dTi|={dti_mag:.3f} should be ~0 (centrosymmetric)")
        else:
            if not sg.startswith("P4mm"):
                fails.append(f"spacegroup {sg} != P4mm")
            if spec.disp_scale == 1 and abs(dti_mag - REF.target_delta_Ti_A) > REF.tol_delta_A:
                fails.append(f"|dTi|={dti_mag:.3f} off target {REF.target_delta_Ti_A} A")
        # orientation gate: the polar displacement must sit on the intended axis vs the beam
        if spec.disp_scale > 0:
            if spec.polar_axis == "z" and along < 0.9 * dti_mag:
                fails.append(f"c||z but only {along:.3f}/{dti_mag:.3f} A is along beam")
            if spec.polar_axis == "x" and perp < 0.9 * dti_mag:
                fails.append(f"c||x but only {perp:.3f}/{dti_mag:.3f} A is perpendicular")
        if scan_record is not None:
            scan_record[spec.disp_scale] = along        # ladder for linearity check

    status = "OK  " if not fails else "FAIL"
    pb = f"|dPb|={dpb_mag:.3f}A" if has_pb else "(no Pb)"
    print(f"  [{status}] {name:10s} sg={sg:8s} a={a:.3f} c={c:.3f} c/a={max(a,b,c)/min(a,b,c):.4f} "
          f"|dTi|={dti_mag:.3f}A (beam {along:.3f} / perp {perp:.3f})  {pb}")
    for f in fails:
        print(f"          -> {f}")
    return fails


# ---------------------------------------------------------------- realistic vortex tilt
def extract_real_tilt(vasp_path: str, min_delta_A: float = 0.15) -> tuple[float, float]:
    """Measure the along-beam polarisation distribution from the labyrinth PTO/STO vortex.

    Beam = raw-POSCAR x (the sim does rotate(-90 deg, y): prepared z <- raw x). For every
    polar (PTO-like, |delta|>min_delta) Ti we take the angle of its off-centering delta from
    the beam. Returns (representative tilt in deg, representative |delta| in A) so real_cell
    reproduces a genuine vortex column's along-beam component (cos tilt). Approximate: uses
    the pure rotation, not the residual orthogonalise axis-relabel (a robust-statistic use)."""
    import ase.io
    big = ase.io.read(vasp_path)
    dti = offcentering(big, Z_TI, 6)                    # (Nti,3), RAW frame
    mag = np.linalg.norm(dti, axis=1)
    sel = mag > min_delta_A                             # polar Ti = PbTiO3 regions
    along = np.abs(dti[sel, 0])                         # |delta . beam|, beam = raw x
    tilt = np.degrees(np.arccos(np.clip(along / mag[sel], 0, 1)))   # angle from beam [0,90]
    rep_tilt, rep_mag = float(np.median(tilt)), float(np.median(mag[sel]))
    print(f"  labyrinth: {sel.sum()}/{len(mag)} polar Ti (|d|>{min_delta_A}A) | "
          f"|d| med {rep_mag:.3f} max {mag.max():.3f}A | angle-from-beam med {rep_tilt:.0f} "
          f"deg (IQR {np.percentile(tilt,25):.0f}-{np.percentile(tilt,75):.0f}) | "
          f"along-beam frac med {np.median(np.cos(np.radians(tilt))):.2f}")
    return rep_tilt, rep_mag


# ---------------------------------------------------------------- core-hole site selection
def _cart_center(atoms: Atoms) -> np.ndarray:
    return np.array([0.5, 0.5, 0.5]) @ atoms.cell.array


def excited_index(atoms: Atoms, site_value: str) -> int | None:
    """Index of the atom with the given site_label closest to the supercell centre (the atom
    whose periodic core-hole images are furthest away -> least self-interaction)."""
    labels = atoms.get_array("site_label")
    idx = np.where(labels == site_value)[0]
    if len(idx) == 0:
        return None
    d = np.linalg.norm(atoms.get_positions()[idx] - _cart_center(atoms), axis=1)
    return int(idx[int(np.argmin(d))])


def corehole_sites(spec: C.CellSpec) -> dict[str, str]:
    """Excited sites to emit per structure (tag -> site_label). Cubic O are all equivalent."""
    if spec.kind in C.BENCH:                            # benchmark: Ti + O (site_label=element)
        return {"Ti": "Ti", "O": "O"}
    if spec.kind == "cubic":
        return {"Ti": "Ti", "O": "O_ap", "Pb": "Pb"}
    return {"Ti": "Ti", "Oap": "O_ap", "Oeq": "O_eq", "Pb": "Pb"}


# ---------------------------------------------------------------- CASTEP .cell writer
def write_castep_cell(atoms: Atoms, path: str, kgrid: tuple, exc_index: int | None = None) -> None:
    """Minimal, explicit CASTEP .cell. If exc_index is set, that atom gets the core-hole
    species label (e.g. 'O:exc') and SYMMETRY_GENERATE is omitted (the hole breaks symmetry)."""
    cell = atoms.cell.array
    symbols = atoms.get_chemical_symbols()
    frac = atoms.get_scaled_positions()
    L = ["%BLOCK LATTICE_CART", "  ANG"]
    for v in cell:
        L.append(f"  {v[0]:18.12f} {v[1]:18.12f} {v[2]:18.12f}")
    L += ["%ENDBLOCK LATTICE_CART", "", "%BLOCK POSITIONS_FRAC"]
    for i, (sym, f) in enumerate(zip(symbols, frac)):
        label = f"{sym}:{DFT_CFG.corehole_label}" if i == exc_index else sym
        L.append(f"  {label:8s} {f[0]:18.12f} {f[1]:18.12f} {f[2]:18.12f}")
    L += ["%ENDBLOCK POSITIONS_FRAC", "",
          f"KPOINT_MP_GRID {kgrid[0]} {kgrid[1]} {kgrid[2]}"]
    if exc_index is None:
        L.append("SYMMETRY_GENERATE")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


# ---------------------------------------------------------------- driver
def cells_to_build(preset: str) -> dict[str, tuple[C.CellSpec, Atoms]]:
    out: dict[str, tuple[C.CellSpec, Atoms]] = {}
    names = list(C.CELLS) if preset in ("all", "scan") else [preset]
    if preset != "scan":
        for n in names:
            spec = C.CELLS[n]
            if spec.kind == "real":
                try:                            # ground the vortex tilt in the real labyrinth
                    tilt, _ = extract_real_tilt(C.find_vasp())
                    spec = dataclasses.replace(spec, tilt_deg=tilt)
                except Exception as e:          # missing/odd .vasp -> fall back to the 45 deg default
                    print(f"  [warn] real_cell tilt from .vasp failed ({e}); using {spec.tilt_deg} deg")
            out[n] = (spec, make_cell(spec))
    if preset in ("all", "scan"):
        for s in C.SCAN_SCALES:
            spec = C.CellSpec(f"scan_{s:.2f}", "tetragonal", polar_axis="z", disp_scale=s,
                              desc=f"M5 calibration point, disp_scale={s}")
            out[spec.name] = (spec, make_cell(spec))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Build + validate PbTiO3 ELNES test cells (M1).")
    ap.add_argument("--preset", default="all",
                    help="all | scan | " + " | ".join(C.CELLS))
    ap.add_argument("--check-only", action="store_true", help="validate but do not write files")
    ap.add_argument("--corehole", action="store_true",
                    help="also emit per-site core-hole .cell files (excited atom labelled X:exc)")
    args = ap.parse_args()

    built = cells_to_build(args.preset)
    print(f"== M1 build + validation ({len(built)} cells, beam=+z) ==")
    n_fail = 0
    scan_ladder: dict[float, float] = {}
    for name, (spec, atoms) in built.items():
        fails = validate(name, spec, atoms, scan_record=scan_ladder if name.startswith("scan_") else None)
        n_fail += len(fails)
        if not args.check_only:
            super_atoms = atoms.repeat(spec.supercell)
            write_castep_cell(atoms, os.path.join(C.OUT_DIR, f"{name}_1cell.cell"), DFT_CFG.kgrid_1cell)
            write_castep_cell(super_atoms, os.path.join(C.OUT_DIR, f"{name}_222.cell"), DFT_CFG.kgrid_2x2x2)
            if args.corehole and (spec.kind in ("cubic", "tetragonal", "real") or spec.kind in C.BENCH):
                for tag, site in corehole_sites(spec).items():
                    ei = excited_index(super_atoms, site)
                    if ei is not None:
                        write_castep_cell(super_atoms,
                                          os.path.join(C.OUT_DIR, "corehole", f"{name}_{tag}.cell"),
                                          DFT_CFG.kgrid_2x2x2, exc_index=ei)

    # scan ladder must be a clean linear calibration axis (dichroism vs |P| needs this at M5)
    if len(scan_ladder) >= 3:
        s = np.array(sorted(scan_ladder))
        d = np.array([scan_ladder[k] for k in s])
        r = np.corrcoef(s, d)[0, 1]
        ok = r > 0.9999 and d[0] < REF.tol_delta_A and np.all(np.diff(d) > 0)
        print(f"  [{'OK  ' if ok else 'FAIL'}] scan ladder along-beam |dTi| = "
              f"{np.array2string(d, precision=3)} A  (linearity r={r:.5f}, monotonic)")
        if not ok:
            n_fail += 1

    print(f"== {'ALL PASSED' if n_fail == 0 else str(n_fail) + ' CHECK(S) FAILED'} ==")
    if not args.check_only:
        print(f"wrote .cell files -> {C.OUT_DIR}")
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
