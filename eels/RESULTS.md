# RESULTS — PbTiO₃ ELNES polarisation study

Milestone-by-milestone evidence. Append to this as HPC gates clear. (See README.md for the
full staircase and pass criteria.)

---

## M1 — structure build + validation — ✅ PASSED

`~/hyperspy-bundle/bin/python build_cells.py --corehole` (beam = +z). All 9 cells passed
every crystallographic gate; 16 base `.cell` + 35 core-hole `.cell` written.

| cell | space group | c/a | \|δ_Ti\| (Å) | along-beam / perp | \|δ_Pb\| (Å) |
|------|-------------|-----|-------------|-------------------|-------------|
| cubic | Pm-3m | 1.0000 | 0.000 | 0.000 / 0.000 | 0.000 |
| tet_Pz | P4mm | 1.0635 | 0.331 | **0.331 / 0.000** | 0.487 |
| tet_Px | P4mm | 1.0635 | 0.331 | 0.000 / 0.331 | 0.487 |
| real | P4mm | 1.0635 | 0.331 | 0.046 / 0.328 | 0.487 |
| scan_0.00 | P4/mmm | 1.0635 | 0.000 | — | 0.000 |
| scan_0.25→1.00 | P4mm | 1.0635 | 0.083→0.331 | along z | 0.122→0.487 |

- c/a = 1.0635 (target 1.0635 ✓); \|δ_Ti\| = 0.331 Å (target ~0.30, tol 0.08 ✓); \|δ_Pb\| =
  0.487 Å (lit. ~0.47 ✓).
- Orientation correct: `tet_Pz` puts **all** of δ_Ti along the beam, `tet_Px` **all**
  perpendicular → the two limiting cases for the M4 dichroism.
- Scan ladder along-beam δ_Ti = [0, 0.083, 0.165, 0.248, 0.331] Å, **linear r = 1.00000**,
  monotonic → clean M5 calibration axis.
- Core-hole cells: exactly one atom per file relabelled `X:exc`, chosen nearest the 2×2×2
  centre (max separation from its periodic images).

### Early scientific finding (from the real labyrinth, sim orientation)
Over **2560 / 3888** polar Ti (\|δ\|>0.15 Å) in `PTO6_STO6_18_18_labyrinthPoscar.vasp`:
median \|δ\| 0.227 Å (max 0.397); **angle-from-beam median 82°** (IQR 73–87°) → **median
along-beam fraction only ~0.14**. Viewed side-on (the pipeline's zone axis) the vortex
polarisation is *mostly in-plane*; the along-beam component we're chasing is intrinsically
small here. `real` cell built at the median 82° tilt captures this (0.046 Å along beam).
**Implication:** `tet_Pz` (full along-beam) is the upper bound; expect the real measurable
signal to be a small fraction of it. A different zone axis would raise the along-beam share.

---

## M6 — detectability geometry model — ✅ VALIDATED (local, synthetic)

`~/hyperspy-bundle/bin/python analyze_elnes.py --selftest`. The collection-aperture
averaging model is validated before any real spectrum:

- **θ_E(O K, 300 keV) = 1.09 mrad**; **magic angle = 4.32 mrad = 3.98 θ_E** (textbook ≈3.97 ✓).
- Surviving intrinsic dichroism vs collection semi-angle β:

  | β (mrad) | 1 | 2 | 5 | 10 | 25 | 100 |
  |---|---|---|---|---|---|---|
  | fraction surviving | 0.62 | 0.28 | −0.04 | −0.17 | −0.26 | **−0.33** |

- **Your 100 mrad ptychography probe → −0.33**: the anisotropy is averaged away / inverted.
  To measure along-beam polarisation, use a *small* EELS collection aperture, ideally
  ≲ 2 mrad (≪ the 4.3 mrad magic angle).
- SNR: a 2% fractional anisotropy needs ~4.5×10⁴ counts/channel at SNR 3 — sets the dose target.

---

## M6b — dynamical forward model (abtem multislice EELS) — ⚙️ API validated, EELS backend gated on gpaw

`simulate_eels.py`. Confirmed the exact abtem 1.0.5 core-loss call path works end-to-end
(`SubshellTransitions(Z,n,l).get_transition_potentials` → `Probe.transition_potential_scan`
→ detectors), on the same abtem engine as `sim/simulate_4dstem.py`.
- **Elastic multislice path validated without gpaw** (`--selftest` passes): build cells →
  potential → probe → scan produces finite channelling output.
- **Core-loss EELS needs `gpaw`**: abtem 1.0.5's atomic radial wavefunctions are gpaw-only
  (no hydrogenic fallback). gpaw is the *accurate* (all-electron) backend, so it fits the
  accuracy goal — but it's a real dependency to add to the bundle:
  `~/hyperspy-bundle/bin/pip install gpaw && ~/hyperspy-bundle/bin/gpaw install-data ...`
  (macOS may need libxc + a C toolchain). Once present, `--edge O_K --thickness-series`
  gives the Channel-B along-beam contrast (tet_Pz vs tet_Px) vs thickness.
- **Route comparison (the ask):** abtem (isolated-atom potentials) = Channel B only; OptaDOS =
  Channel A only → comparing them shows which dominates the measurable along-beam signal.
  BiP-PRISM/scatterem is an *approximation to* exact abtem multislice, so on unit cells exact
  abtem is the accuracy reference and BiP-PRISM adds nothing until scale.

## Pre-flight (CASTEP-independent) — ✅ ALL GREEN (Mac + Blythe HPC)

Also validated **on Blythe** (SLURM job 1240612, COMPLETED/0:0, in the `abtem` env, Python 3.11,
NumPy 2.x, no spglib): 10/10 tests, ALL PASSED build, and the abtem-1.0.9 elastic multislice
self-test passed unchanged — HPC env + code + `.cell` generation + the channelling engine all
confirmed. (`run_pyeels.slurm`; graceful `sg=n/a` without spglib; NumPy-2 `trapezoid` shim.)

`~/hyperspy-bundle/bin/python test_eels.py` → **10/10 passed**. Covers structure crystallography,
**CASTEP `.cell` writer validity (ASE round-trips the files, incl. the `O:exc` custom species)**,
the exactly-one-excited-atom invariant, the benchmark crystals, the M6 geometry model, the
OptaDOS parser (synthetic), and a submit-readiness check that all HPC templates exist. Re-run
this any time the code changes.

**M2(b) benchmark inputs READY:** `srtio3` (Pm-3m, a=3.905) and `tio2` (rutile P4₂/mnm,
a=4.593 c=2.959) build + validate; core-hole cells `srtio3_{Ti,O}` / `tio2_{Ti,O}` emitted.
So the first CASTEP calculation has its structures waiting.

## M0, M2–M5 — pending (require CASTEP/OptaDOS on Blythe)

Next actions when CASTEP lands: M0 env check → M2(a) convergence + double-well → M2(b)
TiO₂/SrTiO₃ benchmark to lock the cutoff, k-points, OTFG core-hole string, and final-state
treatment. Structures + templates are prepared; only the version-specific keyword/OTFG values
remain to fill. Record numbers here as they land.
