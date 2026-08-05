# PbTiO₃ core-loss EELS: does along-beam polarisation show up?

**Question.** STEM ptychography of the PTO/STO vortices recovers only the *in-plane*
ferroelectric polarisation (the projected Ti off-centering, `analysis/figures/pol_vortex.py`).
The component *along the beam* is invisible to projected imaging. Can core-loss EELS (ELNES
fine structure) recover it, giving the full 3-D polarisation of a vortex?

**Observable.** The along-beam signal is an **ELNES dichroism**: the spectrum for momentum
transfer **q ∥ polar axis c** differs from **q ⊥ c**, because the uniaxial polar distortion
makes the projected unoccupied DOS anisotropic. We compute this with **CASTEP** (core-hole
ground state) + **OptaDOS** (`core_geom : polarized`, `core_qdir` = the q knob).

### Two facts that frame every result
1. **Magnitude yes, sign no.** In the dipole limit ELNES sees the *magnitude* of along-beam
   polarisation but not its *sign* (|q·r|² is even in q). Distinguishing +P/−P needs
   beyond-dipole / off-axis large-q terms — a flagged stretch goal.
2. **Collection angle matters more than anything.** The intrinsic dichroism is an upper
   bound; a real EELS aperture averages over q-directions and the signal vanishes at the
   **magic angle ≈ 4·θ_E** (≈4.3 mrad for O K at 300 keV). Your 100 mrad ptychography probe
   sits far past it — a core deliverable is *what collection aperture preserves the signal*.

Orientation matches `sim/simulate_4dstem.py` (`rotate(-90°, y)` → **beam = +z**), so unit-cell
results map onto the real acquisition: "P along beam" ≡ polar axis **c ∥ z**.

---

## Milestone staircase (each gate must pass before the next)

| # | What | Where | Pass gate | State |
|---|------|-------|-----------|-------|
| M0 | CASTEP + OptaDOS available, core-hole OTFG generates | Blythe | binaries resolve; OTFG report shows reduced core occ. | ☐ HPC |
| M1 | build + validate the test cells | local | space group, c/a≈1.06, \|δ_Ti\|≈0.3 Å, linear scan | ✅ **passed** |
| M2a | DFT convergence + ferroelectric double-well | Blythe | E converged <1 meV/atom; tetragonal lower; c/a right | ☐ HPC |
| M2b | benchmark ELNES vs TiO₂/SrTiO₃ | Blythe | onset+peaks match to ~1–2 eV; locks OTFG + final-state | ☐ HPC |
| M3 | **cubic null test** | Blythe | q∥z spectrum == q⊥ spectrum (must, by symmetry) | ☐ HPC |
| M4 | tetragonal dichroism, all edges | Blythe | nonzero symmetry-consistent Δ(E); tet_Pz q⊥ == tet_Px q∥ | ☐ HPC |
| M5 | calibration: Δ vs \|P\| (scan series) | Blythe+local | monotonic Δ-metric vs δ_Ti | ☐ HPC |
| M6 | detectability (analytic magic-angle model) | local | surviving-dichroism vs β, SNR, aperture | ✅ **model validated** |
| M6b | dynamical forward sim (abtem multislice EELS) | local | tet_Pz vs tet_Px map/qEELS vs thickness; thin-limit = M6 | ☐ abtem |
| M7 | scale-up: BiP-PRISM full vortex + CASTEP→TP coupling | local/GPU | full-labyrinth EELS map; intrinsic anisotropy through the scope | ☐ stretch |

M3 (null test) + M4 (rotational-invariance cross-check) are the load-bearing proofs that any
reported anisotropy is physics, not numerics.

**Two channels (why there are two forward-model milestones).** *Channel A — spectroscopic*:
the anisotropic unoccupied DOS (CASTEP+OptaDOS, M2–M5) — "does the atom's spectrum change".
*Channel B — dynamical/channeling*: thickness + convergent-probe channeling + atom depth
(multislice EELS, M6b) — "what the microscope records". M6's analytic magic-angle model is a
kinematic stand-in for B; M6b (abtem `inelastic.core_loss`, the same engine as `sim/`) is the
real thing. abtem's isolated-atom transition potentials capture B but not A's fine structure;
carrying A through the microscope needs CASTEP→transition-potential injection (M7). BiP-PRISM /
scatterem (Pelz 2026, arXiv:2607.00756, code on publication) is the scale path for the full
19,440-atom labyrinth vortex EELS map — comparable to your existing ptychography sim of it.

---

## Files

| File | Role |
|------|------|
| `config.py` | reference crystallography, structure presets, edges, q-dirs, DFT + optics knobs |
| `build_cells.py` | M1: build/validate cells → CASTEP `.cell` (+ `--corehole` variants) |
| `analyze_elnes.py` | M4–M6: dichroism, calibration, geometry-averaging, SNR (`--selftest`) |
| `simulate_eels.py` | M6b: dynamical abtem multislice STEM-EELS forward model (Channel B); `--selftest` runs without gpaw, core-loss needs gpaw |
| `test_eels.py` | unit tests for everything not needing CASTEP (`python test_eels.py`, 10/10) |
| `templates/groundstate.param`, `geomopt.param` | M2a SCF / relaxation |
| `templates/coreloss.param` | M3+ core-hole spectral task |
| `templates/species_pot_corehole.cellblock` | the core-hole OTFG block (**lock at M2**) |
| `templates/coreloss_qc.odi`, `coreloss_qperp.odi` | the two q-directions to difference |
| `submit_castep.slurm` | Blythe job: CASTEP once, OptaDOS twice (both q) |
| `structures/*.cell` | generated geometries (`_1cell` for M2a, `_222` for core hole) |
| `structures/corehole/*.cell` | one excited site per file (O_ap/O_eq/Ti/Pb labelled `X:exc`) |
| `RESULTS.md` | milestone-by-milestone evidence log |

## Running it

```bash
# M1 — build + validate all cells (writes structures/, exits nonzero on any failed gate)
~/hyperspy-bundle/bin/python build_cells.py --corehole

# M6 — validate the detectability model now (synthetic; reproduces the magic angle)
~/hyperspy-bundle/bin/python analyze_elnes.py --selftest

# --- on Blythe (M0/M2+), per (structure, excited site) ---
# append templates/species_pot_corehole.cellblock to each structures/corehole/<seed>.cell,
# copy templates/coreloss_q{c,perp}.odi next to it, then:
sbatch --export=ALL,SEED=tet_Pz_Oap submit_castep.slurm

# --- back on the Mac, once <seed>.qc.* / <seed>.qperp.* are in eels/runs/ ---
~/hyperspy-bundle/bin/python analyze_elnes.py --seed tet_Pz_Oap --edge O_K
```

## Honest caveats
- **Ti L₂,₃ absolute lineshape** is multiplet/spin-orbit-limited in single-particle DFT — use
  it for the *anisotropy* of the e_g/t₂g DOS, lean on **O K** for the quantitative result.
- **Core-hole final-state** (full vs excited-electron) shifts onsets; the *dichroism* (a
  difference) is largely insensitive — locked by the M2 benchmark.
- **Real vortex geometry is unfavourable**: viewed side-on, the labyrinth's P is mostly
  in-plane (median ~14% along beam, see RESULTS) — the along-beam signal is intrinsically
  small in this zone axis. M6 quantifies what's achievable; a different zone axis would help.
