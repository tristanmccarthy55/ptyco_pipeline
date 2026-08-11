# Handoff — M4 & M5 results (for the agent that built M0–M3)

You built and validated M0–M3 (setup, recipe, cubic Ti-L null). **M4 and M5 are now done** — this
is the delta. Full numbers in `RESULTS.md`; how-to-run in `HANDOVER.md`.

## M4 — orientation dichroism
`tet_Pz` O-K, q∥c vs q⊥c = **78% peak dichroism** (σ*/π*): beam∥P (q∥c) → **σ\*** (higher-E
features ~535 & ~540 eV); beam⊥P (q⊥c) → **π\*** (sharp ~530 eV peak).
Data: `runs/tetPz_OK/`; extracted `runs/exc/tet_Pz_Oap.{qc,qperp}_core_edge.exc.txt`.

## M5 — displacement scan + the key decomposition
Ran `scan_0.00 → 0.75` (+ `tet_Pz` = s=1), q∥c & q⊥c, all finished. Lattice held fixed (a, c) —
only the internal Ti/O displacement s varies, so the tetragonal **strain** contribution is
constant and the s-dependence is **purely the off-centering**.

| s | dichroism q∥c−q⊥ (peak) | along-beam (q∥c) change vs s=0 |
|---|---|---|
| 0.00 (zero-P) | **70.1%** | 0 |
| 0.25 | 70.2% | 1.2% |
| 0.50 | 70.5% | 4.3% |
| 0.75 | 71.4% | 10.7% |
| 1.00 (full P) | **78.2%** | **18.7%** |

**Result:** ~70 of the 78% is the **tetragonal backbone** — the O's Ti–O–Ti σ/π split, present at
*zero* polarization (so q∥c ≠ q⊥c even for the non-polar cell; that is *not* a bug — the O site
has a defined backbone axis). The ferroelectric off-centering's **own** signature is a
**monotonic ~19%** change in the along-beam-probed spectrum (and +8% to the dichroism). Reproduces
Bugnet et al. "O-K responds to Ti off-centering" from first principles.

**Honest headline:** the O-K edge is 78% dichroic to the polar-**axis** orientation vs the beam
(what maps a vortex, since P and the tetragonal axis co-rotate); of that, ~19% is the
polarization's own displacement effect. Dipole limit: axis + magnitude, **not sign**.

Data: `runs/scan_*_Oap/` + `runs/exc/scan_*.exc.txt`. PI figures + `SLIDE_NOTES.md` on the Mac at
`~/Desktop/eels_figs/` (fig1 validation, fig2 result+object+zero-P, fig3 decomposition,
fig3a along-beam scan).

## Method note that emerged running these
OptaDOS's `<seed>_core_edge.dat` holds **one block per atom/edge**; the physical ELNES is the
**excited atom's `:exc` block only** — select it, never sum (`analyze_elnes.load_optados_core`).
The broadened spectrum is column 3 (col 2 = raw). OptaDOS is **serial and slow** for PbTiO₃
(~90 min/pass — it recomputes every atom/edge, incl. Pb's many core states) → submit
`run_coreloss.slurm` with `-n 1`, `--time=12:00:00`, and **`SKIP_CASTEP=1`** to reuse an existing
`.elnes_bin` when only OptaDOS needs re-running.

## Left to do (yours to pick up)
1. **M4 cross-checks** — `tet_Px_Oap` (rotational invariance: `tet_Px` q∥c must equal `tet_Pz`
   q⊥c) and `tet_Pz_Oeq` (equatorial O → full O-K = 1×apical + 2×equatorial). Submitted with the
   M5 batch — **check status** (`squeue -u $USER`; `ls runs/tetPx_OK runs/tetPz_Oeq`) and resubmit
   if not COMPLETED. Cells already in `structures/corehole/`.
2. **M6 detectability** — fold the intrinsic dichroism through the 300 keV collection geometry
   (magic angle ≈ 4·θ_E ≈ 4.3 mrad; the 100 mrad ptychography probe averages it → a small EELS
   collection aperture is needed). Model already in `analyze_elnes.py`.
3. **Extend the scan** to `s = 1.25, 1.50` — the labyrinth vortex reaches Ti off-centering ≈ 0.40 Å
   ≈ **s ≈ 1.2**, above bulk (enhanced/"supertetragonal" regions). Add to `config.SCAN_SCALES`,
   rerun `build_cells.py --corehole`, submit. Effect accelerates → bigger signal there.
