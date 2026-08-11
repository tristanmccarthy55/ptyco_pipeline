# STEM-EELS forward simulator — handover

How to run (and finish building) the full STEM-mode EELS simulator: a scanned probe through a
thick, phonon-averaged specimen producing **simultaneous HAADF + a background-included core-loss
EELS spectrum**, with the CASTEP/OptaDOS ELNES injected for the edge shape. Companions:
**HANDOVER.md** (Blythe ops + the CASTEP/OptaDOS recipe), **RESULTS.md** (the M4/M5 ELNES results
this injects), **LITERATURE.md**.

---

## 0. Status (2026-08-11)

| Piece | State |
|---|---|
| STEM-EELS acquisition config (`config.STEMEELS`) | ✅ done |
| Spectrum core: inject ELNES (aperture-avg over β) + power-law background | ✅ built, **validated on real M4 O-K data** (`--selftest`) |
| Reuse `sim/simulate_4dstem.py` for scattering (probe/potential/phonons/scan) | ✅ `simulate_stem_eels.probe_channelling` imports it — no duplication |
| EELS detector-hole geometry (β = det/2) + simultaneous HAADF | ✅ wired |
| abtem core-loss coupling (`transition_potential_scan`) | ⏳ needs **gpaw** in the abtem env + a Blythe run (untested) |
| Per-atom ELNES for the labyrinth vortex | ☐ design (the payoff — see §7) |
| Plasmon low-loss background (Phase 2) | ☐ (power-law is in now) |

**Immediate next step:** install gpaw (§2), then run `probe_channelling` on Blythe and wire its
output weights into `assemble_spectrum`.

---

## 1. Architecture — two codes, each doing what only it can

CASTEP/OptaDOS is a **bulk, kinematic** ELNES calculator (edge shape + q-anisotropy); it has no
STEM probe, thickness, channelling, phonons, or HAADF. abtem is a **dynamical multislice** code
(all of those) but its atomic transition potentials give a *smooth* edge, not the solid-state
ELNES. So the full spectrum combines them, in the **channelling (local) approximation**:

```
I(scan, E) = Σ_species [ channelling-summed core-loss coupling ]  ×  σ_species(E; β)  +  background(E)
             └─────────── abtem, reusing simulate_4dstem ──────┘    └── OptaDOS ──┘   └ power-law ┘
```

- **channelling coupling** = `simulate_4dstem` probe (overfocused 100 mrad) + Lobato potential +
  frozen-phonon TDS + multislice, with the EELS **detector-hole**: the inner half of the pixelated
  detector radius (β = `DETECTOR_MAX_ANGLE_MRAD/2` = **100 mrad**) → the spectrometer; the outer
  `[100, 200]` mrad annulus → **simultaneous HAADF**.
- **σ_species(E; β)** = the OptaDOS q-resolved ELNES (M4/M5), **aperture-averaged over β** with the
  magic-angle model (`analyze_elnes.parallel_weight`), placed at the tabulated edge onset.
- **background** = power-law A·E^-r now; abtem plasmon low-loss later (Phase 2).

**gpaw's only role:** abtem's core-loss transition potential (the *atomic* wavefunction that sets
which atoms the probe ionises — the EELS *spatial* map). It is NOT the spectroscopy (CASTEP does
that). A gpaw-free alternative exists (sample the elastic probe intensity at atoms by hand) but is
more custom code; with Blythe + a small gpaw install, `transition_potential_scan` is the clean route.

---

## 2. Set up gpaw in the abtem env  ← DO THIS

gpaw must live in the SAME env abtem imports it from: `$SHARE/phucrh/envs/abtem` (per
`sim/run_sim.slurm`). It needs compiled deps (libxc, BLAS) best handled by conda-forge. Add it to
the existing env, protecting the validated pip abtem/cupy with `--freeze-installed`:

Run these **one line at a time** (a scrambled multi-line paste is what breaks it). Two gotchas
baked in below: `unset PYTHONPATH` (a stray PYTHONPATH gives `conda: No module named 'conda'`),
and `conda shell.bash hook` (robust init — don't rely on `conda info --base`, chicken-and-egg).

```bash
# Blythe login node -- ONE LINE AT A TIME
unset PYTHONPATH                                   # stray PYTHONPATH breaks conda's own python
module purge && module load Miniconda3/24.7.1-0
eval "$(conda shell.bash hook)"
conda --version                                    # MUST print a version before continuing

export CONDA_PKGS_DIRS="$SHARE/phucrh/conda/pkgs"
# install BY PREFIX (-p) -> no activation needed; --freeze-installed protects pip abtem/cupy:
conda install -y -p "$SHARE/phucrh/envs/abtem" -c conda-forge --freeze-installed gpaw libxc
#  ^ if the solver insists on changing abtem/cupy/numpy, Ctrl-C -> use the dedicated env below.

# PAW datasets + verify (the last line runs the exact atomic solver abtem calls):
"$SHARE/phucrh/envs/abtem/bin/gpaw" install-data "$SHARE/phucrh/envs/abtem/share/gpaw-setups"
export GPAW_SETUP_PATH="$(ls -d $SHARE/phucrh/envs/abtem/share/gpaw-setups/gpaw-setups-* | head -1)"
"$SHARE/phucrh/envs/abtem/bin/python" -c "import abtem, gpaw, ase; print('abtem', abtem.__version__, 'gpaw', gpaw.__version__)"
"$SHARE/phucrh/envs/abtem/bin/python" -c "from gpaw.atom.all_electron import AllElectron; AllElectron('O').run(); print('gpaw atomic O OK')"
```

**Fallback — a dedicated env** (use only if the `--freeze-installed` solve fails; keeps the
ptychography env pristine but re-installs abtem):
```bash
conda create -p "$SHARE/phucrh/envs/stemeels" -c conda-forge python=3.11 \
    abtem gpaw libxc ase dask h5py scipy numpy matplotlib -y
# GPU: add cupy matching CUDA 12 ->  "$SHARE/phucrh/envs/stemeels/bin/pip" install cupy-cuda12x
```
Add `export GPAW_SETUP_PATH=...` to `sim/run_sim.slurm` (near `module load CUDA/12.6.0`) so batch
jobs find the datasets. gpaw installs can need a round of iteration — capture any error and debug.

---

## 3. Files

| File | Role |
|---|---|
| `config.py` → `STEMEELS` dataclass | the **parameter home**: energy, convergence α, collection β (=det/2), HAADF inner/outer, thickness, slice, sampling, scan step/window, frozen-phonon count + per-species B, e-loss axis (min/max/dispersion), edges + onsets + `elnes_source`, background model/exponent/fraction, dose, `inelastic_model`, device |
| `simulate_stem_eels.py` | **spectrum core** (gpaw-free, tested): `aperture_averaged_elnes`, `powerlaw_background`, `assemble_spectrum`, `eloss_axis`; **STEM dynamics** (Blythe): `probe_channelling` (imports `simulate_4dstem`, EELS hole + HAADF), `eels_collection_mrad` |
| `sim/simulate_4dstem.py` | REUSED (not copied) for `load_and_prepare_atoms`, `build_potential`, `build_probe`, frozen-phonon scan, `make_scan`, `DETECTOR_MAX_ANGLE_MRAD` |

Inject ELNES from either the raw OptaDOS `_core_edge.dat` (`:exc` block auto-selected) or a
pre-extracted `.exc.txt` (`_load_elnes` handles both). `elnes_source[edge] = {"qc": ..., "qperp": ...}`.

---

## 4. The EELS-hole geometry + an honest caveat

The simultaneous ptychography+EELS setup (ONE pixelated detector, hole at half radius) gives
**β = 100 mrad**. That is *enormous* for EELS — far past the O-K magic angle (~4 mrad). Per the M6
model the along-beam anisotropy **surviving fraction at 100 mrad is ≈ −0.33**: the σ/π dichroism
(the whole along-beam-polarisation signal from M4/M5) is **largely averaged out** at this
collection. So this geometry yields a clean O-K edge + HAADF but is *anisotropy-washing*. To
actually measure the along-beam polarisation you'd want a dedicated **small** EELS aperture
(≲2 mrad); the simulator takes β as a parameter, so both can be modelled. This is a real output,
not a bug — flag it in any writeup of the combined-detector design.

---

## 5. How to run

```bash
# spectrum core, no abtem/gpaw (works now, on the Mac or Blythe):
$SHARE/phucrh/envs/abtem/bin/python simulate_stem_eels.py --selftest    # full O-K spectrum vs β

# full sim incl. abtem channelling + HAADF (Blythe GPU, AFTER gpaw is installed):
#   run via sim/run_sim.sh conventions (GPU partition, CUDA/12.6.0); simulate_stem_eels imports
#   simulate_4dstem, so the same env + GPU apply. Wire a run_stem_eels.slurm mirroring run_sim.slurm.
$SHARE/phucrh/envs/abtem/bin/python simulate_stem_eels.py --run --cell tet_Pz
```
(`--run` calls `probe_channelling`; finish wiring its `weights` into `assemble_spectrum` and add a
`run_stem_eels.slurm` GPU submit script modelled on `sim/run_sim.slurm`.)

---

## 6. Open items

1. **gpaw install** (§2) → unblocks `probe_channelling`.
2. **Blythe test** of `probe_channelling` (HAADF annulus + `transition_potential_scan` within β);
   confirm the frozen-phonon wrap (mirror `simulate_4dstem.run_scan_binned`'s one-config-at-a-time
   loop to bound memory).
3. **Wire weights → spectrum**: `assemble_spectrum(eloss_axis(cfg), weights, cfg)` per scan pixel →
   the spectrum-image; save alongside the HAADF (+ optional Poisson via `dose_e_per_A2`).
4. **`run_stem_eels.slurm`** GPU submit script (copy `sim/run_sim.slurm`, add `GPAW_SETUP_PATH`).
5. **Phase-2 background**: abtem plasmon low-loss + multiple-scattering convolution.
6. **Per-atom ELNES** (§7) — the physics-faithful vortex sim.

## 7. The payoff — a physics-faithful vortex STEM-EELS

In the labyrinth each O sees a *local* polar-axis orientation (the vortex). M4/M5 showed the O-K
ELNES depends strongly on that orientation vs the beam (78% dichroism; ~19% is the off-centering).
So a faithful vortex STEM-EELS injects a **per-atom** σ(E) keyed to each atom's local displacement
(magnitude + angle-to-beam, exactly the quantities `analysis/figures/pol_vortex.py` / `build_cells`
already compute). Then the simulated EELS map/spectrum-image encodes the 3-D polarisation — closing
the loop from the whole ELNES study back onto the ptychography sample. Start with a single injected
ELNES (uniform), then interpolate σ(E; angle) from the M5 tilt series per atom.
