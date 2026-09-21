# Next phase — relaxing the assumptions toward publication

Written 2026-09-16 at the close of the round-α campaign. Read `HANDOVER.md` first for the state
this builds on, and `PSF_KERNELS.md` for the kernel rules every step below must keep.

## Why this phase exists

The four-α result (`HANDOVER.md` § Results) rests on six idealisations, and each one is a referee's
objection:

| assumption now | what a real experiment has |
|---|---|
| probe known exactly | C3/C5 measured by the corrector, **defocus uncertain** |
| noiseless patterns | finite dose — Poisson counting noise |
| fully coherent source | chromatic focal spread + finite source size |
| static lattice | thermal vibration — frozen-phonon TDS |
| exactly on zone axis | a few degrees off axis |
| 5-cell (19.5 Å) slab | the full 18-cell (70 Å) labyrinth |

Relax them one at a time, measure what each costs, then run them all together. A result that
survives the combined configuration is publishable; one that survives only the idealised one is not.

## Ground rules for every step

1. **One relaxation at a time**, compared against the *previous* step as well as the original
   baseline. Otherwise two effects that cancel look like no effect.
2. **Probe α points: 70 and 90 mrad.** 70 is where depth first resolves and where the Ronchigram
   stops being flat; 90 is the best point. Only the final combined configuration re-runs the full
   50/70/90/100 sweep. This keeps each step to ~6 recons instead of ~12.
3. **Kernels stay byte-identical to their lab leg.** Every relaxation applied to the labyrinth is
   applied identically to its Pb/Ti grid legs — same dose, same phonons, same coherence, same tilt,
   same thickness, and **the same fitted probe** (not the true one). This rule is what made the
   kernels matched; breaking it confounds every comparison that follows.
4. **Engine constants do not move**: `REGLAYER=0` on every leg (it low-passes the depth axis), one
   `BETA_LSQ` for all legs (driver-pinned at 0.05).
5. **Scan field larger than the probe** (scan / d90 ≳ 1.5). Below ~1 the reconstruction has no
   positional diversity and returns speckle — that is what ended the sweep at 110 mrad.
6. **Everything on Blythe lands under `$SHARE/phucrh`** (`/springbrook/share/physics/phucrh/`) —
   never the group root `$SHARE`, which holds other users' folders. The drivers pack to
   `$SHARE/$USER` since `5af792b`; any new script, `sbatch --output`, tarball, scratch or sim
   directory must resolve there too. After each submission, check `ls /springbrook/share/physics/`
   shows nothing new of ours. (Until 2026-09-15 every tarball leaked into the group root.)
7. **Record the same numbers every time**, appended to one table,
   `aberration_experiment/results/<week>/relaxation_ladder.csv`: precision; bulk and whole-slab
   recall per species; xy-RMS; z-RMS; species confusion; false positives by species; and the
   recovered C1 whenever the probe is fitted. `analysis/relaxation_ladder.py` writes the row (it
   reuses `collate_atomfind_depth.load_run`; false positives by species are the confusion matrix's
   `<species>->none` entries). Step 0 is seeded in `results/2026-W38/`.

**Baseline to beat (step 0)** — matched kernels, corrected registration, noiseless, known probe:

| α | precision | bulk recall Pb / Ti / O | xy-RMS | z-RMS | confusion | false positives |
|---|---|---|---|---|---|---|
| 70 | 0.98 | 98 / 82 / 75% | 0.04 Å | 0.56 Å | 1.3% | 10 of 486 (O 6) |
| 90 | 0.94 | 95 / 94 / 95% | 0.03 Å | 0.37 Å | 0.0% | 35 of 581 (O 33) |

## The order, and why

Cheap and reconstruction-side first; then simulation-side knobs that already exist; then the ones
that need new code; the most expensive last. The probe fit leads because it is the assumption most
specific to the claim, it reuses the existing sims, and every later step reuses its probe model.

### Step 1 — Fix C3 and C5, fit C1 only

**Why constrained.** Free complex-probe fitting failed five times on this thin weak-phase slab
(grid junk → noise → NaN → NaN → presolve junk). A one-parameter family is far better posed, and
it is the operator's real situation: C3 and C5 come from the corrector's aberration tableau;
defocus is the knob that is never quite where you think.

**What exists.** `PROBE_START` for free probe updates (the route that failed). The PtychoShelves
GPU engine has **no parametric defocus refinement** — only `probe_fourier_shift_search`,
`estimate_NF_distance` and detector scale/rotation searches. The reconstruction h5 stores **no
error history** (only `reconstruction/object`, `probes`, two flags). *Correction (2026-09-17):* the
error is not in the slurm logs either — 0 lines in every packed log. It lived only inside MATLAB
(`out.error_metric`) and in the unpacked `Niter*.mat`.

**What to build.**
1. **Save the error trace** from `run_synthetic_recon_ML.m` (`fdb.score` / the LSQML
   `fourier_error`) into the output, so an objective can be read without parsing logs.
2. **A probe writer** — `make_probe.py` or a `--probe-only` mode of `sim/simulate_4dstem.py` —
   that writes `probe_initial.mat` for (α, C3, C5, C1, BIN) through the sim's own
   `build_initial_probe`, so the probe grid is guaranteed identical to the data's.
3. **An outer-loop search over C1**, probe fixed per trial: a coarse grid around the start point at
   reduced `NITER` (~50), golden-section refinement, then one full 200-iteration recon at the
   winner. Validate first that the objective actually has its minimum at the true C1 on a
   known case — don't assume it.

**Start point.** The planned C1 offset by a realistic focus error, ±20–40 Å, trying both signs.

**Accept if** the recovered C1 lands within the width of the objective's minimum of the truth
(define that width from the curve, not in advance), and the atomfind numbers sit within noise of
step 0.

**As built (2026-09-17)** — all three pieces exist and are tested locally; nothing submitted yet.
- `ptycho/run_synthetic_recon_ML.m` writes `<run>_error_trace.csv` and `<run>_layer_stats.csv` beside
  the h5 (after it, in `try/catch`). MATLAB is not available locally, so the smoke test is its test.
- `sim/make_probe.py` reproduces both packed true probes (overlap 1.000000, max relative difference
  2.0e-7 at a70 and 2.8e-6 at a90) and self-checks against the truth in every job at the true C1.
- `campaign/run_c1_search.sh` + `analysis/c1_objective.py` (usage in `campaign/README.md`).

**Two facts that shape the search.**
- *The probe changes fast with C1.* Overlap of the probe with itself shifted in C1: a70 0.90 / 0.48 /
  0.18 at 2 / 5 / 10 Å; a90 0.75 / 0.03 at 2 / 5 Å. A 10 Å grid alone can step over the minimum, so
  the grid is 10 Å coarse plus 2 Å fine across ±10 Å.
- *Defocus is degenerate with the object's depth.* A focus error and an equal depth shift of the object
  give the same far-field data, and the full-box recon has 4 Å of vacuum each side to absorb it. Expect a
  **flat bottom about ±4 Å wide**, with the reconstructed slab moving one-for-one with the trial C1
  inside it. So C1 is reported as the **centre of the flat bottom with its half-width**, not an argmin.
  `c1_objective.py` panel (c) tests the one-for-one shift directly.

**Estimator.** Bottom = the run of grid C1 around the minimum within 2σ of the bottom's level (two passes,
because the minimum of noisy points is biased low). σ comes from repeats at the same C1 (the GPU engine
reseeds its RNG from the clock, so repeats differ in batch order and layer starts); with fewer than 4
repeat degrees of freedom it takes the larger of that and the second-difference scatter. Monte Carlo on
the planned grid, truth off-grid: ±4 Å bottom → rms error 1.4–1.9 Å, coverage 0.82–0.85; ±2 Å → 0.8–1.1 Å,
0.92–0.93; a sharp minimum → 0.5–0.6 Å, 0.98–0.99. `--blind-start` replays an operator walking downhill in
10 Å steps from a ±20/±40 Å start, then scanning 2 Å inside the bracket.

**Run order** (agreed 2026-09-17): smoke test at a70 → a70 grid (25 trials at NITER 50, plus 7 at NITER
200 to check that 50 iterations rank the trials as 200 do) → review → a90 grid (h5s packed: the depth-shift
test needs them, see below) plus one
a90 known-C1 NITER 200 run, which also re-baselines a90 at `BETA_LSQ` 0.05 → final lab + Pb + Ti at the
fitted C1 with the same fitted probe → atomfind → `analysis/relaxation_ladder.py --step 1`.

**Smoke test (2026-09-17, job 1284049, a70 at the true C1, NITER 50) — green.** Writer self-check passed on
Blythe's abTEM 1.0.9; 1600 positions; no NaNs; 20 trace rows as the solver's schedule predicts; error falls
41.6 → 25.3 and is still falling ~4% per 10 iterations at 50 (hence the NITER 200 ranking check). Against the
200-iteration known-probe a70 recon on the same data and probe: whole-volume phase correlation 0.93 (0.90–0.96
in the atomic layers, high-passed, central 18 Å), atomic contrast 75–85% of the 200-iteration value, slab at
4.9–24.0 Å (true atomic band 4.0–23.5). Two findings: (1) the sidecar's raw per-layer phase std is dominated by
the phase-ramp gauge (vacuum layer 2: 0.095 raw, 0.025 deramped) and so biases the slab position —
`c1_objective.py` now measures the slab from the h5 (deramped, 1 Å high-pass) and falls back to the sidecar
only when no h5 was packed; (2) therefore every C1 grid packs its h5s.

**If it fails on the thin slab**, record that and retry after Step 6. The campaign notes flag the
thin weak-phase slab as intrinsically under-constraining the probe; a thicker sample may be what
makes it work.

### Stage 2.5 — probe update from the corrector's C3/C5 and an offset C1 (the experimentalist's route)

**Why.** An operator would not scan fixed probes: they build the start probe from the corrector tableau
(C3, C5) and their defocus estimate, then let the probe update refine it. Defocus is the usual free
parameter. We don't need a C1 number out — the refined probe can feed the kernels wholesale — only
convergence from a valid start.

**The earlier "failed five times" was mostly a bug (found 2026-09-17).** Every run with the TEM aperture
constraint (`PROBE_SUPPORT_FFT=1`, rounds v3–v5) NaN'd at iteration 1 of the full engine, including the
`[40 Inf]` runs that froze the probe there; every run without it finished. Cause: `load_from_p.m` builds the
aperture mask from `abs(fft2(probe_initial))` in unshifted order, and `rescale_inputs.m` resizes it for the
presolve with `crop_pad`, which keeps the array centre — pure high frequency. Ported and run on our probes, the
presolve mask was all zeros at a70 and a90, so the probe was erased at `probe_change_start`, and the full
engine's initial probe rescaling turned the zero probe into 0·Inf = NaN. Fixed by storing the mask centred and
shifting it back in `init_solver.m` (both lines inside the TEM branch; nothing else changes): the full-engine
mask is bit-identical to before and the presolve mask now keeps all of the probe's power. So a constrained
probe update has never actually run on this data; v1 (no constraint, blind 4 Å start) and v2 (a scalar
`probe_support_fft` bug) were the other two.

**Runs.** `campaign/run_c1_search.sh` with `PSTART` (presolve release) and optionally `PSTART2` (full engine):
starts at ±10% and ±20% of the true C1 and at the truth, both variants, plus fixed-probe controls at the same
starts, NITER 200, α = 70 and 90. `analysis/probe_refit_check.py` reports per run whether the probe survived,
its shift-invariant overlap with the true probe at the start and the end, the C1 it now looks like, and the
slab position.

**Result (2026-09-18, 24 runs).** The full-pixel update does not recover defocus: started at the truth it
stays (overlap 0.999); started 6–32 Å away it moves the probe by at most ~1 Å or destroys it (overlap
0.02–0.04). Four runs diverged after iteration 100 of the full engine. And it can *lower* the error while
the object collapses into the vacuum (a70, −6 Å start: 14.0 vs 23.7 fixed, slab at 1.9–4.7 Å).
`results/2026-W38/probe_refit.csv`.

### Stage 2.5b — a one-parameter (defocus-only) probe update inside the engine

The user's point stands: an operator pins C3/C5 and lets the solver fit C1 *with* the object. Built
2026-09-21 as an engine option, `p.probe_defocus_only` (`PROBE_DEFOCUS_ONLY=1` with `PROBE_START`;
driver `PDFO=1`, campaign `c1dfo`, dirs `_ps<P>[x<P2>]dfo`). In `engines.GPU.LSQML` the full-pixel probe
update direction is projected onto dP/dC1 of the current probe, the existing LSQ step machinery gives the
scalar step, and the step is applied as the exact Fresnel propagation `fft2(P)·exp(dz·iπλk²)` — checked
against abTEM to 1e-7 including the sign — so the aperture and C3/C5 cannot drift by construction. The
cumulative shift is printed every iteration (`defocus-only probe update: cumulative C1 shift`), carried
across the two engines, and written to the error-trace header (`probe_defocus_shift_A`);
`probe_refit_check.py` prints it beside the C1 it measures from the refined probe. Not testable without
MATLAB: its first run is its test. Launched at a70 (±6, ±12, 0) and a90 (±16, ±32, 0), NITER 200.

### Step 2 — Shot noise

**What exists.** `sim/add_poisson_noise.py --in-dir <sim>/.. --dose <e/Å²> [--seed N]` — a
post-process, **no re-simulation**. It writes `<sim_dir>_dose<tag>/01/data_dp.hdf5` with positions,
probe and metadata symlinked from the noiseless sim. Electrons per pattern = dose × step².
The PX915 report's four-dose series (`analysis/atomfind/dose_series.py`) is precedent.

**Plan.** A dose ladder at α = 70 and 90: 10⁷ (near-noiseless sanity check) → 10⁶ → 10⁵ → 10⁴ e/Å².

**Built (2026-09-17).** `DOSES="..."` on `campaign/run_thin_atomfind.sh`: per α, dose and leg, a CPU job
writes `sim_out_af_a<A>_<mode>_dose<D>` (independent seeds for lab / Pb / Ti; it also links the true probe
and `aberrations.json`, and records `poisson_noise.json`), and the known-probe recon runs on it after that
job succeeds. Same recon settings as step 0, so the 10⁷ rung doubles as the a90 `BETA_LSQ` 0.05 check.

**Watch.** Kernels get the same dose — 25-blob averaging now fights noise, so report kernel
peak/background at each dose. Check whether the recon's likelihood model suits Poisson data at
low dose. The `_dose<tag>` directory naming means `RECON_ONLY` must be pointed at the noisy copy.

### Step 3 — Partial coherence of the source

**Why here.** It changes the probe model that Step 1 fits, so validate it before the expensive steps.

**What exists.** On the reconstruction side, `PROBE_MODES` (mixed-state ptychography) and
`VARIABLE_PROBE`. On the simulation side, **nothing**.

**What to build (sim).**
- **Temporal first** (cheaper): an incoherent sum of diffraction patterns over a Gaussian spread of
  defocus, Δ = Cc·ΔE/E. For a cold FEG at 300 kV (ΔE ≈ 0.3–0.4 eV, Cc ≈ 1.2–1.7 mm) that is
  Δ ≈ 12–23 Å at 1σ; 5–7 quadrature points. Its effect grows with α, so expect 90 to feel it more.
- **Spatial** second: an incoherent sum over Gaussian-distributed source offsets (demagnified source
  σ ≈ 0.2–0.4 Å).

**Recon.** `PROBE_MODES` 3–5 to absorb the incoherence, combined with Step 1's C1 fit.

### Non-round aberrations (pulled forward from "after this phase", 2026-09-21)

`campaign/nonround_sweep.tsv` re-based on the current a070 balance (C3 −5 µm, C1 −60), so its `nr0`
row *is* the ladder's step-0 a70 leg; the C56/C34 magnitudes are unchanged (0.6 / 1.2 / 2.5 waves of
six-fold astigmatism, then 1.2 w six-fold + 1.2 w four-fold). Run through the known-probe pipeline by
label: `TSV=campaign/nonround_sweep.tsv LABELS="..." bash campaign/run_thin_atomfind.sh`; the row's JSON
reaches the sim as `ABERRATIONS_JSON`, `aberrations.json` carries the full tableau, and `make_probe.py`
passes it through unchanged, so the C1 search works on a non-round probe as it is.

### Step 4 — Phonons (thermal diffuse scattering)

**What exists.** `sim/simulate_4dstem.py --phonons N --phonon-sigma σ --per-species-sigma
--phonon-seed`; `run_sim.slurm` reads `PHONONS` / `PHONON_SIGMA`. **Not yet forwarded by
`campaign/run_thin_atomfind.sh`** — add them to `sim_job`'s export.

**Plan.** 8–16 configurations with `--per-species-sigma` (room-temperature RMS for Pb/Sr/Ti/O).
Simulation cost scales ×N; fine at BIN=4/2. **Built 2026-09-21:** `PHONONS=16 PER_SPECIES=1` on
`run_thin_atomfind.sh` forwards them to every leg (dirs `_ph16`); the noiseless sims took 1–4 min, so ×16
is under an hour. Launched at α = 70 and 90.

**Watch.** atomfind scores against mean positions, so xy- and z-RMS now include thermal smear.
Report the per-species σ alongside as the floor. Kernels get the same phonons: thermal blur is part
of the matched PSF.

### Step 5 — Specimen tilt, a few degrees off axis

**Magnitude (confirmed with the user, 2026-09-17): a few degrees off the zone axis** — not a small
accidental mistilt. Suggested ladder: 1°, 2°, 3°, 5° (17, 35, 52, 87 mrad). Note that this is the
same order as the convergence semi-angle itself.

**Why it may help.** Tilt maps depth onto in-plane position: a column's atoms advance laterally by
z·tanθ, so the finder gets a second, lateral handle on depth.

**How big the shear is** — lateral offset between a column's entrance and exit atoms, t·tanθ:

| tilt | 5-cell slab (19.5 Å) | 18-cell slab (70 Å) |
|---|---|---|
| 1° | 0.34 Å | 1.22 Å |
| 2° | 0.68 Å | 2.44 Å |
| 3° | 1.02 Å | 3.67 Å |
| 5° | 1.71 Å | 6.12 Å |

The nearest distinct columns in [001] projection (A-site to O column) are **~1.95 Å** apart. On the
thin slab the shear stays below that at every tilt on the ladder, so tilted columns remain separable.
On the 70 Å slab it passes 1.95 Å at about 2°, and **tilted columns cross their neighbours** in
projection — the column picture itself breaks down. So validate tilt on the thin slab first, and
expect Step 5 × Step 6 to need a finder that handles crossing tubes, not just sheared ones.

**What exists.** Nothing in the simulation. atomfind's CLEAN runs in **beam-parallel column tubes**
(`find.py` `clean_tube`, `find_atoms_v3`), and its lattice-aware species typing and guided
re-detection assume vertical columns.

**What to build.**
- **sim**: rotate the crystal about an in-plane axis — a physical tilt, surfaces tilting with it —
  rather than tilting the probe. Grow the box to hold it: the rotated slab needs ~W·sinθ of extra
  height (70 Å × sin 5° ≈ 6.1 Å) plus in-plane padding, and the exit surface shifts laterally by
  t·tanθ, so re-check that the scan field still sits over material at the exit side. Build the
  ground truth from the *same* rotated atoms.
- **GT cache**: regenerate from the rotated structure.
- **atomfind**: tilted tubes with their axis along the tilt vector; config keys for tilt magnitude
  and azimuth; audit species typing and guided re-detection, which assume vertical columns; and,
  for the thick slab past ~2°, a way to separate crossing tubes (e.g. fit neighbouring tubes jointly).
- **Kernels**: tilt the grid identically — the PSF itself tilts.

### Step 6 — Thickness: the full 70 Å labyrinth

**What exists.** `sim/PTO6_STO6_18_18_labyrinthPoscar.vasp` — 19,440 atoms, **70.01 Å along the beam
after orientation = 18 cells**. `build_thin_sample(n_cells)` behind `THIN_CELLS`, so `THIN=18`
reuses this driver unchanged. The PX915 NL70 pipeline (`NL70_coherent` preset, 70-layer recons)
proves 70-layer reconstructions run.

**Settings.** `THIN=18 ZVAC=4` → box ≈ 78.3 Å. Nyquist slices over that box: **NL ≈ 20 / 39 / 64 / 80**
at 50 / 70 / 90 / 100 mrad. *Measured 2026-09-21 by building the slab:* 19,139 atoms, atom span 69.94 Å,
actual box **77.94 Å**, 87 slices at 0.9 Å; the driver's formula gives 78.29 (CELL_Z 3.905) or 78.00
(`CELL_Z=3.889`), either way NL = 39 / 64 at 70 / 90. Launched with `THIN=18 CELL_Z=3.889 GROUPING=16
RTIME=20:00:00` (dirs `_thin18`). Still to do on the analysis side before its atomfind run: a thick GT cache
(`make_gt_cache --thin-cells 18 --z-vacuum 4`), a `thick` preset (trims, `zmax_show_A`, `clean_max_atoms`
derived from the 77.94 Å box), and `--dz` from `sim_meta.beam_thickness_A / NL`, not the driver's box.

**Watch.**
- **Walltime**: the thin a100 lab recon (NL 28) took ~7.5 h; NL 80 is ~3× that, so the driver's
  10 h BIN=2 walltime is too short. Raise `rtime_for`.
- **Cell spacing**: the driver hard-codes `CELL_Z=3.905` Å, but the structure's is 70.01/18 =
  3.889 Å — invisible at 5 cells, 0.3 Å of box error at 18. Derive `BOXZ` from the structure.
- **Ground truth**: regenerate with `make_gt_cache --thin-cells 18 --z-vacuum 4`. Do **not** reuse the
  packaged NL70 cache — its box and vacuum differ.
- **Typed constants**: `trim_z_A`, `zmax_show_A` and `clean_max_atoms` (~40 atoms per tube plus
  margin) are thickness-coupled. Three of this campaign's bugs came from constants typed for one
  geometry and silently wrong for the next — derive them from the GT and box instead.
- **Physics**: a thick specimen brings dynamical scattering and channelling, so thin-slab behaviour
  may not transfer. This is also where to retry Step 1 if it failed on the thin slab.

### Final — the combined configuration, full sweep

Every relaxation on together at 50 / 70 / 90 / 100 mrad (110 is out: its probe outgrows any scan
field the 70 Å box allows), matched kernels, atomfind with its calibrated error bars. Regenerate the
publication figures with the same scripts — `collate_atomfind_depth.py`,
`make_mep_volumes_fig.py`, `make_ronchigram_fig.py` — so the figures and the numbers cannot drift.

## Gaps worth closing along the way

These bit this campaign and will bite again on thicker, tilted, noisier data:

- **atomfind should fail loudly.** Add a self-check for the fraction of found atoms mapped outside
  the slab — it would have caught the one-unit-cell registration alias instantly. Derive the trim
  bands, `zmax_show_A` and `clean_max_atoms` from the GT rather than typing them.
- **`report.json` should store the full `Alignment`** (`mZ`, `X0`, `Y0`, `dx`, `dz`) so figure
  scripts can map ground truth without re-registering.
- **The sim's `[geom]` halo check ignores the aberrated probe size** (it uses thickness·tanα) and
  reports comfortable margins where there are none. Use the probe's d99.
- **The a90/a100 `BETA_LSQ` inconsistency** (kernels 0.05, labs 0.1) closes itself: every step here
  re-runs the labs at the pinned 0.05.

Later, after this phase: non-round aberrations (`campaign/nonround_sweep.tsv`, 70 mrad with
six-fold astigmatism escalated), then a C5-corrected comparison instrument.
