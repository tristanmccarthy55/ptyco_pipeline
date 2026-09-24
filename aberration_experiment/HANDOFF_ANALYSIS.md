# Handoff — the CEOS-approx sweep: an older corrector, opened wide, probe known

Written 2026-09-24. **Read this whole file before running anything; it is short on purpose.**
Then `HANDOVER.md` for the experiment, `PSF_KERNELS.md` for the kernel rules.

---

## Where it stands

- **Round-α campaign: complete.** Known probe, 50–100 mrad, depth error falls with aperture; 110 mrad is
  not reconstructible (its probe outgrows the scan field). `HANDOVER.md` has the table.
- **Relaxation ladder** — one cumulative table, `results/relaxation_ladder.csv`, 18 rows:

  | step | assumption removed | outcome |
  |---|---|---|
  | 0 | — (known probe, noiseless) | a70 98/82/75 % Pb/Ti/O, depth error 0.56 Å; a90 95/94/95 %, 0.37 Å |
  | 1 | focus unknown | no cost, by an outer search over fixed trial probes; the in-solver fit is closed |
  | 2 | shot noise | holds to 10⁵ e/Å²; at 10⁴ no single-atom reference can be measured |
  | 3 | thermal vibration | costly: depth error roughly doubles, oxygen recall 10–17 % |
  | 4 | known non-round probe (six-fold) | **no cost** — see below |
  | 5 | the 70 Å sample | fails at both solver step sizes |

- **A known non-round probe costs nothing** (2026-09-24, fixed engine, lab + fresh matched kernels):

  | leg | residual | Pb / Ti / O bulk recall | depth error |
  |---|---|---|---|
  | 70 mrad round | 22.6 | 98 / 82 / 75 % | 0.56 Å |
  | 70 mrad six-fold 0.1 waves | 24.5 | 98 / 82 / 76 % | 0.55 Å |
  | 70 mrad six-fold 0.45 waves | 25.5 | 98 / 82 / 71 % | 0.60 Å |
  | 90 mrad round | 5.8 | 95 / 94 / 95 % | 0.37 Å |
  | 90 mrad six-fold 0.1 / 0.3 / 0.45 waves | 5.8 / 5.7 / 5.6 | 93–96 / 94 / 94–95 % | 0.38–0.39 Å |

  Objects correlate with the round control layer by layer at 0.98 (70 mrad, 0.45 waves) and ≥ 0.99 (90 mrad).
  The small leftovers (oxygen 71 vs 75 %, oxygen false positives 35 → 49 at 90 mrad) are single solves; the
  engine reseeds every run, so they are not separated from run-to-run spread.

## The next job — the CEOS-approx sweep (built 2026-09-24, first submission pending)

The question: can an older, widespread hexapole-corrected column still give depth-resolved ptychography by
**tanking** its aberrations — probe known, aperture opened well past the corrector's 30 mrad design?
The user chose a **CEOS approximation** (CESCOR-like, the pre-DELTA generation; the target era's GrandARM-class
columns are 200 kV, the simulation stays at 300 kV). Defined in `campaign/aberration_waves.py` `ceos_tableau`:

- round part exactly as the round sweep: C5 = 1 mm fixed, C3/C1 retuned per α by `plan_probe.py`;
- non-round, all **fixed** (each grows as α^(n+1)), split per CEOS: tuned A1 B2 A2 S3 A3 D4, hardware B4 A4 A5;
- **A5 = 1.0 mm is CEOS's figure; every other term is unsourced** and set to 0.1 waves at 30 mrad (inside the π/4
  flatness a corrector is tuned to) — A1 0.44 nm, B2/A2 22 nm, S3/A3 1 µm, D4/B4/A4 40 µm. S5/R5 left out
  (CEOS names A5, A4, B4 as limiting). One seeded draw of orientations; probe sizes vary < 1 Å across draws.

**The operator fights for probe size, not flatness** (the user's call: only a small, simulable probe matters to a
known-probe reconstruction). A knob can only cancel an aberration of its own symmetry, so: the tunable A1 A2 S3 A3 D4
are retuned at each aperture (to the same 0.1-wave accuracy), coma B2 is set against fourth-order coma B4, and C1/C3
are rebalanced — `plan_probe.py --ceos`, Nelder-Mead on d90. **A5 and A4 have no partner knob and set the floor**:
A5 alone gives 25 Å at 70 mrad. Fighting shrinks the probe 20–35 % and moves the wall from ~60 to ~65 mrad.

Rows in `campaign/ceos_sweep.tsv`; `ceosopt_a<A>` is the sweep:

| α | round control | `ceosopt` probe d90 / d99 | reconstruction | as built (`ceos_a<A>`, C1/C3 only) |
|---|---|---|---|---|
| 40 | BIN 4 | 4.0 / 10.3 Å (4 Å by defocus, as the round rule) | BIN 4 | 5.4 / 12.5 Å |
| 50 | BIN 4 | 6.4 / 12.8 Å | BIN 4 | 9.8 / 19 Å |
| 60 | BIN 4 | 14.6 / 25 Å | BIN 2, **WIN=22 STEP=0.55** | 18 / 35 Å |
| 65 | BIN 4 | 21.1 / 36.5 Å | BIN 1, **WIN=32 STEP=0.8** | — |
| 70 | — | 29.8 / 50 Å | **not run**: needs a 45 Å scan field, the box allows ~34 | 35 / 63 Å |

**The wall is the result, not a failure of the sweep.** Past ~65 mrad this instrument's probe outgrows anything
this sample can host, and the same numbers bound a real experiment: the reconstruction window is λ / detector
pixel angle, so a 50 Å probe at 70 mrad needs ~1000 detector pixels across against ~200 for the 6 Å probe at 50.

**Watch on the first logs:** 40 mrad is new ground — NL 4 (6.9 Å slices), `extract_psf --zdrop 1`; the round
sweep skipped 30 mrad because the bright-field disc filled too little of the detector, and 40 may be marginal.
And every leg is the first reconstruction on point-sampled data (below).

## Changed 2026-09-23/24 — if the next run breaks, look here first

1. **Engine, `save_to_p.m` (`40d28a3`).** The probe now leaves each engine in the frame of `probe_initial.mat`,
   so each engine transposes it exactly once (`custom_data_flip`). Before, the full engine re-flipped the
   presolve's already-flipped probe. Round probes are unaffected (a round probe is its own transpose).
   Signature of the old bug: the full engine's first-iteration error far above the round leg's (98 vs 41).
2. **Simulator, detector sampling.** `sim/simulate_4dstem.py` now point-samples the fine detector at the
   reconstruction's own k-grid (every BIN-th pixel from the zero-angle pixel) instead of summing 4×4 blocks.
   The sum is a pixel-integrating detector the engine does not model; on its own it produced a residual floor
   of 20.2 at 70 mrad / BIN 4 and 5.3 at 90 / BIN 2 (recorded floors 22.6 and 5.8 — the "factor of four"), and
   because its blocks start on the zero-angle pixel it displaced every pattern by (BIN−1)/2 fine pixels (checked:
   −0.375 coarse px), which the object absorbed as the diagonal phase ramp analysis has been subtracting.
   **Expect** lower residual floors, `extract_psf` reporting almost no ramp (was ~0.28 rad at 70 mrad, ~0.055 at
   90), and `detector_sampling: point` in `sim_meta.mat`. **To reproduce old data:** `DETECTOR_SAMPLING=sum`.
   BIN 1 is identical either way. Never reconstructed yet — the first run on it is the test.
3. **`PRESOLVE_NDP`** stays as a diagnostic knob only. The presolve crops the *detector* to ±100 mrad, which
   keeps the whole probe aperture; forcing it to full width gave 73.4 against 74.0 and no better object.

## What is solid — do not re-litigate

- **Focus**: an outer search over fixed trial probes recovers it at no cost. The in-solver fit is closed.
- **Shot noise**: 10⁷ and 10⁶ e/Å² hold, 10⁵ works with the light atoms going, 10⁴ fails.
- **Known non-round probe**: no cost at 0.45 waves of six-fold at 70 and 90 mrad (table above).
- **The 70 Å sample fails** at solver steps 0.05 and 0.02. Separate problem, picked up later.
- **Aperture scaling** (`campaign/aberration_waves.py`) is analytic.

## How this project works

- **The user runs every cluster job.** No SSH. Give exact blocks; assume nothing until you see a log.
- **Everything lands under `$SHARE/phucrh`**, never the group root. Ship a matching `scp -O` line with every
  block, globbing the stable part of the tarball name — or the exact name when an older tarball shares its stem —
  with `tar` on its own line into a fresh directory. The user's `tar` step has twice not run; check before triage.
- **The share is 3.9 TiB and shared with the department.** Submit with `CLEANDATA=1`. Check quota with
  `mmlsquota --block-size auto`, never `df`.
- Local Python is `~/hyperspy-bundle/bin/python`. Commit code, figures and results; leave the user's report
  and `analysis/atomfind/*.md` alone.
- The published record is one page, https://claude.ai/artifact/7ve93UM6yqCmiRcbfNuiJM, rebuilt from
  `page/logbook.html` by `page/build_page.py` and republished to that same URL.

## Traps

| symptom | cause | rule |
|---|---|---|
| a fix to how the probe is loaded changes recall but not the final residual | the fix reached only the presolve: `save_to_p` handed the full engine the flipped probe and it flipped it back (fixed `40d28a3`). The final residual is the full engine's | check the saved probe and the full engine's first-iteration error, not just recall. Once the probe was right in both engines the residual fell 74 → 25.5 |
| a probe or focus change lowers the residual while the object degrades | a free object can absorb some probe error (the in-solver focus fit ended at a lower error with the object collapsed) | judge a probe change on the object as well as the residual |
| a leg passes triage and is still junk | it reported COMPLETED, wrote an h5 and logged no NaN, but the phase is saturated | `analysis/triage_recon.py`: genuine legs wrap ≤ 5×10⁻⁵ of pixels (edge of the illuminated field), saturated ones ~3×10⁻² |
| GPU jobs die with "Disk quota exceeded" naming a path with no quota | `$HOME` is 2 GB and the CUDA JIT cache `~/.nv` fills it | both slurm scripts redirect every cache off `$HOME` |
| grepping a recon's slurm log reads the wrong run | recon dirs keep logs from every attempt | check the job id, not the glob order |
| raising NL to chase a residual makes it worse | NL is Nyquist; more layers with `REGLAYER=0` destabilises the solve (NL 28 gave 739 against 74 at NL 14) | leave NL alone |
| six jobs pending on `DependencyNeverSatisfied` | a simulation whose output exists exits 1 in zero seconds | `RECON_ONLY=1` when only the solver changed, `OVERWRITE=1` to re-simulate |
| `scp` finds nothing although jobs finished | the tarball name carries the submission timestamp | glob the stable part, never the date |

**Two rules that are physics, not plumbing:** `REGLAYER=0` on every leg (it low-passes the depth axis, which
is the measurement), and one `BETA_LSQ` for every leg in a comparison.

## State

`origin/main`; Blythe needs `git pull`. Raw simulation data is deleted after packing (`CLEANDATA=1`), so any leg
not on disk needs re-simulating. The fixed-engine six-fold analysis is in `~/Desktop/sixfold_0923_analysis`
(kernels in `psf/`, one `atomfind_<label>/` per leg), tarballs in `~/Desktop/sixfold_0923_a70` and `_a90`.
