# Handoff — next job: a representative non-round microscope on the round sweep

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

## The next job

Put a **representative non-round corrector tableau on the round sweep**: the same round balance per aperture
(`campaign/round_sweep.tsv`: C5 = 1 mm fixed, C3 and C1 retuned), plus every non-negligible non-round residual,
held fixed so that each grows as α^(n+1). The result wanted is the round sweep's — atomfind improving with α —
for a realistic probe.

**Blocked on one decision: which corrector.** The magnitudes decide whether the sweep can exist at all:

| α (mrad) | probe d90 / d99, round only | + the campaign's `ASSUMED_TABLEAU` | largest terms (waves at the edge) |
|---|---|---|---|
| 50 | 3.9 / 8.2 Å | 11 / 21 Å | three-lobe 3.2, five-fold 1.6, six-fold 1.3 |
| 70 | 4.0 / 10.3 Å | 45 / 79 Å | three-lobe 17, six-fold 10, five-fold 8.5 |
| 90 | 6.5 / 16 Å | 120 / 168 Å | three-lobe 60, six-fold 45, five-fold 30 |
| 100 | 11 / 26 Å | 133 / 173 Å | three-lobe 100, six-fold 85, five-fold 51 |

The box is 70 Å and the scan field at most ~34 Å, so that tableau is unsimulable from 70 mrad up. What is
sourced, and what is not:

- **CEOS, CESCOR hexapole probe corrector** (ceos-gmbh.de): operator-adjustable C1, A1, B2, A2, C3, S3, A3,
  **D4**, C5; limiting **A5 (intrinsic) and A4, B4 (parasitic)**. The tableau had D4 (three-lobe) filed as
  hardware at 100 µm; it is adjustable. Corrected in `campaign/aberration_waves.py`.
- **A5 = 1.0 mm** — CEOS's figure for hexapole-type correctors. 10 waves at 70 mrad on its own.
- **DCOR/ASCOR**: A5 0.03–0.4 mm (Ultramicroscopy 2019, "On the residual six-fold astigmatism in DCOR/ASCOR").
  **JEOL DELTA** (fifth-order corrector): aberration-free illumination to ~70 mrad at 300 kV.
- **A4, B4**: no sourced typical value yet. Every other value in `ASSUMED_TABLEAU` is an unsourced order of
  magnitude. A measured tableau from the target instrument would replace all of them.

So a Cs-only hexapole instrument cannot be opened past ~50 mrad even with the probe known — a real result in
itself — and a sweep to 100 mrad needs an A5-corrected class of corrector. That choice is the user's.

**Whatever is chosen, re-simulate the round controls in the same submission**: the detector model changed
(below), so every existing ladder row is on the old summed detector and is not a like-for-like baseline.

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
