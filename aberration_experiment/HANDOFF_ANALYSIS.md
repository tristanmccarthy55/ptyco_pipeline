# Handoff — the CEOS-approx sweep: an older corrector opened to 80 mrad, probe known

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

## The current job — the CEOS-approx sweep to 80 mrad (running since 2026-09-24 18:39)

The question: can an older, widespread hexapole-corrected column still give depth-resolved ptychography by
**tanking** its aberrations — probe known, aperture opened far past the corrector's 30 mrad design? The user chose
a **CEOS approximation** (CESCOR-like, pre-DELTA; the target era's columns are 200 kV, the simulation stays 300 kV).

**The instrument** (`campaign/aberration_waves.py` `ceos_tableau`): A5 = 1.0 mm (CEOS; consistent with those
correctors' ~35 mrad usable range, and DCOR/ASCOR, built to beat it, hold < 0.2 mm); parasitic A4 = B4 = 10 µm
(unsourced central value — A5 must stay the limiting term, which caps them near 27 µm; the first choice, 40 µm,
broke that); tuned terms at 0.1 waves at the aperture they are tuned at. C5 = 1 mm fixed, as the round sweep.

**The operator fights for probe size, not flatness** (`plan_probe.py --ceos`): A1 A2 S3 A3 D4 retuned at each
aperture, coma B2 set against B4, C1/C3 rebalanced by Nelder-Mead on d90. A5 and A4 have no same-symmetry knob and
set the floor (six-fold alone ≈ the whole probe from 70 mrad).

**The geometry that makes 70–80 mrad simulable** (per-row columns `side win step detmax`, read by the driver):
every leg is a **210 Å cut of the tiled, periodic labyrinth** (`simulate_4dstem --region-side`), centred on the
same material the 70 Å box put under its scan (all 1171 atoms within 15 Å identical to 1e-14 Å) — bulk crystal all
round, where the 70 Å box padded the 47.9 Å x-period with 22 Å of vacuum ~10 Å from the scan. Recon window =
side/BIN from d99; the detector is recorded to ±200 mrad except at 80 mrad, cropped to ±133 (`--detector-max-angle`;
the potential stays sampled for 200 mrad) so N stays ≤ 1424 px; scan field = 1.5 × d90 at 1600 positions.

| α | fought d90 / d99 | window (BIN) | detector | scan / step | resource class |
|---|---|---|---|---|---|
| 40 | 4.2 / 11.0 Å | 17.5 (12) | ±200 | 20 / 0.5 | BIN-4 |
| 50 | 4.8 / 10.8 Å | 17.5 (12) | ±200 | 20 / 0.5 | BIN-4 |
| 60 | 11.6 / 19.6 Å | 35 (6) | ±200 | 20 / 0.5 | BIN-2 |
| 65 | 17.5 / 27.4 Å | 35 (6) | ±200 | 27 / 0.675 | BIN-2 |
| 70 | 25.2 / 38.7 Å | 70 (3) | ±200 | 38 / 0.95 | BIN-1 |
| 75 | 35.4 / 55.1 Å | 70 (3) | ±200 | 54 / 1.35 | BIN-1 |
| 80 | 49.1 / 75.7 Å | 105 (2) | ±133 | 74 / 1.85 | BIN-1 |

Round controls `round_a040`…`round_a080` (4 Å probes) share the region geometry at BIN 12.

**Preflight (2026-09-24).** Driver dry-run with a stub `sbatch`: every submission carries its region, window, step
and detector angle, with resources sized from the pattern size. Found and fixed: at 80 mrad the half-width presolve
of a ±133 mrad detector would keep ±66 mrad and clip the probe's own 80 mrad aperture; `run_synthetic_recon_ML.m` now
widens any presolve that would clip the aperture to 1.2α (80 mrad: 1024 px, ±96 mrad; no earlier leg affected).
Budget, from the 110 mrad BIN-1 logs (~7 h at 1426 px, 34 layers): the ~1420 px legs take 3–4.5 h against a 24 h
request (48 h cap); raw data peaks ~160 GB and `CLEANDATA` removes it; the heavy block runs `PACK_H5=0` and is
analysed on Blythe by `campaign/run_analysis.sh` (`analysis/analyse_sweep.py`, validated locally: it reproduces the
six-fold numbers exactly).

**Caveats to carry into the write-up.** 80 mrad needs a detector of ~1400 px across its ±133 mrad — a modern
large-format camera on an old column. A4/B4 are unsourced (40 µm would give 30 / 56 Å at 70 / 80 mrad). 40 mrad is
new ground (NL 4). The analysis needs the region GT (`make_gt_cache --region-side 210`, already built at
`~/Desktop/ceos_region_gt`), `--set scan_center_xy=105,105`, and each leg's `dx` from its h5.

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

## State (2026-09-26) — first CEOS results are in, and mostly FAILED: diagnose before anything else

**Results** (`~/Desktop/ceos80_analysis/analysis_round_a040-…_20260924_183951/summary.csv`, run on Blythe by
`campaign/run_analysis.sh`; per-leg logs under `logs/`):

| outcome | legs | what the logs say |
|---|---|---|
| **worked** | CEOS 60, 65 (BIN 6, 35 Å window) | Pb/Ti/O 92/81/78 % and 88/65/74 %, z-RMS 0.85 / 0.86 Å; kernels the CLEANEST yet (peak/bg 899–1193 vs 379 before; **phase ramp 0.001 rad, was 0.28 — the point-sampling fix works**) |
| **kernels: 0 grid atoms** | every BIN 12 leg (17.5 Å window): round 40/50/60/70, CEOS 40/50 | `extract_psf`: "only 0 grid atoms found in the inner 13 A", for Pb AND Ti (summary.csv showed only Ti — fixed) |
| **saturated** | round 75, 80 (BIN 12); CEOS 70 (BIN 3), 80 (BIN 2); round 65's Pb leg | triage phase std 0.2–0.93, wrapped up to 1e-2 |
| **junk** | CEOS 75 (BIN 3) | kernels are speckle ("156 atoms at 0.69 A spacing", peak/bg 24); Pb recall 7 % |

**Hypothesis 1 — windows too small for point sampling (my change).** Point sampling is exact only if the whole exit wave
fits the window: probe tails plus electrons scattered sideways through the 27.5 Å slab (±~5 Å at 200 mrad). Whatever
leaves the window aliases back in. The old 4×4 sum filtered those components, which is why 17.5 Å windows used to work;
`plan_probe.region_geometry`'s thresholds (d99 < 15 → 17.5 Å) were calibrated on summed data. Evidence: every BIN 12
leg fails at every α (even 40 mrad, 4 Å probe); every BIN 6 leg is clean. **Test (rows built, not run):**
`round_a070_w35` and `ceosopt_a050_w35` = the same legs at a 35 Å window (BIN 6), ~40 min recons:
```bash
CLEANDATA=1 TSV=campaign/ceos_sweep.tsv LABELS="round_a070_w35 ceosopt_a050_w35" bash campaign/run_thin_atomfind.sh
LABELS="round_a070_w35 ceosopt_a050_w35" GT_REGION=210 DEP=<its pack job> bash campaign/run_analysis.sh
```
If both give clean kernels: raise the minimum window for point-sampled runs to 35 Å (region_geometry) and rerun the
BIN 12 legs. (The alternative, the engine's own detector-upsampling model of a summed detector, costs a 16x bigger model.)

**Hypothesis 2 — large probes break the solve.** CEOS 70/75/80 (d90 25/35/49 Å) fail with 70–105 Å windows, so not
aliasing; the round 110 mrad leg (d90 24.5 Å) failed the same way in the old campaign. Could be the solver (beta 0.05,
200 iterations, step 0.95–1.85 Å) or a real limit. First look at the recon logs and `*_error_trace.csv` (diverging?
stalling?) before changing anything; candidates: smaller BETA_LSQ, more iterations, smaller step / more positions.

**Data to pull for the diagnosis** (the recon logs + error traces are in the sweep tarballs; the two light ones also
carry their h5s, ~1–2 GB each; the heavy one has logs/sidecars only):
```bash
mkdir -p ~/Desktop/ceos80_sweeps && cd ~/Desktop/ceos80_sweeps
scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:/springbrook/share/physics/phucrh/atomfind_results_round_a040-*_20260924_183859.tgz' .
scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:/springbrook/share/physics/phucrh/atomfind_results_ceosopt_a040-*_20260924_183900.tgz' .
scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:/springbrook/share/physics/phucrh/atomfind_results_ceosopt_a070-*_20260924_183900.tgz' .
for t in *.tgz; do mkdir -p "${t%.tgz}" && tar xzf "$t" -C "${t%.tgz}"; done
```
The heavy legs' h5s (~1 GB each) stay on Blythe at `recon_af_ceosopt_a0{70,75,80}_*/analysis/`.

**Nothing from this sweep goes in the ladder, figures or logbook until both hypotheses are settled.**

Older data: the fixed-engine six-fold analysis is in `~/Desktop/sixfold_0923_analysis`; the region GT at
`~/Desktop/ceos_region_gt` and on Blythe `$SHARE/phucrh/gt_region210`. Raw sim data is deleted after packing.
