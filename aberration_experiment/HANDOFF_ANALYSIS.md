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

### 2026-09-28 — read this first; it overturns parts of what follows

**An analysis bug made several "failures".** `extract_psf.load_vol` looked for engine checkpoints (`Niter*.mat`) before
the h5. On Blythe a recon dir holds one per engine — the half-resolution presolve too — and one set per earlier run, all
`Niter200.mat`, so `run_analysis.sh` read whichever the filesystem listed (the logs name the file: CEOS 60/65 got the
full solve by luck, CEOS 75 the presolve). Local tarballs carry only the h5, so local analyses were right. Fixed: the h5
first, never `analysis.prev_*`, and `analyse_sweep` passes the h5 itself. What changes:
- **The 17.5 Å window works.** Every BIN 12 kernel is clean from its h5 (25 atoms, peak/bg 101–799). Re-scored locally
  (`~/Desktop/ceos80_sweeps/analysis_local_fixed`): round 70 **98/83/81 %, z 0.56 Å** (= old a70 row; its 35 Å run gave
  91/74/66, 0.62), round 50 60/60/51 (1.04), round 60 91/82/71 (0.88), CEOS 50 69/58/61 (1.12; 35 Å: 68/53/52). Round 65's
  Pb kernel leg genuinely half-failed (ε 0.11 → 0.066). The 40 mrad legs crash atomfind (`find.py` `np.gradient` on a
  one-layer kernel at NL 4) at either window. **The 35 Å minimum is withdrawn** (`region_geometry` restored); the
  09-27 rerun at 35 Å stands as the window comparison once re-analysed.
- **CEOS 75 was never judged**: lab ε 0.075, vacuum phase std 0.006 — healthy. Re-analysis handed over.

**Round 75/80, answered by the hia_* tests** (`~/Desktop/ceos_hia`): not the chain — `hia_a090_w35` gives 95/85/87 %,
z 0.42 Å (old a90 95/94/95, 0.37) and `hia_a090_std_sum` converges (ε 0.021). Round 80's probe **diverges in the
presolve from iteration 2** at NL 18/23 (presolve ε 0.13 → 0.65 at NL 23) and converges at NL 14 (ε 0.073); a90's knobs
at 80 mrad converge at NL 18 (ε 0.027). Milder presolve rises appear in many legs and the full engine usually recovers.
**Candidate rule** (fits every compact-probe leg we have; not yet tested ahead of a run): the probe's intensity pattern
must change ≳ 6 % per recon slice through the box (worked 6.2–7.7 %, diverged ≤ 4.2 %; round 40 at NL 4 is the
exception). Large CEOS probes change 1–2 % and still split, so it is not universal. Metric: 1 − mean Pearson r of
|P(z)|² vs |P(z+dz)|², z over the 27.5 Å box, abTEM probe (session script; rebuild if needed).

**Submitted 09-28** (rows at the end of `ceos_sweep.tsv`): round controls `round_a075_b8` (C3 −8 µm, C1 −130: d90 4.1,
d99 5.2 Å, 5.6 %/slice — also tests the rule) and `round_a080_b9` (= hia_a080_k90 with kernels); large CEOS, one change
each against ceosopt_a070/080: `_f20` (scan 20 Å at 0.5 Å, drops Rule 4 — which was never shown to help: widening
a110's field did not rescue it), `_s05` (CEOS 70 scan 38 Å at 0.5 Å, 5776 positions — step only), `_nl10` (NL only).
`_s05`'s sim failed (job 1294737): its float64 region buffer needs ~140 GB against the sim's 128 GB (now float32,
bit-identical output), and its recon would need ~4× the 47 GB of data, above a GPU node. Replaced by `_s06`: 38 Å at
0.594 Å, 4096 positions (~33 GB data, ~135 GB recon peak). **Position ceiling at 1422 px: ~5500 per recon.**

**Results** (`~/Desktop/ceos80_analysis/analysis_round_a040-…_20260924_183951/summary.csv`, run on Blythe by
`campaign/run_analysis.sh`; per-leg logs under `logs/`):

| outcome | legs | what the logs say |
|---|---|---|
| **worked** | CEOS 60, 65 (BIN 6, 35 Å window) | Pb/Ti/O 92/81/78 % and 88/65/74 %, z-RMS 0.85 / 0.86 Å; kernels the CLEANEST yet (peak/bg 899–1193 vs 379 before; **phase ramp 0.001 rad, was 0.28 — the point-sampling fix works**) |
| **kernels: 0 grid atoms** | every BIN 12 leg (17.5 Å window): round 40/50/60/70, CEOS 40/50 | `extract_psf`: "only 0 grid atoms found in the inner 13 A", for Pb AND Ti (summary.csv showed only Ti — fixed) |
| **saturated** | round 75, 80 (BIN 12); CEOS 70 (BIN 3), 80 (BIN 2); round 65's Pb leg | triage phase std 0.2–0.93, wrapped up to 1e-2 |
| **junk** | CEOS 75 (BIN 3) | kernels are speckle ("156 atoms at 0.69 A spacing", peak/bg 24); Pb recall 7 % |

**Error traces** (pulled 2026-09-26 to `~/Desktop/ceos80_sweeps/`). Relative amplitude error ε ≈ residual × N / 2e5
(the residual itself scales ~1/N): old summed a70 at 17.5 Å 0.040; **round 70 point-sampled at 17.5 Å 0.050 (worse)**;
round 75/80 at 17.5 Å **0.47** (first full-engine iteration ~370 vs ~37 at 70); CEOS 60/65 at 35 Å **0.033** (best);
CEOS 75 at 70 Å 0.075 (converges, yet object is speckle); CEOS 70 / 80 0.27 / 0.45. Logs clean: no NaN, 80 mrad
`presolve widened to 1026 px`, runtimes 3.9–5.2 h, 1600 positions.

**Hypothesis 1 — the 17.5 Å window fails with point sampling.** Fact: every BIN 12 leg is bad, every BIN 6 leg good.
The mechanism is NOT the beam spreading through the slab: measured probe spill outside the window stays 0.3–0.5 % from
entrance to exit and is the same for round 70 (half-works) as round 75/80 (saturate); the old 4x4 sum filtering the
aliased part is still a candidate, unproven. **Test (rows built, not run)** — the same legs at a 35 Å window (BIN 6):
```bash
CLEANDATA=1 TSV=campaign/ceos_sweep.tsv LABELS="round_a070_w35 round_a080_w35 ceosopt_a050_w35" bash campaign/run_thin_atomfind.sh
LABELS="round_a070_w35 round_a080_w35 ceosopt_a050_w35" GT_REGION=210 DEP=<its pack job> bash campaign/run_analysis.sh
```
If they come back clean: minimum window 35 Å for point-sampled runs (`plan_probe.region_geometry`), rerun the BIN 12 legs.

**Window test RESULT (2026-09-27, `~/Desktop/ceos_w35`, analysed locally in `analysis_local/`).** At a 35 Å window the
50–70 mrad legs recover: kernels clean (25 atoms, peak/bg 494–1066, ramp 0.000); round 70: ε 0.050 → 0.030, Pb/Ti/O
91/74/66 %, z-RMS 0.62 Å, precision 0.989 (old summed 70 Å-box a70: 98/82/75, 0.56 — the small gap is open); CEOS 50:
68/53/52 %, z 1.11 Å (≈ old round a50 67/43/54, 1.20 — the CEOS probe costs nothing at 50). **So the minimum window for
point-sampled runs must be 35 Å** — change `plan_probe.region_geometry` (17.5 Å only when...: never) and rerun the BIN 12 legs.
**But round 80 still saturates at 35 Å (ε 0.47): a third problem at ≥75 mrad**, round and CEOS alike, although the old
pipeline reconstructed a90/a100 at the same 35 Å window and 712 px. **Test (rows built, not run)** — round 80 in the OLD
70 Å box (lazy scan path, no region) at BIN 2, once summed, once point-sampled:
```bash
CLEANDATA=1 DETECTOR_SAMPLING=sum TSV=campaign/ceos_sweep.tsv LABELS=round_a080_std_sum bash campaign/run_thin_atomfind.sh
CLEANDATA=1 TSV=campaign/ceos_sweep.tsv LABELS=round_a080_std_pt bash campaign/run_thin_atomfind.sh
```
sum works + pt fails → point sampling is the cause at high α; both work → the region box / batched sim; both fail →
the round a080 balance (C3 −7 µm, C1 −100 Å) or the solve at NL 18. Judge on the error trace + triage (fast), then kernels.

**High-α test RESULT (2026-09-27, `~/Desktop/ceos_a80std`).** round 80 in the OLD 70 Å box at BIN 2 fails with the
summed detector (ε 0.46) AND point sampling (ε 0.44) — exactly as in the region box. So **neither point sampling nor the
region box / batched sim causes the ≥75 mrad failure.** Remaining suspects: the planner's round 75/80 rows themselves
(C3 −6/−7 µm, C1 −80/−100 Å: the 4 Å target reached by heavy defocus, "free"/"floor; ok" regime) or the solve at NL
16/18. The old pipeline reconstructed a90 (C3 −9 µm, C1 −160 Å, the smallest probe, NL 23, BIN 2, summed) cleanly
(step-0 ladder row). **Next tests, cheapest first:** (1) the OLD a90 row through the new chain (region, point, BIN 6)
— if it works the chain is fine at high α and the 75/80 rows are the problem; (2) round 80 at the smallest-probe
("floor") balance instead of 4 Å-by-defocus; (3) only then the slice count. The CEOS 75/80 legs share whatever this is.

**High-α diagnosis (2026-09-27 evening) — facts checked locally, tests BUILT, NOT RUN.**
- *Test (2) as written is void*: the round 80 row (C3 −7 µm, C1 −100 Å, d90 4.02 Å) already sits at the smallest-probe
  balance (C3 −7 µm, C1 −96 Å, d90 3.77 Å; searched C3 −12…−2 µm with the planner's own probe builder). "4 Å by heavy
  defocus" is not what round 80 is. (At 75 it is: floor 2.87 Å at C1 −68, row −80.)
- *The probe is sampled fine*: the packed round 80 probe (old box, 35 Å) has d90 4.01 / d99 10.8 Å and 3e-4 of its
  intensity beyond 0.45 W — better contained than the a90 probe that worked (6.5 / 15.3 Å, 5e-4).
- *The object is lost before the full engine*: failed legs have phase std 1.0–1.3 rad in EVERY layer, vacuum included,
  |O| 0–1.6, near-Nyquist noise (good round 70 w35: 0.003 rad in vacuum, 0.05–0.08 in the slab). The full engine starts
  at ε 0.46–0.62 (round 70: 0.08), jumps at iteration 12 (201, 241) and never recovers. The presolve's trace was never
  recorded — `core/ptycho_recons.m` now prints every engine's trace to the slurm log (`engine 1 error trace`).
- *Nothing in the recon driver depends on NL or α* except the presolve guard, which did not fire.

The tests (rows `hia_*` in `ceos_sweep.tsv`; each changes ONE thing against a leg we have):

| label | against | changes | if it works |
|---|---|---|---|
| `hia_a090_w35` (lab+Pb+Ti) | old a90, 95/94/95 %, z 0.37 | region + point + today's engine | the chain is fine at 90; compare recall |
| `hia_a090_std_sum` (lab, `DETECTOR_SAMPLING=sum`) | old a90 | today's engine only | only matters if the first fails: engine vs new sim |
| `hia_a080_nl14`, `hia_a080_nl23` (lab) | round_a080_w35 (NL 18) | slice count only (col 15 `nl_force`) | NL 16/18 is the cause |
| `hia_a080_k90` (lab) | round_a080_w35 | a90's C3/C1 at 80 mrad: d50/d90/d99 4.5/5.6/6.9 vs 2.2/4.0/11.3 Å | the round 80 probe is the cause |

Judge on the error trace, the layer stats (phase std in vacuum layers) and, for `hia_a090_w35`, recall.
NL 14/16/18 are also what CEOS 70/75/80 ran: if the slice count is the cause, job 3 (large CEOS probes) may share it.

**Window rule (2026-09-27) — WITHDRAWN 09-28, see the top of this section**: `plan_probe.region_geometry` minimum window 35 Å; the BIN 12 rows of
`ceos_sweep.tsv` are BIN 6. Rerun: round 40–70 + CEOS 40/50 (70 and CEOS 50 repeat their w35 runs: a free run-to-run
spread). Round 75/80 wait on the diagnosis.
**Packing (2026-09-27)**: `pack_results.sh` now packs only the sims this submission's recons link to — a one-label
tarball was 940 MB, 475 MB of it 169 other runs' probes.

**Hypothesis 2 — large probes break the solve.** CEOS 70/75/80 spill only 0.03–0.1 % outside their 70–105 Å windows, so
not geometry. Round 110 (d90 24.5 Å) failed the same way in the old campaign. No crash; the fit is poor from the first
full-engine iteration (70, 80) or fits while the object is speckle (75). Next: look at the objects (h5s on Blythe,
`recon_af_ceosopt_a0{70,75,80}_*/analysis/`) and try one lever at a time on CEOS 70 — BETA_LSQ, NITER, positions.

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
