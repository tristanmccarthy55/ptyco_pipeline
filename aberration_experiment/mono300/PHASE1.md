# Phase 1: focal spread and source size, with probe modes (designed 2026-10-12; NOT submitted as designed)

**THE PAPER AND THE PLAN (the user, 2026-10-16).** One story in short sections, each exploring one limiting factor:
a bad (pessimistic, ARM) set of aberrations and what happens when the aperture opens; partial coherence (the focal-spread
and source walls above); then cooling; then the chosen instrument (the Spectra 300 X-FEG **Mono**, not UltiMono) and its
outputs at doses x-y-z. Runs, in order (cheap "possible at all?" first, realism last):
- **A. ceiling, coherence only** at the realistic Mono settings: fs 2.5 A (~90 meV through Cc ~1.6 mm) with a 0.1 and a
  0.2 A source, mode legs from the mixed state (`ceiling_submit.sh`, `_mixF`), beside step 2's fixed controls.
- **B. ceiling + phonons at 100 K (LN2) and 20 K (LHe; Debye model at its zero-point floor)**, fs 2.5 + src 0.1, fixed
  control + mode legs (`ceiling_submit.sh`, `_mixF` / `_mixF_fixed`). Room temperature is a no-go (user).
- **C. realistic start** (`search_submit.sh`, built 2026-10-16): C1, C3, A1 errors DRAWN from the corrector's quoted
  precision (ARM run 3 ci95: C1 39.5 A, C3 0.804 um, A1 59 A; sigma = ci95/1.96, A1 per vector component; seed 1 -- the
  first; it draws C1 +7 A, C3 +0.34 um, A1 (+9.9, -39.2) A, so the start A1 is 36.8 A @ -0.82 rad against the true
  12.4 A @ 1.47, and the truth lies inside every grid, as ~81 % of draws would), a staged outer search bounded by the same
  intervals (stage 1 C1 x C3 7 x 7 about the start, stage 2 the A1 vector 5 x 5 at stage 1's answer; 74 fixed-probe trials
  of 50 it; round and two-fold errors are orthogonal over the aperture, C1 and C3 are not, hence the pairing), the mixed
  state built on the SEARCHED probe (`make_mixed_probe.py --base`), then refine with 8 modes. Higher orders held at their
  true values (stated). At fs 2.5 + src 0.1 on the kept joint sims; ~185 GPU-h. Pieces: `campaign/psearch_plan.py`
  (draw + grids), `campaign/select_grid.py` (one stage's minimum, the C1/C3 search's own fit, in grid-normalised
  coordinates), `campaign/search_trial.sh` (columns by name; A1 absolute), `sim/make_probe.py --set` (absolute terms, last),
  the slurm's PROBE_FROM "C1 C3 C12 phi12" and PROBE_MIXED on the searched probe, driver `PSEARCH=1` (legs tagged `_ps`).
  Why A1 matters: the old 5 % kick moved C12 by 0.6 A, the corrector's precision is ~30 A per component (~5 waves at the
  80 mrad edge); a 10 A C1 error alone drops the start's overlap with the truth to 0.18.
- **D. dose series 1e5-1e8 with the search**, at the final survivor and temperature. Last (most expensive).
**SUBMITTED 2026-10-16 17:14** (logs `logs/ceiling_20261016_171444.txt`, `logs/search_20261016_171458.txt` on Blythe; both
dry runs passed; the seed-1 draw as planned). A: fs2p5_src0p1 `_mixF` legs 1329288/89/90, pack 1329291, analysis 1329292 ·
fs2p5_src0p2 1329293/94/95, pack 1329296, analysis 1329297. B, 100 K: sims 1329298/301/304, fixed 1329299/302/305, modes
1329300/303/306, pack 1329307, analyses 1329308 (`_mixF`) / 1329309 (`_mixF_fixed`) · 20 K: sims 1329310/313/316, fixed
1329311/314/317, modes 1329312/315/318, pack 1329319, analyses 1329320 / 1329321. C: stage-1 array 1329323, select 1329324,
stage-2 array 1329325, select 1329326, refine legs 1329327/28/29, pack 1329330 (`atomfind_results_armf_a080_fs2p5_src0p1_2026*`,
no `_ps` in the name), analysis 1329331 (`analysis_armf_a080_fs2p5_src0p1_ps_*`).
Kept as limitations, with numbers: no measured S-CORR tableau (ARM pessimistic), Cc assumed, the single-Wien elliptical /
energy-position-correlated source bounded by the round scan, drift through the scan-time budget (D's minimum dose ->
dwell -> scan time against 0.8-3.3 A/s at LHe), true-probe ceiling for A/B.

**Status 2026-10-13 (the user's choice): two steps.** Step 1 = the controls only, `phase1_fixed_submit.sh`: the same ten
points, the same `_<point>_fixed` dirs, the sims KEPT -- DONE, results below ("Step 1 results"). Step 2 = the mode legs below, RECON_ONLY on those sims, at the points
whose control fails, started from the exact mixed state (to build). Why not the design as it stands: the engine seeds each
extra mode at `mode_start_pow` 0.02 (`ptycho/run_synthetic_recon_ML.m`), so an 8-mode leg starts with 0.86 in mode 1 against
a true 0.30 at 7.29 A, and every mode run so far ended with mode 1 at 0.84-0.90 (`_true4` 0.94 -> 0.84 in 500 + 500 it;
`_long_m3` / `_long_m6` 0.90); and its packs made the control wait for the 34 h legs. `phase1_submit.sh` is kept, unrun.

Supersedes the fixed-probe-only proposal in `HANDOFF_MONO300.md` ("Phase 1 as proposed"). The user's rule: probe update
with modes is the method, the fixed probe is the side control. Stub dry run of `phase1_submit.sh`: 120 sbatch calls,
`campaign/test_coh_names.sh` 26/26 first.

## A, B, C results (2026-10-20; phase images viewed: `figs/mono300/phase1_ceiling.png`; `phase1_look.py`)

All at focal spread 2.5 A, ARM tableau, 80 mrad, noiseless; right place (0.5 A) vs guessing; mode 1 = the sim's share.

| run | Pb right / guess | Ti right / guess | vacuum / slab | verdict |
|---|---|---|---|---|
| A: src 0.1, fixed (step 2) | 91 / 26 | 60 / 18 | 0.41 | good |
| **A: src 0.1, modes** | 81 / 22 | 47 / 20 | 0.36 | **good** |
| A: src 0.2, fixed (step 2) | 45 / 29 | 38 / 20 | 0.47 | off |
| **A: src 0.2, modes** | 55 / 27 | 38 / 23 | 0.35 | **okay** |
| B: src 0.1 + phonons 100 K, fixed | 29 / 30 | 22 / 18 | 0.47 | off |
| **B: src 0.1 + phonons 100 K, modes** | 51 / 27 | 19 / 17 | 0.38 | **okay** |
| B: src 0.1 + phonons 20 K, fixed | 41 / 31 | 26 / 21 | 0.46 | off |
| **B: src 0.1 + phonons 20 K, modes** | 70 / 30 | 24 / 19 | 0.36 | **okay** |
| C: src 0.1, searched start, modes | no kernels | | 0.99 | failed |

- **A**: the Mono near its best (0.1 A source) is good with modes; a realistic 0.2 A source is okay with modes where the
  fixed probe was off. Here the mode leg sits a little below its fixed control (81/47 vs 91/60) -- likely the engine's
  object step divided by the mode count (slower object convergence), not the physics. Recovered state drifts toward mode 1
  (0.45 -> 0.65 at src 0.1, 0.26 -> 0.53 at src 0.2).
- **B**: with phonons, the fixed probe has no depth at either temperature; with modes Pb comes back to okay (51 % at 100 K,
  **70 % at 20 K**) but **Ti stays at the guessing level even at 20 K** (Debye zero-point floor, Ti 0.044 A per axis).
  Images: modes turn the fixed probe's leaning streaks back into discrete Pb sites; the B-site columns stay continuous.
- **C failed, and the cause is the SEARCH ORDER, not the method.** Stage 1 searched C1 x C3 while A1 was still ~40 A
  wrong: its error surface is a shallow diagonal C1-C3 valley whose minimum the A1 error pulled to the grid corner (C1
  -991 A, true -958.6: 32.5 A off, ~5 waves at the edge; C3 690 A off). Stage 2 then found A1 well (15.8 A @ 1.51 rad vs
  12.4 @ 1.47: ~3.6 A) despite the wrong C1 -- A1 is robust to a C1 error, C1 is not to an A1 error. The refine raised the
  probe overlap from 0.038 to 0.74 but the object moved into the vacuum (the C1/depth trade), and the Pb/Ti grid legs
  found 2 and 0 grid atoms, so no kernels. Fix: a third stage, C1 x C3 again at the searched A1 (or A1 first).

## Step 1 results (2026-10-14; controls: true probe fixed, 1 mode, noiseless, own matched kernels)

Submitted 2026-10-13 ~08:08 (tarballs `*_fixed_20261013_0808*`, job ids not recorded), pulled to `~/Desktop/mono_phase1`.
Scores and sheets: `phase1_look.py` -> `phase1_scores.json`, `figs/mono300/phase1_focal.png`, `phase1_source.png` (phase
images viewed). Right place (right element within 0.5 A in depth) vs guessing; mode 1 = the sim's first-mode share.
Both references reproduce the page's numbers (coherent 90/58, `_coh` 23/7).

| point | mode 1 | Pb right / guess | Ti right / guess | wrong el. | vacuum / slab rms | verdict |
|---|---|---|---|---|---|---|
| coherent | 1.00 | 90 / 26 | 58 / 18 | 1 % | 0.31 | good |
| fs 0.75 / 1.5 / 2.5 A | 0.95 / 0.84 / 0.68 | 92 / 90 / 89 | 57 / 60 / 58 | 1-2 % | 0.31 / 0.32 / 0.39 | good |
| fs 3.5 A | 0.55 | 84 / 27 | 58 / 18 | 2 % | 0.50 | good |
| fs 5 A | 0.42 | 39 / 26 | 28 / 16 | 11 % | 0.74 | off |
| fs 7.29 A | 0.30 | 31 / 24 | 15 / 15 | 12 % | 1.03 | off |
| src 0.05 / 0.1 A | 0.87 / 0.61 | 90 / 96 | 56 / 54 | 1-2 % | 0.30 / 0.30 | good |
| src 0.2 A | 0.26 | 82 / 28 | 39 / 21 | 4 % | 0.35 | okay |
| src 0.4 A | 0.17 | 21 / 15 | 1 / 1 | 47 % | 0.48 | off |
| `_coh` (fs 7.29 + src 0.4) | 0.13 | 23 / 19 | 7 / 2 | 34 % | 1.19 | off |

- **Focal-spread wall: between 3.5 A (good) and 5 A (off), sharp.** 0.75-2.5 A look like the coherent run; 3.5 A starts to
  load the vacuum layers; from 5 A the vacuum fills and the columns become leaning streaks (the coherence artefact).
- **Source-size wall: between 0.1 A (good) and 0.4 A (off); 0.2 A keeps Pb (82 %) but loses Ti (39 %).** A different
  failure: the vacuum stays clean (0.48 vs 0.31) while the slab's contrast fades into smooth streaks, like phonons (0.39).
  The sim's 16-configuration source is the optimistic case (a real 0.2 A source needs 17 modes for 99 %, the sim's 9).
- **The old `_coh` failure was both at once**: each half fails alone. Its vacuum loading is the focal spread's.
- The two series do not collapse on the mode-1 share: src 0.2 (0.26) is okay where fs 5 (0.42) is off, so focal spread
  costs depth more per unit of lost coherence (a likely reason: defocus is what multislice trades against depth).
- **Instruments (focal spread = 17.4 A per eV x mm, before high-tension / lens ripple):** good needs dE <= 0.14 eV at
  Cc 1.4 mm (<= 0.08 eV at 2.4 mm); off by 0.21 eV (0.12 eV). ARM cold FEG 0.3 eV (7.3 A) and a Titan Wien at 0.2 eV
  (8.3 A at Cc 2.4) are past it; a Wien at its best 0.1 eV (4.2 A) or the Spectra Mono at its 0.2 eV spec with Cc 1.4
  (4.9 A) sit on it; the Spectra Mono at 90 meV (2.2-3.8 A) and UltiMono (1.2-2.1 A) are inside. Real sources at the
  specimen are 0.24-0.6 A (CEOS 2006): at or past the source wall, so the source image must be shrunk (current ~ size^2).

## Step 2 (handed over 2026-10-14): joint points + probe-mode legs from the mixed state

**SUBMITTED 2026-10-14 20:58** (log `logs/phase1_step2_20261014_205843.txt` on Blythe; both dry runs passed; all four
mode points reused their step-1 sims). A joint, fixed: fs2p5_src0p1 sims 1323755/57/59, recons 1323756/58/60, pack
1323761, analysis 1323762 · fs2p5_src0p2 1323763/65/67, 1323764/66/68, 1323769, 1323770 · fs3p5_src0p1 1323771/73/75,
1323772/74/76, 1323777, 1323778 · fs3p5_src0p2 1323779/81/83, 1323780/82/84, 1323785, 1323786. B mode legs (lab/Pb/Ti),
pack, analysis: fs5 1323787/88/89, 1323790, 1323791 · fs7p29 1323792/93/94, 1323795, 1323796 · src0p2 1323797/98/99,
1323800, 1323801 · src0p4 1323802/03/04, 1323805, 1323806. Pull to `~/Desktop/mono_phase1/{results,analysis}`.

**2026-10-15: 4 of the 12 mode legs died in their first 1-3 presolve iterations** (Pb + Ti fs5, Ti src0p2, lab src0p4:
"Convergence failed, error contains NaNs", LSQML line 224); 8 run (past the probe release and into the full engine). Not
VRAM (6.2-7.0 of 47.7 GB used) and not a node (each failed node also runs a healthy leg). The start is right: the
simulator's own defocused probe equals the writer's construction to overlap 1.000000. Likely cause, in the engine: LSQML
takes the step length from MODE 1 ALONE (`if ll == 1` -> gradient_projection_solver on chi{1}, probe{1}) and divides
the object step by the number of modes (`beta_object / par.probe_modes`). Harmless with 2 % seeds (mode 1 held 86-98 %);
with a real mixed state (mode 1 at 9-43 %) the step is fitted on a weak signal -> ill-conditioned for sparse grid objects
in small batches, stochastically (batches are shuffled). It also slows the object before the probe is released: lab fs5
mode leg 30.8 -> 18.4 over iterations 1-36 (probe fixed) against its fixed control's 14.2 -> 5.2; after the release it
falls to 2.88 at the end of the presolve, below the control's final 3.49 -- the modes do carry the coherence.
**Interpretation caveat for every mode leg (this campaign's earlier ones too, e.g. `_true4`): the object step is 1/M of
the single-mode one.** Retry: `phase1_step2_retry.sh` resubmits only the dead legs with their point's settings (one
operator per point), cancels the three original analyses and chains new ones on both packs. If a retried leg dies again,
rerun that point's three legs at a larger grouping (e.g. `16;8`: 4;2 used only ~6-7 GB of 47.7).
**Retry SUBMITTED 2026-10-15 13:35** (log `logs/phase1_step2_retry_20261015_133450.txt`): fs5 Pb/Ti 1325696/97, pack 1325698, analysis 1325699 (waits on 1325698 + 1323790) · src0p2 Ti 1325700, pack 1325701, analysis 1325702 (+ 1323800) · src0p4 lab 1325703, pack 1325704, analysis 1325705 (+ 1323805). Original analyses 1323791 / 1323801 / 1323806 cancelled. Current analyses for step 2: 1323762/70/78/86 (joint), 1323796 (fs7p29), 1325699/702/705.

One script, `phase1_step2_submit.sh` (stub dry run: 55 sbatch calls; guards: exported job variables, both dry runs, half-
deleted sims). **A, joint points** (true probe fixed, new sims kept, `_<point>_fixed`): focal spread 2.5 and 3.5 A, each
with a 0.1 and a 0.2 A source -- a real instrument has both at once. Submitted first. **B, mode legs** (`_<point>_mix`):
fs 5 A (7 modes), fs 7.29, src 0.2, src 0.4 (8 modes), probe released (PSTART 40 / 20, TEM aperture constraint), 250 it,
`4;2`, 96G, 48 h, lab + Pb + Ti, on the step-1 sims (RECON_ONLY; re-simulated in the same submission if deleted).

The start (`KICK_START=mixed` -> `PROBE_MIXED` in the job -> `sim/make_mixed_probe.py`): the top M eigenmodes of a
continuous Gaussian focal spread / source with the sim's own values (sim_meta.mat) on the true probe, N x N x 1 x M
(the live loader keeps a 4-D probe as modes; its 3-D reshape is commented out). Checked on the packed probes:
occupancies match mode_spectrum.txt's real-gun rows within its sampling noise (7.29 A: 0.314/0.258/0.186; 0.4 A: 0.090);
modes orthogonal to 1e-15, powers to 1e-13; their summed intensity matches a brute-force defocus average to its 2 %
noise (the coherent probe is 48 % off). The start holds 100 / 99.8 / 90.7 / 50.6 % of the real-gun state and 100 /
99.9 / 91.2 / 53.0 % of the sim's own 16-configuration state (fs 5 / 7.29 / src 0.2 / 0.4).
Engine changes, for multi-mode starts only (a single probe takes the old lines): the TEM aperture mask is built from the
modes' summed far-field intensity (`load_from_p.m`; per-mode masks would cut holes and a 4-D mask crashes the
`custom_data_flip` permute) -- same aperture as the true probe's to 12 edge pixels of 573 297. `kick_probe.py` reads a
multi-mode start (its strongest mode for ov_start) and records the start's mode powers. Not testable without MATLAB:
the engine run itself -- a load / init failure shows within ~30 min of a leg starting (check line in the handover).

Cost: A ~120 GPU-h; B ~400 GPU-h (1.31 x 250 x (29 + 43 M) s: 30 h at 7 modes, 34 h at 8), + 14 GPU-h per re-simulated
point; ~2-3 days wall (12 mode legs at 96G take ~156 of the 160 CPUs). Disk: the joint sims, kept (4 x 39 GB); the four
mode points' step-1 sims must stay until their legs finish (4 x 39 GB); the other six step-1 points' raw data (fs 0.75,
1.5, 2.5, 3.5, src 0.05, 0.1: 6 x 39 GB) are no longer needed.
What to judge: the mode legs against their fixed controls (step 1 above) on honest scores and the phase images; the
recovered mode powers against the start's (kick/*.json `start_mode_power`): if they collapsed back toward one mode, the
reconstruction left the mixed state -- say so before blaming the physics.

## Step 2 results (2026-10-16; joint points complete, modes partial; phase images viewed)

`phase1_look.py` (now step 1 + 2) -> `phase1_scores.json`, `figs/mono300/phase1_joint.png`, `phase1_modes.png`; mode 1 =
the sim's 16-configuration first-mode share (`mode_spectrum.txt`, joint rows added; the earlier rows unchanged).

| run | mode 1 | Pb right / guess | Ti right / guess | vacuum / slab | verdict |
|---|---|---|---|---|---|
| fs 2.5 + src 0.1, fixed | 0.43 | 91 / 26 | 60 / 18 | 0.41 | good |
| fs 2.5 + src 0.2, fixed | 0.26 | 45 / 29 | 38 / 20 | 0.47 | off |
| fs 3.5 + src 0.1, fixed | 0.36 | 72 / 27 | 42 / 19 | 0.54 | okay |
| fs 3.5 + src 0.2, fixed | 0.23 | 43 / 31 | 40 / 20 | 0.61 | off |
| fs 7.29, fixed (step 1) | 0.30 | 31 / 24 | 15 / 15 | 1.03 | off |
| **fs 7.29, probe modes** | 0.30 | **70 / 28** | 32 / 20 | **0.49** | **okay** |
| src 0.2, probe modes (lab only) | 0.26 | no Ti kernel | | 0.34 | not scored |

- **The two effects compound.** fs 2.5 A is good and src 0.2 A okay on their own; together off (Pb 45 %). A 0.1 A source
  keeps fs 2.5 good and drops fs 3.5 to okay. Joint points fail like the source alone (smooth leaning streaks, vacuum
  fairly clean), with the focal spread's vacuum loading growing on top at 3.5 A.
- **Probe modes move the focal-spread wall.** At 7.29 A the mode leg is okay (Pb 70 vs 28 guessing) where its fixed control
  is off (31 vs 24); the vacuum loading halves (1.03 -> 0.49) and the columns show discrete Pb sites again. Ti 32 % is
  still short of "good" (45 %). The recovered state stayed mixed: mode powers 0.32/0.26/0.19/0.12 at the start ->
  0.43/0.19/0.15/0.11 (truth 0.30/0.25/0.19/0.12). eps 0.0109 against the fixed control's 0.0286.
- src 0.2 A with modes: the lab leg looks cleaner than its control (vacuum 0.34, eps 0.0097 vs 0.0197) and drifted
  further toward one mode (0.31 -> 0.58, truth 0.26); its Ti kernel leg died, so no atoms yet.
- **The retry repeated the failures**: Pb fs5 (3 it), Ti src0p2 (6 it), lab src0p4 (0 it) NaN'd again, so those legs are
  not stochastic but close to deterministic; Ti fs5 (1325697) runs. The shuffled-batch explanation was wrong for them; the
  ill-conditioned mode-1 step stands.
- **Next, handed over 2026-10-16: `phase1_step2_rerun16.sh`** -- all three legs of fs 5, src 0.2, src 0.4 at grouping
  `16;8` (RTAG `_mix16`, new dirs; the 4;2 legs stay), cancels the 4;2 Ti fs5 retry 1325697, and reads the finished 4;2
  lab legs (fs 5, src 0.2) with their step-1 fixed-probe kernels as a labelled, provisional fallback
  (`_mix_kern-armf_a080_<pt>_fixed`). Stub dry run 17 calls; ~280 GPU-h. A NaN at 16;8 would show within minutes.
- **16;8 was WORSE: 7 of 9 legs NaN'd at their start** (only lab fs5 and Ti src0p4 ran). Bigger batches -> bigger
  summed updates, so the cause is a step that is too LARGE, not ill-conditioned. **Root cause found (an engine bug):**
  LSQML builds the object's illumination -- the preconditioner the object update is divided by -- from probe mode 1 ALONE
  (`init_solver.m`, `ptycho_solver.m`: `aprobe2 = |probe{1}|^2`), while the update sums all modes. A real mixed state
  (mode 1 at 9-43 %) therefore got steps ~1/w1 too large on average and far more where mode 1 is weak: lab src0p4
  (w1 0.09) died at iteration 0 every time. **Fixed 2026-10-16**: the illumination is the sum over all modes (a one-mode
  run is unchanged to the bit; the loop does not run). The object step stays divided by the mode count (an engine
  heuristic, roughly right for strongly mixed states), so mode legs converge their object more slowly than single-mode
  ones -- read their error traces with that in mind. Every earlier mode leg (incl. `_mix` fs 7.29 and `_true4`) ran the
  old normalisation.
- **Handed over: `phase1_step2_rerunF.sh`** -- all four mode points (fs 7.29, 5, src 0.2, 0.4), every leg, the original
  `4;2`, 250 it, RTAG `_mixF`, on the kept sims; refuses unless the fix is pulled; cancels queued `_mix16` legs.
  Stub dry run 20 calls; ~400 GPU-h. The `_mix` fs 7.29 result (okay) is kept as the old-engine comparison.

## Design

armf_a080 (ARM tableau, 80 mrad, 210 A region, BIN 2, NL 18), noiseless, phonons off, 16 coherence configurations.
Each point: lab, Pb and Ti sims with the same coherence setting, then per sim two legs:

- `_fixed`: true probe fixed, 1 mode, 200 it, `32;16` (the control).
- mode leg (no tag): TRUE start (`KICK=1.0`, C1 -958.6 A, C3 -40.66 um), probe released (PSTART 40/20, TEM aperture
  constraint), M modes, 250 it, grouping G, 48 h.

| point | dirs | M | G | 1 mode holds (sim / real gun) | M modes hold (sim / real gun) |
|---|---|---|---|---|---|
| focal spread 7.29 A rms | `_fs7p29` | 8 | 4;2 | 30 / 32 % | 100 / 99.6 % |
| 5 | `_fs5` | 7 | 4;2 | 42 / 44 % | 100 / 99.9 % |
| 3.5 | `_fs3p5` | 5 | 4;2 | 55 / 57 % | 100 / 99.9 % |
| 2.5 | `_fs2p5` | 4 | 8;4 | 68 / 70 % | 100 / 99.9 % |
| 1.5 | `_fs1p5` | 4 | 8;4 | 84 / 85 % | 100 / 100 % |
| 0.75 | `_fs0p75` | 3 | 8;4 | 95 / 96 % | 100 / 100 % |
| source 0.05 A FWHM | `_src0p05` | 4 | 8;4 | 87 / 87 % | 100 / 100 % |
| 0.1 | `_src0p1` | 8 | 4;2 | 61 / 63 % | 100 / 99.6 % |
| 0.2 | `_src0p2` | 8 | 4;2 | 26 / 30 % | 97.6 / 90.2 % |
| 0.4 | `_src0p4` | 8 | 4;2 | 16.5 / 9.9 % | 78 / 53 % |

M = the modes holding 99 % of the real-gun state, plus one, capped at 8 (8 modes have never run; at `4;2` they
put the same waves in flight as `_true4`'s 4 modes at `8;4`, so the GPU memory should fit). "sim" = the sim's own 16 configurations, "real gun" = a continuous Gaussian (1000
samples). Source: `mode_spectrum.py` -> `mode_spectrum.txt` (exact eigen-spectrum of the mixture; the common aberration
phase cannot change it, so only aperture, wavelength and offsets matter). The source scan dropped 0.8 A (8 modes hold 62 % /
20 %) and gained 0.05 and 0.1 A: at 80 mrad the probe's finest detail is ~lambda/alpha = 0.25 A (the probe itself is
~32 A across, d50), and a shift of that size stops it overlapping itself, so even a 0.2 A source leaves a quarter of the
state in one mode. Real sources at the specimen are 0.24-0.6 A (CEOS 2006), i.e. the 0.2 and 0.4 A points.

## What to watch

- **Does the reconstruction find the mixed state?** Compare its mode occupancies with `mode_spectrum.txt`. `_true4`
  (full physics, 4 modes) kept 0.84 in mode 1, while the true probe state there (the old `_coh`, 7.29 A + 0.4 A) has only
  0.13 in mode 1 and 4 modes can hold at most 42 %: it did not find it.
- **Source size in the sim is optimistic.** 16 configurations make a rank-16 state; a real 0.4 A source needs 35 modes for
  95 % (the sim's needs 12). Focal spread is faithful (stratified quantiles; sim and real agree to a few %).
- `_fs7p29` (fixed) and `_src0p4` (fixed) split the failed `_coh` run exactly into its halves.

## Cost (estimates)

Run time of a mode leg ~ 1.31 x N x (29 + 43 M) s, fitted to `_true4` (4 modes, 500 + 500 it: 8.9 h presolve + 27.8 h
full engine, logs in ~/Desktop/arm_true4) against the fixed legs (5.3 h at 200 it); the 6-mode 250 it `4;2` long kick
also fitted in 48 h. 8 modes ~34 h, 7 ~30, 5 ~22, 4 ~18, 3 ~14. Total ~1070 GPU-h (sims 138, fixed 159, modes 770);
~4-5 days wall (CPU quota: a recon's 73G ~ 10 of 160 CPUs). Raw data 39 GB a point, kept until the pack, which waits on
the slowest leg: up to ~390 GB on the share at once.

## Outputs

Per point: pack `atomfind_results_armf_a080_<point>_<ts>.tgz`; analyses `analysis_armf_a080_<point>_<ts>.tgz` (mode legs)
and `analysis_armf_a080_<point>_fixed_<ts>.tgz` (control). Pull to `~/Desktop/mono_phase1/{results,analysis}`.
