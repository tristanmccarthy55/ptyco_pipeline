# Handoff — the true-to-life run: an ARM200F-class corrector, at 300 kV

Written 2026-09-29. **Read this whole file first**, then `HANDOFF_ANALYSIS.md` (its "State" section is the CEOS sweep
that led here), `HANDOVER.md` (the experiment), `PSF_KERNELS.md` (the kernel rules). Where they disagree, this file is
newer.

---

## State (2026-09-30) — job 1 DONE; BLOCK 1 handed over (70 mrad and up in the 80 mrad geometry)

**Read this block first; it overrides the numbers further down where they differ.**

- **The user's rule (2026-09-30, [memory: dense-scan-not-big-window])**: large probes need dense scan points, not bigger
  windows, regions or smaller fields. 20 Å at 0.5 Å (1600 positions) everywhere, never cut to 10–15 Å; denser (0.4 Å,
  2500 positions) fits at 1419 px. "Loses ≤ 0.31 %" is the most any run has lost, not a limit — test it by cropping.
  Everything from 70 mrad runs in the ceosopt_a080_f20 geometry (210 Å region + its GT, BIN 2 = 105 Å, ±133 mrad,
  N 1419 px); the window loss is reported as a number. If 105 Å fails: detector crop first (add FULL_NDP beside
  PRESOLVE_NDP in run_synthetic_recon_ML.m, one sim kept with CLEANDATA=0, RECON_ONLY), region and scan last.
- **Block 1** (rows in `ceos_sweep.tsv`; plan in `results/2026-W40/block1_plan.json`), lab + Pb + Ti, CLEANDATA=1
  PACK_H5=0, analysed on Blythe with GT_REGION=210. Loss outside the 105 Å window / outside the 210 Å region:
  ARM 70 0.77 / 0.01 %, 75 3.18 / 0.01, 80 6.09 / 0.21, 90 17.35 / 4.04 (the sim itself wraps 4 %); CEOS 90 3.43 / 0.01,
  100 18.29 / 1.08; **calibration `ceosopt_a080_f20_b4`** (the working 80 mrad probe at BIN 4 = 52.5 Å, N 709) loses 5.54 %
  — about ARM 80's loss, on a probe that reconstructs at 105 Å (100/67/75 %, 0.48 Å). Est. 1–8 h per recon.
  Judge on the object and the residual, and look at the phase images (`<out>/phase/`).
- **RUNNING ON BLYTHE (2026-10-02 12:00)** -- block 2 (armf_a040..a090, kicks 40/60/80): pack 1301450, analysis 1301451
  (pulled the analysis fix first). Hail mary: sims 1301883/86/89 (16 configs x ~17.5 min), noise 1301884/87/90, kicked
  recons 1301885/88/91, pack 1301892, analysis 1301893. Block-1 pack 1296796 was HELD -- check it was released.
  (12:05 commit 15203bc came from a SECOND Claude session still open since 2026-10-01 23:36.) Status 12:15: 1296796
  COMPLETED (released); hail-mary sims still PENDING (Priority) -- no sim log yet, so ~17.5 min/config is an estimate;
  block 2 only armf_a060_kick running, the rest of a060-a090 queued. **armf_a100 NOT submitted** (Blythe had not pulled
  its row); its block re-dry-run with a stub sbatch 12:17 (sims 12 h, recons NL 28 / 175 G / 36 h / "32;16", kick C1 -1819
  C3 -58.17 um A1 x1.05 / 3 modes / "8;4" / 48 h, pack + CLEANDATA, analysis GT_REGION=210 on the captured pack id).
  armf_a100 PLANNED 2026-10-02 (d90 208, d99 374 A; 23 % outside 105 A, 7.3 % outside the 210 A region -- the sim wraps it) and HANDED OVER (block in chat; analysis tarball analysis_armf_a100_2*.tgz). Was: its plan ran locally, output aberration_experiment/results/2026-W40/plan_runs/armf_a100.*
  (the scratchpad dies with a session; write plans here). Next: append its row to ceos_sweep.tsv under the BLOCK 2 header
  (the row line in armf_a100.tsv), dry-run, submit LABELS=armf_a100 KICK_LABELS=armf_a100 GROUPING="32;16" RTIME=36:00:00.
- **HAIL MARY (2026-10-02): armf_a080 with everything real** -- PHONONS=16 PER_SPECIES=1, partial coherence (ARM200F
  cold FEG 0.3 eV FWHM -- 0.26-0.4 eV measured for this gun; STEM objective Cc 1.4 mm (JEM-ARM200cF spec) -> 7.29 A rms
  focal spread at 300 kV = 1.19 waves rms at the 80 mrad edge; source 40 pm FWHM, the cold-FEG mid estimate of Quigley et
  al. 2021), DOSES=1e7 e/A^2, and KICK_ALL=1 (lab, Pb, Ti all from C1/C3/A1 x 1.05, probe released, 4 modes, "8;4", 48 h).
  Coherence is sampled jointly with the phonons (simulate_4dstem `_scan_partially_coherent`). Dirs
  `*_ph16_coh_dose1e7`, analysis SUFFIX=_ph16_coh_dose1e7. Risk: probe update on the single-plane Pb/Ti grids may not
  hold -- if the kernels fail, block 2's known-probe armf_a080 kernels are the fallback for atomfind:
  `PSF_DIR=$(ls -d $SHARE/phucrh/analysis_armf_a040-armf_a050-armf_a060-*/psf | tail -1) PSF_TAG=armf_a080 LABELS=armf_a080
  SUFFIX=_ph16_coh_dose1e7 GT_REGION=210 bash campaign/run_analysis.sh` (2cf4204; tarball `..._kern-armf_a080_<ts>.tgz`;
  those kernels are NOT from the lab leg's operator -- say so beside any number).
  **Code check of the hail-mary chain (2026-10-02 afternoon, all local CPU tests):** the configuration average equals a
  direct abTEM reference (each phonon config with its own defocus and its own shifted CustomScan) to 1.1e-6 on the
  batched region path (batch boundaries inside the scan) and 4e-8 on the lazy path; `build_probe(pot, ddf)` differs from
  the central probe by exactly pi lambda ddf k^2 with the full armf_a080 tableau on; 0.3 eV + Cc 1.4 mm -> 7.294 A in
  sim_meta; recorded positions identical to a coherent run; the true probe written is the central coherent one (shape
  overlap 1 - 2e-13); the DOSES chain (sim -> noise afterok -> kick afterok), CLEANDATA (noiseless AND noisy dirs),
  PROBE_SCALE in the environment, make_probe's C1 convention (abTEM defocus, x 1.05) and the noise step's links
  (aberrations.json, probe_initial_true.mat) all read right. Fixed: GPU.m's 600 s wait clock counted solver time
  (503309c); sim_meta coherence_samples now 0 when coherent. **Reading the kick result with coherence on:** kick_probe
  compares mode 1 with the CENTRAL coherent probe, which is not the mixed state's dominant mode (1.19 waves rms of focal
  spread at the edge), so ov_final < 1 even for a perfect recovery -- compare with the block-2 coherent armf_a080 kick
  and read the mode powers, do not read ov_final against 1. (A mixed-state reference was not computed: blocked.)
  The user, on species: below ~70 mrad species labels are not expected to be reliable -- that is the point of the sweep.
  First kicked result: armf_a040_kick start overlap 0.985 -> recovered 0.998, mode 1 = 97 % of the power.
- **BLOCK 2 = THE ARM SET WE MOVE FORWARD WITH (the user, 2026-10-01; replaces the CEOS-approx data going forward).** ARM200F
  tableau at 300 kV in lengths, the knobs the user is SURE the ARM's CEOS tuning moves set for the smallest probe: C1, C3,
  B2 (against B4) and A1 (against S3) -- `plan_probe.py --arm --fight B2 A1` (90/100 --quick). A2, S3, A3 held at run3
  (editable on the ARM not confirmed; A2 against D4 would shrink 80 mrad further: 55 vs ~70 A). Rows `armf_a040`..`armf_a100`,
  105 A-capped geometry, 20 Å at 0.5 Å. Plus KICK legs at 40/60/80/100 (`KICK_LABELS`): probe released from C1, C3, A1 x 1.05,
  3 modes, presolve/full probe start 40/20, aperture constraint, GROUPING "8;4", 48 h -- does the update recover the true
  probe (`analysis/kick_probe.py`, `<out>/kick/`)? This opens the "probe not perfectly known" front the user wants answered
  next: the real approach is free variables with an outer search + probe update from a good-enough start.
  User on knob precision: roughly the reference's standard deviations or a bit finer (ARM 95 % intervals: C1 4 nm,
  A1 6 nm, B2 19 nm, C3 0.8 µm) -- a 5 % kick on C1/C3 is that order.
  Engine: GPU.m now gives up after 3 OOM retries / 600 s waiting (was unbounded). Main legs GROUPING "32;16" (100 mrad's
  presolve OOMs at 64), RTIME 36 h.
- **Block 1 results so far (2026-10-01, `~/Desktop/arm_block1`)**: residual ε (lab) and first look at the phase images, against
  ceosopt_a080_f20 (ε 0.013, 100/67/75 %, 0.48 Å): ARM 75 ε 0.013, clean, 100/67/71 %, 0.56 Å · ARM 80 0.013, clean (corner
  streaks), 100/67/73 %, 0.47 Å · ARM 90 0.030, lattice over 2/3, right third fringes, 100/74/53 %, 0.42 Å · CEOS 90 0.020,
  clean, 100/75/81 %, 0.38 Å · calibration (52.5 Å, 5.5 % lost) 0.039, reconstructs noisier, 100/67/71 %, 0.54 Å.
  **ARM 70 re-run (2026-10-01, fixed engine): works** -- presolve 898 px ran, eps 0.0135 (= reference), clean lattice and
  slices, 100/67/66 %, z 0.48 A, precision 0.97 (`~/Desktop/arm_block1/analysis_a070`). CEOS 100 re-run with GROUPING
  "32;16" iterates normally (full engine ~102 s/iteration at NL 28, ~6 h a leg) -- the presolve batch was the hang.
  **ARM 70 crashed** (first run): engine bug, `rescale_inputs.m` ceil(1420*(898/1420)) = 899 (fixed by round, 32c0305). **CEOS 100 hung
  17 h** in the GPU engine's unbounded OOM retry loop (~8000 "Failed due to GPU issue", presolve batch 64 at 1282 px x 28
  slices). Re-runs queued 2026-10-01 (RECON_ONLY, own packs): ARM 70 recons 1299664-66 / pack 1299667 / analysis 1299675;
  CEOS 100 with GROUPING="32;16" RTIME=36 h, recons 1299668-70 / pack 1299671 / analysis 1299676. Block-1 pack 1296796 is
  HELD (it would delete the raw data the re-runs read): `scontrol release 1296796` once both are back. Pull: the four
  tarballs `{analysis,atomfind_results}_{arm_a070,ceosopt_a100_f20}_2*.tgz` (packs carry the error traces).
  Offered, not built: a retry cap in `engines/GPU.m` so an OOM fails fast instead of burning the walltime.
- `plan_cost.py`'s big geometries (280–420 Å regions, N up to 4600, 400–900 positions) are NOT the plan; the tool stays
  for a later window discussion. `plan_probe.region_geometry` now caps the window at 105 Å; `plan_arm` reports
  `window_loss` and the region loss, never gates. The driver now exports each row's own C5 in LABELS mode.

- **Convention, resolved and tested** (`campaign/aberration_waves.py` `HAIDER`, `haider_to_abtem`, `haider_chi`;
  `campaign/test_haider.py`, 4 checks pass). The CEOS software writes the Haider complex notation (Uhlemann & Haider 1998;
  same form in arXiv 2510.01493 eq. A1 and arXiv 2603.23958 Table 1): the displayed angle is the complex coefficient's
  argument = m × azimuth (read off the data: every term's angles fill ±180°), with the sign flipped for the mixed terms,
  and **B2, S3, B4, D4 carry no 1/(n+1): abTEM C21 = 3 B2, C32 = 4 S3, C41 = 5 B4, C43 = 5 D4.** The factor is what matters:
  read naively the 60 mrad ARM probe comes out half size (d90 14 vs 28 Å); the sign alone changes it by ~8 %.
- **Consequence for the published CEOS page and the numbers below**: the page's "How the model compares with a real
  ARM200F" table had put the simulation's coefficients beside the tableau's. Corrected in page v10 (every value in CEOS
  notation): the CEOS approximation is MORE optimistic than shown on B2 (7.3 nm), S3 (0.24 µm), B4 (2 µm), D4 (8.1 µm).
  Its runs are unchanged. The "waves at 80 mrad" further down in this file were computed without the factors — the
  measurement uncertainty at 80 mrad (300 kV) is A1 ±9.6, A2 ±2.9, B2 ±4.9, S3 ±8.1, A3 ±5.0, A4 ±10, B4 ±50, D4 ±33,
  A5 ±44, C5 ±89 waves (the explorer lists every term at every aperture).
- **Planner**: `plan_probe.py --arm` (`plan_arm`): measured terms held, C5 fixed, C1/C3 by a 7×7 grid walked until its best
  point is interior + Nelder-Mead from the 3 best; `--set A1=0nm B4=10um ...`, `--scale 0.785`, `--c5`, `--quick` (90/100:
  5×5, one start). Scan 20 Å at 0.5 Å. Rows `arm_a040…arm_a090` live in `ceos_sweep.tsv`, `arm_a100` planned only.
- **The ARM plan** (d90 / d99, window): 40: 7.1 / 18 Å, 35 Å · 50: 14.2 / 33, 70 · 60: 28.0 / 65, 105 · 65: 40.2 / 90, 105 ·
  70: 45.8 / 96 · 75: 80.1 / 168 · 80: 95.6 / 206 · 90: 153 / 326 · 100: 213 / 435 Å (CEOS approx. 4.2/4.8/11.6/17.5/25/35/49).
  From 70 mrad all run in the 105 Å window (block 1). The probes are comets (B4, fourth-order axial coma, ~12 waves at 60 and 50 at 80, dominates, then D4) —
  both from the 1998 reference, not measured on this ARM.
- **What the open decisions change** (`results/2026-W40/arm_plan_variants.json`): A1 nulled: 40 mrad 7.1 → 5.0 Å, 60 mrad −5 %.
  C5 3 / 5 mm at 60: within the planner's ±1 Å. Waves-matched (×0.785): 60 −6 %, **80 −34 %** (63 vs 96 Å). **B4/D4 at CEOS's
  achievable 10 / 4 µm** (Zeltmann, Cornell probe-corrector guide §2.3.1, CEOS notation): 60 −43 % (16 Å), 70 33 Å and
  runnable, 80 56 Å.
- **Probe-size caveat**: `sizes()` measures d90 about the brightest pixel (every plan did). For these speckled probes that
  point moves with the grid: 75 mrad d90 68 Å at 1.1 α sampling vs 80 Å at 1.55 α. Centroid/axis measures and the window
  loss (`plan_probe.window_loss`, about the beam axis) do not move. `plan_cost.py` and the explorer therefore take the
  sizes from the planner's note — one source.
- **Page v10** (same URL): explorer third column "ARM200F tableau" + 90/100 tabs; §2 notation callout, three generated tables
  (`make_ceos_figs.py --figs 10` → `results/2026-W40/arm_plan_tables.html`), corrected comparison table.
- **Job-2 driver dry-run** (stub sbatch) passes for `LABELS="arm_a040 arm_a050 arm_a060 arm_a065"`; pass `C5=4e7` (the JSON
  wins anyway, but the export would log 1e7). 40 mrad has NL 4: atomfind cannot run there (phase images only).
- **Job 4 prerequisites** (not built): keep the lab sims' raw data (CLEANDATA=0 for the lab legs job 4 reuses);
  `sim/make_probe.py` full-tableau override; `run_c1_search.sh` by label/region (today: round_sweep.tsv by alpha).
  Open design point: the tableau's ± were measured at a ~26–40 mrad tilt; at 80 mrad they are tens of waves.

**Waiting on the user** (for 40–65 mrad): A1 (measured / nulled), C5 bracket, coefficients vs waves, B4/D4 source
(reference / CEOS achievable / both), and whether the job-4 lab sims are kept (CLEANDATA=0). ARM 100 (d90 213, d99 435 Å,
29 % outside 105 Å) is planned only.

---

## Why this run

The campaign asks whether depth-resolved multislice ptychography works when an older, widespread corrected column is
opened far past its design aperture with the probe known. The CEOS sweep (`ceosopt_*`) answered "yes" — but on a
tableau that assumed the operator re-tunes every adjustable term to 0.1 waves at *each* aperture. The user's own
microscope's measurements show that is unreachable: at 80 mrad it needs A2 at 1.2 nm, S3/A3 at 0.02 µm, D4 at 0.30 µm,
20–70× below what the corrector can even measure. So the `ceosopt` results are an idealised bound.

**This run models a real instrument instead**: the measured tableau of a JEOL ARM200F with its 2009 CEOS CESCOR-type
probe corrector — the user's lab microscope, and in the user's words "a very realistic, not particularly difficult to
achieve in many labs setup; a Titan would be better than this, so it is the right set of values". The simulation stays
at **300 kV** (the ARM runs at 200 kV; that is accepted). **Never** model a JEOL GRAND ARM: it is flat to π/4 out to
80–100 mrad and is not the question.

## The instrument to model

Source data: `campaign/arm200f_tableau.tsv` — hand-transcribed from the user's photographs (check the originals before
quoting a number):
- **run1–run3**: three CEOS STEM CsCorrector tableaus from the ARM200F, 200 kV, 23 Sep 2026, taken during C1/A1 tuning.
  Each value has the software's 95 % interval. They measure only up to A4.
- **manual**: Table 3 of a photographed reference, "wave aberration coefficients of the corrected microscope, measured at
  a tilt of |τ| = 26 mrad", value ± standard deviation, to fifth order. It looks like Uhlemann & Haider,
  *Ultramicroscopy* 72 (1998) 109, Table 3 — confirm with the user. It is the only source here for B4, D4, C5, A5.

**The proposed tableau** (confirm with the user before running; see "Decisions"):

| term | value | source |
|---|---|---|
| A1 | 4.15 nm / −70.5° | run3 (the last tableau) — *or a knob*, see Decisions |
| A2 | 16.4 nm / +116.5° | run3 |
| B2 | 12.8 nm / −85.9° | run3 |
| A3 | 0.633 µm / −168.9° | run3 |
| S3 | 0.436 µm / +11.8° | run3 |
| A4 | 4.93 µm / −90.9° | run3 |
| B4 | 30 µm / −141° | manual |
| D4 | 19 µm / −10° | manual |
| A5 | 1 mm / 150° | manual |
| C5 | **+4 mm** | **fixed — factory set, NOT an operator knob.** The manual measures 6 ± 4 mm (± = measurement standard deviation, not a tuning range); the user: take 3–5 mm "and say it's a well tuned microscope, to be a little generous" — i.e. a generous choice of the instrument's fixed value |
| C1, C3 | set for the smallest probe at each aperture | the knobs an operator always has (as `ceosbuilt_a080`) |

**Who can change what (the user, 2026-09-29 — be careful here):** the operator tunes C1, A1, B2, A2, C3, S3, A3 in
routine use (what the CEOS tuning measures and corrects). Everything else is fixed on the day: C5 factory set, A4 and B4
parasitic, A5 intrinsic, D4 not a routine knob. "Adjustable" in a vendor page is not "tuned to zero", and a ± in a
measured tableau is measurement uncertainty, never a range anyone can tune over. So in the planner only C1 and C3 move.

Why run3 as one snapshot rather than a mean: run-to-run the angles wander (A2: +159°, −104°, +116°) because several terms
sit at their measurement noise floor; one real, internally consistent tableau is more honest than an average of vectors.

**Coefficients or waves at 300 kV?** Recommended: keep the measured coefficients in length units. At 300 kV the same
length gives 1.27× more waves than at 200 kV (λ 1.97 vs 2.51 pm), so this is slightly pessimistic; the alternative
(scale by 0.785 to keep the waves at the tuning aperture) is defensible too. Say which you used.

**What the numbers mean at 80 mrad** (from the last session; recompute, don't trust): a coefficient's phase grows as
(α/26 mrad)^(n+1) — ×9, 29, 90, 276, 849 for orders 1–5. The ARM's own 95 % intervals become A2 ±2.3, B2 ±1.3, S3 ±1.6,
C3 ±3.5, A3 ±4, A4 ±8 waves at 80 mrad; the manual's standard deviations B4 ±8, D4 ±5, A5 ±35, C5 ±70 waves. A corrector
tuned at ~26–30 mrad leaves an 80 mrad probe with several waves of each residual, and cannot tell you that probe better
than several waves per term.

## Collected 2026-09-29 (all in `~/Desktop/ceos_figdata`; the CEOS page v9 shows them)

| run | result |
|---|---|
| `ceosbuilt_a080` — CEOS at 80 mrad, tuned terms kept at their 30 mrad values, only C1 and C3 re-set | **100 / 67 / 74 %, depth error 0.47 Å** — the same as the idealised re-tuned run (100 / 67 / 75 %, 0.48 Å). With the probe known, tuning precision made no measurable difference. |
| `ceosopt_a075_f20` | 95 / 86 / 90 %, 0.55 Å (the old 54 Å scan gave 96 / 64 / 52 %, 0.73 Å) |
| `round_a075_b7` — predicted stable before running (7.2 %/slice) | converged: 95 / 75 / 78 %, 0.60 Å |
| `ceosopt_a070_s06` — old 38 Å field, 0.594 Å step, 4096 positions | reconstructs (ε 0.085 → 0.038, clean lattice and depth sections; some field-corner artefacts): **the step, not the field**, was the failure |
| figure data for the page's figs 4–5 | in; fig 5 complete |
| depth sections, checked against the raw arrays (09-29) | the blocks down a section are the reconstruction's slices (7 of 3.9 Å at 50 mrad … 18 of 1.5 Å at 80), not a drawing or convergence fault: every lab recon reaches its error plateau by iteration ~20–30 and then wanders ±3–10 % with no trend. A real plotting bug was fixed: the image sat half a slice above the atom markers (see Traps). |
| Pb rows vs B-site rows (09-29) | a Pb column (Pb every 3.9 Å, nothing between) splits into separate atoms from ~70 mrad: the phase between two Pb atoms is ~90 % of the on-atom phase at 50 mrad, 30–55 % at 70–80 (every Pb column in the field). A B-site column (Ti and O alternating every 1.95 Å) stays a bar at 80; λ/α² reaches ~2 Å only near 100 mrad. The page's tabs now show both cuts. |

## Your jobs, in order

1. ~~Collect the waiting runs and update the CEOS page~~ — DONE 2026-09-29 (page v9,
   **https://claude.ai/artifact/VqmcbqPVGYkgqdExRxdFnb**; to change it: read it first with the Artifact tool — it was
   published from another conversation — then rebuild with `page/build_page.py --template page/ceos_2026-09-28.html`
   and publish with `url=`).
2. **Build the ARM200F-class tableau and plan it.** Add a planner mode to `campaign/plan_probe.py` that reads
   `arm200f_tableau.tsv` (run3 + manual B4/D4/A5, C5 = +4 mm), holds every measured term fixed, and optimises C1 and C3
   for the smallest d90 at each aperture (as the `ceosbuilt_a080` search did: Nelder-Mead from 3 starts). First resolve the
   **angle convention**: CEOS may report the phase of the complex coefficient (= m × azimuth) where abTEM's `phi_nm` is the
   azimuth — check against a CEOS/Uhlemann–Haider definition and against the probe abTEM builds; the round terms cannot
   tell you (see Traps). Compute probes, windows (`region_geometry`: the rule allows d99 < 93 Å for the 105 Å window at 210 Å; `ceosbuilt_a080`
   ran at d99 97.7 Å because only 0.25 % of its probe fell outside, less than round 70 loses from its 17.5 Å window),
   and write rows `arm_a0XX` into `ceos_sweep.tsv`. Put the new setups into the page's explorer (`--figs 8`) as a third
   column so the user sees the tableau and the probes **before** anything runs.
3. **Run it**: 40–80 mrad as the user chooses, the 20 Å scan at 0.5 Å (1600 positions), lab + Pb + Ti, `CLEANDATA=1
   PACK_H5=0`, analysis on Blythe with `GT_REGION=210`. Dry-run the driver with a stub `sbatch` first.
   **3b. Stretch: 90 and 100 mrad** (the user, 2026-09-29: "in case it magically works"), on the ARM tableau AND the
   CEOS-approx one (`ceosopt_a090` / `ceosopt_a100`, planned and sitting commented-out at the end of `ceos_sweep.tsv`).
   They need engineering first, so plan them, show the user the probes, and say what each costs before building:
   - CEOS-approx, re-tuned (optimistic): d90/d99 **88/136 Å at 90 mrad, 149/223 Å at 100** (80 mrad: 49/76). The ARM
     tableau (C5 4 mm, B4 30 µm, D4 19 µm fixed) will be at least as large — compute it.
   - `region_geometry` has no window between 105 Å and the whole 210 Å box, and caps N at 1424 px (`MAX_NDP`), which at
     210 Å crops the detector to ±66 mrad — inside the aperture. That is why the planner commented both rows out.
   - 90 mrad: a ~140 Å window, detector 1.2–1.66 α (108–150 mrad; 80 ran at 1.66) → N ≈ 1540–2130 px, beyond the
     engine's proven 1424. 1600 positions = 15–29 GB of data (fits the node); full-engine GROUPING ~7–14 on an L40
     (16 fits at 1426). Presolve half-N must be even. Region 210 Å is still enough (d99 + scan + 20 Å margin ≈ 176 Å).
   - 100 mrad: a ~230 Å window → a region side ≥ ~280 Å (new sim region, new GT cache via `make_gt_cache
     --region-side`), and N ≈ 2800 px at 1.2 α. 1600 positions = 50 GB of data → ~200 GB peak, over the 192 GB node:
     cut to ~900 positions (15 Å at 0.5 Å — the step is what matters, not the field) or fewer. Slow: GROUPING ~4.
   - Controls exist: round 90 through the region chain (`hia_a090_w35`: 95/85/87 %, z 0.42 Å) and the old round 100
     (`~/Desktop/thin_ab_af_final/recon_af_a100_lab_NL28`). Nyquist NL is 23 at 90 and 28 at 100.
4. **The key test for a real instrument: the probe known only as well as the corrector measures it.** Simulate with the
   true tableau; reconstruct with probes built from the tableau perturbed by its measurement uncertainty (run3's 95 %
   intervals for A1–A4, the manual's standard deviations for B4, D4, C5, A5; a few random draws). Mechanism exists:
   `campaign/run_c1_search.sh` reconstructs an existing sim with a fixed trial probe written in the job by
   `sim/make_probe.py` — which today overrides only C1/C3/C5, so add a full-tableau override (a JSON). This decides whether
   "probe known" is achievable on this microscope at 80 mrad, which is the real question behind the whole campaign.
5. **Report** on a page — ask the user whether it extends the CEOS page or is a new one.

Nothing from this goes into the relaxation ladder, the figures of record or the logbook without asking.

## Decisions to raise with the user

- **A1**: take run3's 4.15 nm (≈ 7 waves at 80 mrad at 300 kV), or treat A1 as a knob the operator sets on the day with the
  stigmator (it can be measured at any aperture)? The 95 % interval is ±6 nm, so the tableau value is noise-dominated.
- **C5**: +4 mm proposed (user said 3–5). A 3 / 5 mm bracket costs two more runs.
- **Coefficients vs waves** at 300 kV (above).
- **Which apertures** (40–80 mrad as before, or fewer) — and whether the 90/100 stretch (job 3b) runs with the first
  block or after 80 mrad has worked on the ARM tableau.

## What carries over from the CEOS sweep (do not re-litigate)

- **Large probes need a fine scan step (≤ ~0.6 Å); the field size does not matter**: a 20 Å field at 0.5 Å reconstructs
  CEOS 70 (98/86/86 %, 0.57 Å), 75 (95/86/90 %, 0.55 Å) and 80 mrad (100/67/75 %, 0.48 Å), and so does CEOS 70 at its old
  38 Å field with a 0.59 Å step; 1–2 Å steps fail. Fewer slices do not help. The 20 Å / 0.5 Å scan is the default.
- **Tuning precision did not matter at 80 mrad with the probe known**: tuned terms kept from a 30 mrad tune, only C1 and C3 re-set (`ceosbuilt_a080`), gives
  100/67/74 %, 0.47 Å — the same as the idealised re-tuned run.
- **Compact round probes at 75–80 mrad** make the presolve diverge at Nyquist slicing; a probe whose intensity changes
  ≥ 6.2 % per slice converges, ≤ 5.6 % fails (candidate rule, compact probes only).
- **17.5 Å windows work** (the 35 Å minimum was an analysis artifact and is withdrawn); window from d99.
- Numbers: `results/2026-W40/ceos_sweep_results.csv` (each row names its source file).

## How this project works

- **The user runs every cluster job.** No SSH. Give exact blocks; wait for logs. Every block ships its `scp -O` line
  (glob the stable part of the tarball name), and `tar` goes on its own line into a fresh directory.
- Everything lands under `$SHARE/phucrh`. `CLEANDATA=1`. Max walltime 48 h. Half a TB is fine. Big outputs are analysed
  **on Blythe** (`campaign/run_analysis.sh`): it writes kernels, atomfind, `summary.csv`, **phase images for every leg**
  (`phase/`) and figure data (`figdata/`); pull only those small tarballs.
- Local Python: `~/hyperspy-bundle/bin/python`. No MATLAB locally.
- **Verify against code or data before asserting anything.** Test first, then explain. Two conclusions last session were
  analysis artifacts and one explanation was wrong; the user caught the last one by looking at the images.
- **Judge on the object AND the residual.** ε = residual × N / 2e5 (N = pattern width) compares runs of different sizes.
  **Look at the phase images** (`analysis/render_phase.py`) before calling a run good or failed, and check which file an
  analysis read (`extract_psf` logs print `src :`).
- Dry-run the driver with a stub `sbatch` before handing over a block.
- Budget: a 710 px recon takes 0.3–1.4 h, a 1420 px one about 5 h, plus simulation and queue — a block is an overnight
  turnaround. One 1422 px recon holds at most ~5500 positions (it peaks at ~4× its data in host RAM).
- Commit and push (`origin/main`); leave `report/` and `analysis/atomfind/*.md` alone.
- **Figures** (the user's standing preferences): build up from the instrument (wavefront wrapped at one colour cycle per
  wave, Ronchigram, probe) to the reconstruction to the numbers; reuse the existing figure code; every figure gets a
  "how to read it" naming each line, marker and colour; image sets across conditions go in a tabbed viewer (one tab per
  aperture, big images, small found-atom rings — `make_ceos_figs.py --figs 9`); page text runs at the figures' width;
  comparison tables are judged good / okay / off with a colour AND a word; no before/after "first scan vs second scan" figures — write failure
  modes in words; one colour per species everywhere (Pb `#2a78d6`, Ti `#eb6834`, O `#1baf7a`) and never those colours
  for anything else; no number typed into a figure script. Pages inline their figures (base64) via `page/build_page.py`.
- Plain language for a professor-level reader; no buzzwords.

## Traps from the last session

| symptom | cause | rule |
|---|---|---|
| kernels "0 grid atoms" / speckle on Blythe but clean locally | `extract_psf` read a solver checkpoint (`Niter*.mat`, presolve or older run) | fixed; if it recurs, read the `src :` line |
| "50 atoms at 1.33 Å" in a kernel | noise between grid atoms counted as sites | `--min-sep` = grid spacing (analyse_sweep passes 3 Å) |
| a non-round wavefront disagrees with abTEM's probe while round terms agree | abTEM arrays are `[x, y]` (first axis x) | `make_ceos_figs.chi_waves` returns abTEM orientation; its `check_chi` guards it |
| a coefficient 10× off | the tsv/JSON is in Å (1 µm = 10⁴ Å, 1 mm = 10⁷ Å) | convert with code, never by eye |
| "tunable" read as "zero" | CEOS lists terms as adjustable without saying to what | always ask "tuned to what"; the measured tableau answers it |
| a 5776-position sim killed | float64 pattern buffer | fixed (float32); keep recons ≤ ~5500 positions at 1422 px |
| `DependencyNeverSatisfied` | the simulation failed | read `logs/af_sim_<id>.err`; cancel the recon and its pack together |
| depth-section image half a slice off its atom dots; a blank strip at the bottom; phase above the vacuum line | `imshow` extent put slice i at i·dz; `site_to_index` (dots, rings) puts it at (i+0.5)·dz | fixed 09-29 in every x–z figure: `extent=[0, nx*dx, nL*dz, 0]`, markers at `(layer+0.5)*dz`. atomfind's `z_A` is `layer*dz` (its own index frame) — never plot it against a physical depth axis |
| "the columns are bars, not converged?" | at ≤ 80 mrad a B-site column (Ti, O every 1.95 Å) cannot split; a Pb column can from ~70 | show a Pb row next to a B-site row (`--figs 9` does) and check the error trace plateau before blaming convergence |

## Where things live

`campaign/arm200f_tableau.tsv` (measured) · `campaign/ceos_sweep.tsv` (every row, with comments on why) ·
`campaign/aberration_waves.py` (`ceos_tableau`, waves per term; the tuning-assumption note) · `campaign/plan_probe.py`
(`plan_ceos`, `region_geometry`) · `campaign/run_thin_atomfind.sh` (driver; cols 11–15 per-row geometry and `nl_force`) ·
`campaign/run_analysis.sh` → `analysis/analyse_sweep.py` · `analysis/render_phase.py` · `analysis/figdata.py` ·
`analysis/make_ceos_figs.py` (figures 1–9: 8 = setup explorer, 9 = tabbed reconstructions; `RESULTS` manifest) · `aberration_experiment/page/` (page templates, builder) ·
`sim/make_probe.py` · `campaign/run_c1_search.sh`. Region GT: `~/Desktop/ceos_region_gt` (Blythe: `$SHARE/phucrh/gt_region210`).
