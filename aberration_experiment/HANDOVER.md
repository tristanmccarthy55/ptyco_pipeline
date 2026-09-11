# Aberration experiment — handover (2026-W37, updated 2026-09-11)

Continuation notes; pairs with the `aberration-campaign` auto-memory. Overview in `README.md`.
PSF kernel rules live in `PSF_KERNELS.md` — read that before touching the grid legs.

## Goal
Thin PTO/STO slab, 300 keV, multislice electron ptychography (PtychoShelves GPU LSQ-ML). Model a
**Cs-corrected-at-30-mrad scope opened up**: C5 = 1 mm is a FIXED corrector residual; retune C3(Cs)+
C1(defocus) — the `campaign/round_sweep.tsv` balance — to keep the probe compact as α opens. Show
(a) depth resolution improving with α, (b) **atomfind** localising atoms in depth. ROUND only for now.

## Infra / conventions
- **Repo**: this `ptychoshelves-clean` == `origin/main` (github `tristanmccarthy55/ptyco_pipeline`).
  Blythe: `$SHARE/phucrh/ptyco_baseline/ptyco_pipeline` (= `/springbrook/share/physics/phucrh/...`),
  pulls `origin/main`. **Runs are sbatch-only — give the user commands, they run them** (no SSH access
  from here). Analysis env: `~/hyperspy-bundle/bin/python` (abtem/skimage/h5py).
- **Outputs → git, week-ordered**: figures to `aberration_experiment/figs/<ISO-week>/`, small products
  to `results/<week>/`. Big raw data (`*.tgz`, recons, `*_vol.npy`) is gitignored — Blythe `$SHARE` /
  Desktop scratch (`~/Desktop/thin_ab_af`).
- **scp** (2FA — AVOID nested `ssh $(ssh ...)`): one connection with a remote glob, then `tar xzf`
  SEPARATELY (it gets skipped when chained).
- Only uncommitted files are the user's report/atomfind **docs** (`report_v2.tex`, atomfind `*.md`,
  `*.bak_*`) — leave them; commit only code/figures, push to `origin/main`.

## What WORKS (banked)
- **Known/calibrated probe**: `ab_known ≈ perfect` across α. Depth sharpens with α →
  `figs/2026-W37/round_depth.png`. Sweep tool: `analysis/analyze_thin_campaign.py <dir>`.
- **Ronchigram/probe evolution**: `figs/2026-W37/ronchigram_evolution.png`. Probe d90 holds 4 Å to
  70 mrad, then grows (6.5/11/21/26 Å at 90/100/110/120) — probe compactness (not a flat Scherzer χ)
  is the C3-selection criterion.
- **atomfind runs end-to-end on the thin recons** and the depth trend is real (see Results).

## What's DEAD
- **Blind probe fitting**: 5 failures. Intractable on this thin weak-phase slab. Known-probe is the
  deliverable. Next goal after atomfind: blind from a BEST-GUESS probe (needs a vacuum/edge region
  in the scan, or a thicker sample — not a solver flag).

## 2026-09-10/11: bugs found and fixed (all pushed)

1. **`probe_initial_true.mat` never written** (`92d1b4b`). `run_thin_atomfind.sh` passed
   `PROBE_INITIAL=true`, but the sim writes `probe_initial_true.mat` only for a *nominal* initial
   probe (`simulate_4dstem.py:584`); with `=true` it writes the true probe AS `probe_initial.mat`
   and skips the other. `recon_job` symlinks `probe_initial.mat -> probe_initial_true.mat`, so the
   link dangled and all 12 recons span forever on "File corrupt". Now `PROBE_INITIAL=nominal`,
   matching `run_campaign.sh`'s proven ab_known leg.
   *Recovery without re-simulating*: the old run's `probe_initial.mat` already IS the true probe —
   `cp` it to `probe_initial_true.mat`, then `RECON_ONLY=1`.
2. **Wrong in-plane origin** (`889ce2c`) — the big one. The `thin` preset set `dx` for the full-field
   object but kept `X0=30/Y0=10`, the SCAN-WINDOW corner, valid only when the object covers the scan.
   These recons put a 20 Å scan on a 37 Å (753 px, BIN=4) / 54.6 Å (1109 px, BIN=2) object → off by
   8.75 / 17.25 Å, not a lattice vector, so every atom mapped to the wrong site. `X0/Y0 = None` now
   → `align.resolve_origin()` derives them from the object size (verified against `p/illum_sum`,
   whose centroid sits at the object centre).
3. **Wrong ground truth** (`889ce2c`). The cached GT is the full 18-cell/70 Å box; the campaign sims
   a 5-cell/27.5 Å slab. `make_gt_cache --thin-cells 5 --z-vacuum 4` builds the matching GT
   (`align._prepare_gt_thin`, mirrors the sim's `build_thin_sample`). Thin GT spans z 4.00–23.52 Å,
   independently confirming `trim_z_A=(4.0,23.5)`. Lives at `~/Desktop/thin_ab_af/gtdata/` — pass it
   with `--data-dir` (it shadows the packaged 70 Å cache by basename).
   *Combined effect of 2+3 at a50*: Pb recall 2→55%, xy-RMS 0.46→0.17 Å, depth corr 0.30→0.71,
   species confusion 67→17.5%.
4. **Halo included in the analysis** (`fea29ea`). The finder ran over the whole 37–55 Å object though
   only the ~20 Å scan is data-constrained. `cfg.fov_A` (thin preset: 20.0) crops to the scan field
   and shifts X0/Y0. Precision a50 0.25→0.75, a90 0.11→0.73, a100 0.11→0.77, recall held.
5. **One-blob kernel extractor** (`014d4bf`). Even from the clean S1 grids it mis-picked 3/8 kernels:
   it searched the whole object (junk scan edge), took one blob, never removed the phase-ramp
   gauge. Now it deramps, locates the interior grid atoms, and AVERAGES 25 of them. See
   `PSF_KERNELS.md` § Extraction. `--zdrop` is required now (old default 12 erased NL=7 volumes).
6. **Relative `SIM_BASE` in a hand-written sbatch** (`014d4bf`). MATLAB starts in `ptycho/` (`-sd`),
   so a relative SIM_BASE resolved there, and job 1272524 (a70 lab re-run) died in 2 min with its
   reason in the unpacked `.err`. The wrapper now anchors relative paths to the submit dir and
   **preflights** sim_meta / data_dp / data_position / probe_initial (`-e` follows links), failing to
   stdout before MATLAB. This also stops the endless "File corrupt ... Retrying" loop of bug 1 from
   ever burning walltime again. **Prefer the driver** (`run_thin_atomfind.sh`, which always passes
   absolute paths) over hand-written sbatch lines.

## PSF kernels — the rule (`PSF_KERNELS.md`, `35c6acb`)
A kernel is the matched system PSF only if the grid leg and the lab leg agree on **everything except
the object**: same per-α `NL` (7/14/23/28 at 50/70/90/100) and `dz` (27.525/NL), same box, same
probe, same engine settings. **`REGLAYER` stays 0 on every leg** — it "symmetrizes information
content between layers" (Odstrčil), i.e. low-passes the depth axis, which is both the measurement
and the entire content of a kernel; `run_fusion_hollow.m` warns anything >0.05 "destroys the signal"
and prints `regularize_layers <- MUST be 0`. An earlier `REGLAYER=0.5` proposal here was WRONG and is
recorded as rejected in the note. **Fix the object/scan, never the operator.**

Why the kernels broke at high α: NOT blob overlap (the reconstructed blob is **0.20 Å FWHM** against
a 4.01 Å probe — ptycho resolves ~20× below the probe), but an under-constrained solve — a single
atomic plane, little total scattering, 23–28 depth layers, and a 14 Å scan window where the lab legs
use 20 Å. Ladder: S1 scan identity + density → S2 denser → S3 sparse multi-plane (α≥90 only,
spacing ≳4·δz). Thinning the box is rejected: it changes NL.

**S1 worked** (`WIN=20, GRIDSP=3`, all four α, `atomfind_results_20260910_2109`): with the averaging
extractor, all 8 kernels average 25 atoms, argmax at the crop centre, peak/bg 44–1694. S2/S3 not
needed. Kernels live at `~/Desktop/thin_ab_af2/psf/psf_{Pb,Ti}_a<A>_vol.npy`.

## Results — matched byte-identical kernels (2026-09-11)
Bulk recall (interior band, the honest metric), all fixes in, S1 kernels,
`~/Desktop/thin_ab_af2/out/atomfind_a*`:

| α | dz (Å) | prec | Pb | Ti | O | z-RMS | species confusion |
|---|---|---|---|---|---|---|---|
| 50 | 3.93 | 0.75 | 58% | 78% | 57% | 0.97 Å | 10.9% ⚠ unreliable |
| 70 | 1.97 | — | — | — | — | — | lab recon pending (see In flight) |
| 90 | 1.20 | 0.76 | **98%** | **98%** | **95%** | **0.41 Å** | 0.5% |
| 100 | 0.98 | 0.79 | **100%** | **98%** | 86% | **0.46 Å** | 2.9% |

Figure: `figs/2026-W37/atomfind_depth_vs_alpha.png` (log depth axis; hollow markers = species
labels unreliable); CSV in `results/2026-W37/`; collator `analysis/collate_atomfind_depth.py`.

**How to frame it.**
- **Depth error falls with α**: z-RMS 0.97 → 0.41–0.46 Å, always well below δz = λ/α² (7.9 → 2.0 Å)
  because atomfind fits atom centres. ⚠ CORRECTION: the 2026-09-10 version of this file said z-RMS
  was flat (0.6–1.0 Å) and the α dependence showed only in recall. That was an artefact of the
  data-derived stand-in PSF; with matched kernels the depth error halves. Don't reuse the old framing.
- **Matched kernels are worth ~2× in depth error** at a90/a100 (0.72→0.41, 0.84→0.46 Å vs the
  data-derived PSF) — the byte-identity rule is not pedantry.
- **a50: atoms are found in-plane** (precision 0.75, xy-RMS 0.15 Å) **but species are not
  resolvable in depth**: at 3.93 Å slices the AO/BO₂ planes (1.95 Å apart) merge, so labels swap
  (confusion 10.9% > atomfind's 5% health threshold). Its per-species split (Ti > Pb) is label
  swapping, not physics. This is the expected "mediocre at low α" end.
- **a100 slightly behind a90** on O recall (86 vs 95%) and z-RMS (0.46 vs 0.41 Å), as the probe grows
  (d90 6.6 → 11 Å). Possibly the start of the break — ONE point, don't claim it; a110/a120 test it.

## IN FLIGHT / TO RUN
- **a70 lab recon** — the only missing point. Degenerate original (`|obj|` collapses to 8e-4 in
  patches → inf/NaN in the v1-spike baseline); the damped re-run 1272524 died on bug 6. Re-run via
  the driver (absolute paths, reuses the existing sim):
  `ALPHAS=70 MODES=lab RECON_ONLY=1 BETA_LSQ=0.05 bash campaign/run_thin_atomfind.sh`
  then pull just its h5 (`recon_af_a70_lab_NL14/analysis/S00000-00999/S00001/*_recons.h5`).
  The a70 S1 kernels are already extracted and clean.

## Next steps
1. a70: pull the lab h5 → `run_atomfind.py --preset thin --recon <h5> --dz 1.966 --data-dir
   ~/Desktop/thin_ab_af/gtdata --single-atom-vol .../psf_Pb_a70_vol.npy --ti-kernel-vol
   .../psf_Ti_a70_vol.npy --out ~/Desktop/thin_ab_af2/out/atomfind_a70`, then re-collate all four.
2. **a110** (BIN=1, NL=34; heavy — 175G, GROUPING=16, 24 h): tests the break the a100 dip hints at.
   `ALPHAS=110 bash campaign/run_thin_atomfind.sh` (lab + kernels in one go; the driver now defaults
   to the validated S1 grid, `WIN=20 GRIDSP=3`). **a120 is probably not simulable**: the planner put
   its probe at d90 ≈ 47 Å, and the campaign memory records it exceeding even the 70 Å BIN=1 box.
   Only try it after checking the sim's `[geom]` halo margins on a110.
3. MEP phase-volume figure for the report (x–z + in-plane per α, from the lab h5s).
4. Then: non-round aberrations (`nonround_sweep.tsv`, 70 mrad + C56), then a C5-corrected comparison
   ("Grand ARM 3") — and the blind-from-a-guess probe goal.

## Key facts
- atomfind run: `run_atomfind.py --preset thin --recon <lab h5> --dz 27.525/NL --data-dir <gtdata>
  --single-atom-vol <Pb.npy> --ti-kernel-vol <Ti.npy> --out <dir>`. `align.load_object` reads
  `*_recons.h5` (`reconstruction/object`), `.mat` or `.npy`. Writes `report.json` (per-species
  recall/`recall_bulk`/`z_rms_A`, precision, xy-RMS) — the collator's input.
- Driver: `run_thin_atomfind.sh` env — `ALPHAS`, `MODES` (lab/Pb/Ti), `GRIDSP`, `WIN`, `THIN`, `ZVAC`,
  `RECON_ONLY`, `OVERWRITE`, `BETA_LSQ`, `NITER`. `nl_full()` = Nyquist over the 27.525 Å box.
- Recon `ptycho/run_synthetic_recon_ML.m` env: NLAYERS, PROBE_START(+2), PROBE_SUPPORT_FFT, GROUPING,
  NITER, REGLAYER, BETA_LSQ, SIM_BASE, RESTART_DIR.
- h5: `reconstruction/object` (NL,1,Ny,Nx) complex; `.../probes`; `p/probe_initial`; `p/dx_spec`×1e10
  = 0.049 Å; `p/positions` (scan-box relative, NOT object frame); `p/illum_sum/illum_sum_0` (object
  frame — use THIS for the origin). `measurement/` is a broken external link — skip it when walking.
- Blythe: GPU 3×L40 48GB / 192GB, `--mem` 175G for BIN=1, GROUPING 16@BIN1 / 32@BIN2, none@BIN4.
