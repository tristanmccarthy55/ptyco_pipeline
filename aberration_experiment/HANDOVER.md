# Aberration experiment — handover (2026-W37, updated 2026-09-10)

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

## 2026-09-10: four bugs found and fixed (all pushed)

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

## Results so far — PROVISIONAL
Bulk recall (interior band, the honest metric), after fixes 2–4, `~/Desktop/thin_ab_af/out/atomfind_a*`:

| α | dz (Å) | prec | Pb | Ti | O | z-RMS | kernel |
|---|---|---|---|---|---|---|---|
| 50 | 3.93 | 0.75 | 68% | 71% | 65% | 0.97 Å | matched grid |
| 70 | 1.97 | 0.58 | — | — | — | — | recon degenerate, see below |
| 90 | 1.20 | 0.73 | 92% | 90% | 88% | 0.72 Å | **data-derived** |
| 100 | 0.98 | 0.77 | **96%** | 88% | 84% | 0.84 Å | **data-derived** |

Figure: `figs/2026-W37/atomfind_depth_vs_alpha.png`; CSV in `results/2026-W37/`; collator
`analysis/collate_atomfind_depth.py <run dirs>` (reads each run's `report.json`).

**Caveats to carry forward.** (a) a90/a100 used the DATA-DERIVED PSF because the grid kernels were
junk — label them as such until S1 lands. (b) a70's lab recon is degenerate (`|obj|` collapses to
8e-4 in patches → inf/NaN in the v1-spike baseline; it crashes there *after* the v3 finder succeeds),
so a70 is excluded from the figure. (c) Measured z-RMS is 0.6–1.0 Å across the sweep, well BELOW
δz = λ/α² (7.9→2.0 Å): fitting atom centres beats the resolution limit, so the α dependence shows up
in **recall**, not in localisation error. Frame the result that way.

## IN FLIGHT (submitted 2026-09-10 ~21:09)
- **1272524** — a70 lab recon re-run damped (`BETA_LSQ=0.05`), fixes caveat (b).
- **1272533–1272548** — S1 kernel rebuild: 8 grid sims + 8 recons,
  `ALPHAS="50 70 90 100" MODES="Pb Ti" GRIDSP=3 WIN=20 OVERWRITE=1 BETA_LSQ=0.05`.
  All four α rebuilt so the sweep is internally uniform.
- **1272549** — pack → `$SHARE/atomfind_results_20260910_2109.tgz`. NOTE: it depends only on the
  kernel jobs, so if 1272524 hasn't finished the a70 lab h5 won't be in the tarball — check, and
  re-run `campaign/pack_results.sh af <out.tgz> campaign/round_sweep.tsv` if needed.
- Cancelled 1272511 (old-geometry a100_Ti) — it collided with S1's `recon_af_a100_Ti_NL28`.

## Next steps
1. scp the tarball to `~/Desktop/thin_ab_af/` (fresh dir), `tar xzf` separately.
2. Check every leg produced `analysis/*_recons.h5`; a COMPLETED slurm job with no h5 = MATLAB
   crashed internally (grep the `slurm_*.out` for `contains NaNs`).
3. Re-extract all kernels: `extract_psf.py recon_af_a<A>_<El>_NL<NL> <El>_a<A> --out psf --zdrop R`
   with **`R = round(4.0/dz)` = 1/2/3/4** at a50/70/90/100 (the default 12 erases a NL=7 volume).
   Judge them by `PSF_KERNELS.md` § "How to tell a kernel is good".
4. If S1 kernels are clean → re-run all four α with matched `--single-atom-vol` + `--ti-kernel-vol`;
   if still speckled → S2 (`GRIDSP=2.5`), then S3 (multi-plane, α≥90).
5. Rebuild the depth-vs-α figure with a70 restored and matched kernels; drop the data-derived flag.
6. Then: non-round aberrations (`nonround_sweep.tsv`, 70 mrad + C56), then a C5-corrected comparison
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
