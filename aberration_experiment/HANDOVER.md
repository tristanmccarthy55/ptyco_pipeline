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
- **Ronchigram/probe evolution**: `figs/2026-W37/ronchigram_evolution.png` (rebuilt 2026-09-11:
  simulated Ronchigrams, wrapped χ, probe, phase-ramp profiles Δx = ∂W/∂θ, probe size vs α). Probe
  d90 held at the 4 Å TARGET to 70 mrad, then 6.6 / 11.0 / 24.5 / ~60 Å at 90/100/110/120 (converged
  140 Å box — the old figure's 21/26 Å at 110/120 came from a 30 Å box). **a120 is not simulable**:
  d99 ≈ 105 Å > the 70 Å BIN=1 window; a110 (d99 51 Å) is tight but valid. The figure reads its probes
  from `round_sweep.tsv` and its flatness metric from `plan_probe.py`, so it cannot drift from the plan.
- **Probe plan REVISED 2026-09-11 (user): the traditional recipe at low α.** Objective: d90 = 4 Å
  first, flattest Ronchigram second. Flatness = P-V of the NON-defocus aberration (C3θ⁴ + C5θ⁶ minus
  its best-fit θ²) across the aperture; flat ≤ λ/4 (π/2 rad). Where 4 Å is reachable (30–70 mrad,
  the "free" regime; smallest reachable 1.2/0.8/2.2 Å) the probe is deliberately ENLARGED by
  **defocus**, with Cs kept flat: the realistic +1 µm corrector residual if already flat, else the
  flattest 1 µm step against C5. Where it is not (≥90, "floor"), C3 and C1 both fight C5 as before.

  | α | old (Cs, C1) | new (Cs, C1) | non-defocus P-V |
  |---|---|---|---|
  | 30 | +1 µm, −50 Å | +1 µm, −50 Å (unchanged — already traditional) | 0.33 rad, flat |
  | 50 | −4 µm, 0 Å (Cs enlarged it) | **−2 µm, +38 Å** | 2.34 → **1.14 rad, flat** |
  | 70 | −4 µm, +2 Å | **−5 µm, −60 Å** | 8.23 → 6.23 rad (flattest possible; C5 dominates) |
  | ≥90 | unchanged | unchanged (floor) | 30 / 55 / 94 / 163 rad |

  a50/a70 therefore need re-simulating (all three legs); a90–a120 rows are byte-identical, so the
  running a110 stays valid. The old plan's tie-break ("least defocus, then least Cs") is what had
  made Cs do the enlarging.
- **`round_depth.png` rows a110/a120 are a RECONSTRUCTION failure, not the aberration break** — do not
  cite them. Their "diamonds" have a depth period of exactly 2 slices (1.67 Å at NL14, the layer
  grid's Nyquist; the physical AO/BO₂ period 1.95 Å is a different Fourier bin): the odd/even-layer
  ambiguity of unregularised multislice, drawn bilinearly. The aberration-FREE "perfect" leg also
  fails there, which a C5 break cannot explain. Watch for the same mode in the new a110.
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
7. **Depth registration aliased by one unit cell** (a90/a100 OFF −3.86/−3.92 Å = −c). atomfind's comb
   registration groups GT Pb atoms into columns by rounding (x, y) to a 0.5 Å grid; the labyrinth's
   wandering columns (polar displacements) fragment into ONE-atom "columns", so each comb has a
   single tooth and fits any atom along the column. Diagnosed by mapping the found atoms with the
   fitted map: 11% landed OUTSIDE the physical slab (0% at OFF = 0). Fix: `depth_register="atoms"`
   (thin preset only) scores every GT Pb atom at its own (x, y, z+OFF) inside the atomic band — no
   column grouping; OFF now +0.19/+0.12/+0.01 Å. Default stays "comb" so the published presets are
   bit-identical. The PX915 report's NL70 numbers use the comb method on the same labyrinth and
   were CHECKED: not aliased (comb +0.36 vs per-atom +0.22 Å; see Next steps 0).

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

## Results — matched kernels + corrected depth registration (2026-09-11)
All fixes in (bugs 1–7), S1 kernels, `~/Desktop/thin_ab_af2/out/atomfind_a*`
(the pre-bug-7 runs are kept in `out_comb_aliased/` for comparison):

| α | dz (Å) | prec | recall bulk Pb / Ti / O | xy-RMS | z-RMS | species confusion |
|---|---|---|---|---|---|---|
| 50 | 3.93 | 0.90 | 57 / 40 / 58% | 0.12 Å | 1.17 Å | 16.6% ⚠ — OLD Cs-enlarged probe; re-sim pending |
| 70 | 1.97 | — | — | — | — | lab recon pending (see In flight) |
| 90 | 1.20 | **0.94** | **95 / 94 / 95%** | **0.03 Å** | **0.37 Å** | 0.0% |
| 100 | 0.98 | **0.99** | **96 / 94 / 88%** | **0.03 Å** | **0.45 Å** | 3.1% |

Whole-slab recall now matches bulk (a90 94/94/95%), i.e. no surface plane is lost. Figure:
`figs/2026-W37/atomfind_depth_vs_alpha.png` (log depth axis; hollow = species labels unreliable);
CSV in `results/2026-W37/`; collator `analysis/collate_atomfind_depth.py`.

**How to frame it.**
- **Depth error falls with α**: z-RMS 1.17 → 0.37–0.45 Å, always well below δz = λ/α² (7.9 → 2.0 Å)
  because atomfind fits atom centres. ⚠ Two earlier framings are WRONG and must not be reused:
  (a) 2026-09-10: "z-RMS flat, α shows only in recall" — an artefact of the data-derived PSF;
  (b) 2026-09-11 morning: precision 0.76/0.79, xy-RMS 0.15 Å, whole-slab recall ~75% — all
  artefacts of the one-cell depth-registration alias (bug 7).
- **The phase volumes resolve depth-dependent polar displacements** (`mep_volumes.png`): at a90/a100
  each column breaks into one blob per atom and zig-zags ~0.5 Å laterally with depth. That is REAL
  labyrinth structure — in the GT, atoms within one column wander in-plane by a median 0.24 Å (Pb) /
  0.18 Å (Ti), up to 0.7–0.8 Å — and each found atom sits **0.03 Å** (under a pixel) from its own GT
  atom. So the displacements are tracked essentially exactly. At a50 the columns are unbroken
  streaks. Ties to the px915 found-atom polarisation result.
- **Matched kernels are worth ~2× in depth error** vs the data-derived PSF (like-for-like, both under
  the old registration: 0.72→0.41, 0.84→0.46 Å) — the byte-identity rule is not pedantry.
- **a50: atoms are found in-plane** (precision 0.90, xy-RMS 0.12 Å) **but species are not
  resolvable in depth**: at 3.93 Å slices the AO/BO₂ planes (1.95 Å apart) merge, so labels swap
  (confusion 16.6% > atomfind's 5% health threshold). The per-species split is label swapping,
  not physics. This is the expected "mediocre at low α" end.
- **a100 behind a90** on O (88 vs 95% bulk) and z-RMS (0.45 vs 0.37 Å), and it SURVIVES the
  registration fix, so it is not an alias artefact. It matches the probe physics: the aperture
  area whose rays land within ±2 Å halves from 45% to 25% between 90 and 100 mrad
  (`ronchigram_evolution.png`). Still one point — a110 tests the trend.
- **Not yet strictly byte-identical at a50/a90/a100**: kernels at `BETA_LSQ=0.05`, those lab recons
  at the default 0.1. A step size, not a penalty (far milder than REGLAYER), but it can move a
  200-iteration result. The driver now pins one BETA_LSQ (0.05) for all legs; close the gap with
  `ALPHAS="50 90 100" MODES=lab RECON_ONLY=1 bash campaign/run_thin_atomfind.sh`. a70/a110 are clean.

## IN FLIGHT / TO RUN
- **a110** (1273379–1273385): lab + kernels + pack, BIN=1, plan unchanged by the revision.
- **a50 + a70 — re-simulate with the revised probes** (all legs; supersedes the old a70 lab re-run
  1273358, which must be cancelled first if still running — `OVERWRITE=1` rewrites the sim it reads):
  `ALPHAS="50 70" OVERWRITE=1 bash campaign/run_thin_atomfind.sh` (driver defaults = S1 grid,
  BETA_LSQ 0.05 on every leg, so a50/a70 come out strictly byte-identical). Then re-extract their
  kernels (`--zdrop 1` / `2`) and re-run atomfind.

## Next steps
0. ~~Check the PX915 report's NL70 registration for the bug-7 alias~~ CHECKED 2026-09-11: NOT
   affected. On `NL70_coherent` (`atomfind/paper/inputs/NL70_new_vol.npy`) comb OFF = +0.36 Å vs
   per-atom OFF = +0.22 Å — agree to 0.14 Å, no one-cell alias. The 18-cell columns keep several
   teeth per fragment and that recon has no in-volume vacuum band; the alias is specific to the
   thin 5-cell full-box geometry.
1. a50/a70 (revised probes): pull the tarball, extract kernels (`extract_psf.py ... --zdrop 1` at
   a50, `2` at a70), run `run_atomfind.py --preset thin --recon <lab h5> --dz 3.932|1.966 --data-dir
   ~/Desktop/thin_ab_af/gtdata --single-atom-vol .../psf_Pb_a<A>_vol.npy --ti-kernel-vol
   .../psf_Ti_a<A>_vol.npy --out ~/Desktop/thin_ab_af2/out/atomfind_a<A>`, then re-collate and
   re-render (`collate_atomfind_depth.py`, `make_mep_volumes_fig.py`).
2. **a110** (BIN=1, NL=34; heavy — 175G, GROUPING=16, 24 h): tests the break the a100 dip hints at.
   `ALPHAS=110 bash campaign/run_thin_atomfind.sh` (lab + kernels in one go; the driver now defaults
   to the validated S1 grid, `WIN=20 GRIDSP=3`). **a120 is probably not simulable**: the planner put
   its probe at d90 ≈ 47 Å, and the campaign memory records it exceeding even the 70 Å BIN=1 box.
   Only try it after checking the sim's `[geom]` halo margins on a110.
3. ~~MEP phase-volume figure~~ DONE: `figs/2026-W37/mep_volumes.png` via
   `analysis/make_mep_volumes_fig.py --recons <dir> --atomfind <dir>/out --gt <gtdata>` (in-plane +
   A-site-row and B-site-row x–z per α, equal aspect, GT planes + atomfind atoms). Re-run it when
   a70/a110 land — α are auto-discovered from `recon_af_a*_lab_NL*`.
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
