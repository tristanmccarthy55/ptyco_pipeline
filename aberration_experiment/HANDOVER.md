# Aberration experiment — handover (2026-W37)

Continuation notes; pairs with the `aberration-campaign` auto-memory. Overview in `README.md`.

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
  to `results/<week>/`. `analyze_thin_campaign.py` and `make_ronchigram_fig.py` default there. Big raw
  data (`*.tgz`, recons, `*_vol.npy`) is gitignored — lives on Blythe `$SHARE` / Desktop scratch.
- **scp** (2FA — AVOID nested `ssh $(ssh ...)`, it prompts twice / hangs): one connection with a remote
  glob: `scp -O 'phucrh@blythe.scrtp.warwick.ac.uk:/springbrook/share/physics/<name>_results_*.tgz' <dir>/`
  then `tar xzf` separately (it gets skipped when chained).
- Only uncommitted files are the user's report/atomfind **docs** (`report_v2.tex`, atomfind `*.md`,
  `*.bak_*`) — leave them; commit only code/figures, push to `origin/main`.

## What WORKS (banked)
- **Known/calibrated probe**: `ab_known ≈ perfect` across α — aberrated wide-aperture data reconstructs
  fully if the probe is known. Depth sharpens with α → `figs/2026-W37/round_depth.png` (from
  `round_results_20260907`). Sweep tool: `analysis/analyze_thin_campaign.py <extracted-dir>`.
- **Ronchigram/phase evolution**: `figs/2026-W37/ronchigram_evolution.png` (`analysis/make_ronchigram_fig.py`).
  Probe d90 holds 4 Å to 70 mrad, then grows (6.5/11/21/26 Å at 90/100/110/120) — probe compactness (not
  a flat Scherzer χ) is the C3-selection criterion (ptycho recovers the phase).

## What's DEAD
- **Blind probe fitting**: 5 failures (grid-junk → noise → NaN → NaN → presolve-junk-NaN). Even the
  presolve-only fix (`Nst_probe=[40 Inf]`, commit `1c2d140`) NaN'd — the presolve recovers a junk probe
  that destabilises the full engine. Intractable on this thin weak-phase slab. **USER'S NEXT GOAL after
  known-probe atomfind works: blind from a BEST-GUESS probe** (the real-experiment case). Likely needs a
  vacuum/edge region in the scan (direct probe constraint) or a thicker sample — not a solver flag.

## CURRENT TASK — atomfind campaign (`campaign/run_thin_atomfind.sh`)
Per α (default 50 70 90 100): full-box KNOWN-probe labyrinth recon (`THIN=5` cells + `Z_VACUUM=4` +
`--recon-full-box` → surface artifacts dump into vacuum, edges clean = the interface fix) + matched Pb/Ti
GRID PSFs through the SAME aberrated probe/box. NL = Nyquist over the full 27.5 Å box. Sim flags verified.

**JAM (2026-09-09):** launched twice → 2nd batch (`1271xxx`) `DependencyNeverSatisfied` (sims hit the
overwrite-guard on existing `sim_out_af_*`). One `1270xxx` recon was running (7.5 h — full-box NL≈28 is
HEAVY; watch walltime). Fixed script (commit `ee0dc7d`): `sim_job` passes `OVERWRITE`; `RECON_ONLY=1`
reuses sims with no dependency. **Recovery** (see the message that shipped this file): diagnose the real
sim-fail cause from `logs/af_sim_*.err`, `scancel` the dead/dup batch, `rm -rf sim_out_af_* recon_af_*`,
`git pull`, re-run once. If full-box NL≈28 recons are too slow, reduce NL (coarser depth) or the box.

## Next steps (ordered)
1. Recover + finish the atomfind campaign → `atomfind_results_<ts>.tgz` on `$SHARE`.
2. scp down (glob form above) → e.g. `~/Desktop/thin_ab_atomfind/`, `tar xzf`.
3. `python analysis/atomfind/extract_psf.py recon_af_a<A>_Pb_NL<NL> Pb_a<A>` (and Ti) → `psf_{Pb,Ti}_a<A>_vol.npy`.
4. Point atomfind `thin` preset `single_atom_vol`/`ti_kernel_vol` at them; run
   `python analysis/atomfind/run_atomfind.py --preset thin --recon recon_af_a<A>_lab_NL<NL>/.../*_recons.h5 --dz <27.5/NL>`.
5. Build depth-localisation-vs-α figure → `figs/<week>/` (expect real results ≥~100 mrad; BO₂/AO past
   depth-res below that).
6. Then: blind-from-a-guess.

## Key facts
- Sim `sim/simulate_4dstem.py` flags: `--recon-full-box --z-vacuum --grid-box-z --thin-cells --aberrated
  --cs --c5 --defocus --probe-defocus --single-atom --grid-spacing --scan-window --convergence --bin-factor
  --aberrations-json`. `round_sweep.tsv` cols: label alpha c5 c3 c1 df_perf bin nl aber_json (a30 commented
  — low-conv fails even known-probe).
- Recon `ptycho/run_synthetic_recon_ML.m`: env NLAYERS, PROBE_START(+PROBE_START2), PROBE_SUPPORT_FFT
  (→p.probe_support_tem, electron aperture mask), GROUPING, NITER, REGLAYER, SIM_BASE. beam_thickness from
  sim_meta; dz = beam_thickness/NL.
- h5 (`*_recons.h5`): `reconstruction/object` (NL,1,Ny,Nx) complex; `.../probes` final probe; `p/probe_initial`;
  `p/dx_spec`×1e10 = 0.049 Å; `p/illum_sum/illum_sum_0`. store_images=0 (no tiffs — render from h5).
- Blythe: GPU 3×L40 48GB / 192GB, `--mem` 175G for BIN=1 (>32 cores fails), GROUPING 16@BIN1 / 32@BIN2.
