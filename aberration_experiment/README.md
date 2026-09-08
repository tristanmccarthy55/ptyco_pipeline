# Aberration experiment

Can multislice electron **ptychography** recover a real microscope's residual aberrations well
enough to keep depth resolution as you open the aperture far past a Cs-corrector's spec — and can
**atomfind** then localise the atoms in depth? Thin PTO/STO slab, 300 keV, PtychoShelves GPU
LSQ-ML. This folder is the experiment's **home for outputs**; the code lives with the pipeline.

## Where things live
- **Sweep driver / planner**: [`../campaign/`](../campaign/) — `plan_probe.py` → `round_sweep.tsv`
  (the C3-coarse + C1-fine balance to a 4 Å probe per α), `run_campaign.sh` (perfect/known/fitprobe
  legs), `run_thin_atomfind.sh` (full-box recons + matched aberrated PSFs for atomfind), `README.md`.
- **Analysis scripts**: [`../analysis/`](../analysis/) — `analyze_thin_campaign.py` (depth x–z,
  in-plane, probe overlap), `make_ronchigram_fig.py` (the χ / Ronchigram evolution figure),
  `atomfind/` (the finder + its `thin` preset).
- **Sim + recon**: `../sim/simulate_4dstem.py`, `../ptycho/run_synthetic_recon_ML.m`.

## Outputs, ordered by week (ISO `YYYY-Www`)
- `figs/<week>/` — figures committed to git (small PNG/PDF).
- `results/<week>/` — small analysis products (json/csv/npy summaries, extracted PSF kernels).
- **Big raw data stays out of git** (see `.gitignore`): recon `*.tgz` tarballs and multi-GB
  volumes live on Blythe (`$SHARE`) and are pulled to a scratch dir when analysed, not committed.

## Status / key findings so far
- **Known/calibrated probe works**: `ab_known ≈ perfect` across α (a50–a100) — aberrated wide-aperture
  data is fully usable if the probe is known; depth structure sharpens with α (δz = λ/α²).
- **Blind probe retrieval is hard** on this thin weak-phase slab (grid-junk → noise → NaN across
  fixes; presolve-only fit under test). Real fix is likely a vacuum/edge scan.
- **Where it breaks** (`figs/2026-W37/ronchigram_evolution.png`): C5 = 1 mm is a fixed corrector
  residual; retuning Cs(C3)+C1 holds the probe ~4 Å to **70 mrad**, then it grows (6.5/11/21/26 Å at
  90/100/110/120). Probe compactness — not a flat Scherzer χ — is the C3-selection criterion (ptycho
  recovers the phase).
- **atomfind**: interface fix (full-box recon) + matched aberrated PSFs built; sweep pending →
  depth-localisation-vs-α (expected real only ≥~100 mrad).
