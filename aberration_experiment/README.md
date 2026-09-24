# Aberration experiment

Can multislice electron **ptychography** recover a real microscope's residual aberrations well
enough to keep depth resolution as you open the aperture far past a Cs-corrector's spec — and can
**atomfind** then localise the atoms in depth? Thin PTO/STO slab, 300 keV, PtychoShelves GPU
LSQ-ML. This folder is the experiment's **home for outputs**; the code lives with the pipeline.

**Start here if you are new**: the logbook, https://claude.ai/artifact/7ve93UM6yqCmiRcbfNuiJM — the permanent checkpoint record
(what was done, the numbers, the rules, the figure conventions). Source: `page/logbook.html`.

## Where things live
- **Sweep driver / planner**: [`../campaign/`](../campaign/) — `plan_probe.py` → `round_sweep.tsv`
  (the C3-coarse + C1-fine balance to a 4 Å probe per α), `run_campaign.sh` (perfect/known/fitprobe
  legs), `run_thin_atomfind.sh` (full-box recons + matched aberrated PSFs for atomfind), `README.md`.
- **This page's source**: [`page/`](page/) — `logbook.html` + `build_page.py` (figures in as data URIs).
- **Analysis scripts**: [`../analysis/`](../analysis/) — `analyze_thin_campaign.py` (depth x–z,
  in-plane, probe overlap), `make_ronchigram_fig.py` (the χ / Ronchigram evolution figure),
  `atomfind/` (the finder + its `thin` preset).
- **Sim + recon**: `../sim/simulate_4dstem.py`, `../ptycho/run_synthetic_recon_ML.m`.

## Outputs, ordered by week (ISO `YYYY-Www`)
- `figs/<week>/` — figures committed to git (small PNG/PDF).
- `results/<week>/` — small analysis products (json/csv/npy summaries, extracted PSF kernels).
- **Big raw data stays out of git** (see `.gitignore`): recon `*.tgz` tarballs and multi-GB
  volumes live on Blythe (`$SHARE`) and are pulled to a scratch dir when analysed, not committed.

## Status / key findings (2026-09-24)
- **Known probe, round aberrations: works.** Depth error falls as the aperture opens (0.56 → 0.37 Å from 70 to
  90 mrad with matched kernels); 110 mrad is out of reach because its probe outgrows the scan field.
- **Known probe, non-round aberrations: costs nothing.** Six-fold astigmatism at 0.45 waves reconstructs like the
  round leg at 70 and 90 mrad (ladder step 4), once two probe-orientation bugs in the engine were fixed.
- **Focus** is recovered at no cost by an outer search; **shot noise** holds to 10⁵ e/Å²; **phonons** are costly;
  **the 70 Å sample** fails. Blind probe retrieval (probe updated by the solver) does not work on this thin slab.
- **Next**: a representative non-round corrector tableau on the round sweep. With a hexapole corrector's
  six-fold A5 ≈ 1 mm the probe outgrows the simulation box from 70 mrad, so the corrector class is the open choice.
- Front door with the details and what changed recently: `HANDOFF_ANALYSIS.md`.
