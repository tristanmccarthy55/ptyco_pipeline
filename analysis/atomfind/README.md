# atomfind — atom finding & oxygen detection for 3-D electron ptychography

Finds atoms (Pb, Ti, and **oxygen**) in the reconstructed 3-D object of the PTO/STO
labyrinth and localises them **along the beam** — the hard axis — validated against the
ground-truth structure via the recon↔model map (`HANDOVER.md` §3.5/§11.4). Runs
**unchanged** across volumes (0.15 Å coherent NL70, 0.1 Å-binned 16-phonon dose series) —
switch with `--preset`.

**Docs:** this file is the entry point (what / run / modules). The method is in
[METHODS.md](METHODS.md); measured numbers and open items are in [RESULTS.md](RESULTS.md);
kernel provenance is archived in [docs/history/](docs/history/).

## The problem, and the idea
In-plane is easy: the 100 mrad probe gives ~0.1 Å in-plane resolution (2 object pixels), so
columns are crisp and xy is trivially precise. **Depth is the crutch**: dz ≈ 1 Å/layer and
atoms ~4 layers apart overlap through the axial PSF into a weakly-modulated streak. Each
reconstructed atom *is* the measured **system PSF** — tight in-plane, elongated along z
(missing cone) — so `recon ≈ atoms ⊛ PSF`. We measured that PSF empirically and validated it
as the forward model (comb-at-GT-z reproduces a real column's z-profile at corr 0.94), so the
blur is just kernel overlap and inverting it recovers the atoms. Two complementary machineries:

- **Blind finder (`find.py`, `find_atoms_v3`)** — the deliverable: raw phase → an ASE object
  of typed atoms with 3-axis error bars, using **no ground truth**. Background subtraction →
  column detection → per-column 3-D matching-pursuit CLEAN with the measured kernel + joint
  Gauss-Newton refinement → lattice-aware species → guided re-detection of empty lattice slots
  (tagged `guided=1`). Full walk-through in [METHODS.md](METHODS.md).
- **GT-seeded amplitude detector (`fit.py` + `validate.py`)** — the calibrated O-contrast
  harness: place the PSF at every *known* site, solve per-site amplitude β, score vs an
  off-lattice null (ROC/AUC). The honest oxygen-contrast measurement.

## Run it
Needs abtem + skimage + scipy — the **hyperspy-bundle** Python, not system `python3`:
```bash
cd ptychoshelves-clean/analysis
~/hyperspy-bundle/bin/python atomfind/run_atomfind.py                 # NL70 + gold Pb PSF (default)
~/hyperspy-bundle/bin/python atomfind/run_atomfind.py --psf all       # empirical vs data vs synthetic
~/hyperspy-bundle/bin/python atomfind/run_atomfind.py --preset reviewer2 --dose 1e8   # when it lands
```
The gold Pb + Ti kernels and preset paths are wired into `config.preset()`. Outputs →
`~/Desktop/atomfind_out/`: **`found_atoms.extxyz`** (ASE object: element + xyz + per-atom
σ/amplitude/quality, in the prepared-cell frame — overlays the GT VASP directly in
OVITO/VESTA), **`found_atoms.csv`** (same + error bars), `found_atoms.npy`, plus figures
(`detection_overlay.png`, `z_accuracy.png`, `psf_compare.png`, `amplitude_vs_Z.png`,
`roc_oxygen.png`), `deconvolved_vol.npy`, and `report.json`.

## Modules
| file | role |
|---|---|
| `config.py` | every per-volume knob (path, dx/dz, dose, PSF, finder, window). `preset()` selects a volume + its dose-matched gold kernel. |
| `align.py`  | recon↔GT map: depth registration + per-axis affine in-plane calibration + `refine_with_atoms` fiducial refinement (calibration only; the finder never sees the map). |
| `psf.py`    | 3-D PSF: empirical (gold sim kernel, default), data-derived, synthetic. `species_kernels()` = {Pb, Ti}; `axial_kernel()` = the 1-D z-response. |
| `deconv.py` | Richardson–Lucy (+ MEM) 3-D deconvolution: clip-to-zero, interior-trim + re-embed, compact kernel, dose-capped iters, divergence guard. |
| `find.py`   | **blind finder**: 3-D tube CLEAN + Gauss-Newton refinement + error bars + species + `export_atoms` (CSV/ASE). No GT. |
| `fit.py`    | GT-seeded per-site NNLS (the O amplitude detector), Gram-solved, with off-lattice null sites. |
| `validate.py` | greedy GT matching (recall / precision / z-RMS) for the finder; ROC/AUC/O-split for the detector; `health_warnings`. |
| `uncertainty.py` | model σ (joint Cramér–Rao) + split-conformal per-stratum calibration of the error bars. |
| `run_atomfind.py` | end-to-end driver + figures + `report.json` + printed verdict. |
| `fig_check.py` | visual sanity checks (cross-sections, single-column overlay, 3-D atoms). |
| `dose_series.py` | portability harness: NL105 `.mat` loader + per-dose run. |

## Results (headline)
On NL70 (0.15 Å, coherent, noiseless), gold Pb+Ti kernels: bulk recall (z 10–56 Å)
**Pb 96 / Ti 96 / O 96 %** at 95 % precision, xy-RMS 0.03 Å, z-RMS 0.37 Å, species labels
~99 % correct. The finder's edge over a competent peak-picker is **oxygen**, species ID,
precision control, and calibrated per-atom σ. Full tables, the detector-comparison ladder,
noise/dose robustness, and open items: [RESULTS.md](RESULTS.md).

## Caveats (the short list — detail in [METHODS.md](METHODS.md))
- The finder uses **no GT**; every score is a post-hoc match to the model via the §11 map.
- One kernel for all species: the **Pb shape** is used for Pb/Ti/O (system response is
  element-independent; the element sets amplitude). Isolated O/Ti single-atom kernels are
  noise and must not be used — see [docs/history/PSF_SIM_RESPONSE.md](docs/history/PSF_SIM_RESPONSE.md).
- RL/MEM deconvolution is the image-space *view*; the **CLEAN model fit is the detector**.
  RL amplifies noise (iters capped, interior-only, divergence-guarded).
- Error bars combine the Gauss-Newton covariance with resolution floors, then are conformally
  calibrated to coverage; junk is rejected by **fit quality**, not amplitude. Per-atom quality
  is exported so you can re-threshold downstream.

## Working on this code
Documentation convention (Doxygen; a `Doxyfile` at the repo root renders HTML to `docs/html/`):
- **Module docstring:** `@file` + `@brief` + a 2–4 line summary. Deep rationale goes in
  METHODS.md, measured numbers in RESULTS.md — not in inline comments.
- **Public function/class:** one-line `@brief`; add `@param`/`@return` only where the
  signature isn't self-evident. Keep the `# ---- section ----` banners.
Build the API docs: `cd ptychoshelves-clean && doxygen Doxyfile`.
