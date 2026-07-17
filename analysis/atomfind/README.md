# atomfind — atom finding & oxygen detection for 3-D electron ptychography

Finds atoms (Pb, Ti, and **oxygen**) in the reconstructed 3-D object of the PTO/STO
labyrinth and localises them **along the beam** — the hard axis — validated against the
ground-truth structure via the recon↔model map (`HANDOVER.md` §3.5/§11.4). Runs
**unchanged** on the current NL70 volume (0.15 Å, coherent) and on the better data coming
(0.1 Å-binned, 16-phonon TDS, Poisson dose) — switch with `--preset`.

## The problem, and the idea
In-plane is easy: the 100 mrad probe gives ~0.1 Å in-plane resolution (2 object pixels), so
columns are crisp and xy is trivially precise. **Depth is the crutch**: dz ≈ 1 Å/layer and
atoms ~4 layers apart overlap through the axial PSF into a weakly-modulated streak (the
"3–5-pixel blur"). Each reconstructed atom is the measured **system PSF** — tight in-plane,
elongated along z (missing cone) — so `recon ≈ atoms ⊛ PSF`.

We measured that PSF empirically (single-atom sim, `psf_Pb_NL70_vol.npy`, "gold" per
`PSF_SIM_RESPONSE.md`) and **validated it as the forward model**: a comb of that kernel at
the true atom depths reproduces a real column's z-profile at **corr 0.94**. So the blur is
just kernel overlap, and inverting it recovers the atoms. Two complementary machineries:

- **Blind finder (`find.py`, v3 = `find_atoms_v3`)** — the deliverable: **raw phase → an ASE
  object** of typed atoms with 3-axis error bars. Pipeline: smooth per-layer **background
  subtraction** → column detection → per-column 3-D tube **matching-pursuit CLEAN** with the
  measured system kernel + **joint Gauss-Newton refinement** (sub-voxel xyz, 1-σ errors) →
  **lattice-aware species**: columns self-classify as A / B–O / pure-O from their own
  amplitude stats (measured zero-overlap bands), B–O atoms split Ti/O by per-column amplitude
  with a dead zone + **translation-invariant local parity** (immune to comb-phase drift) →
  **guided re-detection**: empty lattice slots (predicted from the column's OWN atoms — still
  no GT) are re-fit with a species-appropriate kernel, a lower evidence bar for O only, and a
  lattice-consistency amplitude gate; guided atoms are **tagged** (`guided=1`).
  (v2 `find_atoms_v2` and v1 1-D spike kept as baselines.)
- **GT-seeded amplitude detector (`fit.py` + `validate.py`)** — the calibrated O-contrast
  harness: place the PSF at every *known* site, solve for per-site amplitude β, score β vs
  an off-lattice null (ROC/AUC). This is the honest oxygen-contrast measurement.

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
σ/amplitude/quality arrays, in the prepared-cell frame — overlays the GT VASP directly in
OVITO/VESTA), **`found_atoms.csv`** (same + error bars), `found_atoms.npy`, plus figures
`detection_overlay.png`, `z_accuracy.png` (raw vs v1-spike vs v2-CLEAN), `psf_compare.png`,
`amplitude_vs_Z.png`, `roc_oxygen.png`, `deconvolved_vol.npy`, `report.json`.

## Modules
| file | role |
|---|---|
| `config.py` | every per-volume knob (path, dx/dz, dose, PSF, finder, window). `preset()` selects a volume + its dose-matched gold kernel. |
| `align.py`  | recon↔GT map: depth registration + **per-axis AFFINE in-plane calibration** (a constant offset cannot register this data — dx = window/N differs from the physical pixel by ~0.6%, i.e. ~2 px across the field) + `refine_with_atoms` fiducial refinement on matched heavy atoms (calibration only; the finder never sees the map). |
| `psf.py`    | 3-D PSF: **empirical** (gold sim kernel, default), data-derived, synthetic. `species_kernels()` = {Pb, Ti} measured kernels; `axial_kernel()` = the v1 1-D z-response. |
| `deconv.py` | Richardson–Lucy 3-D: clip-to-zero (no pedestal), interior-trim + re-embed, compact kernel, dose-capped iters, max-growth divergence guard. |
| `find.py`   | **blind finder**: v2 = 3-D tube CLEAN + Gauss-Newton refinement + error bars + species + `export_atoms` (CSV/ASE); v1 = 1-D spike baseline. No GT. |
| `fit.py`    | GT-seeded per-site NNLS (the O amplitude detector), Gram-solved, with off-lattice null sites. |
| `validate.py` | greedy GT matching (recall / precision / **z-error RMS**) for the finder; ROC/AUC/amplitude-vs-Z/O-split for the detector. |
| `run_atomfind.py` | end-to-end driver + figures + `report.json` + printed verdict. |
| `fig_check.py` | visual sanity checks: `fig_cross_sections.png` (Pb-row + Ti-row cuts, raw AND RL-deconvolved, blind atoms + GT + yellow mis-ID rings — the O-on-blobs audit), `fig_column_overlay.png` (single sketchiest B-O column), `fig_atoms_3d.png` (pyvista, parallel-projection side view so each column reads as its own row + end-on; crop is display-only — atoms are found on the full field first). |

## Results on NL70 (0.15 Å, coherent, noiseless), gold Pb+Ti kernels

**Blind finder v3 — recall (all / BULK z 10–56 Å) + accuracy:**
| species | raw-peak | v1 spike | v2 | **v3** | **v3 bulk** | z-RMS |
|---|---|---|---|---|---|---|
| Pb (n=446) | 34 % | 91 % | 90 % | 89 % | **96 %** | 0.38 Å |
| Ti (n=371) | 24 % | 65 % | 87 % | **91 %** | **97 %** | 0.40 Å |
| O (n=1182) | 15 % | 49 % | 52 % | **83 %** | **92 %** | 0.38 Å |

Overall **precision 97 %** (up from 94), **xy-RMS 0.035 Å** (after the affine map
refinement exposed it — a ~1.6 px map bias + 0.6 % scale residual had been inflating it to
0.13 Å), z-RMS 0.38 Å. The v3 jumps:
**O bulk 58 %→92 %** (guided re-detection recovers the O absorbed under Ti spikes),
**species confusion 3.1 %→1.2 %** with Ti→Pb = 0 (lattice parity replaced global amplitude
bands; guided re-detection restricted to O — the few guided Pb/Ti were mislocated with
overconfident bars, guided-Ti z-coverage measured at 9 %, so they were removed rather than
papered over), and **error bars calibrated to 68 % COVERAGE** — the honest metric (fraction
of atoms whose true error is within ±1σ). Median-ratio calibration was tried and rejected:
it reported "perfect" while actual coverage was ~50 % (heavy error tails). Final floors:
σ_xy = 0.015 Å (post-affine-refinement); σ_z per species = 0.24 (Pb) / 0.27 (Ti) / 0.31 (O)
Å; guided ×1.4; exit band (z>56 Å) ×1.4 — audited per species, axis, depth band, and
blind/guided: every subset ≥ 68 % at 1σ and ≥ 89 % at 2σ. Guided atoms (~23 %) are tagged
in the export; blind-only numbers are in `report.json: finder.v3_blind_only` so the lattice
prior is never silently laundered into the headline.

**GT-seeded O amplitude detector (calibrated contrast, vs off-lattice null):**
Pb 100 %/AUC 1.00, Ti 93 %/0.97 · **O all AUC 0.87** (46 % @5% FPR) · O in-plane isolated
AUC 0.83 · O axial-overlap AUC 0.95. Amplitude-vs-Z is monotonic (null→O→Ti→Pb) with O on
the light-atom Z-scaling line — the fitted O amplitude is physically calibrated to Z.

## What the current data CAN and CAN'T show
**CAN:** blind atom finding → **an ASE object of typed atoms with calibrated 3-axis error
bars**; bulk recall Pb 96 / Ti 97 / **O 92 %** at 97 % precision, sub-Å in every axis;
species confusion 2 %. Per-atom oxygen IS now recovered in bulk — via the lattice prior
(guided slots), with guided atoms honestly tagged. **CAN'T (yet):** the last few % —
domain-wall/dim atoms (blind Pb misses at kinks), edge layers (z<10, >56 Å: recall 66–77 %,
entrance/exit recon artifacts; atoms at z>66 Å are inside the trimmed exit region and
structurally unreachable), and the residual 2 % Ti↔O label swaps at column ends. The
**in-situ vacancy-difference kernels** (PSF_SIM_REQUEST.md, request 2) are the next lever:
a measured in-crystal O matched filter should lift guided-O acceptance and the edge bands.

## What the better data should unlock
Higher dose → SNR (lifts contrast-limited pure-O columns); finer sampling/overlap → tighter,
better-conditioned axial PSF (less Ti→O leakage); 16-phonon TDS → more O contrast. The
pipeline + the same harness re-report the tables above, so the gain is **measured, not
assumed**. Use the dose-matched gold kernel (`psf_Pb_rev2_d1e8/d1e10`) and the capped RL
iters (keyed to dose in `config.effective_rl_iters`).

## Honesty / caveats
- The finder uses **no GT**; every score is a post-hoc match to the model via the §11 map.
- The gold kernel is validated as the forward model (comb-at-GT-z reproduces real profiles,
  corr 0.94). We use the **Pb shape for all species** (system response is element-independent;
  the element sets amplitude) — per `PSF_SIM_RESPONSE.md`; the O/Ti single-atom kernels are
  noise/artifact (isolated light atoms fall below the recon floor) and must not be used.
- RL deconvolution amplifies noise: iterations capped (15 noiseless, dose-scaled otherwise),
  interior-only, with a max-growth divergence guard. RL is the image-space view; the **CLEAN
  fit** is the quantitative detector.
- Error bars are the formal Gauss-Newton covariance combined with physical resolution floors
  (registration ~0.10 Å in-plane; axial kernel-mismatch ~0.30 Å) — on noiseless data the
  formal CRB alone is ~40× too optimistic (systematic-, not counting-, limited). The floors
  are validated: |GT error|/σ ≈ 1 per axis.
- Junk is rejected by **fit quality** (`quality_min_corr`), not amplitude — an amplitude cut
  alone traded Ti recall against precision. At the default cut, precision 94 %; raise it for a
  cleaner list (0.6 → 97 %) at some recall cost. `.csv`/`.extxyz` carry per-atom quality so
  you can re-threshold downstream.

Dependency history / kernel provenance: `PSF_SIM_REQUEST.md` (what was asked) and
`PSF_SIM_RESPONSE.md` (what was delivered + why only the Pb kernel is usable).
