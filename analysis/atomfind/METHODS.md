# Raw MEP phase → atom coordinates with error bars: the algorithm

*(Assumes you know what a PSF, multislice electron ptychography, and Richardson–Lucy are.)*

> Entry point: [README.md](README.md) (what / run / modules). Measured numbers and open
> items: [RESULTS.md](RESULTS.md). Kernel provenance: [docs/history/](docs/history/).

## The model we invert
The reconstructed phase volume is, to a good approximation, a sparse sum of identical blobs:

```
phase(r) ≈ Σ_i  β_i · K(r − r_i)  +  background
```

where **K is the system PSF — measured, not assumed**: a single atom simulated and
reconstructed through the byte-identical sim→PtychoShelves pipeline. K is violently
anisotropic: FWHM ≈ 0.1 Å in-plane (the 100 mrad probe) but 2–4 Å along the beam (the
missing cone), with an hourglass cross-section. We verified K is the right forward model:
a comb of K's placed at the ground-truth depths reproduces a real B–O column's z-profile
at corr 0.94 — the "3–5 px axial blur" is nothing but K-overlap of atoms 1.9–3.9 Å apart.

Atom finding = estimating {r_i, β_i}. Free-form deconvolution (RL) sharpens the volume but
spreads the ill-posedness everywhere; **fitting K directly** (sparse deconvolution — the
CLEAN / SMLM branch of the deconvolution literature) is the statistically efficient use of
the same information. We do both: RL for the human eye, model fitting for the numbers.
**We measured the field-standard alternative** rather than asserting it. The established
STEM depth-sectioning workflow is 3-D deconvolution of the volume by RL **and MEM**
followed by peak analysis (Ishizuka, Ishizuka, Ishikawa, Shibata, Ikuhara, Hashiguchi &
Sagawa, *Microscopy* **70**, 241 (2021) — where MEM was the stronger routine). We ran
both with our measured PSF, each at its best-effort threshold: RL+peaks → Pb 93 / Ti 87 /
**O 58 %**; MEM (Gull–Daniell iteration, converged χ²) + peaks → Pb 93 / Ti 88 /
**O 59 %** — and peak-picking the *raw* volume does slightly better still (**O 71 %**).
Deconvolution compresses dynamic range into the bright species (MEM most aggressively:
O only appears at a 0.001 relative floor), and the missing-cone information it cannot
restore is exactly what the direct model fit exploits. None of these baselines yields
species labels or error bars. The full machinery's edge on the same volume: **O 83 %
(92 % bulk) vs 58–71 %**, precision 97 % vs 84–90 %, plus species and calibrated σ. The
deconvolve-then-detect family was built for sparse features (isolated dopants, adatoms);
on a dense lattice its limit is structural, not a tuning artifact.

## The five algorithmic stages

**1. Preprocess.** Per-layer median subtraction (vacuum → 0), subtraction of a heavy
Gaussian blur of each layer (kills the slowly-varying channeling haze that biases weak-O
amplitudes), clip ≥ 0, trim the entrance/exit artifact layers (z < 2, > 66 Å).

**2. Column seeding.** In-plane the problem is easy (K is 2 px wide): local maxima of the
depth-mean image over a ~1.4 Å non-max window → ~98 column seeds.

**3. 3-D CLEAN + joint refinement (the deconvolution).** Around each seed, take a fixed
tube (±0.5 Å in-plane, full interior depth) and run matching pursuit with unit-norm K:
correlate residual with K → take the global max as an atom → subtract β·K → non-max
suppress ±1.5 layers × ±0.2 Å → repeat to a low amplitude floor. This is CLEAN from radio
astronomy = greedy sparse deconvolution; because each atom carries its own (x, y, z), the
ferroelectric column lean costs nothing. Then 2 sweeps of **joint Gauss–Newton
refinement**: for atom i, subtract all others' model and fit (β, x, y, z) on its local
patch by linearizing K (∂K/∂x etc. by central differences) — sub-voxel positions. Junk is
rejected by **normalized correlation of the patch with β·K** (fit quality ≥ 0.5), *not* by
amplitude — an amplitude cut cannot separate faint O from junk, a shape test can.

**4. Lattice-aware species + guided re-detection.** Columns self-classify by k-means on
their own amplitude statistics (per-column p75: Pb ≈ 4.1 ≫ Ti–O ≈ 1.5 ≫ pure-O ≈ 0.66,
zero band overlap). On B–O columns, Ti and O alternate every ~1.95 Å; species is assigned
by per-column amplitude with a dead zone, ambiguous atoms resolved by **local parity** —
distance of z_i − z_j (mod the fitted period) to confident neighbours; deliberately
translation-invariant because a global comb phase accumulates period error toward column
ends and flips labels there.
**Guided re-detection**, precisely: enumerate expected sites by extending ±k·p from each
found same-species atom (and, for O, ±p/2 from each Ti) — the prior comes from the
column's *own* atoms, never the ground truth; where a site is empty (> 0.8 Å from any
same-species atom), do a matched-filter search on the residual volume within tight gates
(±0.7 Å z, ±0.35 Å xy of the prediction) + GN refinement; accept only if (i) fit quality
clears an evidence bar, (ii) the amplitude is lattice-consistent (0.45–1.6× the column's
same-species median — this is what stops a guided-O slot from swallowing a dim Ti). The
lower bar (0.35 vs 0.5) is **allowed only for O**: the position prior substitutes for the
missing contrast. We tried guided Pb/Ti: the handful found were mislocated with
overconfident errors (z-coverage 9%), so heavy atoms are blind-only. Guided atoms carry a
`guided=1` flag through to the output — the lattice prior is never silently laundered.

**5. Errors.** Each GN fit yields the standard nonlinear-LS covariance
σ² = s²·diag[(JᵀJ)⁻¹]. On (near-)noiseless data this is systematic-limited and far too
optimistic, so it is combined in quadrature with **empirical floors calibrated to 68%
coverage against ground truth**: the floor per axis/species is the measured 68th
percentile of |true error| (σ_xy 0.015 Å; σ_z 0.24 / 0.27 / 0.31 Å for Pb / Ti / O),
with ×1.4 for guided atoms and ×1.4 in the exit band (z > 56 Å) where the recon is
artifact-prone. We explicitly audited **coverage** — the fraction of atoms whose true
error lies within ±1σ — per species × axis × depth band × blind/guided: every subset
sits at 68–96% at 1σ and ≥ 89% at 2σ. (We first calibrated on the median-ratio statistic;
it read 1.00 while actual coverage was ~50% — heavy tails. Coverage is the honest metric.)

## Why the xy error bars are ~10× smaller than the blobs (this is correct)
An error bar measures the uncertainty of the blob's **fitted centre**, not the blob's
**width** — the same reason single-molecule localization microscopy beats the diffraction
limit: σ_centre ≈ FWHM / SNR_eff, so a high-SNR blob localizes to a small fraction of its
own width. Expecting bars ≈ half the blob would make them the *resolution*, which is a
different quantity (the ability to separate two nearby atoms — that is what limits O under
Ti, and what CLEAN fights). Our measured ratios are consistent on both axes:
in-plane σ ≈ 0.01–0.015 Å vs K-FWHM 0.1 Å (≈ /10); axial σ ≈ 0.24–0.31 Å vs axial
FWHM 2–3 Å (≈ /10). The bars only *look* asymmetric because the PSF is ~30× wider along z.
And they are not taken on faith: the 68%-coverage audit against ground truth is exactly
the statement "the true error really is this small, verified per subset."

## Calibration infrastructure (validation only — the finder never sees it)
Scoring found atoms against the model needs the recon↔model coordinate map. It must be
**affine per axis**, not a constant shift: the working pixel size (window/N = 0.0495 Å)
differs from the physical object pixel (0.0492 Å) by 0.6%, i.e. ~2 px of drift across the
404-px field; and blob-peak calibration differs from the kernel-fit convention by ~1 px,
so the map is refined fiducial-style on matched heavy atoms. This exposed the true
in-plane accuracy (xy-RMS 0.13 → 0.035 Å) and lifted the GT-seeded O detector to AUC 0.94.

## The field-native story (verified citations)
Our method is not an import from another field — it is the merger of two things the EM
community already does, applied at scale:

**(a) The field already extracts 3-D atom positions from MEP volumes by fitting the
axial response.** Chen et al., *"Electron ptychography achieves atomic-resolution limits
set by lattice vibrations"*, **Science 372, 826 (2021)** located embedded Pr dopants in
all three dimensions from a single projection — with the dopant **depth obtained by
Gaussian fits to the phase depth-profile**. Follow-ups pushed the same idea: 3-D
localization of interstitials by MEP (arXiv:2407.18063), sub-nm depth resolution +
single-dopant visualization by tilt-coupled MEP (**Nat. Commun. 16, 2025**), and
oxygen-site quantification with multislice ptychography (Dong et al., *"Visualization of
oxygen vacancies and self-doped ligand holes in La₃Ni₂O₇₋δ"*, **Nature, 2024**). The STEM
depth-sectioning lineage did it first: Ishikawa et al., *"Three-Dimensional Location of a
Single Dopant with Atomic Precision by Aberration-Corrected STEM"*, **Nano Lett. 14, 1903
(2014)**. So "fit a localized peak model to the depth response to get z" is established
practice — for one dopant at a time.

**(a′) The field-standard "PSF deconvolution" workflow is our measured baseline, not our
method.** For depth sectioning, the established route is 3-D deconvolution (RL/MEM) of
the volume with the probe PSF — Ishizuka et al., *Microscopy* **70**, 241 (2021);
foundations: Behan et al., Phil. Trans. R. Soc. A **367**, 3825 (2009); Xin & Muller,
J. Electron Microsc. **58**, 157 (2009). We run both RL and MEM with our measured PSF and
report them in the benchmark ladder; on a dense lattice they sharpen the heavies but
cannot separate the axially-buried O (see measured numbers above) — consistent with that
literature's own scope, which targets sparse dopants/adatoms.

**(b) For many overlapping columns, the field's standard is model-based least-squares
fitting of a peak superposition.** De Backer, van den Bos, Van den Broek, Sijbers &
Van Aert, *"StatSTEM…"*, **Ultramicroscopy 171, 104 (2016)** — the accepted quantitative
route to positions/intensities in atomic-resolution images, whose defining feature is
*explicitly handling overlap between neighbouring peaks*, in 2-D. (The same school has
combined depth sectioning with model-based atom counting toward 3-D structure retrieval:
*"Depth sectioning combined with atom-counting in HAADF STEM to retrieve the 3D atomic
structure"*, Ultramicroscopy — check authors when citing.)

**What we did = (a) × (b), in 3-D, for every atom at once**: fit the whole reconstructed
volume as a superposition of single-atom responses — the response *measured* by
reconstructing an isolated atom through the identical sim+recon pipeline (essential
because the axial response is asymmetric, and its overlap is precisely what buries O
under Ti) — with candidate atoms proposed by iterative residual peak detection and
refined by local least-squares with the neighbours' contributions subtracted (the 3-D
analogue of how 2-D fitting handles overlapping columns). Per-atom uncertainties are the
least-squares covariance with floors calibrated to 68 % empirical coverage on the
simulated validation volume — precision quantification in exactly the Van Aert-school
sense, here validated against known ground truth.

Suggested methods sentences: *"Atomic coordinates were extracted by model-based fitting
of the reconstructed phase volume as a superposition of single-atom response functions,
with the anisotropic 3-D response measured by reconstructing an isolated atom through the
identical simulation and reconstruction pipeline. Candidate atoms were proposed by
iterative residual peak detection and refined by local least-squares fits with
neighbouring atoms' contributions subtracted, extending model-based quantification of
atomic-resolution images [De Backer 2016] to three dimensions; depth positions thereby
follow from fitting the axial response, as established for single-dopant depth
measurements in multislice ptychography [Chen 2021] and ADF depth sectioning [Ishikawa
2014]. Position uncertainties combine the fit covariance with floors calibrated to 68 %
empirical coverage on a simulated volume with known ground truth."*

*(Cross-field context, footnote-level only: the initialization loop is algorithmically
equivalent to matched-filter subtraction — CLEAN [Högbom 1974] — and the σ ≪ FWHM
behaviour of fitted centres is the localization-microscopy result [Thompson et al. 2002].
Richardson 1972 / Lucy 1974 for the RL comparison panels.)*

## Inputs / outputs / assumptions
**In:** the complex recon volume; measured single-atom PSFs (Pb, Ti; O's is unusable —
isolated O reconstructs below the noise floor, hence the pending in-situ vacancy-difference
kernels); pixel calibrations. Ground truth is used only to score and to calibrate floors.
**Out:** `found_atoms.extxyz`/`.csv` — element, xyz (model frame), σx/σy/σz, amplitude,
fit quality, column id, guided flag. NL70: ~1,760 atoms, precision 97%, bulk recall
Pb 96 / Ti 97 / O 92%, xy-RMS 0.035 Å, z-RMS 0.37 Å, species labels 98.8% correct.
**Assumptions:** one K for all species (shape is the imaging system's; amplitude carries
Z) — verified at corr 0.94; perovskite lattice knowledge used to *type and complete*
columns, never to place atoms (every position is fit to data behind an evidence bar);
floors and map are sim-calibrated constants — on experimental data they transfer as stated
tolerances; the entrance/exit ~4 Å is trimmed/de-weighted.
