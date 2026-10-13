# Handoff — a monochromated 300 kV microscope: where is the wall in temperature and energy spread?

Written 2026-10-12 at the close of the ARM200F campaign. **Read this whole file first**, then the newest bullets of the
State block in `HANDOFF_REALISTIC.md` (blocks 3–5: every number behind the summary below), `PSF_KERNELS.md` (kernel
rules), and `page/ARM_PAGE.md` (the results page and the honest atom scoring). Where they disagree, this file is newer.

## The goal (the user, 2026-10-12)

- **Target: a widespread, last-generation 300 kV probe-corrected microscope with a monochromator** — not a next-generation
  instrument flat to 80–100 mrad (never a JEOL GRAND ARM 2). The candidate class is Thermo Fisher's Spectra 300 with an
  X-FEG monochromator (e.g. Cambridge's, below); the older Titan Themis is the comparison (see the facts).
- **Find the wall**: the specimen temperature and the energy spread (and anything else that turns out to matter) at
  which depth-resolved multislice ptychography stops working, with the probe found from the data (outer search, then probe
  update with modes; probe update always on — memory `probe-update-always-on`).
- **Aberrations**: a measured tableau for the target if one can be found; otherwise the ARM200F tableau as a
  **pessimistic** baseline (its 2009 CEOS corrector does not correct A5; the Spectra's S-CORR does).
- **Hoped-for headline**: it works at liquid-nitrogen temperature with a monochromator.
- **Expectations set by the user**: seeing oxygen is close to impossible; judge Pb and Ti. Phonon blur along the beam is
  acceptable (from the images); the partial-coherence artefact (vacuum loading, one blob per column) is not. A slight
  specimen tilt may help with the phonon blur (not built).
- **Ruled out (the user, 2026-10-12): the Nion.** SuperSTEM's Nions run at 40–100 kV; at 100 kV the 3.1 Å depth resolution
  of 300 kV / 80 mrad needs 110 mrad (126 at 60 kV), so a low-voltage instrument is no use however good its monochromator.

## What is established (ARM200F tableau, 300 kV, 5-cell PTO/STO slab, 20 Å scan at 0.5 Å)

Atom numbers are "right place" (right element within 0.5 Å in depth) against "guessing" (the same found atoms slid along
their columns): `analysis/make_arm_page.py atom_scores`. Never judge on summary.csv recall (it accepts any element within
2 Å in depth, so a Pb column passes at any depth).

| result | numbers (Pb / Ti right, guessing) |
|---|---|
| probe known, coherent: works 70–80 mrad | 80 mrad 90 / 58 vs 26 / 18 % |
| 90 mrad works once the window holds the probe | 140 Å window 92 / 50 vs 23 / 13 %; 105 Å window 76 / 32 % |
| 90/100 mrad failures were the 105 Å window | 80 mrad probe in a 52.5 Å window fails like 100 mrad |
| outer search finds C1, C3 from realistic data | to 2.3 Å and 0.14 µm (corrector: ±40 Å, ±0.8 µm); 4× closer than a 5 % kick |
| shot noise not the cause | 1e7 = 1e8 e/Å² |
| RT phonons alone remove depth (80 mrad, noiseless, probe fixed) | 24 / 11 vs 18 / 14 %, wrong element 33 %; clean recon, smooth column streaks |
| partial coherence alone removes depth (0.3 eV through Cc 1.4 mm = 7.3 Å rms focal spread, plus a 0.4 Å source) | 23 / 7 vs 19 / 2 %, wrong element 34 %; vacuum layers loaded, one blob per column |
| phonons at 100 K (Debye model) recover part | 49 / 22 vs 29 / 19 %, wrong element 8 % |
| 70 mrad full physics is no better | 12 / 0–2 vs 7–9 / 2–3 % |
| refined probe (searched start, 4 modes, 500 it) under full physics | 19–20 / 12–22 % vs 11–13 / 8–10 %: off at both doses |

Results page (generated, republished to one link): https://claude.ai/artifact/MUkJ7yYrsek3nTLXcZgCrk (v9).

## Facts checked 2026-10-12 (sources at the end)

| instrument | generation | energy spread at 300 kV | corrector | UK example found |
|---|---|---|---|---|
| Thermo Fisher Spectra 300, X-FEG Mono | current mainstream (30–300 kV) | tunable 1 eV to <0.2 eV (datasheet); **Cambridge quotes 90 meV at 300 kV** | S-CORR, 5th order, corrects A5 | **Cambridge Wolfson EM Suite**: probe-corrected, Merlin direct-electron camera, Panther STEM detector |
| Spectra 300, X-FEG UltiMono | current, high end | <55 meV at 300 kV (<25 meV at 60 kV) | S-CORR | none found |
| Spectra 300, X-CFEG | current | 0.31–0.39 eV, no monochromator | S-CORR | — |
| FEI Titan Themis / Titan 80-300, X-FEG + Wien monochromator | previous generation, very widespread | ~0.2–0.3 eV (Leeds: "approx. 0.25 eV") — barely better than the ARM200F's cold FEG | CEOS DCOR (Themis) | Leeds Titan3 Themis 300 (probe correction not stated on its listing); Imperial monochromated Titan 80/300 |
| JEOL GRAND ARM 300F (first generation) | cold FEG, no monochromator found | ~0.3 eV | double-corrected (type not checked) | ePSIC (Diamond), Oxford |

- A search summary reports a Spectra 300 X-CFEG + EMPAD used for multislice ptychography at 300 kV (30 mrad, 0.42 Å steps);
  the paper was not identified — confirm before citing.
- **Not found**: a measured residual-aberration tableau for a Spectra 300 / S-CORR; the Cc of the Spectra objective (the
  runs so far used the ARM200cF's 1.4 mm); the probe current at 90 meV. Ask the user whether Cambridge can share a tableau
  and these numbers.
- The Wien/X-FEG sources are Schottky emitters: larger than a cold FEG, so the **source-size** term (0.4 Å FWHM so far) may
  be worse on these machines. It has never been tested on its own.

## FINAL CHECK (2026-10-12, the user: "2019 is erring on too new") — target the Titan generation, scan the focal spread

| generation | introduced | monochromator at 300 kV | probe corrector | how far it opens |
|---|---|---|---|---|
| FEI Titan 80-300 | 2005–06 | Wien filter: ~0.3 eV spec, ~0.1 eV measured at some sites | CEOS CESCOR (hexapole): intrinsic A5 ≈ 3.0 mm; total Cc 2.4 mm (objective 2.0 mm, standard 5 mm gap) | optimum ~25 mrad, up to ~35 mrad |
| FEI Titan Themis | Aug 2013 | X-FEG + Wien monochromator, ~0.2–0.3 eV (Leeds ~0.25) | CEOS DCOR: all axial aberrations to 5th order (incl. B4, D4, A4, C5), A5 < 0.2 mm | "beyond 40 mrad", 30–300 kV |
| Thermo Spectra 300 | Aug 2019 | X-FEG Mono <0.2 eV (Cambridge 90 meV); UltiMono <55 meV | S-CORR: corrects A5 and D6 | no published flat-phase angle found |

Sources: Müller, Uhlemann, Hartel, Haider, Microsc. Microanal. 12, 442 (2006) (CEOS; the Titan 80-300 CESCOR numbers; the
improved design needs A5 < 0.5 mm for 37.5 mrad at 300 kV; its monochromator gives 0.3 eV at half the current); CEOS DCOR page;
Titan Themis launch (FEI, 1 Aug 2013); Titan 80-300 facility pages.
- **None of these is flat to 80 mrad** (that is GRAND ARM 2 territory, excluded). At 80 mrad every one needs the C1/C3
  balancing used for the ARM, against its own residuals.
- **The ARM tableau is not pessimistic on every term**: against a CESCOR Titan 80-300 its A5 (1 mm, 1998 reference) is smaller
  than the Titan's 3 mm and its Cc (1.4 mm) smaller than 2.4 mm. Against a DCOR Titan Themis it is pessimistic (DCOR corrects
  B4, D4, C5; A5 < 0.2 mm).
- **Cc matters as much as the energy spread.** Focal spread at 300 kV = 17.4 Å per (eV × mm) of ΔE·Cc:
  ARM cold FEG 0.3 eV × 1.4 mm = 7.3 Å (fails) · Titan Wien mono 0.3 / 0.2 / 0.1 eV × 2.4 mm = 12.5 / 8.3 / 4.2 Å ·
  Spectra Mono 90 meV × 1.4–2.4 mm = 2.2–3.8 Å · UltiMono 50 meV = 1.2–2.1 Å · unmonochromated X-FEG ~0.7 eV × 2.4 = 29 Å.
- **So scan the FOCAL SPREAD in Å, not the energy spread** (e.g. 7.3, 5, 3.5, 2.5, 1.5, 0.75 Å rms), find the wall, and map it
  onto instruments with the table above. That needs no Cc for any particular machine and cannot cherry-pick one.
- Source size: the CEOS paper uses 45–60 pm (CESCOR example) and 24–26 pm (improved design, max current) at the specimen;
  our 40 pm sits inside that range. Monochromated probe currents are ~pA (8 pA in the CEOS example).
- Recommended target class for the write-up: **Titan Themis (2013) with DCOR and a Wien monochromator**: widespread, old
  enough, 300 kV; aberrations: the ARM tableau (pessimistic for DCOR) or a DCOR-like variant (A5 0.2 mm, B4/D4 at measurement
  noise) if no measured tableau is found.

## Is that cherry-picking? (checked 2026-10-12)

- The Spectra 300 was unveiled in August 2019 (M&M 2019), so the fleet is 2020–2025. Installs found: Cambridge
  (X-FEG Mono, 90 meV at 300 kV), Imperial (X-FEG UltiMono, <30 meV, K3 camera), Jülich ER-C (monochromated X-FEG, <0.03 eV
  at 60 kV, double-corrected, EMPAD), LSU, Notre Dame, Stanford, PolyU, Twente, Stuttgart, several in China. No public sales
  figures; some sites take the X-CFEG (no monochromator) instead.
- **90 meV is Cambridge's best quoted figure, not the guaranteed spec** (X-FEG Mono: <0.2 eV). The widespread previous
  generation (Titan Themis, Wien monochromator) gives ~0.2–0.3 eV. So do not adopt 90 meV as "typical": scan the energy spread
  and place instruments on the result afterwards — 0.3 eV (cold FEG / Titan Wien), 0.2 eV (Mono spec), ~0.1 eV (Mono at its
  best, Cambridge), ~0.05 eV (UltiMono at 300 kV: Imperial, Jülich-class), 0.025 eV.

## The group's own study: Lei & Wang, arXiv 2506.22352 (Warwick; Yu Lei wrote the reconstruction code used here)

Multislice hollow ptychography; the depth part is simulation only: abTEM, 300 kV, **80 mrad, aberration-free probe (defocus
−0.7 nm only), fully coherent (no energy spread or source size anywhere in the paper)**, 8 nm PrScO3, 256² detector to
160 mrad, Poisson 1e5–1e8 e/Å², frozen phonons at 300 / 90 / 4 K, LSQML with the probe updated from the first slice (no modes
mentioned). Single-atom depth FWHM 9.5 / 3.5 / 2.8 Å at 20 / 40 / 80 mrad. Result: at 300 K only layers 8–9 Å apart are
resolved (Pr, one O set); **at 90 K and 4 K every layer, oxygen included, is resolved at 1e8 e/Å²**; at 1e7 the 3–4 Å-spaced Sc
and O layers are not, even at 90 K. Thermal model: per-axis σ = √(tr U / 3) at RT (Pr 0.097, Sc 0.087, O 0.103 Å), scaled by
√(T/300) to 90 K (0.053 / 0.048 / 0.057 Å); 4 K taken as 23 % of RT (0.022–0.024 Å).
How this campaign differs, and questions for Yu Lei:
- Their microscope is ideal (no aberrations, coherent); ours adds the ARM tableau and a cold FEG. Their LN2 conclusion holds
  only if the illumination is coherent enough — which is exactly the energy-spread question here.
- Displacements: our 100 K (Debye fit with zero-point motion) per-axis u = Pb 0.064, Ti 0.051, O 0.072 Å, larger than their
  90 K values; their 4 K (23 % of RT) is below the zero-point motion our Debye fit gives for O (~0.065 Å) — likely optimistic.
- Our simulator's own check found abTEM 1.0.5's `sigmas` is the 3-D rms magnitude (per-axis = sigmas/√3; see the comment above
  `_B_BY_SPECIES`). If their per-axis σ went into abTEM directly, their displacements were √3 smaller than intended. Ask.
- Their metric is line profiles of layers; ours is atoms in the right place against guessing. Spacings differ (their Pr 8 Å,
  our Pb 3.9 Å).

## What this means (computed; recheck)

Focal spread = Cc (σE/E)(1+E/E0)/(1+E/2E0), σE = FWHM/2.355; at 300 kV with Cc 1.4 mm it is 7.29 Å rms × (ΔE / 0.3 eV):
0.3 eV 7.3 Å (fails) · 0.25 eV 6.1 Å · 0.2 eV 4.9 Å · **90 meV 2.2 Å** (5.2 Å FWHM, against 3.1 Å depth resolution at
80 mrad) · 55 meV 1.3 Å · 25 meV 0.6 Å. So a Titan Themis monochromator hardly helps; a Spectra Mono at 90 meV cuts the
focal spread 3.3×; whether that is past the wall is the first question. **Everything stays at 300 kV**, so the existing
pipeline runs as it is: no beam-energy change is needed (it is hard-coded at 300 kV in the simulator, planner, driver,
make_probe and the atom finder's config — only matters if a 200 kV target comes back).

## The experiment (the user's decisions, 2026-10-12)

- **Leave the corrector be**: the ARM200F tableau, rows `armf_a080` (and `armf_a090_w140` later), unchanged. No DCOR variant.
- **Phase 1 — coherence**: scan the **focal spread** (Å rms; it lumps Cc × energy spread, and also high-tension ripple and lens-
  current noise, which add in quadrature and set a floor once the beam is monochromated) and the **source size** (Å FWHM),
  everything else as the coherent known-probe reference: 80 mrad, 210 Å region / 105 Å window, noiseless, true probe fixed
  (cheap first). Suggested grid: focal spread 7.3, 5, 3.5, 2.5, 1.5, 0.75 Å with no source size, source size 0.2, 0.4, 0.8 Å
  with no focal spread, then a few joint points near the wall. Then map instruments onto the wall (table in FINAL CHECK).
- **Phase 2 — temperature**: add phonons (RT, 100 K, colder) at the coherence values that survived, to find the temperature
  wall. Then the realistic pipeline (outer search + refine with modes, dose 1e8) on what survives both.
- **Kernels carry the same physics — verified 2026-10-12**: the driver passes PHONONS / PER_SPECIES / PHONON_T and every
  coherence setting to the Pb and Ti grid simulations exactly as to the lab one (sim_meta of the split, 100 K, 70 mrad and
  true-start runs agree leg for leg), each grid using its own element's displacement. Keep it that way: never reuse kernels
  from different physics except as a labelled fallback (`PSF_DIR`).
- **Cc answered**: the 1.4 mm used for the ARM is the whole-instrument value. CEOS (Müller et al. 2006) report a TOTAL Cc of
  1.42 mm, objective plus hexapole corrector, for a CEOS-corrected JEOL with the small-gap URP pole piece (A5 1.2 mm), the
  ARM200F's class. Instruments are mapped by their total Cc (Titan standard gap: 2.4 mm total, 2.0 mm objective alone).

**Built 2026-10-12 — sim dirs now carry the coherence values** (the bare `_coh` put every setting in one sim dir):
`_fs<rms Å>` (FOCAL_SPREAD), `_dE<eV>_cc<mm>`, `_src<FWHM Å>`, `_cs<n>`, with `.` written `p` (`_fs3p5`: the engine's
MATLAB `fileparts` would read a `.` as an extension). The cold-FEG setting every earlier run used (ENERGY_SPREAD_EV=0.3
CC_MM=1.4 SOURCE_FWHM=0.4) keeps `_coh`. Settings the sim would ignore are refused (ENERGY_SPREAD_EV without CC_MM,
FOCAL_SPREAD with it, COHERENCE_SAMPLES with phonons, zero or non-numbers). Dry run: `bash campaign/test_coh_names.sh`
(26 cases); the 15 launch lines in the page manifest give byte-identical sbatch calls under the old and new drivers.

**SUPERSEDED — current plan and status: `mono300/PHASE1.md` (top).** **Phase 1 as proposed 2026-10-12**: armf_a080, noiseless, true probe fixed, 1 mode, 200 it,
`32;16`, PACK_H5=0, CLEANDATA=1, one submission per point with its analysis chained on its pack (GT_REGION=210).
Focal spread alone 7.29 (= the old 0.3 eV × 1.4 mm, 7.294 Å, so with source 0.4 alone it splits the failed `_coh` run
exactly), 5, 3.5, 2.5, 1.5, 0.75 Å; source alone 0.2, 0.4, 0.8 Å. Already in hand: coherent (0, 0) 90/58 vs 26/18 and
`_coh` (7.29 + 0.4) 23/7 vs 19/2. Cost per point: 3 sims × 4.6 h (16 configurations, the hail-mary sims) + 3 recons ×
5.3 h (measured on the coherence-only split: 1.25 h presolve + 4.0 h full) ≈ 30 GPU-h; 9 points ≈ 270 GPU-h. The CPU
quota binds (sims ask 128G ≈ 18 CPUs, recons 73G ≈ 10, against 160): ≥ ~23 h wall for all nine. Raw data 39 GB a point,
up to ~350 GB at once before the packs clean it.

## The tools, as they stand

- Driver `campaign/run_thin_atomfind.sh` (header documents every option): `LABELS`/`TSV` (rows in `campaign/ceos_sweep.tsv`),
  `PHONONS=16 PER_SPECIES=1 [PHONON_T=<K>]` (Debye model fitted per species to the RT B, `--phonon-temperature`),
  `ENERGY_SPREAD_EV CC_MM SOURCE_FWHM` (or `FOCAL_SPREAD`), `DOSES`, `FIXED_TWIN`, `KICK_ALL KICK KICK_MODES KICK_NITER
  KICK_GROUPING KICK_PSTART`, `SEARCH_C1F SEARCH_C3F` (outer search → `_rk` refined legs), `RTAG`, `RECON_ONLY`, `CLEANDATA`,
  `PACK_H5`, `GROUPING`, `NITER_PRE`. Analysis: `campaign/run_analysis.sh` (`LABELS SUFFIX DEP GT_REGION [PSF_DIR PSF_TAG]`).
- Partial coherence is our own Monte Carlo, not abTEM's: `simulate_4dstem._scan_partially_coherent` gives each of 16
  configurations its own phonon displacements, a defocus offset (stratified Gaussian, exactly the rms focal spread) and a
  probe-position offset (the source size); patterns are the incoherent average. Verified against direct abTEM calls (≤1e-6).
  Focal spread only: set ENERGY_SPREAD_EV + CC_MM (or FOCAL_SPREAD) without SOURCE_FWHM; source only: SOURCE_FWHM alone.
- Results page: `analysis/make_arm_page.py` + `page/arm_results.json` + `page/ARM_PAGE.md`. A new page is probably a new
  template on the same builder — ask the user.

## Traps (each cost a round trip)

| symptom | cause / rule |
|---|---|
| array jobs on `compute` nodes, `[CPU-MLs]`, OOM | an sbatch'd wrapper without the `#SBATCH --gres` header; slurm reads only the submitted script's header. A stub-sbatch dry run cannot catch it: grep every new wrapper for `--gres` |
| large-aperture failure | the reconstruction window, not the aperture: keep < ~5 % of the probe outside it |
| region sim OOM in seconds | grid with a large prime factor (now rounded to 7-smooth) and batch size (2 above 7000 px) |
| `error contains NaNs` in the full engine | probe released at full resolution with tiny GROUPING (90 mrad, "2;2"); release in the presolve only or use larger groups |
| full-physics runs look "okay" on summary.csv | the 2 Å any-element match; use atom_scores |
| tarballs not extracted | the user's tar step often does not run; extract before triage |
| `DependencyNeverSatisfied` | the simulation failed; read its log |
| scratch scripts gone | macOS wipes `/tmp` after ~3 days; keep tests in the repo |

## Working rules

The user runs every cluster job; give exact blocks with their `scp -O` lines (memory `blythe-scp-with-every-script`);
dry-run every driver change with a stub sbatch; look at the phase images before judging (`look-at-the-object`); honest
scores (`atomfind-honest-scoring`); plain language for a professor-level reader (`viva-explanation-level`); new
experiments may need their own repo (`experiments-own-repo`) — ask whether this work stays in ptychoshelves-clean (it
reuses sim/, campaign/, analysis/ unchanged at 300 kV).

## Sources

- Cambridge Spectra 300: https://www.wems.msm.cam.ac.uk/electron-microscopy-equipment/tf-spectra-300
- Spectra 300 datasheet (Thermo Fisher DS0297): https://www.utwente.nl/en/mesaplus/infrastructure/nanolab/analysis/facilities/tem/doc-1-spectra-300-stem-datasheet-ds0297.pdf
- Spectra 300 product page: https://www.thermofisher.com/uk/en/home/electron-microscopy/products/transmission-electron-microscopes/spectra-300-tem.html
- Spectra 300 launch (Aug 2019): https://www.technologynetworks.com/tn/product-news/thermo-fisher-scientific-introduces-revolutionary-stem-platform-323058
- Jülich ER-C Spectra 300: https://er-c.org/index.php/facilities-2/tfs-spectra-300/ ; Imperial Spectra 300: https://www.imperial.ac.uk/materials/facilities/em/titan/
- Lei & Wang, Multislice Hollow Ptychography: https://arxiv.org/abs/2506.22352
- Müller et al. 2006, Advancing the hexapole Cs-corrector for the STEM: https://doi.org/10.1017/S1431927606060600 (CEOS PDF: https://www.ceos-gmbh.de/en/produkte/stemkorrektoren/dcor-ascor/Mueller2006a.pdf)
- CEOS DCOR/ASCOR: https://www.ceos-gmbh.de/en/produkte/stemkorrektoren/dcor-ascor
- Titan Themis launch (2013): https://www.globenewswire.com/news-release/2013/08/01/564108/11292/en/FEI-Revolutionizes-the-World-of-TEM-With-the-Introduction-of-Three-New-Systems.html
- Titan 80-300 (Jülich, Chalmers): https://er-c.org/index.php/facilities-2/fei-titan-80-300-stem/ ; https://www.chalmers.se/en/infrastructure/cmal/instruments/electron-microscopy/tem-titan-80-300/
- Leeds LEMAS (Titan3 Themis 300, ~0.25 eV): https://www.rms.org.uk/community/facilities-database/facilities-database-details.html?slug=lemas
- Titan Themis Wien monochromator (~0.2–0.3 eV at 300 kV): https://afmm.iisc.ac.in/?p=1150
- Imperial (monochromated Titan 80/300): https://www.imperial.ac.uk/engineering/departments/materials/facilities/em/
- ePSIC (ARM300F): https://www.diamond.ac.uk/Instruments/Imaging-and-Microscopy/ePSIC/Capabilities.html ; Oxford: https://eia.web.ox.ac.uk/facilities
- SuperSTEM (the Nions, 40–100 kV): https://www.superstem.org/facility/instrumentation
