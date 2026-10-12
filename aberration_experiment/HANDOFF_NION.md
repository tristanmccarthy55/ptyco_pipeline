# Handoff — the Nion direction: where is the wall in temperature and energy spread?

Written 2026-10-12 at the close of the ARM200F campaign. **Read this whole file first**, then the newest bullets of the
State block in `HANDOFF_REALISTIC.md` (blocks 3–5: every number behind the summary below), `PSF_KERNELS.md` (kernel
rules), and `page/ARM_PAGE.md` (the results page and the honest atom scoring). Where they disagree, this file is newer.

## The goal (the user, 2026-10-12)

- **Target instrument: a monochromated Nion** (UltraSTEM 100MC "HERMES" class; the UK one is SuperSTEM 3 at Daresbury).
  The user: a relatively common microscope with some stability issues, not the state of the art.
- **Find the wall**: the specimen temperature and the energy spread (and anything else that turns out to matter) at
  which depth-resolved multislice ptychography stops working, with realistic aberrations and the probe found from the data
  (outer search, then probe update with modes; probe update always on — memory `probe-update-always-on`).
- **Aberrations**: a measured Nion tableau if one can be found; otherwise the ARM200F tableau as a pessimistic baseline.
- **Hoped-for headline**: it works at liquid-nitrogen temperature with a monochromator.
- **Expectations set by the user**: seeing oxygen is close to impossible; judge Pb and Ti. Phonon blur along the beam is
  acceptable (from the images); the partial-coherence artefact (vacuum loading, one blob per column) is not. A slight
  specimen tilt may help with the phonon blur (not built).

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
| partial coherence alone removes depth (0.3 eV through Cc 1.4 mm = 7.3 Å rms focal spread; 0.4 Å source) | 23 / 7 vs 19 / 2 %, wrong element 34 %; vacuum layers loaded, one blob per column |
| phonons at 100 K (Debye model) recover part | 49 / 22 vs 29 / 19 %, wrong element 8 % |
| 70 mrad full physics is no better | 12 / 0–2 vs 7–9 / 2–3 % |
| refined probe (searched start, 4 modes, 500 it) under full physics | 19–20 / 12–22 % vs 11–13 / 8–10 %: off at both doses |

Results page (generated, republished to one link): https://claude.ai/artifact/MUkJ7yYrsek3nTLXcZgCrk (v9).

## Facts checked 2026-10-12 (sources at the end)

- **UK Nions are at SuperSTEM** (Daresbury): SuperSTEM 3 = Nion UltraSTEM 100MC "HERMES", **40–100 kV**, monochromator
  "for ultimate 5 meV energy resolution", C5 Nion QO corrector, Dectris ELA (EELS); SuperSTEM 2 = UltraSTEM 100,
  40–100 kV, C5 QO corrector, ±35° tilt. Their instrumentation page lists **no 4D-STEM/pixelated detector** — find out what
  they use for ptychography before planning dose and scan time.
- **Monochromator**: 5.7 meV at 60 kV and 4.2 meV at 30 kV (HERMES; UCI's page), long-term energy stability of the order
  of 10 meV/min; SuperSTEM "routinely sub 20 meV". The best numbers are EELS resolutions at very low current.
- **A 200 kV HERMES exists** (UCI "Nion200"); none found in the UK.
- **Cc** of the UltraSTEM 100 objective: **1.1 mm at 100 keV** (4 mm gap) — from a search summary whose source page could not
  be opened: CONFIRM before relying on it.
- **No published Nion residual-aberration tableau found yet.** The first place to look: Krivanek et al., Ultramicroscopy
  108 (2008) 179 ("An electron microscope for the aberration-corrected era"); later Nion corrector papers (e.g. "Tuning
  fifth-order aberrations in a quadrupole-octupole corrector", ORNL). Ask the user whether SuperSTEM can share a measured
  tableau.

## What this changes (computed 2026-10-12; recheck)

- **Voltage sets the aperture needed.** Depth resolution λ/α²: 3.1 Å (our 300 kV / 80 mrad) needs **110 mrad at 100 kV**,
  126 at 60 kV, 90 at 200 kV. At 100 kV and 80 mrad it is 5.8 Å. Larger apertures mean larger probes and windows (the
  campaign's window/compute limit), unless the Nion corrector keeps the phase flatter than the ARM's.
- **Focal spread** = Cc (σE/E)(1+E/E0)/(1+E/2E0), σE = FWHM/2.355 (simulate_4dstem's formula). With Cc 1.1 mm at 100 kV:
  350 meV (no monochromator) 17.8 Å rms — as bad as our failing 300 kV case; 50 meV 2.5 Å; 20 meV 1.0 Å; 10 meV 0.51 Å;
  5 meV 0.25 Å. The monochromator fixes the temporal part.
- **It does not fix the spatial part**: the source-size term (0.4 Å FWHM in the runs so far) was never tested alone. The
  coherence-only run had both on, so the first new test is focal spread only vs source size only.
- **Dose at monochromated currents**: 1e8 e/Å² over the 20 × 20 Å field takes 21 min at 5 pA, 1.8 h at 1 pA; 1e7 takes
  ten times less. Stability (drift, energy drift) then matters.

## What must be built first

1. **Beam energy as a parameter.** Today 300 kV is hard-coded: `sim/simulate_4dstem.py` `ENERGY_EV = 300e3` (no CLI
   option); `campaign/plan_probe.py` `LAM` and `abtem.Probe(energy=300e3, ...)` (three places); `campaign/run_thin_atomfind.sh`
   `LAM=0.0196877` (slice count); `sim/make_probe.py` refuses a sim whose energy differs from the module's. The
   reconstruction already reads `meta.energy_kev`. Also 300 kV in: `analysis/atomfind/config.py` (`wavelength_A = 0.01969`,
   the atom finder's own), `analysis/c1_objective.py`, `campaign/run_c1_search.sh`, and the figure scripts
   (`make_meeting_figs.py`, `make_ceos_figs.py`, `make_ronchigram_fig.py`, `campaign/test_haider.py`); `make_arm_page.py`
   reads the simulator's. The region grid sampling depends on λ and the detector angle: re-check memory per batch.
2. **A tableau**: a measured Nion one, or the ARM200F one re-planned at the new voltage (`plan_probe.py --arm --fight B2 A1`
   needs the energy option) — say which, and that the ARM is the pessimistic case.
3. **The wall scan** (cheap first): one voltage and aperture, true probe fixed, noiseless; energy spread × temperature
   (e.g. 350 / 100 / 50 / 20 / 10 meV × RT / 100 K / colder), plus focal-spread-only and source-size-only legs. Then the
   realistic pipeline (search + refine with modes, dose) only on the survivors.
4. **Tilt** (a few degrees; not built in the simulator).

## The tools, as they stand

- Driver `campaign/run_thin_atomfind.sh` (header documents every option): `LABELS`/`TSV` (rows in `campaign/ceos_sweep.tsv`),
  `PHONONS=16 PER_SPECIES=1 [PHONON_T=<K>]` (Debye model fitted per species to the RT B, `--phonon-temperature`),
  `ENERGY_SPREAD_EV CC_MM SOURCE_FWHM` (or `FOCAL_SPREAD`), `DOSES`, `FIXED_TWIN`, `KICK_ALL KICK KICK_MODES KICK_NITER
  KICK_GROUPING KICK_PSTART`, `SEARCH_C1F SEARCH_C3F` (outer search → `_rk` refined legs), `RTAG`, `RECON_ONLY`, `CLEANDATA`,
  `PACK_H5`, `GROUPING`, `NITER_PRE`. Analysis: `campaign/run_analysis.sh` (`LABELS SUFFIX DEP GT_REGION [PSF_DIR PSF_TAG]`).
- Partial coherence is our own Monte Carlo, not abTEM's: `simulate_4dstem._scan_partially_coherent` gives each of 16
  configurations its own phonon displacements, a defocus offset (stratified Gaussian, exactly the rms focal spread) and a
  probe-position offset (the source size); patterns are the incoherent average. Verified against direct abTEM calls (≤1e-6).
- Results page: `analysis/make_arm_page.py` + `page/arm_results.json` + `page/ARM_PAGE.md`. A Nion page is probably a new
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
scores (`atomfind-honest-scoring`); plain language for a professor-level reader (`viva-explanation-level`); 300 kV was a
deliberate choice for the ARM work (`target-arm200f`); new experiments may need their own repo (`experiments-own-repo`) —
ask the user whether the Nion work stays in ptychoshelves-clean (it reuses sim/, campaign/, analysis/).

## Sources

- SuperSTEM instrumentation: https://www.superstem.org/facility/instrumentation
- UCI Nion200 / HERMES monochromator: https://imri.uci.edu/facilities/tem/tem-instruments/nion200
- Nion products: https://nion.com/products.html ; Nion HERMES (Bruker): https://www.bruker.com/en/products-and-solutions/microscopes/electron-microscopes/nion-hermes.html
- ORNL HERMES: https://ornl.gov/content/nion-hermes-mac-stem
- Krivanek et al. 2008: https://doi.org/10.1016/j.ultramic.2007.07.010
