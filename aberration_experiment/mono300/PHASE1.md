# Phase 1: focal spread and source size, with probe modes (designed 2026-10-12; NOT submitted as designed)

**Status 2026-10-13 (the user's choice): two steps.** Step 1 = the controls only, `phase1_fixed_submit.sh`: the same ten
points, the same `_<point>_fixed` dirs, the sims KEPT -- DONE, results below ("Step 1 results"). Step 2 = the mode legs below, RECON_ONLY on those sims, at the points
whose control fails, started from the exact mixed state (to build). Why not the design as it stands: the engine seeds each
extra mode at `mode_start_pow` 0.02 (`ptycho/run_synthetic_recon_ML.m`), so an 8-mode leg starts with 0.86 in mode 1 against
a true 0.30 at 7.29 A, and every mode run so far ended with mode 1 at 0.84-0.90 (`_true4` 0.94 -> 0.84 in 500 + 500 it;
`_long_m3` / `_long_m6` 0.90); and its packs made the control wait for the 34 h legs. `phase1_submit.sh` is kept, unrun.

Supersedes the fixed-probe-only proposal in `HANDOFF_MONO300.md` ("Phase 1 as proposed"). The user's rule: probe update
with modes is the method, the fixed probe is the side control. Stub dry run of `phase1_submit.sh`: 120 sbatch calls,
`campaign/test_coh_names.sh` 26/26 first.

## Step 1 results (2026-10-14; controls: true probe fixed, 1 mode, noiseless, own matched kernels)

Submitted 2026-10-13 ~08:08 (tarballs `*_fixed_20261013_0808*`, job ids not recorded), pulled to `~/Desktop/mono_phase1`.
Scores and sheets: `phase1_look.py` -> `phase1_scores.json`, `figs/mono300/phase1_focal.png`, `phase1_source.png` (phase
images viewed). Right place (right element within 0.5 A in depth) vs guessing; mode 1 = the sim's first-mode share.
Both references reproduce the page's numbers (coherent 90/58, `_coh` 23/7).

| point | mode 1 | Pb right / guess | Ti right / guess | wrong el. | vacuum / slab rms | verdict |
|---|---|---|---|---|---|---|
| coherent | 1.00 | 90 / 26 | 58 / 18 | 1 % | 0.31 | good |
| fs 0.75 / 1.5 / 2.5 A | 0.95 / 0.84 / 0.68 | 92 / 90 / 89 | 57 / 60 / 58 | 1-2 % | 0.31 / 0.32 / 0.39 | good |
| fs 3.5 A | 0.55 | 84 / 27 | 58 / 18 | 2 % | 0.50 | good |
| fs 5 A | 0.42 | 39 / 26 | 28 / 16 | 11 % | 0.74 | off |
| fs 7.29 A | 0.30 | 31 / 24 | 15 / 15 | 12 % | 1.03 | off |
| src 0.05 / 0.1 A | 0.87 / 0.61 | 90 / 96 | 56 / 54 | 1-2 % | 0.30 / 0.30 | good |
| src 0.2 A | 0.26 | 82 / 28 | 39 / 21 | 4 % | 0.35 | okay |
| src 0.4 A | 0.17 | 21 / 15 | 1 / 1 | 47 % | 0.48 | off |
| `_coh` (fs 7.29 + src 0.4) | 0.13 | 23 / 19 | 7 / 2 | 34 % | 1.19 | off |

- **Focal-spread wall: between 3.5 A (good) and 5 A (off), sharp.** 0.75-2.5 A look like the coherent run; 3.5 A starts to
  load the vacuum layers; from 5 A the vacuum fills and the columns become leaning streaks (the coherence artefact).
- **Source-size wall: between 0.1 A (good) and 0.4 A (off); 0.2 A keeps Pb (82 %) but loses Ti (39 %).** A different
  failure: the vacuum stays clean (0.48 vs 0.31) while the slab's contrast fades into smooth streaks, like phonons (0.39).
  The sim's 16-configuration source is the optimistic case (a real 0.2 A source needs 17 modes for 99 %, the sim's 9).
- **The old `_coh` failure was both at once**: each half fails alone. Its vacuum loading is the focal spread's.
- The two series do not collapse on the mode-1 share: src 0.2 (0.26) is okay where fs 5 (0.42) is off, so focal spread
  costs depth more per unit of lost coherence (a likely reason: defocus is what multislice trades against depth).
- **Instruments (focal spread = 17.4 A per eV x mm, before high-tension / lens ripple):** good needs dE <= 0.14 eV at
  Cc 1.4 mm (<= 0.08 eV at 2.4 mm); off by 0.21 eV (0.12 eV). ARM cold FEG 0.3 eV (7.3 A) and a Titan Wien at 0.2 eV
  (8.3 A at Cc 2.4) are past it; a Wien at its best 0.1 eV (4.2 A) or the Spectra Mono at its 0.2 eV spec with Cc 1.4
  (4.9 A) sit on it; the Spectra Mono at 90 meV (2.2-3.8 A) and UltiMono (1.2-2.1 A) are inside. Real sources at the
  specimen are 0.24-0.6 A (CEOS 2006): at or past the source wall, so the source image must be shrunk (current ~ size^2).

## Design

armf_a080 (ARM tableau, 80 mrad, 210 A region, BIN 2, NL 18), noiseless, phonons off, 16 coherence configurations.
Each point: lab, Pb and Ti sims with the same coherence setting, then per sim two legs:

- `_fixed`: true probe fixed, 1 mode, 200 it, `32;16` (the control).
- mode leg (no tag): TRUE start (`KICK=1.0`, C1 -958.6 A, C3 -40.66 um), probe released (PSTART 40/20, TEM aperture
  constraint), M modes, 250 it, grouping G, 48 h.

| point | dirs | M | G | 1 mode holds (sim / real gun) | M modes hold (sim / real gun) |
|---|---|---|---|---|---|
| focal spread 7.29 A rms | `_fs7p29` | 8 | 4;2 | 30 / 32 % | 100 / 99.6 % |
| 5 | `_fs5` | 7 | 4;2 | 42 / 44 % | 100 / 99.9 % |
| 3.5 | `_fs3p5` | 5 | 4;2 | 55 / 57 % | 100 / 99.9 % |
| 2.5 | `_fs2p5` | 4 | 8;4 | 68 / 70 % | 100 / 99.9 % |
| 1.5 | `_fs1p5` | 4 | 8;4 | 84 / 85 % | 100 / 100 % |
| 0.75 | `_fs0p75` | 3 | 8;4 | 95 / 96 % | 100 / 100 % |
| source 0.05 A FWHM | `_src0p05` | 4 | 8;4 | 87 / 87 % | 100 / 100 % |
| 0.1 | `_src0p1` | 8 | 4;2 | 61 / 63 % | 100 / 99.6 % |
| 0.2 | `_src0p2` | 8 | 4;2 | 26 / 30 % | 97.6 / 90.2 % |
| 0.4 | `_src0p4` | 8 | 4;2 | 16.5 / 9.9 % | 78 / 53 % |

M = the modes holding 99 % of the real-gun state, plus one, capped at 8 (8 modes have never run; at `4;2` they
put the same waves in flight as `_true4`'s 4 modes at `8;4`, so the GPU memory should fit). "sim" = the sim's own 16 configurations, "real gun" = a continuous Gaussian (1000
samples). Source: `mode_spectrum.py` -> `mode_spectrum.txt` (exact eigen-spectrum of the mixture; the common aberration
phase cannot change it, so only aperture, wavelength and offsets matter). The source scan dropped 0.8 A (8 modes hold 62 % /
20 %) and gained 0.05 and 0.1 A: at 80 mrad the probe's finest detail is ~lambda/alpha = 0.25 A (the probe itself is
~32 A across, d50), and a shift of that size stops it overlapping itself, so even a 0.2 A source leaves a quarter of the
state in one mode. Real sources at the specimen are 0.24-0.6 A (CEOS 2006), i.e. the 0.2 and 0.4 A points.

## What to watch

- **Does the reconstruction find the mixed state?** Compare its mode occupancies with `mode_spectrum.txt`. `_true4`
  (full physics, 4 modes) kept 0.84 in mode 1, while the true probe state there (the old `_coh`, 7.29 A + 0.4 A) has only
  0.13 in mode 1 and 4 modes can hold at most 42 %: it did not find it.
- **Source size in the sim is optimistic.** 16 configurations make a rank-16 state; a real 0.4 A source needs 35 modes for
  95 % (the sim's needs 12). Focal spread is faithful (stratified quantiles; sim and real agree to a few %).
- `_fs7p29` (fixed) and `_src0p4` (fixed) split the failed `_coh` run exactly into its halves.

## Cost (estimates)

Run time of a mode leg ~ 1.31 x N x (29 + 43 M) s, fitted to `_true4` (4 modes, 500 + 500 it: 8.9 h presolve + 27.8 h
full engine, logs in ~/Desktop/arm_true4) against the fixed legs (5.3 h at 200 it); the 6-mode 250 it `4;2` long kick
also fitted in 48 h. 8 modes ~34 h, 7 ~30, 5 ~22, 4 ~18, 3 ~14. Total ~1070 GPU-h (sims 138, fixed 159, modes 770);
~4-5 days wall (CPU quota: a recon's 73G ~ 10 of 160 CPUs). Raw data 39 GB a point, kept until the pack, which waits on
the slowest leg: up to ~390 GB on the share at once.

## Outputs

Per point: pack `atomfind_results_armf_a080_<point>_<ts>.tgz`; analyses `analysis_armf_a080_<point>_<ts>.tgz` (mode legs)
and `analysis_armf_a080_<point>_fixed_<ts>.tgz` (control). Pull to `~/Desktop/mono_phase1/{results,analysis}`.
