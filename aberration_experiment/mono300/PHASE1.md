# Phase 1: focal spread and source size, with probe modes (designed 2026-10-12; NOT submitted as designed)

**Status 2026-10-13 (the user's choice): two steps.** Step 1 = the controls only, `phase1_fixed_submit.sh`: the same ten
points, the same `_<point>_fixed` dirs, the sims KEPT. Step 2 = the mode legs below, RECON_ONLY on those sims, at the points
whose control fails, started from the exact mixed state (to build). Why not the design as it stands: the engine seeds each
extra mode at `mode_start_pow` 0.02 (`ptycho/run_synthetic_recon_ML.m`), so an 8-mode leg starts with 0.86 in mode 1 against
a true 0.30 at 7.29 A, and every mode run so far ended with mode 1 at 0.84-0.90 (`_true4` 0.94 -> 0.84 in 500 + 500 it;
`_long_m3` / `_long_m6` 0.90); and its packs made the control wait for the 34 h legs. `phase1_submit.sh` is kept, unrun.

Supersedes the fixed-probe-only proposal in `HANDOFF_MONO300.md` ("Phase 1 as proposed"). The user's rule: probe update
with modes is the method, the fixed probe is the side control. Stub dry run of `phase1_submit.sh`: 120 sbatch calls,
`campaign/test_coh_names.sh` 26/26 first.

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
