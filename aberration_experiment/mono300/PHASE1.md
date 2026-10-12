# Phase 1: focal spread and source size, with probe modes (submitted 2026-10-12)

Supersedes the fixed-probe-only proposal in `HANDOFF_MONO300.md` ("Phase 1 as proposed"). The user's rule: probe update
with modes is the method, the fixed probe is the side control. Submitted with `phase1_submit.sh` (stub dry run: 120 sbatch
calls, `campaign/test_coh_names.sh` 26/26 first).

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

M = the modes holding 99 % of the real-gun state, plus one, capped at 8 (8 modes at `4;2` is the most the GPU has run:
same waves in flight as `_true4`). "sim" = the sim's own 16 configurations, "real gun" = a continuous Gaussian (1000
samples). Source: `mode_spectrum.py` -> `mode_spectrum.txt` (exact eigen-spectrum of the mixture; the common aberration
phase cannot change it, so only aperture, wavelength and offsets matter). The source scan dropped 0.8 A (8 modes hold 62 % /
20 %) and gained 0.05 and 0.1 A: at 80 mrad the probe is ~lambda/alpha = 0.25 A, so even a 0.2 A source leaves a quarter
of the state in one mode.

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
