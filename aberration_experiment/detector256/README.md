# detector256 — is a 256 × 256-pixel camera enough for the 80 mrad reconstructions?

## Status (2026-10-23)

- **R1 + R2 SUBMITTED 2026-10-23 16:11** (`submit.sh`; log `logs/detector256_r1r2_20261023_161114.txt` on Blythe; every
  guard passed there; the kept lab sim confirmed at (1600, 1420, 1420), the reference geometry).
  R1 `_detc96`: rebin jobs 1339556 / 58 / 60 (lab / Pb / Ti), recons 1339557 / 59 / 61, pack 1339562, analysis 1339563.
  R2 `_detc96b4_up2`: rebin jobs 1339564 / 66 / 68, recons 1339565 / 67 / 69, pack 1339570, analysis 1339571.
  Next: a crash check ~30 min into the recons (the MATLAB changes never ran before: `grep -E "upsampl|kept as it
  is|Probe amplitude|Ndp|NaN" recon_af_armf_a080_lab_fs2p5_src0p1_detc96*_NL18/slurm_*.out`); pull to
  `~/Desktop/detector256_r1r2` once 1339563 and 1339571 are COMPLETED. R3 (optional): `... submit.sh r3` once R2's rebin
  jobs are done. R4: `... submit.sh r4`, only if R2 works.
- **The frame question is settled by test, not argument.** `test_upsampling_frame.py` ports the engine's data and
  model path line by line to numpy. The 256-px data equal the engine's own 4 × 4 binning of the full-pixel truth (up to
  the one scale factor the engine fits at iteration 0), in the same array frame (max relative difference 3e-8). On the
  test's random patterns, blocks shifted by 1, 2 or 3 fine pixels, or binned centred, miss by 33–47 %, so the test can
  fail. It also passes on the jobs' own outputs, at the real geometry and with the real ARM probe.
- **Five engine/driver changes were needed** before upsampling could run on our data. Without them it would have
  shrunk the probe 4×, started from an opaque object ring, and crashed at iteration 200 (list below). None of them
  changes a run without `UPSAMPLING`: all 15 recorded launch lines are unchanged and the naming test passes 26/26.
- **Not yet known:** whether the kept sims are still on Blythe. `submit.sh` checks for the lab, Pb and Ti sims and
  refuses if any is missing. It also checks that the checkout carries these changes.
- **Next, when the results are pulled:** `~/hyperspy-bundle/bin/python aberration_experiment/detector256/det_look.py`
  (scores and depth sections beside the reference); look at the phase images (`<analysis>/phase/`); fill in the results
  table below.

## The question

A referee will ask why the 80 mrad runs reconstruct 1420 × 1420-pixel patterns recorded to ±133 mrad, when lab 4D-STEM
cameras are 128–256 pixels across (Lei & Wang's EMPAD: 128). The solver models a real-space window of N λ / (2 θ_max)
(λ = 0.0197 Å at 300 kV), and that window must hold the probe. The ARM-tableau probe at 80 mrad is 72 Å across (d90;
d50 32 Å, d99 120 Å), which is why the campaign uses a 105 Å window. A 256-pixel camera at ±96 mrad gives a 26 Å
window, which leaves 56 % of the probe outside. So the pixel count is set by probe size × detector angle, not by the
aperture itself.

**The idea** (Batey et al., Phys. Rev. A 89, 043812, 2014): record fewer, bigger pixels and let the solver model them.
Each pattern is computed on a 4×-finer angular grid (the 105 Å window) and summed into the camera's pixels before it
is compared with the data. The angular detail inside a pixel is recovered from the overlap of a dense scan.
PtychoShelves has this as `p.detector.upsampling = k`; it was off in every run so far.

## The test (one change at a time, true probe fixed, the same simulation)

Data: the kept joint sim at the target setting (focal spread 2.5 Å + source 0.1 Å, ARM tableau, 80 mrad, noiseless),
`sim_out_af_armf_a080_{lab,Pb,Ti}_fs2p5_src0p1`. Every leg uses the reference control's settings: NITER 200 + 200,
GROUPING 32;16, BETA_LSQ 0.05, REGLAYER 0, NL 18 (dz 1.53 Å). The Pb and Ti kernel legs get the same camera as the lab
leg, so the kernels come from the identical pipeline (PSF_KERNELS.md).

| run | camera (data) | solver's model | window | object pixel | measurements / unknowns* | recon dirs |
|---|---|---|---|---|---|---|
| reference (done) | 1420 px, ±133 mrad | 1420 px | 105 Å | 0.0739 Å | 31 : 1 | `_fs2p5_src0p1_fixed` |
| **R1** | 1024 px, ±96 mrad | 1024 px | 105 Å | 0.1025 Å | 31 : 1 | `_fs2p5_src0p1_detc96` |
| **R2 — the test** | **256 px**, ±96 mrad (4 × 4 sums) | 1024 px, summed 4 × 4 | 105 Å | 0.1025 Å | 2.0 : 1 | `_fs2p5_src0p1_detc96b4_up2` |
| R3 (optional) | 256 px, ±96 mrad | 256 px | 26 Å (56 % of the probe outside) | 0.1025 Å | — (model wrong) | `_fs2p5_src0p1_detc96b4` |
| R4 (if R2 works) | 354 px, ±132.75 mrad (4 × 4 sums) | 1416 px, summed 4 × 4 | 105 Å | 0.0742 Å | 2.0 : 1 | `_fs2p5_src0p1_detc133b4_up2` |

\* 1600 patterns × pixels, against a 125 Å object (20 Å scan + 105 Å window) × 18 slices × 2 real numbers per voxel.
These are rough; I checked the brief's figures with this arithmetic. With a 0.4 Å step (2500 positions) R2 would be
about 3 : 1.

- **R1** isolates the cost of recording less angle; it ties to the zblur G0c question of whether an atom's depth core
  is set by the recorded angle. Its data are exactly the full data cropped the way the engine's own presolve crops (T2).
  The reference's presolve was 1026 px (±96.2 mrad): the driver's 1.2 α guard rounds 512.00003 up.
- **R2** is the test: the same angle as R1, with 16× fewer pixels. Its data are 0.4 GB a scan instead of 6.7 GB
  (13 GB for the reference). The solver's forward model is the same 1024-pixel model as R1's, with the same probe,
  pixel and window.
- **R3** shows what upsampling buys: the same 256-pixel data, with the solver modelling only the camera's own window.
- **R4** asks whether the angle cut was needed at all: the full ±133 mrad on a 354-pixel camera.
- **Presolve.** At ±96 mrad (1.2 α) the half-width presolve has nothing to crop. In R1 and R2 it therefore runs on
  the full grid, so the schedule stays the reference's two engines. In R4 it crops the fine data to ±96 mrad in whole
  4 × 4 blocks; that crop is block-aligned (T5), and its data are exactly R2's.
- **Cost:** about 2.5–3 h a leg (the reference legs took 1.2 h presolve + 3.8 h full engine at 1420 px; R1 and R2 run
  both engines at 1024 px), so about 18 GPU-h for R1 + R2. The rebin jobs are CPU jobs of a few minutes.

## How to judge (in this order)

1. **Phase images** (`<analysis>/phase/recon_af_*.png`; depth sections in `figs/detector256/det_sections.png` from
   `det_look.py`).
2. **Misfit** eps = final Fourier error × N / 2e5, with N the width of the data the run fitted (1420, 1024, 256 or
   354). `det_look.py` reads N from each run's own `sim_meta.mat`. The error-trace file name carries the model's width
   for an upsampled run, so `make_arm_page.trace` would overstate R2's eps 4×. eps compares like with like only on the
   same data grid (R2 with R3).
3. **Vacuum / slab** rms phase ratio: 0.31 for a clean run; the reference is 0.41.
4. **Right place**, from `make_arm_page.atom_scores`: right element within 0.5 Å in depth, against guessing. Good means
   Pb ≥ 80 % and Ti ≥ 45 %; okay means Pb ≥ 50 %.

**Success:** R2 within a few points of R1 and the reference, and as clean in the images. **If R2 fails, try a denser
scan before anything else.** The sims are at 0.5 Å; a 0.4 Å step (2500 positions, about 3 : 1) needs new lab, Pb and
Ti sims at this point, then a rebin and R2 again.

Reference, reproduced by `det_look.py` from the pulled step-2 analysis: Pb / Ti / O right 91 / 60 / 56 % against
guessing 26 / 18 / 19 %, 3 % wrong element, vacuum / slab 0.41, eps 0.0190 (N 1420): **good**.

## Results (pending)

| run | Pb / Ti right (guess) | vacuum / slab | eps (N) | images | verdict |
|---|---|---|---|---|---|
| reference | 91 / 60 (26 / 18) | 0.41 | 0.0190 (1420) | viewed 2026-10-20 | good |
| R1 | | | | | |
| R2 | | | | | |
| R3 | | | | | |
| R4 | | | | | |

## Caveats to state beside any result

- **An idealised camera.** The data are made by summing 4 × 4 point samples and fitted by a model that sums 4 × 4
  samples. Real pixels integrate continuously and blur into their neighbours (MTF, DQE). A real 256-pixel camera's
  data are worse than these by that blur.
- **Where the beam centre falls in a pixel.** Here the zero-angle direction sits exactly on a pixel corner, in the
  engine's convention, by construction. On a real camera it falls anywhere inside a pixel and must be calibrated or
  fitted. A shift of one fine pixel is a 0.19 mrad tilt, a lateral walk of 0.005 Å through the 27.5 Å slab, so it is
  close to a phase-ramp gauge; but the summed pixels themselves change with it, so it cannot simply be ignored.
- **Noiseless.** Summing raises the counts per pixel 16× at the same dose. Shot noise is not tested here.
- **eps** compares like with like only on the same data grid.

## What had to change, and why (2026-10-23)

Read before use, as the brief asked. The upsampling path is: `initial_checks.m` switches the loader; the loader's
`binning.m` unbins the data onto the fine grid (÷ 16); `matlab_aps.m` calls `apply_binning(p, 2^-k)`; `init_solver.m`
bins the data back; `get_reciprocal_model.m` bins |Ψ|² with `utils.binning_2D` (uncentred blocks, block mean; the
probe's iteration-0 rescale absorbs the 16); `modulus_constraint.m` unbins the ratio back onto the same fine blocks.

1. **`ptycho/+core/initial_checks.m`**: upsampling forced the `matlab_ps` loader. That loader reads data through the
   cSAXS compiled reader, not the `h5read('/dp')` our synthetic data need; it is only in `matlab_aps`, whose fast-reader
   block is commented out. Both loaders carry the same `binning.m`, so `matlab_aps` is now kept when set.
2. **`ptycho/+core/ptycho_model_probe.m`** (Latin-1 file, edited byte-wise): a probe file of any other size than the
   detector's was resampled with `interpolateFT` to the detector size (0.41 Å pixels over 105 Å). `apply_binning`
   then zero-padded it to 1024 px at 0.1025 Å, so the probe came out 4× too small (overlap with the true probe 0.000
   in T8). **The engine wants the probe on the fine grid** (detector size × 2^k). The shipped loader assumes the probe
   fits the camera's own window and zero-pads it, which is useless here. A probe already on the fine grid is now kept
   as it is.
3. **`ptycho/+core/apply_binning.m`**: upsampling grows the start object by the extra window and padded it with
   zeros. That is an opaque ring under the probe's tails at the start (the 26 Å core holds 44 % of the probe), a
   different start from R1's. It now pads with the object's mean, as `load_from_p` already does. Binning only crops,
   so it is unaffected.
4. **`ptycho/+engines/+GPU/LSQML.m`**: the exit-wave buffer was allocated at the model's size, but holds |aPsi| on the
   data grid, so the assignment would have failed at iteration 200 (the end of the presolve). It is now the data
   grid's size; without upsampling the size is unchanged.
5. **`ptycho/run_synthetic_recon_ML.m`**: reads `UPSAMPLING=k` into `p.detector.upsampling`. It gives the presolve
   width to the engine in fine pixels: everything else stays in detector pixels, but the engine sees unbinned data,
   and it read the old width as fine pixels, which would have cropped R2 to ±24 mrad.

Plus `campaign/run_thin_atomfind.sh`: `DET_CROP_MRAD`, `DET_BIN`, `UPSAMPLING`. These run a CPU rebin job per leg
(`sim/rebin_detector.py`) into `sim_out_af_<leg>_<mode><SFX>_det[c<crop>][b<bin>]`, link the copy's true probe on the
modelled grid (`probe_initial_true_up<k>.mat`), and tag dirs and tarballs. A finished copy (`rebin_detector.json`,
written last) is reused under `RECON_ONLY`, and the kept sim is never touched. Known-probe legs only: the switches are
refused with doses, kicks or searches. `campaign/pack_results.sh` also ships `rebin_detector.json`.

**Geometry facts measured on the way** (the kept sim's own `sim_meta.mat`): `Ndpx` = **1420** (the plans' "1419" was
an estimate), 0.18749989 mrad per pixel (±133.1 mrad), object pixel 0.07394 Å, window 105.0 Å, 1600 positions. Crop
rule: the fewest whole 4 × 4 blocks that keep ±θ (to 0.001 px), so ±96 mrad gives 1024 px at both 1 × 1 and 4 × 4, and
±133 mrad gives 1416 px at 4 × 4. The ±96 mrad Fourier crop of the true probe keeps all its power (to 1e-8): the
80 mrad aperture fits inside.

## Files

| file | what |
|---|---|
| `submit.sh` | the hand-over: guards (leaked exports, naming 26/26, old launches unchanged, detector dry run 21/21, frame unit test, sims present at the reference geometry), then R1 + R2 (or `r3` / `r4`) and their analyses; log in `logs/detector256_<what>_<ts>.txt` |
| `test_upsampling_frame.py` | the frame unit test (40 checks: 22 that the engine's source still has the lines it ports or the changes above, 18 on the data, the probes and controls that must fail); `~/hyperspy-bundle/bin/python ...` |
| `test_det_driver.sh` | dry run of the driver's detector switches with a stub sbatch (21 cases) |
| `det_look.py` | scores, eps on the data grid, vacuum ratio, depth sections beside the reference |
| `../../sim/rebin_detector.py` | the post-processing step (crop, B × B sums, probes per grid, sim_meta, record) |
| `HANDOFF_PROMPT.md` | the brief this work started from |
