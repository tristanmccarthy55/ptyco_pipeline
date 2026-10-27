# detector256 — is a 256 × 256-pixel camera enough for the 80 mrad reconstructions?

## Status (2026-10-27): the answer

**A 256 × 256 camera is enough at 80 mrad, if the solver models its pixels and is given about twice the iterations.**
With 400 + 400 iterations the 256-pixel camera finds the Pb depths as well as the full 1420-pixel camera (91 %) and the
Ti depths at "good" (47 %, against the reference's 60). The same data fitted without modelling the pixels find nothing.
Noiseless, true probe, idealised pixels: see the caveats below.

Every run on the same simulation (80 mrad ARM tableau, focal spread 2.5 Å + source 0.1 Å, noiseless, true probe fixed),
kernels cropped to 2.22 Å, all centred on their own atom from 16 atoms. "Camera's own window" is N λ / 2θ, the window the
camera would give on its own. The repo's documents took this to be the window the probe has to fit in (next
section). The scores count 78 Pb, 144 Ti and 432 O atoms, so one standard deviation is about ±3 points for Pb and ±4 for
Ti.

| run | camera | iterations | camera's own window (probe outside it) | Pb right (guess) | Ti right (guess) | vacuum / slab | hours a leg | verdict |
|---|---|---|---|---|---|---|---|---|
| reference | 1420 px, ±133 mrad | 200 + 200 | 105 Å (all inside) | 91 (26) | 60 (18) | 0.41 | 5.1 | good |
| R1 | 1024 px, ±96 mrad | 200 + 200 | 105 Å (all inside) | 87 (27) | 57 (17) | 0.42 | 2.4 | good |
| C512 | 512 px, ±96 mrad | 200 + 200 | 52 Å (15 %) | 84 (27) | 46 (19) | 0.40 | 2.4 | good |
| R4 | 354 px, ±133 mrad | 200 + 200 | 26 Å (56 %) | 80 (29) | 46 (19) | 0.38 | 5.6 | good |
| R2 | **256 px**, ±96 mrad | 200 + 200 | 26 Å (56 %) | 86 (29) | 25 (22) | 0.39 | 2.3 | okay |
| **R2_n400** | **256 px**, ±96 mrad | **400 + 400** | 26 Å (56 %) | **91 (27)** | **47 (18)** | 0.40 | 4.7 | **good** |
| R2_n1000 | 256 px, ±96 mrad | 1000 + 1000 | 26 Å (56 %) | 91 (27) | 49 (18) | 0.43 | 11.6 | good |
| C128 | 128 px, ±96 mrad | 200 + 200 | 13 Å (88 %) | 65 (28) | 22 (18) | 0.40 | 2.3 | okay |
| R3 | 256 px, pixels **not** modelled | 200 + 200 | 26 Å (56 %), and the solver uses it | no atoms | | 0.83 | 0.2 | failed |
| S04 | 256 px on a 0.4 Å scan | 200 + 200 | 26 Å (56 %) | not back | | | | |

O right (guess): reference 56 (19), R1 51 (18), C512 40 (20), R4 46 (21), R2 23 (21), R2_n400 31 (20), R2_n1000 44 (21),
C128 31 (20). Every phase image was viewed (`results/phase/`). Each is a clean crystal with the reference's column lean
and quiet vacuum; C128 is a little softer, with fringes at two corners. All logs are clean: upsampling on, probe kept on
the fine grid, no NaN, full traces. Figures: `results/det_summary.png` (one row per run, with the hours) and
`results/det_limits.png`.

- **R2's Ti loss was mostly convergence, not missing information.** On the same 256-px data, more iterations take Ti from
  25 to 47 % (400) to 49 % (1000), and Pb from 86 to 91 %. Most of the gain comes by 400. Going to 1000 adds 2 points of
  Ti (within the ±4) and 13 of O (31 → 44). The misfit is still falling slowly: the lab leg fell 1 % over the last 200
  iterations.
- **Coarser pixels make each iteration less effective; they do not, so far, lose the depths.** At 200 iterations Ti
  falls with the pixel count: 57 (1024 px), 46 (512), 25 (256), 22 (128). The 256-px camera recovers with more
  iterations; the 128-px one was not given more.
- **More angle on the same pixel converges faster, but costs the same.** R4 has R2's 0.75 mrad pixel to ±133 mrad and
  reaches Ti 46 % in 200 iterations, where R2 needed 400. Its solver grid is 1416 px, though, so its legs take 5.6 h,
  about R2_n400's 4.7 h.
- **No coarse camera has yet reached the reference's Ti.** The best is 49 % (R2_n1000) against R1's 57 and the
  reference's 60, about two standard deviations. This does not yet show missing information: R1 was not run longer,
  and might gain too. S04, the denser scan, is the other test.
- **S04 is not back.** On 2026-10-27 at 15:45 its two tarballs were not on Blythe. Its sims (2500 positions) had 12 h
  and its recons 24 h; whether they are still running or failed is not known (check in its submission log,
  `logs/detector256_limits_<ts>.txt`, and `sacct`). The five other runs were pulled into `~/Desktop/detector256_limits`.

### What this overturns

**"The window is set by the detector" is wrong.** Until today the repo took the window to be the camera's own,
N λ / 2θ_max, and concluded that a large probe needs a camera with fine angular pixels: 0.19 mrad, 1420 px to ±133 mrad
for the 105 Å window. The statement, or rules built on it, appears in:

- `report/aberrations_paper/aberrations_paper.tex` (draft 1), Geometry section: "The window is set by the detector: a
  pattern N pixels across, recorded to an angle θ_max, gives a window Nλ/2θ_max. The solver handles patterns up to about
  1420 pixels ..."
- this README's "The question" below, and `HANDOFF_PROMPT.md`: "the pixel count is set by probe size × detector angle".
- the campaign primer (`artifacts/primer_src/primer_template.html`): "The window (region over detector binning, 17.5 to
  140 Å) must hold the probe."
- `campaign/README.md`: the planner's `bin` column, "probe size → real-space window" (`plan_probe.region_geometry` sets
  the simulated detector's pixel from the window).
- the paper's referee-risk list: "ideal point-sampled detector (~1400 px needed at ±133 mrad)".

**What is true instead: the window is set by the grid on which the solver computes each pattern.** The camera's pixels
can be coarser, provided the solver adds its fine pattern up into the camera's pixels before comparing it with the data.
A camera with 4× coarser pixels (0.75 mrad) gives a 26 Å window on its own, which leaves 56 % of the probe outside, and
it still gives the full camera's Pb depths and good Ti (R2_n400). A camera with 8× coarser pixels (1.5 mrad; its own
window 13 Å, 88 % of the probe outside) still gives the Pb depths at "okay" (65 %).

**Still true: the window must hold the probe**, the solver's window, that is. R3 proves it: the same 256-px data, fitted
on the camera's own 26 Å grid, find no atoms, and the misfit is 42× worse. So the rule is not withdrawn, only moved from
the camera to the solver. The "dense scan, not big window" rule and the large-aperture trap ("the reconstruction window,
not the aperture") stand as written; they concern the solver's window.

**What it means for the experiment.** The 1420-pixel patterns of every simulation describe the forward model, not the
camera the experiment needs. A 256-pixel camera at ±96 mrad records what the 80 mrad runs need. The simulation itself
must still be computed on the fine grid, because the physics happens there; only the recorded data can be coarse.

### Speed: what does and does not make the reconstructions faster

- **The coarse camera does not.** Time per iteration is set by the solver's grid, not the data's. R2 (256-px data on a
  1024 grid) took 2.3 h a leg; R1 (1024-px data) took 2.4 h. R2 also needs about twice the iterations for Ti, so it is
  about twice as slow to the same quality: 4.7 h against 2.4 h. R3 ran 12× faster on its 256 grid, and failed.
- **Cutting the recorded angle does: about 2× faster, at no cost we can measure.** Cutting from ±133 to ±96 mrad
  (1.2 α at 80 mrad) makes the solver's grid 1024 px instead of 1420 at the same 105 Å window. Each leg falls from 5.1 h
  to 2.4 h; the reference's full engine alone took 3.8 h at 1420 px against 1.2 h at 1024. It cost 4 points of Pb and 3
  of Ti (R1 against the reference), within the ±3–4 scatter. Only 0.13 % of the electrons land beyond ±96 mrad. Today
  this is `DET_CROP_MRAD=96` on a known-probe leg. Dose, kick and search legs would need that switch extended, because it
  refuses them now.
- **Part of that gain is the FFT size, which could be had without cutting angle (a hint, untested).**
  - The grid's area grows as θ²: 1420² / 1024² = 1.9×. Yet the full engine took 3.1× longer at 1420 px.
  - R4's presolve at 1032 px took 1.6× longer than the reference's at 1026 px, almost the same size.
  - FFTs are slow on sizes with a large prime factor: 1420 = 4·5·71 and 1032 = 8·3·43, against 1024 = 2¹⁰ and
    1026 = 2·27·19.
  - Choosing grid sizes made of small primes, as the region sims already do (7-smooth), may give about 1.5× at ±133 mrad
    with no angle cut. It is one run to check.
  - Limits: tested at 80 mrad only, noiseless, true probe fixed. At 90 mrad, 1.2 α is 108 mrad, a 1152-px grid for 105 Å
    (about 1.5× faster). Whether a searched probe or probe modes miss the outer angles is untested.
- **What the coarse camera does save is data:** 0.42 GB a scan of 1600 patterns at 256 px, against 6.7 GB at 1024 px and
  12.9 GB at 1420. That means less disk, transfer and host memory. GPU memory was not logged, so whether the space freed
  allows larger windows on the same card is untested. On a microscope a small camera also reads out faster, so a scan is
  shorter and drifts less; that is not tested here.

### Still open

- **S04** (0.4 Å scan): pull when back. It asks whether positions buy Ti faster than iterations.
- **R1 at 400 + 400 iterations**: is 49 against 57 % a real information gap or only convergence? (no new sims: R1's
  copies.)
- **C128 at 400 + 400 or more**: does Ti come back at 0.5 measurements per unknown?
- **Shot noise**: summing raises the counts per pixel 16× at the same dose; untested.

### Record

- **The first analysis was faulty (fixed):** the kernel extractor's 30-px crop is 3.08 Å at the 0.1025 Å pixel, past
  the 2.96 Å grid, so R1's kernels centred on a neighbour (scored 0 %). The fix is `PSF_HALF_A=2.22` (= the reference's
  30 px); re-analyses 1348348 / 1348349 were pulled to `~/Desktop/detector256_r1r2_kh`. `det_look.py` now flags such
  kernels.
- **Runs:**
  - R1 + R2 were submitted 2026-10-23 16:11 and pulled to `~/Desktop/detector256_r1r2` (reconstructions and sims) and
    `..._r1r2_kh` (analyses).
  - R3 was submitted 2026-10-26 13:16 (recons 1348350 / 51 / 52, analysis 1348354) and pulled to
    `~/Desktop/detector256_r3`.
  - The limits wave was submitted 2026-10-26 22:37 (`submit_limits.sh`; tarball stamps `20261026_2237*`) and pulled to
    `~/Desktop/detector256_limits` on 2026-10-27, all but S04. `det_look.py` (no arguments) scores everything into
    `results/`.
- **Submission record:**
  - R1 `_detc96`: rebin jobs 1339556 / 58 / 60 (lab / Pb / Ti), recons 1339557 / 59 / 61, pack 1339562, analysis
    1339563.
  - R2 `_detc96b4_up2`: rebin jobs 1339564 / 66 / 68, recons 1339565 / 67 / 69, pack 1339570, analysis 1339571.
  - Log on Blythe: `logs/detector256_r1r2_20261023_161114.txt`; the kept lab sim confirmed at (1600, 1420, 1420).
  - Limits-wave recon slurm ids, read from the pulled logs:

    | run | lab / Pb / Ti recons |
    |---|---|
    | R2_n400 | 1348891 / 92 / 93 |
    | R2_n1000 | 1348896 / 97 / 98 |
    | C512 | 1348902 / 04 / 06 |
    | C128 | 1348910 / 12 / 14 |
    | R4 | 1348918 / 20 / 22 |
- **The frame question is settled by test, not argument.** `test_upsampling_frame.py` ports the engine's data and
  model path line by line to numpy.
  - The 256-px data equal the engine's own 4 × 4 binning of the full-pixel truth, in the same array frame (max relative
    difference 3e-8), apart from the one scale factor the engine fits at iteration 0.
  - On the test's random patterns, blocks shifted by 1, 2 or 3 fine pixels, or binned centred, miss by 33–47 %, so the
    test can fail.
  - It also passes on the jobs' own outputs, at the real geometry and with the real ARM probe.
- **Five engine/driver changes were needed** before upsampling could run on our data. Without them it would have
  shrunk the probe 4×, started from an opaque object ring, and crashed at iteration 200 (list below). None of them
  changes a run without `UPSAMPLING`: all 15 recorded launch lines are unchanged and the naming test passes 26/26.

## The question

A referee will ask why the 80 mrad runs reconstruct 1420 × 1420-pixel patterns recorded to ±133 mrad, when lab 4D-STEM
cameras are 128–256 pixels across (Lei & Wang's EMPAD: 128). The solver models a real-space window of N λ / (2 θ_max)
(λ = 0.0197 Å at 300 kV), and that window must hold the probe. The ARM-tableau probe at 80 mrad is 72 Å across (d90;
d50 32 Å, d99 120 Å), which is why the campaign uses a 105 Å window. A 256-pixel camera at ±96 mrad gives a 26 Å
window, which leaves 56 % of the probe outside. So the pixel count is set by probe size × detector angle, not by the
aperture itself. **(That was the premise, and it was wrong for the camera: it sets the solver's grid only. See "What
this overturns" above.)**

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

1. **Phase images** (`<analysis>/phase/recon_af_*.png`, 75 % copies in `results/phase/`; depth sections in
   `results/det_sections_*.png` from `det_look.py`).
2. **Misfit** eps = final Fourier error × N / 2e5, with N the width of the data the run fitted (1420, 1024, 256 or
   354). `det_look.py` reads N from each run's own `sim_meta.mat` (the error-trace file name carries the data width
   too: `S00001_256x256_...` for R2, checked 2026-10-26). eps compares like with like only on the
   same data grid (R2 with R3).
3. **Vacuum / slab** rms phase ratio: 0.31 for a clean run; the reference is 0.41.
4. **Right place**, from `make_arm_page.atom_scores`: right element within 0.5 Å in depth, against guessing. Good means
   Pb ≥ 80 % and Ti ≥ 45 %; okay means Pb ≥ 50 %.

**Success:** R2 within a few points of R1 and the reference, and as clean in the images. **If R2 fails, try a denser
scan before anything else.** The sims are at 0.5 Å; a 0.4 Å step (2500 positions, about 3 : 1) needs new lab, Pb and
Ti sims at this point, then a rebin and R2 again.

Reference, reproduced by `det_look.py` from the pulled step-2 analysis: Pb / Ti / O right 91 / 60 / 56 % against
guessing 26 / 18 / 19 %, 3 % wrong element, vacuum / slab 0.41, eps 0.0190 (N 1420): **good**.

## Limits wave (2026-10-26): what limits Ti on a 256-px camera, and where Pb goes too

R2 kept Pb (86 %) and lost Ti (25 against 22 guessing). Each run changes ONE thing against R2 (256 px, ±96 mrad, 0.5 Å
step, 200 + 200 iterations); same sim physics, true probe fixed, lab + Pb + Ti alike, analyses with 2.22 Å kernel crops.

| run | changes | measurements / unknown | question | data |
|---|---|---|---|---|
| R2_n400 | 400 + 400 iterations | 2.0 | was R2 short of convergence? (its misfit was still falling) | R2's copies |
| R2_n1000 | 1000 + 1000 iterations | 2.0 | ditto, ~12 h a leg | R2's copies |
| C512 | 512-px camera (2 × 2 sums, upsampling 1) | 7.8 | how many pixels does Ti need? | new copies of the kept sims |
| C128 | 128-px camera (8 × 8 sums, upsampling 3) | 0.5 | an EMPAD-class camera: does Pb survive? | new copies |
| R4 | ±133 mrad on 354 px (4 × 4, upsampling 2) | 2.0 | does more angle help at this pixel size? | new copies |
| S04 | 0.4 Å scan, 2500 positions | 3.1 | does a denser scan bring Ti back? (the brief's rule) | NEW lab / Pb / Ti sims, label `armf_a080_s04` |

Reading it: if Ti returns with iterations, R2 was under-converged; if it returns with C512 or S04 but not iterations,
the limit is information (measurements per unknown), and C512 against S04 (7.8 vs 3.1) says whether pixels or positions
buy it more cheaply; C128 bounds the low end; R4 isolates the angle. Not in this wave (a second wave if needed): 0.3 Å
scans (4489 positions: the coherence-averaged sim at 1420 px would need ~145 GB, past the sims' 128 GB), 512 px on the
dense scan, full pixels on the dense scan.
Tests for the wave: `test_upsampling_frame.py` now also binds 512 px (k = 1) and 128 px (k = 3) to the engine's binning
(scales 4 and 64, max rel diff < 6e-8); `test_det_driver.sh` 25 cases (512 / 128 px, NITER + RTAG on a reused copy, the
dense-scan label sim -> rebin -> recon chain); the whole wave stub-run (31 sbatch calls, every export checked).

## Where the results are

- **`results/`** (committed; rebuilt by `det_look.py` from the pulled folders):
  `det_scores.json` -- every run's numbers AND where each came from (analysis folder, recon leg, data dir, error trace,
  kernel check), with the scoring rules, each camera's own window and the hours per leg; `det_scores.csv` -- the same
  as one flat table; `det_summary.png` -- one row per run (camera, own window, Pb, Ti, hours); `det_limits.png` -- Pb and Ti
  right place against measurements per unknown and against iterations, guessing beside each point;
  `det_sections_cameras.png` (reference, R1, C512, R2, C128, R3) and `det_sections_limits.png` (R2 and its iteration
  runs, R4, S04) -- the depth sections with the found atoms; `phase/<run>_lab.jpg` -- each lab leg's phase image
  (projection, every slice, x-z), a 75 % copy.
- **Pulled tarballs, full size** (phase images, kernels, figure data, found atoms, logs): `~/Desktop/detector256_r1r2`
  (R1 + R2 reconstructions and sims; their first, faulty analyses), `~/Desktop/detector256_r1r2_kh` (R1 + R2 analyses
  with 2.22 Å kernels), `~/Desktop/detector256_r3`, `~/Desktop/detector256_limits`; the same tarballs stay on Blythe in
  `/springbrook/share/physics/phucrh/` (names in each submission log, `logs/detector256_*.txt`).
- **Reconstructions (h5, ~1 GB each) and sims stay on Blythe** (PACK_H5=0, CLEANDATA=0): `recon_af_<label>_<mode>_fs2p5_
  src0p1_<tag>_NL18/` and `sim_out_af_*` under the pipeline checkout.

## Results

See the status table at the top and `results/` (above). `results/det_summary.png` has one row per run: the camera, its
own window and the share of the probe outside it, Pb and Ti right place against guessing, and the hours a leg took.

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
| `submit_limits.sh` | the limits wave (six runs + analyses; guards as `submit.sh`; log `logs/detector256_limits_<ts>.txt`) |
| `det_look.py` | scores, eps on the data grid, vacuum ratio, measurements per unknown, kernel check, camera's own window, hours per leg; writes `results/` |
| `../../sim/rebin_detector.py` | the post-processing step (crop, B × B sums, probes per grid, sim_meta, record) |
| `HANDOFF_PROMPT.md` | the brief this work started from |
