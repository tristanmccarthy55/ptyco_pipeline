# Prompt for the next session -- can a 256-pixel detector do the 80 mrad runs? (written 2026-10-23)

Paste the block below into a new session.

```text
You are testing whether a 256 x 256-pixel detector is enough for the PX915 aperture campaign's 80 mrad
reconstructions, by modelling the detector pixels in the reconstruction (reciprocal-space upsampling).
Repo: ~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PtycoShelves/ptychoshelves-clean (origin/main, push
allowed). Blythe checkout: /springbrook/share/physics/phucrh/ptyco_baseline/ptyco_pipeline (updated by git pull). Keep notes
and results in aberration_experiment/detector256/ (a README.md with the status at the top).

WHY. A referee's first question about the paper: the 80 mrad runs reconstruct 1419 x 1419-pixel patterns recorded to
+-133 mrad, while lab 4D-STEM cameras are 128-256 px across (Lei & Wang's EMPAD was 128). The window the solver models is
N lambda / (2 theta_max), lambda = 0.0197 A at 300 kV, and it must hold the probe. At 256 px that window is 19 A at
+-133 mrad, 26 A at +-96 mrad, 32 A at +-80 mrad; the ARM-tableau probe at 80 mrad is 74 A across (d90) and needs the
105 A window. So the pixel count is set by probe size x detector angle, not by alpha itself. The user's goal: 256 x 256
if it works.

THE IDEA. Record fewer, bigger pixels and let the solver model them: compute each pattern on a 2^k-finer grid and sum it
into the detector pixels before comparing with the data. The missing angular detail is recovered from the dense scan's
overlap (Batey et al., Phys. Rev. A 89, 043812, 2014). PtychoShelves already has this, off in every run so far:
p.detector.upsampling = k. Read the code path before using it:
 - ptycho/+core/initial_checks.m ~l.214: upsampling switches the data loader to 'matlab_ps' (run_synthetic_recon_ML.m
   uses 'matlab_aps'); check our synthetic data load through it.
 - ptycho/+core/ptycho_prepare_scans.m ~l.212: apply_binning(p, 2^-k) rescales the geometry "except data and mask".
 - ptycho/+engines/+GPU/+initialize/init_solver.m ~l.146: the data are binned back to detector size.
 - ptycho/+engines/+GPU/private/get_reciprocal_model.m ~l.96: |Psi|^2 on the fine grid, binned with utils.binning_2D
   before the comparison (the fast amplitude path is switched off).
 - ptycho/+engines/+GPU/private/modulus_constraint.m ~l.88: unbinned back to the fine grid.
 - +utils/binning_2D.m takes a `centered` argument (default false); the engine bins uncentred blocks.
Why the details matter: the campaign once summed 4x4 detector pixels WITHOUT modelling it in the solver and paid with a
residual floor and a phase ramp (the blocks started on the zero-angle pixel and shifted every pattern by (BIN-1)/2 fine
pixels; HANDOFF_ANALYSIS.md "Changed 2026-09-23/24" item 2; DETECTOR_SAMPLING in sim/simulate_4dstem.py). Upsampling
models the summing, so the data must be binned exactly the way the engine bins its model: same block alignment, same
array frame (FFT-shifted or not, after the engine's data flip and transpose). Prove this with a local unit test in Python
against the engine's code before any cluster job.

THE TEST. One change at a time, true probe held fixed, all on the same simulation.
 Data: the kept joint sim at the target setting (focal spread 2.5 A + source 0.1 A, ARM tableau, 80 mrad, noiseless):
 sim_out_af_armf_a080_lab_fs2p5_src0p1 under the Blythe checkout. Confirm it still exists first: a disk cleanup ran,
 and PHASE1.md says the joint sims were kept. Reference, already done: the fixed-probe control at that point,
 Pb / Ti right place 91 / 60 % (good) -- aberration_experiment/mono300/PHASE1.md, "A, B, C results".
 - R1: crop the detector symmetrically to +-96 mrad (1419 -> 1024 px, zero-angle pixel centred in the engine's
   convention), full pixels. Isolates the cost of recording less angle (ties to the zblur G0c question of whether an
   atom's depth core is set by the recorded angle).
 - R2: R1's 1024 px binned 4 x 4 -> 256 px, reconstructed with upsampling k = 2 at the same 105 A window. THE TEST.
 - R3 (cheap, optional): the 256 px data with no upsampling (26 A window): the expected failure, showing what
   upsampling buys.
 - R4 (only if R2 works): +-133 mrad binned 4 (crop 1419 -> 1416, bin -> 354) with upsampling 2: was the angle cut needed?
 To build:
 - a post-processing step in the style of sim/add_poisson_noise.py (reads data_dp.hdf5 in chunks, writes a new tagged
   dataset, links positions / probe / meta), e.g. sim/rebin_detector.py --crop-mrad 96 --bin 4;
 - the true probe rebuilt for the new geometry (the crop changes the object pixel from 0.074 to 0.1025 A); find out at
   which size the engine wants probe_initial when upsampling is on;
 - a switch (e.g. UPSAMPLING=2) through run_synthetic_recon_ML.m and campaign/run_thin_atomfind.sh so the legs run
   through the normal chain;
 - the presolve: it crops the detector to half width and is guarded to 1.2 alpha (run_synthetic_recon_ML.m); at
   +-96 mrad = 1.2 alpha it has nothing left to crop. Decide what it should do with upsampling, or switch it off for this
   test and say so.
 Kernels: the finder's Pb / Ti kernels must come from the identical pipeline (aberration_experiment/PSF_KERNELS.md):
 apply the same crop / bin / upsampling to this point's Pb and Ti grid sims if they were kept. If not, judge R1 / R2 on
 the phase images and the misfit first and label any score made with the reference kernels.

HOW TO JUDGE. Phase images first (analysis/render_phase.py; analyse_sweep writes <out>/phase/). Then the misfit
eps = final residual x N / 2e5 (N = pattern width; compare like with like), the vacuum / slab phase ratio (0.31 for a clean
run), and right-place scores with make_arm_page.atom_scores (right element within 0.5 A in depth, against guessing;
good = Pb >= 80 and Ti >= 45 %, okay = Pb >= 50). Success: R2 within a few points of the reference and as clean in the
images. Rough redundancy (measurements / real unknowns, 125 A object, 18 slices): full pixels ~30:1; 256 px at +-133 mrad
~1:1; at +-96 mrad ~2:1; with a 0.4 A step (2500 positions) ~3:1. If R2 fails, try a denser scan before anything else.
 Caveat to record beside any result: data made by summing 4x4 point samples, fitted by a model that sums 4x4 samples, is an
 idealised detector (real pixels integrate continuously and blur: MTF / DQE).
 Cost: each leg is about one current 80 mrad reconstruction, a few hours (the forward model still runs on the fine
 1024 grid); only the camera and the data shrink (13 GB -> 0.4 GB per scan).

HOW WE WORK. The user runs every cluster job: give exact command blocks, each with its scp -O pull line (glob the stable
part of the tarball name) and tar on its own line into a fresh directory; never submit anything yourself. Dry-run every
driver change with a stub sbatch, plus campaign/test_driver_legacy.sh and campaign/test_coh_names.sh, before handing it
over. No MATLAB locally: test the binning and frame logic in Python against the engine's code. Local Python:
~/hyperspy-bundle/bin/python. Look at phase images before judging; honest numbers only; plain language for a
professor-level reader. Other sessions have uncommitted work in the repo: stage files by name. Read first: the memory
notes target-mono300, aberration-campaign, dense-scan-not-big-window and blythe-hardware; then
aberration_experiment/mono300/PHASE1.md (top), the STATE block of HANDOFF_MONO300.md, and PSF_KERNELS.md.
```
