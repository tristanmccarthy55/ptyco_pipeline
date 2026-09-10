# fusion — hollow-detector data fusion: ptychography for the direction, EELS for the magnitude

Proof of concept for the idea the PX915 report ends on. Projected ptychography measures the
in-plane polarisation superbly (0.005–0.007 Å, 0.9° direction — `analysis/atomfind`) but **cannot
measure the along-beam component at all**; EELS measures the along-beam **magnitude** but is blind
to its **sign**. One 4D-STEM scan through a **hollow** detector produces both at once, and together
they give the 3-D vector neither can.

Everything here is driven by a single simulation. The hole is applied at *reconstruction* time as a
mask, so sweeping the hollow semi-angle costs reconstructions, not simulations.

> Reference for the hollow geometry: **Yu Lei and Peng Wang**, *Multislice hollow ptychography for
> simultaneous atomic-layer-resolved 3D structural imaging and spectroscopy*, arXiv:2506.22352 —
> the MultiHollowPtycho engine in `ptycho/` is theirs (`readme.txt`, `ptychography_exp_30nm.m`).

---

## Status (2026-09-10) — the EELS channel works, the depth readout does not yet

First full run on Blythe: one simulation (6400 positions, noiseless at 1e10 e/pattern) and four
hollow reconstructions at HSA = 0, 0.50, 0.75, 0.95, all completed in ~36 min each.

**What worked.** Everything up to the reconstruction. The measured dose split at 0.75α is
55.5 % / 44.5 %; the simultaneous HAADF carries 1.0 % of the beam and, as it must, cannot separate
the domains at all (0.99 % spread). The reconstructions resolve the lattice beautifully in
projection — Pb, Ti–O and even the weak O columns, lattice coherence 0.98, Pb:Ti amplitude 1.85.

**What did not.** The reconstructions carry **no depth structure whatsoever**: power at the 4.152 Å
lattice period along the beam is 0.5–0.9× the mean band power, i.e. no atomic planes, even after
detrending, even at HSA = 0. So the sign readout had nothing to work with and scored at chance
(29–56 %). `check_recon.py` reports this in one command and should be run before any physics is
read off a reconstruction.

**It is not a dose problem and not an information problem.** The data were noiseless throughout.
`sign_encoding.py` compares two uniform membranes differing only in sign(δ_z) at identical probe
positions: at 5 cells, **4.9 % of the diffraction pattern differs**, 2.3 % of it outside the hole
(65 % of the Poisson information). The null control — δ purely in-plane, so the flip is a no-op —
gives exactly 0.000e+00. The sign is abundantly present in the data; reconstructing 24 free depth
layers from a 21 Å slab is simply a lossy way to get at it.

**Open.** Whether a better-conditioned reconstruction (fewer layers, more iterations, a known-object
start) recovers the depth, or whether the sign should be read directly from the diffraction data as
a two-hypothesis test — EELS supplies |δ_z| and projection δ_xy, leaving exactly two candidates.
Note also that the engine's `fourier_error_out` scores every detector pixel while the hollow mask is
applied only inside `modulus_constraint.m`, so it is **not** a valid convergence monitor for a
hollow run.

---

## The physics, in three facts

**1. Projected ptychography is exactly degenerate in the sign of P_z.** An up cell and a down cell
have the *same* projected potential, atom for atom — only their stacking order along the beam
differs. Breaking that degeneracy needs depth sectioning, which is why this is a **multislice**
experiment and not a single-slice one.

**2. Dipole EELS is exactly degenerate in the sign of P_z.** The O-K anisotropy enters through the
polar displacement *along the Ti–O–Ti chain*, as |δ·ĉ| — even in δ. Measured here as 0.000%
contrast between the up and down domains (`eels_forward.py`), not merely small.

**3. The 100 mrad probe is not the problem — and the hole costs the spectroscopy nothing.**
The report flagged that a convergent probe washes out the ELNES anisotropy (RESULTS.md M6c/M6d).
It does, to about a third of intrinsic and with inverted sign, but:

| geometry | ⟨cos²(q,z)⟩ | AB↔CD contrast | counts/channel (SNR 3) |
|---|---|---|---|
| α=100, hole 0.10α | 0.110 | 4.89 % | 7.5×10³ |
| α=100, hole 0.50α | 0.121 | 4.66 % | 8.3×10³ |
| **α=100, hole 0.75α** | **0.120** | **4.68 %** | **8.2×10³** |
| α=100, hole 0.95α | 0.121 | 4.67 % | 8.3×10³ |
| near-parallel α=β=2 | 0.500 | (limit 5.3 %) | — |

The contrast is **flat to within 8% across the entire range of hole sizes**, because beyond
α ≈ 5 mrad the convergence cone alone already sets the momentum-transfer distribution — closing the
spectrometer aperture cannot undo it. Two consequences:

* **Size the hole from the ptychography side alone.** Lei & Wang: lateral resolution holds to
  ~0.75α and degrades past ~0.90α (sub-Å to 0.95α at ≥10⁵ e/Å²). Hence the 0.75α default.
* **The convergent probe is not a handicap here.** A near-parallel beam gives a *smaller* limiting
  contrast (5.3% vs 6.8%) because |−0.32| > |+0.25|: past the magic angle the anisotropy inverts
  and saturates rather than vanishing. Reproduce with `eels_forward.py --sweep`.

Where the beam actually goes, measured on the simulated patterns (not estimated):

| hole | → EELS | → ptychography |
|---|---|---|
| 0.50α | 24.3 % | 75.7 % |
| **0.75α** | **55.5 %** | **44.5 %** |
| 0.95α | 89.8 % | 10.2 % |

plus a simultaneous virtual HAADF (100–200 mrad) carrying 1.0 % of the beam. Nothing is lost to the
kept pixels: masking discards *information*, not counts, which is why MHP keeps sub-Å resolution.

---

## The specimen

A free-standing PbTiO₃ membrane, 12×12×5 unit cells (46.8 × 46.8 × 21.4 Å), beam along the polar
axis, PbO-terminated on **both** surfaces so there is no composition asymmetry to confound the sign
test. Every column is identical over all 5 cells along the beam, so EELS reads the polarisation
where the probe lands with no depth deconvolution.

The lattice is tetragonal with c ‖ z **everywhere** (an in-plane-clamped film — what a PTO/STO
superlattice membrane is); only the polar displacement rotates. That keeps the box commensurate and
removes every confound (thickness, strain, composition), so any measured contrast is polarisation.

|δ_Ti| = 0.331 Å in all four quadrant domains; only the direction differs:

| domain | δ (Å) | θ from beam | what it tests |
|---|---|---|---|
| **A** | (+0.113, 0, **+0.311**) | 20° | |
| **B** | (+0.113, 0, **−0.311**) | 20° | **identical EELS *and* identical projection to A** — only depth sectioning separates them |
| **C** | (+0.287, 0, **+0.165**) | 60° | separated from A/B by the EELS magnitude |
| **D** | (0, +0.287, **−0.165**) | 60° | identical EELS to C (a round aperture cannot tell an x-chain from a y-chain) |

---

## How the sign is read out

Cations and anions move in **opposite** directions under the polar distortion, so the three distinct
columns of a unit cell — Pb at (0,0), Ti + apical O at (a/2,a/2), equatorial O at (a/2,0) — shift
apart along the beam by up to **0.92 Å**, far more than the 0.33 Å displacement itself. EELS has
already supplied |δ_z| and projected ptychography δ_xy, which leaves exactly **two** candidate
structures; a matched filter on the depth profiles asks the reconstruction which one it matches.
Only the sign of the answer is used, so no amplitude calibration is needed, and the reconstruction's
arbitrary depth origin is fitted once, globally, by a criterion symmetric in the two hypotheses.

Measured robustness (`analyze_fusion.py --robustness`):

| | per-cell sign correct | per-domain sign |
|---|---|---|
| profile noise 0.5 / 0.7 / 1.0 / 1.5 × RMS | 100 / 100 / 99 / 92 % | 100 % throughout |
| true depth FWHM 2 / 4 / 5 / 6 / 8 Å | 100 / 100 / 92 / 76 / 69 % | 100 / 100 / 100 / 95 / 85 % |

λ/α² at 100 mrad is ≈ **2.0 Å**, and Lei & Wang measure 2.8 Å at 80 mrad — so the readout has
roughly a factor of two of margin in the quantity that matters.

---

## Running it

### 1. Build the specimen (local, seconds)

```bash
cd ptychoshelves-clean/fusion
~/hyperspy-bundle/bin/python toy_sample.py --preview
```
Every gate must pass (atom count, both surfaces PbO, δ exact in every domain interior, column
homogeneity along the beam, one undistorted lattice). Writes `sample/toy_membrane.vasp` +
`sample/toy_truth.npz`.

### 2. Cache the CASTEP ELNES ladder (local, once, ~4 s)

```bash
~/hyperspy-bundle/bin/python eels_forward.py --cache      # needs eels/runs/*_core_edge.dat
~/hyperspy-bundle/bin/python eels_forward.py --sweep      # the hole-size answer above
```

### 3. Check feasibility before spending GPU time (local)

```bash
~/hyperspy-bundle/bin/python analyze_fusion.py --predict
~/hyperspy-bundle/bin/python analyze_fusion.py --selftest
~/hyperspy-bundle/bin/python analyze_fusion.py --robustness
~/hyperspy-bundle/bin/python make_figure.py               # headline figure, synthetic recon
```

### 4. Simulate on Blythe (GPU)

```bash
cd $SHARE/phucrh/ptyco_baseline/ptyco_pipeline && git pull
sbatch --output=fusion_sim_%j.out --error=fusion_sim_%j.err fusion/run_fusion_sim.slurm
# heavier: SCAN_STEP=0.2 SCAN_WINDOW=32 PHONONS=8 sbatch ... (add WALLTIME as needed)
```
Writes `fusion/runs/fusion/01/` (PtychoShelves inputs + one mask per hollow angle),
`hollow_budget.json`, `haadf.npy`, `bf.npy`.

Re-derive the split for other hole sizes with **no re-simulation**:
```bash
$SHARE/phucrh/envs/abtem/bin/python fusion/simulate_fusion.py --rebudget \
    --out-dir fusion/runs/fusion --hsa 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.75 0.8 0.9 0.95 1.0
```

### 5. Hollow reconstructions (GPU, MATLAB)

```bash
bash fusion/run_fusion_recon.sh 0.75              # the operating point
bash fusion/run_fusion_recon.sh 0 0.50 0.75 0.95  # the whole sweep, one sim
```

The driver defaults **`regularize_layers = [0 0]`**, unlike its parent. `regulation_multilayers.m`
is a missing-cone low-pass in k_z: it blurs the depth axis, which is precisely where the sign of
P_z lives. Only raise it (≤ 0.05) to rescue a divergent run. `probe_modes = 1` (the simulated probe
is coherent) and the probe and positions are held fixed at their known values. If a deep run NaNs,
drop `BETA_LSQ` 0.1 → 0.05 → 0.02 before touching `REGLAYER`. Each job prints a preflight block
with all of these so the log is self-documenting.

### 6. Read it out and make the figure

```bash
~/hyperspy-bundle/bin/python fusion/make_figure.py \
    --recon fusion/runs/recon_hsa0.75_NL24/01/*/Niter*.mat \
    --budget fusion/runs/fusion/hollow_budget.json \
    --out fusion/fusion_headline.png
```

---

## Files

| file | what it is |
|---|---|
| `toy_sample.py` | the four-domain membrane + ground truth, with the build gates |
| `eels_forward.py` | per-column O-K spectra: chain-resolved tensor × convergent-aperture average × the M5 displacement ladder |
| `simulate_fusion.py` | one abTEM 4D-STEM scan → PtychoShelves inputs, hollow masks, measured dose split, virtual HAADF/BF |
| `../ptycho/run_fusion_hollow.m` | the MHP driver: identical to `run_synthetic_recon_ML.m` but with `mask1` = the hole |
| `analyze_fusion.py` | the matched-filter sign readout, the EELS magnitude inversion, and the fusion |
| `check_recon.py` | did a reconstruction recover DEPTH? run this before reading physics off one |
| `sign_encoding.py` | is the sign in the data at all, and how does it grow with thickness (noiseless) |
| `make_figure.py` | the headline figure |
| `run_fusion_sim.slurm`, `run_fusion_recon.sh` | Blythe launchers |

**Nothing specific to the PTO/STO labyrinth was touched.** `sim/simulate_4dstem.py` is imported and
driven through its own module globals, never edited, so its rotation, scan region and defaults are
exactly as validated; `ptycho/run_synthetic_recon_ML.m` is untouched and `run_fusion_hollow.m` is a
sibling that differs only in the mask block, the depth-layer default and the output tag.

## Caveats

* The EELS channel is a **kinematic** forward model: first-principles unit-cell spectra (CASTEP +
  OptaDOS core hole) combined by the validated uniaxial tensor form and averaged over the α and β
  cones. It carries no channelling and no thickness dependence. The dynamical route
  (`eels/simulate_stem_eels.py`, abTEM transition potentials + gpaw) exists and runs; folding it in
  as the per-column weight is the natural next step.
* The 4.7% contrast is a **domain-scale** measurement at realistic dose, not a per-probe-position
  one; bin over a domain before inverting for θ.
* Frozen phonons are off by default (`--phonons 0`). Turn them on for the dose-realistic run.
* The synthetic-reconstruction fallback in `make_figure.py` is a *model*, clearly labelled in the
  panel; it exists so the figure and the readout are validated before the HPC run, not to stand in
  for one.
