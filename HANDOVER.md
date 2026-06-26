# Handover — abTEM → PtychoShelves multislice electron ptychography pipeline

End-to-end MVP that **simulates** 4D-STEM in Python (abTEM) and **reconstructs** it
in MATLAB (PtychoShelves, Yu Lei's `MultiHollowPtycho` engine) on the Warwick
**Blythe** HPC, then **validates** the reconstruction against the known model.

Status (as of this handover): the pipeline works end-to-end on a coherent,
noiseless PTO/STO labyrinth dataset. In-plane lattice, depth-dependent atomic
**displacement**, and along-beam **atomic-plane depth resolution** are all
recovered and validated against ground truth. Realism (frozen phonons + Poisson
dose) and finer scans are wired and ready but not yet run.

---

## 1. TL;DR — the one-paragraph mental model

abTEM propagates a 100 mrad overfocused 300 keV probe through the orthogonalized
PTO/STO cell (beam along z), records full diffraction patterns on a scan grid,
bins the detector 4×4, and writes them in the **exact PtychoShelves worked-example
format** (no conversion step). A MATLAB driver feeds them to Yu's multislice LSQ-ML
engine (`MLs`) with hollow angle = 0. Geometry is carried in `sim_meta.mat` so the
driver never hardcodes (and never drifts from) the sim. The detector needs a single
**transpose** (`custom_data_flip=[0,0,1]`) and nothing else. Everything is validated
by overlaying the reconstructed object on the abTEM model.

---

## 2. Codebase layout

```
ptyco_pipeline/                 (this repo; on Blythe under $SHARE/phucrh/ptyco_baseline/)
├── ptycho/                     Yu Lei's MultiHollowPtycho engine (DO NOT EDIT the engine)
│   ├── +engines/+GPU/...       LSQ-ML / MLs multislice GPU engine
│   ├── run_synthetic_recon.m       single-slice recon driver (geometry check)
│   ├── run_synthetic_recon_ML.m    MULTISLICE recon driver  <-- the main driver
│   ├── run_orientation_sweep.m     8-DOF detector-orientation test (phantom)
│   ├── run_baseline_PrScO3.m       shipped PrScO3 worked example, hollow=0 baseline
│   ├── ptychography_exp_30nm.m     Yu's original real-data example (reference only)
│   └── exampleData/01/             PrScO3 example data (folder MUST be named '01')
├── sim/
│   ├── simulate_4dstem.py          THE simulator (single file; abTEM 1.0.x new API)
│   ├── add_poisson_noise.py        post-sim Poisson/dose step (sweep dose, no re-sim)
│   ├── run_sim.sh                  sim launcher (descriptive folder + log inside)
│   ├── run_sim.slurm               sim SLURM worker (driven by run_sim.sh)
│   └── PTO6_STO6_18_18_labyrinthPoscar.vasp   the structure
├── analysis/                   reusable validation (run locally on pulled results)
│   ├── depth_resolution.py         kz plane-frequency test + convergence
│   ├── column_cross_section.py     depth cross-section down a column
│   └── column_cross_section_overlay.py   recon-vs-GT model overlay (the money figure)
├── run_recon_multi.sh          recon launcher (descriptive folder + log+outputs together)
├── run_recon_synthetic_ML.slurm    recon SLURM worker
├── DATA_FORMAT.md              the data contract (read this if touching I/O)
└── HANDOVER.md                 this file

Outputs (all gitignored, live on $SHARE, never in git):
  sim_out_<tag>/01/             one sim dataset (data_dp/position .hdf5, probe_initial.mat, sim_meta.mat)
  recon_<tag>/01/...            one recon (Niter*.mat objects) + slurm_<jobid>.out beside it
```

Top-level (NOT in the repo, on the Mac): `HPC_COMMANDS.md` (login/scp cheatsheet),
Desktop figures + `*_vol.npy` from analysis.

---

## 3. Geometry bindings — THE critical reference

Get any of these wrong and the object reconstructs mirrored/rotated/empty. All were
established empirically and are encoded in the drivers; documented here so they're
never re-derived from scratch.

### 3.1 Physics / calibration (coherent baseline)
| Quantity | Value | Notes |
|---|---|---|
| Beam energy | 300 keV | λ = 0.01969 Å |
| Convergence semi-angle | 100 mrad | very large (electron ptycho regime) |
| Defocus | **−20 Å** | = overfocus 20 Å (2 nm). Sign derived from abTEM source, see §3.6 |
| Detector outer angle | 200 mrad | full pattern, no crop |
| Detector binning | 4×4 (sum) | lazy Dask `coarsen`; `N_u≈1424 → N_b=356` |
| `Ndpx` (`p.asize`) | **356** | binned detector size |
| `d_alpha` | 1.1249 mrad/px | = `BIN_FACTOR·λ/box` |
| `rbf` (BF disk radius) | 88.9 px | = conv/d_alpha |
| Object pixel `dx` | **0.0492 Å** | = `1/(N_b·dk_b)`; set by diffraction angle, NOT scan step |
| In-plane box | 70.008 Å (square) | POSCAR rotated −90° y, orthogonalized, padded square |
| Beam thickness | 69.93 Å | sample along z; +2 Å vacuum each side (box ≈ 74 Å) |
| Scan | centre (40, 20) Å, 20 Å window, 0.1 Å step | ≈ 200×200 ≈ 40k positions |

These are **read by the driver from `sim_meta.mat`** — change them in the sim and the
recon follows automatically.

### 3.2 The data files (PtychoShelves contract — see DATA_FORMAT.md)
- `data_dp.hdf5`  `/dp`  — diffraction intensities. Loader (`matlab_aps`) wants MATLAB
  `dp[Ndet, Ndet, Npos]`, y-fastest, DC-centred.
- `data_position.hdf5`  `/probe_positions_0` — positions, MATLAB `[Npos, 2]` cols `[x, y]`.
- `probe_initial.mat` — `probe` (complex `[Ndet,Ndet]`) + `p` struct; consumed by MATLAB
  `load()` so it's a real `.mat` (scipy `savemat`), not h5py.
- `sim_meta.mat` — geometry the driver reads (Ndpx, d_alpha, rbf, dx, energy, thickness,
  scan_step, n_phonons, …).

### 3.3 The C-vs-Fortran landmine (most important single fact)
MATLAB `h5read` of an HDF5 dataset shape `(s0,s1,s2)` returns a `[s2,s1,s0]` array
(**axes reversed**), `M(p,q,r)=H[r-1,q-1,p-1]`. Therefore the sim writes:
- `/dp` as **`(Npos, Ndet_x, Ndet_y)` = `A.transpose(0,2,1)`** so MATLAB sees `[Ndet,Ndet,Npos]`.
- `/probe_positions_0` as **`(2, Npos)` = `pos_xy.T`** so MATLAB sees `[Npos,2]`.

The sim has a **round-trip self-test** (`selftest_ordering`) that re-reads the files and
emulates MATLAB's reversal — it catches this in Python without needing MATLAB.

### 3.4 Detector orientation — `custom_data_flip = [0,0,1]` (transpose)
Determined by the **8-DOF orientation sweep** (`run_orientation_sweep.m`) on an
asymmetric "F" Pb phantom: of the 8 dihedral orientations `[fliplr, flipud, transpose]`,
**only the pure transpose `[0,0,1]`** renders the F correctly. This is fixed in the
driver. (Scan/detector flatten y-fastest; the transpose reconciles abTEM's detector
axes with PtychoShelves'.)

### 3.5 Reconstruction ↔ ground-truth coordinate map (for validation/overlay)
Locked by NCC of the recon depth-sum vs the abTEM projected potential (corr ≈ 0.79):
- **In-plane:** recon `(row r, col c)` → GT physical **`X = 30 + c·0.0492`, `Y = 10 + r·0.0492`**
  (i.e. recon is the transpose of the GT grid; zero shift to integer px, ~1 px sub-pixel
  residual calibrated out in the overlay script).
- **Depth:** `z_recon ≈ z_GT` (entrance at layer 0; fitted offset +0.5 Å ≈ 0, physical).
- These live in `analysis/column_cross_section_overlay.py` — re-use them for any
  recon-vs-model comparison.

### 3.6 Defocus sign (why −20 Å = overfocus)
From abTEM 1.0.x source: χ uses `+½α²C10`, abTEM `defocus = −C10`, transfer
`T=exp(−iχ)=exp(+iπλ·defocus·k²)`; the Fresnel propagator is `exp(−iπλ·dz·k²)`, so
`T = P(dz=−defocus)`. A **negative** defocus builds the entrance wave as the in-focus
probe propagated forward by `|defocus|` → crossover **above** the entrance surface =
**overfocus**. PtychoShelves' ASM propagator has the same sign (verified). So overfocus
magnitude 20 Å ⇒ `defocus = −20`.

### 3.7 Dose / intensity scaling
abTEM flux-normalises each pattern to total ≈ 1; PtychoShelves' `load_from_p` rejects
avg counts < 1e-4. The **coherent** sim scales to `DOSE_E=1e10 e/pattern` (≈ noiseless).
For realism, leave the sim noiseless and apply **`add_poisson_noise.py --dose <e/Å²>`**
afterwards (incident e/pattern = dose·step², relative totals preserved so TDS/ADF
contrast survives).

---

## 4. Tests run + results

| Test | What | Result |
|---|---|---|
| **Baseline PrScO3** | shipped worked example, hollow=0, through Yu's engine | matched the shipped reference (probe modes + lattice) — engine + our SLURM setup validated |
| **8-DOF orientation sweep** | F phantom, all 8 `custom_data_flip` | only `[0,0,1]` (transpose) renders F correctly |
| **In-plane fidelity** | recon depth-sum vs abTEM projected potential | corr **0.84–0.85** |
| **Depth resolution (kz)** | on-column z-power spectrum, plane peak | NL70: clean peak at **0.257 Å⁻¹ (3.9 Å planes), 4.0× prominence**; NL42 marginal (1.8×, at Nyquist) |
| **Pb vs Ti displacement** | per-column lateral wander vs depth | recon recovers it; GT **Pb 0.29 Å > Ti 0.20 Å** (A-site lone pair), matched |
| **Model overlay** | recon column cross-section + GT atoms | Pb/Ti markers land on the reconstructed blobs down the whole column, including the lean/displacement and a domain-wall kink |
| **70-layer stability** | reg-off, single-mode, fixed probe | converged (final error 28.9 < NL42's 31.8); no blow-up |

Key validated claims you can stand behind: in-plane lattice resolved; depth-dependent
**displacement** real (Pb > Ti, matches model); **atomic planes resolved along the beam**
at ~3.9 Å (NL70). Oxygen is **not** cleanly resolved (Z=8; faint signal at the right
place) — that's a stretch goal / its own project.

---

## 5. What went well / key learnings
- Writing straight into the worked-example format (no converter) removed a whole class
  of bugs; the Python round-trip self-test pays for itself.
- Carrying geometry in `sim_meta.mat` means the driver never drifts from the sim.
- The probe is **known** (we simulated it) → reconstruct with **1 probe mode, fixed**.
  Earlier NaN-at-iter-20 blow-ups were under-constraint at probe-update onset, not a
  bad probe — fixing the probe (`PROBE_START=inf`) and `beta_LSQ=0.1` cured it.
- Depth resolution needs **fine slicing** (1 Å / 70 layers) to put the 3.9 Å plane
  frequency below Nyquist; 42 layers sit right at their limit.
- The recon's **last layer is a dumping ground** (exit-surface artifact) — trim it for
  figures; consider extra exit vacuum in future sims.

---

## 6. How to use it — the workflow

Three stages, each producing a **self-describing folder** (params in the name, log
inside). All commands run from the repo root **on a Blythe login node** unless noted.

```bash
# 1. SIMULATE  -> sim_out_step<step>_slice<slice>_<coherent|ph<N>s<sig>>/
bash sim/run_sim.sh                                   # coherent, 0.1 Å, default params
SCAN_STEP=0.05 PHONONS=8 WALLTIME=1-00:00:00 bash sim/run_sim.sh   # finer + TDS

# 2. ADD DOSE (post-sim, cheap, repeatable)  -> ..._dose<D>/
python sim/add_poisson_noise.py --in-dir sim_out_step0.1_slice2_ph8s0.08 --dose 1e5

# 3. RECONSTRUCT  -> recon_<simref>_NL<n>_reg<r>_p<modes>_b<beta>/
SIM_SRC=sim_out_step0.1_slice2_ph8s0.08_dose1e5/01 REGLAYER=0 PROBE_MODES=1 \
    NITER=120 WALLTIME=2-00:00:00 bash run_recon_multi.sh 70

# pull results to the Mac (see §8), then validate locally:
python analysis/depth_resolution.py
python analysis/column_cross_section_overlay.py
```

---

## 7. Key levers (env knobs)

**Sim** (`sim/run_sim.sh`):
| Var | Default | Effect |
|---|---|---|
| `SCAN_STEP` | 0.1 | probe-position spacing (Å). Smaller = more overlap = better, 4× cost at 0.05 |
| `SLICE_THICKNESS` | 2 | abTEM potential slice (Å); accuracy, not depth info |
| `PHONONS` | 0 | frozen-phonon configs (TDS). 0=coherent; 8–16 realistic; N× cost |
| `PHONON_SIGMA` | 0.08 | rms thermal displacement (Å) |
| `WALLTIME` | (12 h) | bump for phonon sims |
| `SIM_TAG`, `OVERWRITE` | — | force a name / allow overwrite |

**Dose** (`sim/add_poisson_noise.py`): `--dose` (e/Å²; ~1e4–1e6 sensible), `--seed`.

**Recon** (`run_recon_multi.sh`, args = Nlayers list):
| Var | Default | Effect |
|---|---|---|
| `SIM_SRC` | sim_out/01 | which sim dataset to reconstruct |
| `REGLAYER` | [1,0.5] | depth regularizer; **0 = off** (use for sharp depth) |
| `PROBE_MODES` | 1 | keep 1 (probe is known/coherent) |
| `BETA_LSQ` | 0.1 | LSQ step; lower = more stable |
| `PROBE_START` | inf | probe-update onset; inf = **fixed probe** |
| `NITER` | [200,200] | iterations per engine (lower for deep/fine runs to fit walltime) |
| `WALLTIME` | (24 h) | bump for 70+ layers / 0.05 Å scans |

**Hollow angle** (the engine's defining feature) is **0** here (standard baseline):
`mask1 = ones(Ndpx)` in the driver. Non-zero would mask the BF disk.

---

## 8. Commands that need a human (HPC — I cannot run these)

I can edit code and analyse pulled results, but **cannot** ssh, submit SLURM jobs,
or scp. You drive these (cheatsheet in top-level `HPC_COMMANDS.md`):

```bash
# login
ssh phucrh@blythe.scrtp.warwick.ac.uk

# get the latest code onto Blythe
cd $SHARE/phucrh/ptyco_baseline/ptyco_pipeline && git pull

# submit (examples in §6).  Monitor:
squeue -u $USER
sacct -u $USER --starttime today --format=JobID,JobName%30,State,Elapsed

# pull results to the Mac (SFTP is disabled on Blythe -> use scp -O):
scp -O -r 'phucrh@blythe.scrtp.warwick.ac.uk:.../ptyco_pipeline/recon_<tag>/01/*step02*' \
    ~/Desktop/recon_new/<tag>/
```
Connection details: user `phucrh`, host `blythe.scrtp.warwick.ac.uk`, repo at
`/springbrook/share/physics/phucrh/ptyco_baseline/ptyco_pipeline`, MATLAB module
`MATLAB/2024b_Update_3`, GPU `gpu:lovelace_l40:1`, conda env `$SHARE/phucrh/envs/abtem`.

When something needs my eyes, paste the SLURM log (`recon_<tag>/slurm_<jobid>.out`) or
scp the `Niter*.mat` and I'll analyse it.

---

## 9. Gotchas (hit + fixed — don't rediscover)
- **conda activate is unreliable in batch** → call `$CONDA_ENV/bin/python` by absolute path.
- **scp plain fails** ("subsystem request failed") → SFTP disabled, use **`scp -O`**.
- **Repo root in SLURM** → use `$SLURM_SUBMIT_DIR`, not `$BASH_SOURCE` (spool dir).
- **Scan folder must be named `01`** (matches `scan.format='%02d'`), not `scan01`.
- **NaN at iter ~20** → probe-update onset under-constraint; fix the probe + low beta.
- **load_from_p photon-count error** → scale dose above 1e-4 avg counts.
- Engine code (`ptycho/+engines`, `ptychography_exp_30nm.m`) is **Yu's — leave it intact**;
  the two-engine presolve (128→256, here 178→356) is his scheme, only the px rescaled.

---

## 10. Recommended next steps
1. **Realistic baseline (paper backbone):** frozen-phonon sim on the *validated* 0.1 Å
   geometry → dose → NL70 recon. Re-run §4's checks: do planes + displacement survive
   TDS + shot noise? (Verify the first phonon run logs `frozen phonons ON` on abTEM 1.0.9.)
2. **0.05 Å "ceiling" run** if compute allows — the oxygen stretch (more overlap → SNR
   for PSF-deconvolution atom-finding). ~4× data/compute (~80 GB binned).
3. Consider **exit vacuum padding** to kill the last-layer artifact at the source.
4. Not done deliberately: scan **jitter** (position-refinement demo) and **partial
   coherence / source size** — easy future reviewer-proofing if needed.
