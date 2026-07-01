# Handover — abTEM → PtychoShelves multislice electron ptychography pipeline

End-to-end MVP that **simulates** 4D-STEM in Python (abTEM) and **reconstructs** it
in MATLAB (PtychoShelves, Yu Lei's `MultiHollowPtycho` engine) on the Warwick
**Blythe** HPC, then **validates** the reconstruction against the known model.

Status (as of this handover): the pipeline works end-to-end and is validated on a
coherent, noiseless dataset — in-plane lattice, depth-dependent atomic
**displacement**, and along-beam **atomic-plane depth resolution** all recovered and
matched to ground truth. Realism is now built and running: **frozen phonons (TDS)**,
**Poisson/dose**, and **scan tiling** for arbitrarily-large sims are implemented and
validated; a **0.05 Å-step, 16-phonon "reviewer-2" dataset** is simulating on Blythe
(30 tiles). The reconstruction of that large dataset is the current open task (it's a
monolithic job that exceeds the 2-day wall — see §10).

> **If you are the atom-finding agent:** jump to **§11 — For the atom-finding agent**.
> It has the local data paths, how to load the reconstructed 3D object, the
> recon↔model coordinate map, the ground-truth atom positions, and the
> PSF-deconvolution plan. §3 (geometry) and §11 are the two you need; the rest is about
> generating/reconstructing the data (a separate thread owns that).

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
│   ├── simulate_4dstem.py          THE simulator (abTEM 1.0.x; phonons processed 1-at-a-time)
│   ├── add_poisson_noise.py        post-sim Poisson/dose step (streamed; sweep dose, no re-sim)
│   ├── run_sim.sh                  single-job sim launcher (descriptive folder + log inside)
│   ├── run_sim_tiled.sh            TILED sim launcher: scan split into N array jobs + merge
│   ├── merge_tiles.py              reassemble scan tiles into one dataset (streamed)
│   ├── run_sim.slurm               sim SLURM worker (driven by both launchers; tile-aware)
│   └── PTO6_STO6_18_18_labyrinthPoscar.vasp   the structure (ground truth)
├── analysis/                   reusable validation (run LOCALLY on pulled results)
│   ├── depth_resolution.py         kz plane-frequency test + the recon-mat loader
│   ├── column_cross_section.py     depth cross-section down a column
│   ├── column_cross_section_overlay.py   recon-vs-GT model overlay (the money figure)
│   └── dose_compare.py             recon cross-section + slice vs electron dose
├── run_recon_multi.sh          recon launcher (descriptive folder + log+outputs together)
├── run_dose_series.sh          make noisy copies at several doses + launch a recon for each
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
| Scan | centre (40, 20) Å, 20 Å window | step varies (see below) |
| Object grid (ROI) | 404 × 404 px | = the reconstructed in-plane field, `object_roi` |

Scan step is the one knob that changes between datasets:
- **Validated recons (NL42/NL70, the cross-section figures): 0.15 Å step → ~18k positions**,
  coherent, 2 Å multislice slices. (This is the data the atom-finding agent has — §11.)
- **Reviewer-2 run (simulating now): 0.05 Å step → 160k positions**, 16 phonons, 0.5 Å slices.

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
avg counts < 1e-4. The sim applies a **fixed** `×DOSE_E=1e10` (not normalised by the
run's own mean — that is deliberate, so independently-simulated scan tiles share one
scale and merge seamlessly), giving a noiseless dataset. For realism, leave the sim
noiseless and apply **`add_poisson_noise.py --dose <e/Å²>`** afterwards (incident
e/pattern = dose·step², relative totals preserved so TDS/ADF contrast survives; the
`1e10` is renormalised away, so the noisy dose is set purely by `--dose`).

### 3.8 Frozen phonons (TDS) — processed one config at a time
`--phonons N` averages the diffraction **intensities** over N randomly-displaced atomic
configs (σ default 0.08 Å). Configs are simulated **sequentially** (build that config's
potential → scan → bin → accumulate), so peak RAM ≈ a single coherent sim (~14 GB)
regardless of N — wrapping `FrozenPhonons` in the `Potential` instead makes abTEM hold
all N at once and OOMs at 16. Same `PHONON_SEED` across tiles ⇒ identical configs ⇒
seamless tiling. TDS shows as a diffuse dark-field background (validated: ~5–12% DF lift
vs a matched coherent sim; BF disk unchanged).

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
# 1a. SIMULATE (single job) -> sim_out_step<step>_slice<slice>_<coherent|ph<N>s<sig>>/
bash sim/run_sim.sh                                              # coherent, defaults
SCAN_STEP=0.1 SLICE_THICKNESS=0.5 PHONONS=8 WALLTIME=1-00:00:00 bash sim/run_sim.sh

# 1b. SIMULATE (TILED — for big/16-phonon runs that exceed the 2-day wall)
#     splits the scan into TILES array jobs (<=MAXPARALLEL at once) + a merge job
SCAN_STEP=0.05 SLICE_THICKNESS=0.5 PHONONS=16 TILES=30 WALLTIME=18:00:00 \
    bash sim/run_sim_tiled.sh
rm -rf sim_out_step0.05_slice0.5_ph16s0.08/tiles                 # after merge: reclaim ~81 GB

# 2. ADD DOSE (post-sim, streamed, cheap, repeatable) -> ..._dose<D>/
python sim/add_poisson_noise.py --in-dir sim_out_step0.05_slice0.5_ph16s0.08 --dose 1e8

# 3. RECONSTRUCT -> recon_<simref>_NL<n>_reg<r>_p<modes>_b<beta>/   (log+outputs together)
SIM_SRC=sim_out_step0.05_slice0.5_ph16s0.08_dose1e8/01 REGLAYER=0 PROBE_MODES=1 \
    NITER=120 WALLTIME=2-00:00:00 bash run_recon_multi.sh 70

#     dose SERIES in one shot (a recon per dose, for the degradation figure):
SIM_SRC=sim_out/01 DOSES="1e10 1e8 1e6 1e4" NL=70 bash run_dose_series.sh

# 4. pull results to the Mac (see §8), then validate LOCALLY:
python analysis/depth_resolution.py
python analysis/column_cross_section_overlay.py
python analysis/dose_compare.py
```

---

## 7. Key levers (env knobs)

**Sim** (`sim/run_sim.sh` or `sim/run_sim_tiled.sh`):
| Var | Default | Effect |
|---|---|---|
| `SCAN_STEP` | 0.1 | probe-position spacing (Å). Smaller = more overlap = better; 0.05→160k pos |
| `SLICE_THICKNESS` | 2 (use **0.5**) | abTEM multislice slice (Å); forward accuracy, not depth resolution. Use 0.5 with phonons (high-angle TDS) |
| `PHONONS` | 0 | frozen-phonon configs (TDS). 0=coherent; 8–16 realistic. Memory-bounded (sequential) so any N is safe; cost ∝ N |
| `PHONON_SIGMA` | 0.08 | rms thermal displacement (Å); per-element would be more correct for O |
| `TILES` | 16 | (tiled only) scan bands = independent array jobs. More = smaller/faster/robuster each |
| `MAXPARALLEL` | 15 | (tiled only) concurrent tiles; ≤ the 15-GPU per-user cap |
| `WALLTIME` | 12 h | per job; each tile ~10 h at 0.05 Å/16-phonon |
| `SIM_TAG`, `OVERWRITE` | — | force a name / allow overwrite |

**Dose** (`sim/add_poisson_noise.py`): `--dose` (e/Å²; ~1e4–1e6 experimental, 1e8 near-noiseless), `--seed`. Streamed (low memory). Or `run_dose_series.sh` (`DOSES`, `NL`).

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
- **`$HOME` is tiny (2 GB) and fills silently** → jobs die with `Disk quota exceeded`.
  Culprits: cupy's CUDA-kernel cache (`~/.cupy`) and MATLAB's **ServiceHost**
  (`~/.MathWorks/ServiceHost`, ~2 GB). Fixes are permanent: `run_sim.slurm` points
  `CUPY_CACHE_DIR`/`XDG_CACHE_HOME`/`MPLCONFIGDIR` at `$SHARE`; `~/.MathWorks` is a
  symlink to `$SHARE`. Keep **everything** off `$HOME`.
- **16-phonon OOM** (Killed / exit 9) → don't wrap `FrozenPhonons` in the `Potential`
  (holds all configs at once, >128 GB). Loop configs one at a time (already done).
- **cupy kernel cache is per-array-size** → a *new* scan size recompiles; if the cache
  can't write (quota) the whole job dies. See the `$HOME` fix above.
- **git pull asks for a password / fails** → GitHub HTTPS needs a **PAT** (not your
  password) or an SSH remote; a cached token had expired.
- **Tiling introduces NO edge effects** (verified: per-position corr 1.0 across tile
  boundaries) — each scan position is an independent multislice through the *full*
  potential; tiling splits only *which positions* are computed.

---

## 10. Recommended next steps
1. **Reconstruct the 0.05 Å / 16-phonon dataset** (currently simulating). The blocker:
   160k positions is a **monolithic** recon (the object is solved jointly — the recon does
   NOT tile) and won't converge in the 2-day wall. Plan: add **position sub-sampling** to
   the recon driver (reconstruct every k-th position → k=2 gives a tractable 0.1 Å-
   equivalent for the headline figure), and keep the full 0.05 Å recon for when wall-time
   allows (it checkpoints `Niter*.mat`, so partial is usable).
2. **Dose series** (`run_dose_series.sh`) on a tractable dataset → the recon-vs-dose figure.
3. **Per-element phonon σ** (Pb/Ti ≈ 0.08, O ≈ 0.11) — more correct and more O signal.
4. **Exit vacuum padding** to kill the last-layer artifact at the source.
5. Deliberately not done: scan **jitter** and **partial coherence / source size** — future
   reviewer-proofing if asked.
6. **Atom finding / oxygen** — handed to a separate agent (§11).

---

## 11. For the atom-finding agent

**Your goal:** locate atoms in the reconstructed 3-D object — in particular **oxygen**,
which is *not* cleanly resolved (light, Z=8, buried under the Pb/Ti blur tails). The
recommended route is **PSF measurement + deconvolution / model-fitting** (§11.5). You
work entirely on the **Mac, locally** — no HPC needed; the data is already pulled. The
sims/reconstructions are owned by a separate thread; you consume its output.

### 11.1 The data (local, on the Mac)
| File | What |
|---|---|
| `~/Desktop/NL70_new_vol.npy` | **Primary.** Reconstructed object, complex64, shape **(70, 404, 404)** = [depth layer z, y, x]. Reg-off, single-mode, from the 0.15 Å coherent baseline. |
| `~/Desktop/NL42_new_vol.npy` | 42-layer version (dz = 1.665 Å). |
| `~/Desktop/recon_new/NL70/.../Niter120.mat` | The raw PtychoShelves output the npy came from (loader: `analysis/depth_resolution.py::load_recon`; object is in `outputs/object_roi`, a per-layer cell array). |
| `~/Desktop/column_cross_section_overlay.png` | The figure to reproduce/extend: recon cross-section down a column with GT Pb/Ti/O overlaid. |
| `~/Desktop/model_overlay.png` | GT projected potential with atom-type labels (which column is which). |

A **better dataset is coming** (0.05 Å step, 16 phonons/TDS, Poisson dose — higher SNR &
overlap, ideal for deconvolution). Same format, bigger; write your code against
`NL70_new_vol.npy` now and it will drop straight onto the new volume.

### 11.2 How to load + interpret the object
```python
import numpy as np
vol = np.load("~/Desktop/NL70_new_vol.npy")        # complex64 (70, 404, 404)
phase = np.angle(vol).astype(float)                 # the structure ∝ projected potential/layer
phase -= np.median(phase, axis=(1, 2), keepdims=True)   # remove per-layer offset
# axis 0 = depth (entrance→exit), axes 1,2 = (y, x)
```
- **Object pixel** `dx = 0.0492 Å` (in-plane, both axes).
- **Depth spacing** `dz = 1.0 Å` (NL70) / `1.665 Å` (NL42); total depth ≈ 70 Å.
- Atoms appear as **bright phase blobs**; brightness ∝ Z, so Pb ≫ Ti ≫ O.

### 11.3 Resolution / what is and isn't resolved (your target)
- **In-plane:** resolved to sub-Å (0.85 correlation to the model). Pb/Ti columns are clear.
- **Depth (along beam):** information-limited to **~2 Å** (λ/NA² optical axial limit).
  Atomic *planes* along the beam are resolved at **~3.9 Å** (NL70; verified by a kz peak).
- **Depth-dependent displacement** is recovered (~0.25 Å column lean; the labyrinth/vortex).
- **Oxygen: not resolved** — there is faint signal where O should be, but it does not
  separate from the Ti/Pb blur. **This is the job.**

### 11.4 Ground truth + recon↔model coordinate map (where atoms SHOULD be)
Structure file: `sim/PTO6_STO6_18_18_labyrinthPoscar.vasp` (PbTiO₃/SrTiO₃ labyrinth).
To get atom positions in the **same frame as the reconstruction**, apply the exact
transform the sim used, then map to recon pixels:
```python
import ase.io, abtem, numpy as np
a = ase.io.read("sim/PTO6_STO6_18_18_labyrinthPoscar.vasp")
a.rotate(-90, "y", rotate_cell=True); a = abtem.orthogonalize_cell(a)
s = max(a.cell.lengths()[:2]); a.cell[0,0]=s; a.cell[1,1]=s
a.center(axis=0); a.center(axis=1); a.center(axis=2, vacuum=2.0)
pos, Z = a.get_positions(), a.get_atomic_numbers()     # Pb=82, Ti=22, O=8, Sr=38
```
**Map (locked by NCC, corr 0.79; the recon is the transpose of the GT grid):**
- recon **(row r, col c)** → GT physical **X = 30 + c·dx**, **Y = 10 + r·dx**
- inverse (GT→recon pixel): **col = (X−30)/dx + CAL**, **row = (Y−10)/dx**, `CAL ≈ +1.8 px`
  (a measured sub-pixel registration; `analysis/column_cross_section_overlay.py` derives it
  data-drivenly — reuse that, don't hand-tune).
- **depth:** recon layer `i` ↔ GT z ≈ `(i+0.5)·dz` (entrance at layer 0; fitted offset ≈ 0).
- The scan window is X∈[30,50], Y∈[10,30] Å; only atoms in that box are in the field.

`analysis/column_cross_section_overlay.py` already implements the **full** alignment
(transpose + NCC shift + depth registration + sub-pixel CAL) and overlays Pb/Ti/O markers
on a column cross-section — **start from that script**; it's the tested path.

### 11.5 Recommended approach — PSF deconvolution / model fitting
Every atom reconstructs as the system's **point-spread function**: a 3-D blob, **tight
in-plane (~0.5–1 Å)** but **stretched along z (~2 Å+, the missing cone)**. So
`reconstruction ≈ true_atoms ⊛ PSF`.

1. **Measure the PSF empirically (the unfair advantage of a sim pipeline):** run
   `sim/simulate_4dstem.py` on a *single isolated atom* (a 1-atom structure in the same
   box) and reconstruct it through the identical pipeline → the reconstructed blob **is**
   your 3-D PSF (anisotropy, missing-cone and all). Defensible in a paper: measured, not
   assumed. (This needs one small sim+recon — coordinate with the sims thread.)
2. **Deconvolve** with **Richardson–Lucy** (non-negative, Poisson-appropriate — the right
   choice for electron counts; `skimage.restoration.richardson_lucy`), OR — better for
   weak/overlapping O — **model-based fitting**: place a PSF at each candidate lattice site
   (you know the lattice and the PSF) and solve for per-site amplitudes.
3. **Caveats:** deconvolution amplifies noise — use the **highest-dose** recon and cap RL
   iterations / regularise. The 3-D PSF is **anisotropic** (don't use a symmetric one — the
   empirical PSF handles this for free). O is a *contrast/SNR* problem, not a *resolution*
   one, so it works best where O is "resolved-but-blurred," not noise-buried.

### 11.6 Physics parameters (for computing/checking the ideal PSF)
300 keV (λ = 0.01969 Å) · convergence 100 mrad · defocus −20 Å (= 20 Å overfocus) ·
detector 200 mrad · 4×4 binned → Ndet 356, d_alpha 1.125 mrad/px, BF-disk radius 88.9 px ·
object pixel 0.0492 Å. (Full table + derivations in §3.)

### 11.7 Scripts to build on (`analysis/`, run locally)
- `column_cross_section_overlay.py` — recon vs GT overlay down a column (the alignment + figure).
- `column_cross_section.py` — depth cross-section down a column.
- `depth_resolution.py` — the recon-mat loader + the kz plane-frequency test.
- `dose_compare.py` — recon vs dose (useful once you pick a working dose for deconvolution).
