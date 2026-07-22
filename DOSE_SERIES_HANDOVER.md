# Dose-series handover — for the figures agent

**Purpose.** You are making publication figures from a **4-point electron-dose series**: the
same simulated 4D-STEM dataset reconstructed at four doses (10¹⁰, 10⁸, 10⁶, 10⁴ e/Å²) to
show how ptychographic fidelity — especially **depth resolution** — degrades as dose drops.
Everything below is ground truth read back from the reconstructed `.mat` files and the
sim/recon source. You do **not** need HPC access; all data is on the local Mac.

---

## 1. Where the data is (local Mac)

```
~/Desktop/dose_series/
  dose1e10/   dose1e8/   dose1e6/   dose1e4/     # one folder per dose (high -> low)
```

Each dose folder contains:

| Item | What it is |
|---|---|
| `Niter200.mat` | **Final reconstruction** (Niter 200). The 3-D object is under `outputs/object_roi`. Use this for the hero figures. (~1.1 GB — it also carries the exit wave.) |
| `Niter25.mat … Niter175.mat` | Intermediate checkpoints every 25 iterations (convergence panels, if wanted). |
| `O_phase_roi/O_phase_roi_Niter100_Layer{1..105}.tiff` | Per-depth-layer **object-phase images** already rendered to TIFF (Niter 100). Handy for quick per-slice figures without touching the `.mat`. |
| `probe_mag/` | Reconstructed probe magnitude images. |
| `exit_wave_200.mat` | Exit wave at Niter 200. |
| `run_synthetic_recon_ML.m` | **A copy of the exact recon script that produced this folder** — the definitive recon provenance. |

Already-made summary figures (regenerate from `analysis/figures/fig1_depth_resolution.py` + `fig4_species_columns.py`):
- `~/Desktop/dose_compare.png` — 2×4 grid: depth cross-section (top) + in-plane slice (bottom) per dose.
- `~/Desktop/dose_metrics.png` — quantitative degradation vs dose.

---

## 2. How the data was collected (full provenance)

Three stages: **abTEM sim → Poisson dose → PtychoShelves recon.** Same sim for all four
doses; only the noise level and the reconstruction differ.

### (a) Simulation — `sim/simulate_4dstem.py`
Multislice 4D-STEM of a **PbTiO₃ / SrTiO₃ labyrinth** (the ground-truth model).

| Parameter | Value |
|---|---|
| Beam energy | **300 keV** (λ = 1.969 pm) |
| Convergence semi-angle | **100 mrad** |
| Defocus | overfocus 20 Å (abTEM `defocus = −20 Å`) |
| Detector | 200 mrad outer, **4× binned → 356×356 px** diffraction patterns |
| Multislice slice thickness | 2.0 Å |
| Frozen phonons | none (coherent baseline) |
| Structure | `sim/PTO6_STO6_18_18_labyrinthPoscar.vasp`, cell 70.008×70.008×~48 Å, rotated so the **beam axis ≈ 70 Å** |
| Scan | centre **(40, 20) Å**, window **20 Å**, step **0.1 Å** → **200×200 = 40 000 positions** |

### (b) Dose — `sim/add_poisson_noise.py`
Post-sim Poisson noise applied to the noiseless patterns. Incident electrons per pattern =
`dose × step²` with `step² = 0.01 Å²`, i.e. **dose 10¹⁰ → 10⁸ e/pattern**, 10⁸ → 10⁶, etc.
The four doses reconstructed are **1e10, 1e8, 1e6, 1e4 e/Å²**.

### (c) Reconstruction — `ptycho/run_synthetic_recon_ML.m` (copy saved in each dose folder)
PtychoShelves, **Yu Lei's `GPU_MS` MultiHollowPtycho engine (LSQ-ML/MLs, hollow-angle = 0)**,
`custom_data_flip = [0,0,1]` (a transpose — see the coordinate map below).

| Parameter | Value |
|---|---|
| Depth layers | **105** layers, **dz = 0.666 Å**, total ≈ 70 Å |
| Object ROI | **405 × 405 px**, pixel **0.04916 Å/px** (≈ 20 Å field) |
| Diffraction size | 356 × 356 |
| Probe modes | **1** (fixed) |
| Layer regularisation | **off** (`REGLAYER = 0`) — this is the knob that buys the depth resolution |
| `beta_LSQ` | 0.1 |
| Iterations | **200** (run as a restart chain; the pulled result is the full-resolution engine only) |
| Presolve | two-engine schedule (grouping [64,32], Ndp 178→356) on the fresh run; collapsed on restart |

---

## 3. Geometry cheat-sheet (numbers you'll need in captions)

```
dose values      : 1e10, 1e8, 1e6, 1e4  e/Å²
object pixel dx  : 0.04916 Å            (in-plane, both axes)
layer spacing dz : 0.666 Å              (105 layers over ~70 Å)   [dose series]
scan step        : 0.1 Å                (200×200 = 40 000 positions)
scan window      : 20 Å, centred at GT (X,Y) = (40, 20) Å
beam energy      : 300 keV
DP size          : 356 × 356
```

> Note the two z-samplings in this project: the **dose series is NL105 (dz = 0.666 Å)**.
> The separate clean NL70 run (`~/Desktop/NL70_new_vol.npy`) is dz ≈ 1.0 Å — don't mix them.

---

## 4. Recon ↔ ground-truth model map (for overlays)

The reconstruction is **transposed** relative to the model. In-plane map:

```
recon (row r, col c)  ->  GT  X = 30 + c·0.04916 Å ,  Y = 10 + r·0.04916 Å
depth  layer l        ->  GT  z ≈ l·dz            (z_recon ≈ z_GT)
```

(30, 10) is the scan-window origin = centre (40,20) − half-window (10,10). A sub-pixel
in-plane calibration `CAL_X ≈ +0.087 Å` and a depth registration were fit for the NL70
overlay; re-fit per dose since each is an independent reconstruction.

**Canonical overlay code (already written):** `analysis/column_cross_section_overlay.py` —
it loads the `.vasp`, reproduces the sim orientation, does NCC + depth registration, and
draws atom markers on the recon. Reuse its model-loading/orientation rather than
re-deriving the rotation. Ground-truth structure: `sim/PTO6_STO6_18_18_labyrinthPoscar.vasp`.

---

## 5. Relevant code (all in `analysis/`)

| File | What it does — reuse for |
|---|---|
| `figures/dose_fig_common.py` | **Start here.** The canonical `.mat` reader (`load_dose`, `.npy`-cached), reference-column picker, per-dose registration, and kz spectrum — shared by all the `fig*` scripts. |
| `column_cross_section.py` | Depth cross-section down a Pb and a Ti column + phase-vs-depth profile with plane-period ticks (written for NL70 `.npy`; adapt the loader for dose `.mat`). |
| `column_cross_section_overlay.py` | Overlays the GT atom positions on the recon (model registration). |
| `figures/fig1_depth_resolution.py` | Depth resolution vs dose: cross-sections + kz plane-frequency peak — the numbers behind the dose-degradation story. |

---

## 6. Ready-to-use snippets

**Load a dose volume → phase array `[nL, Ny, Nx]` + `dz`** (the canonical reader, `figures/dose_fig_common.load_dose`):

```python
import glob, os, numpy as np, h5py

def load_dose(dose, root="~/Desktop/dose_series"):
    m = sorted(glob.glob(os.path.join(os.path.expanduser(root), f"dose{dose}", "Niter*.mat")),
               key=lambda p: int("".join(filter(str.isdigit, os.path.basename(p)))))[-1]
    with h5py.File(m, "r") as f:
        g = f["outputs"]
        layers = []
        for r in g["object_roi"][:, 0]:            # object_roi is an array of HDF5 refs
            a = f[r][:]
            a = (a["real"] + 1j*a["imag"]) if a.dtype.names else a
            layers.append(a.T)                     # .T undoes MATLAB's column-major order
        V = np.angle(np.array(layers)).astype(float)
        V -= np.median(V, (1, 2), keepdims=True)   # per-layer background subtraction
        dz = float(np.median(g["z_distance"][:, 0][np.isfinite(g["z_distance"][:, 0])])) * 1e10
    return V, dz                                    # V: [105, 405, 405] phase ; dz≈0.666 Å

V, dz = load_dose("1e10")
DX = 0.04916                                        # in-plane Å/px
```

**Auto-pick a strong Pb column (interior), reuse it across doses:**

```python
from scipy.ndimage import maximum_filter
ref, _ = load_dose("1e10")                          # cleanest recon chooses the column
dm = ref.mean(0); dmn = dm - dm.min()
mask = np.zeros_like(dmn, bool); mask[40:-40, 40:-40] = True
pk = (dmn == maximum_filter(dmn, 25)) & (dmn > np.percentile(dmn[mask], 97)) & mask
ys, xs = np.where(pk); j = int(np.argmax(dmn[ys, xs])); yc, xc = ys[j], xs[j]
```

**Depth cross-section down that column (z vs in-plane):**

```python
W = 22                                              # half-strip px (~1.1 Å)
cs = V[:, yc-1:yc+2, xc-W:xc+W].mean(1)             # [nL, 2W]
extent = [-W*DX, W*DX, (V.shape[0]-0.5)*dz, 0.5*dz] # x (Å) ; z entrance->exit
# imshow(cs, extent=extent, aspect="auto", cmap="inferno",
#        vmin=np.percentile(cs,5), vmax=np.percentile(cs,99.3))
```

**A mid-depth in-plane slice:**

```python
sl = V[V.shape[0] // 2]                              # [405, 405], imshow with inferno
```

---

## 7. The story the figures should tell (from our analysis)

Depth resolution is the **first casualty of dose**:

- **1e10 vs 1e8** — essentially identical (in-plane correlation ≈ 0.998; depth plane-period
  peak strong at both). 1e8 e/Å² is a realistic experimental dose → the method survives it.
- **1e6** — still usable in-plane (corr ≈ 0.85) but depth contrast roughly halves.
- **1e4** — broken: in-plane correlation collapses (≈ 0.31) and **depth resolution is gone**
  (the periodic plane signal down a column disappears).

So the headline figure is the **depth cross-section down a Pb column across the four doses**
(top row of `dose_compare.png`): the periodic bright plane spots persist at 1e10/1e8, fade
at 1e6, and vanish at 1e4.

---

## 8. Gotchas

- `.mat` files are **HDF5 v7.3** → read with `h5py`, not `scipy.io.loadmat`. `object_roi` is an
  array of object references; transpose each layer (`.T`) to undo MATLAB column-major order.
- Take **`np.angle`** of the complex object (it's a phase object) and **subtract the per-layer
  median** before displaying — raw phase has a per-layer offset.
- Each dose is an **independent reconstruction**: don't assume identical pixel registration.
  Auto-find the column on the reference dose and reuse those indices, or re-register per dose.
- Don't confuse the **dose series (NL105, dz 0.666 Å, `.mat`)** with the clean **NL70 run
  (dz ≈ 1.0 Å, `~/Desktop/NL70_new_vol.npy`)** — different depth sampling.
- `O_phase_roi/` TIFFs are at **Niter 100**; the `.mat` volumes go to **Niter 200**.
```
