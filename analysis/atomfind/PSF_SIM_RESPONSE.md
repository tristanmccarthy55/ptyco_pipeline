# PSFs delivered + data pointers — response to `PSF_SIM_REQUEST.md`

The empirical PSFs are on `~/Desktop`. **Read §1 first — not every kernel is usable, and the
recommended path is deliberately narrow.** The sims thread verified each kernel against your
own `empirical_psf()` (which does `np.angle → argmax → crop`), so a kernel whose brightest
voxel is off-centre or on a surface artifact will silently give you a garbage PSF.

---

## 1. The recommended run (highest chance of success — do this first)

Your **current, best-conditioned target is the NL70 noiseless volume**, and its matched,
rock-solid kernel is the **Pb NL70** PSF. Coherent + high-count + strongest scatterer = the
cleanest point response we have (tight in-plane, textbook missing-cone z-elongation, 0.01 %
of the volume above |phase| 0.5).

```python
from config import preset
from dataclasses import replace
cfg = replace(preset("NL70_coherent"),
              single_atom_vol="~/Desktop/psf_Pb_NL70_vol.npy",   # <-- the gold kernel
              single_atom_species=82)                            # Pb
```
```
~/hyperspy-bundle/bin/python run_atomfind.py --preset NL70_coherent
```

**Use the Pb kernel for every species, not just Pb.** The reconstruction blur is a property of
the *imaging system* (probe × missing cone × recon), which is essentially element-independent
in shape — the element mostly sets amplitude. The per-element O/Ti kernels are noisier (O) or
broken (see §2), so the clean Pb PSF is both the safest and the most defensible system PSF.
Your `config` already defaults `single_atom_species=82`, so this is the natural choice.

---

## 2. PSF inventory (all on `~/Desktop`, `psf_<tag>_vol.npy`, complex64, same object format)

Verified with the exact `empirical_psf()` logic (argmax of `np.angle − per-layer median`):

| kernel | geometry | argmax centred? | verdict |
|---|---|---|---|
| **`psf_Pb_NL70_vol.npy`** | NL70, coherent | ✅ (24,30,30) | **GOLD — use this** |
| `psf_Pb_vol.npy` | NL70, coherent | ✅ | duplicate of the above (earlier extract) |
| `psf_Ti_NL70_vol.npy` | NL70, coherent | ✅ (24,30,30) | clean — usable if you want an element-specific Ti kernel |
| `psf_O_NL70_vol.npy` | NL70, coherent | ❌ atom sub-noise (SNR 0.39) | **not a usable kernel** — isolated O is ~1 % of Pb and *below* the recon noise floor here (peak 0.007 vs 99%-noise 0.017); argmax grabs a noise spike 1.15 Å off. Keep only as a detectability number (see note). Detect O with the Pb shape. |
| `psf_Pb_rev2_d1e10_vol.npy` | NL105, 16-ph, dose 1e10 | ✅ (43,30,30) | clean — for the reviewer-2 data (§4) |
| `psf_Pb_rev2_d1e8_vol.npy` | NL105, 16-ph, dose 1e8 | ✅ (43,30,30) | clean — for the reviewer-2 data (§4) |
| `psf_O_rev2_d1e10/d1e8` | NL105, 16-ph, dosed | ❌ corner artifact | **broken — do not use.** O does not localise as an isolated atom under TDS+dose. |
| `psf_Ti_rev2_d1e10/d1e8` | NL105, 16-ph, dosed | ❌ corner artifact | **broken — do not use.** Same reason. |
| `psf_Pb_NL70_z10_vol.npy` | NL70 depth series (entrance) | ✅ (z=6 → atom z=10 Å) | **FIXED** (re-extracted `--zdrop 4`) — clean, your `hz` crop fits. |
| `psf_Pb_NL70_z64_vol.npy` | NL70 depth series (exit) | ✅ (z=60 → atom z=64 Å) | **FIXED** (`--zdrop 4`) — but **edge-limited**: the atom is only 6 Å from the exit surface, so your `hz=4` crop clips ~3 layers. That *is* the exit band — which is exactly why recall is weak there, not an extraction fault. |

Every kernel also has a `psf_<tag>_check.png` next to it (in-plane max-projection + axial
slice) — glance at these; the good ones look like a single compact blob, the broken ones are
obviously noise/edge junk.

**Physics takeaway (real, not a bug):** only the **heavy Pb** column reconstructs as an
*isolated* atom. Even coherent/noiseless, isolated O measures at **SNR ≈ 0.39** (peak ~1 % of
Pb, ≈ the median noise) — below the reconstruction noise floor; under 16-phonon TDS + dose,
O and Ti (Debye–Waller-attenuated, buried in the incoherent background) go the same way. This
is not an extraction bug: **isolating the atom removes the crystalline support that makes O
visible in the first place.** O *is* present in the real crystal recon (`NL70_new_vol`, between
the Ti columns) because there it's held up by its neighbours + the periodic lattice. So detect
O/Ti with the **Pb system-PSF shape** (element-independent; it's the shape that matters for
deconvolution/matched-filtering) plus a low expected amplitude — exactly what your model-based
O detector (windowed NNLS) is built for. Don't feed an O/Ti kernel to `empirical_psf()`.

---

## 3. Datasets on `~/Desktop`

| volume | what it is | use |
|---|---|---|
| **`NL70_new_vol.npy`** (70,404,404) | **primary target.** Coherent, noiseless, 0.15 Å scan, **dz = 0.999 Å**, object pixel 0.0492 Å, 70 depth layers. Best depth resolution we have. | `--preset NL70_coherent` |
| `NL42_new_vol.npy` | same data, 42 layers (**dz = 1.665 Å**) — coarser depth. | `--preset NL42_coherent` |
| `NL42_vol.npy`, `NL21_vol.npy` | older/coarser layerings — ignore unless you want a depth-sampling sweep. | — |
| `dose_series/dose{1e10,1e8,1e6,1e4}/Niter200.mat` | a **separate** experiment: 0.1 Å / NL105, coherent, 4 doses. HDF5-v7.3 `.mat` (read with h5py; `outputs/object_roi`). Good extra targets to test detection-vs-dose. | optional |

The reviewer-2 production volume (0.1 Å-combined, 16-phonon, mixed-state-probe recon) is **not
on the Desktop yet** — its PSFs (§2, the `rev2` Pb ones) are staged ahead of it.

---

## 4. Reviewer-2 preset (fix these when that volume lands)

The `reviewer2` preset in `config.py` is a template with placeholder values. When the 16-phonon
production recon arrives, set:

```python
"reviewer2": Config(name="reviewer2",
                    recon_vol="~/Desktop/<the reviewer2 vol>.npy",
                    dz=0.666,                 # NL105, NOT the template's 1.0
                    dose_e_per_A2=1e8,        # match the recon you point at (1e8 or 1e10)
                    single_atom_vol="~/Desktop/psf_Pb_rev2_d1e8_vol.npy"),  # dose-matched Pb
```
Notes: dose-**match** the kernel to the volume (`d1e8` kernel ↔ 1e8 recon, `d1e10` ↔ 1e10).
It's 0.1 Å step (the production data was binned to 0.1 Å from the 0.05 Å sim), and the recon
used released mixed-state probes — but you still consume the object volume exactly as before.

---

## 5. Experimental / calibration details (already baked into `config.py` — confirming)

- **Beam:** 300 keV (λ = 0.01969 Å), 100 mrad convergence, 20 Å overfocus (defocus −20 Å).
- **Detector:** 200 mrad, 4× binned → 356×356 patterns.
- **Structure (ground truth):** `sim/PTO6_STO6_18_18_labyrinthPoscar.vasp` — `config.vasp()`
  resolves it. PbTiO₃/SrTiO₃ labyrinth, ~70 Å along the beam.
- **Recon ↔ GT map (validated in `column_cross_section_overlay.py`, already in `config`):**
  `recon(row r, col c) → GT X = 30 + c·dx, Y = 10 + r·dx`, `dx ≈ 0.0492 Å`, `z_recon ≈ z_GT`;
  a sub-pixel `CAL_X ≈ +0.087 Å` and a data-driven depth registration (`depth_branches`) refine it.
- **Object convention:** take `np.angle` of the complex volume and subtract the per-layer median
  (atoms are positive-phase peaks) — this is what `empirical_psf()` and the overlay both assume.
- **Interpreter:** `~/hyperspy-bundle/bin/python` (abtem + skimage + scipy), per your config header.

---

## 6. How every volume was generated (full pipeline — so you handle each one correctly)

Both the target data AND the PSF grids come from the **same two stages** — an abTEM multislice
4D-STEM simulation, then a PtychoShelves multislice reconstruction — differing ONLY in (a) the
object placed in the box and (b) the noise level. Knowing this fixes how to read the volume.

### 6.1 Stage 1 — abTEM multislice 4D-STEM (`sim/simulate_4dstem.py`, the forward model)
Common to every dataset:
- **Beam** 300 keV; **probe** 100 mrad convergence, 20 Å overfocus (defocus −20 Å).
- **Potential**: abTEM `Potential`, Lobato parametrization, infinite projection, on the FULL
  square box (no cropping — a cropped box aliases the broadened exit wave). Slice thickness
  **2 Å** (NL70/NL42 data) or **0.5 Å** (reviewer-2).
- **Scan**: `GridScan` over a window centred at prepared-cell **(40, 20) Å**, positions
  flattened **y-fastest**. Step **0.15 Å** (NL70), **0.1 Å** (dose series + reviewer-2).
- **Detector**: full ~200 mrad pixelated → binned **4×4 → 356×356** patterns, DC-centred.
- **Frozen phonons (reviewer-2 only)**: N displaced configs, diffraction **intensities**
  averaged (incoherent TDS). Per-species room-T displacement from tabulated isotropic B
  (Pb 0.90, Sr 0.55, Ti 0.45, O 0.80 Å²) → σ = √(3B/8π²); Pb & O vibrate ~2× Ti.
- Patterns flux-normalised then scaled to a fixed electron count (noiseless base).

### 6.2 The object in the box — the ONLY thing that differs
- **Target data** (`NL70_new_vol`, `NL42`, dose series): the real PbTiO₃/SrTiO₃ labyrinth from
  `PTO6_STO6_18_18_labyrinthPoscar.vasp` — read → **rotate −90° about y** (beam→z) →
  `orthogonalize_cell` → pad to a **square 70.008 Å** in-plane box → `center` + 2 Å vacuum.
  ~70 Å of material along the beam.
- **PSF grids** (`psf_*` volumes): a sparse **2-D grid of ONE element at ONE depth** (z = 37 Å)
  via `build_atom_grid` — atoms on a `spacing`-Å grid (**4 Å** NL70, **2.5 Å** rev2) across the
  scan window, placed directly in the identical 70.008 × 70.008 × 74 Å box. A grid, not a lone
  atom, because a lone atom is too sparse for the reg-off recon (it fills the volume with
  noise); the grid gives enough in-plane density to CONVERGE while each atom stays isolated
  in-plane (spacing ≫ the ~1 Å blob) and axially (single plane). **Everything else — optics,
  detector, slice thickness, phonons, dose, and the recon below — is byte-for-byte the
  production pipeline, so each reconstructed grid blob IS the system point-spread response.**
  (This is exactly why light/isolated O collapses — see §2 — while the recon of the *full
  crystal* still shows O: the grid removes the neighbour support.)

### 6.3 Stage 2 — dose (optional, `sim/add_poisson_noise.py`)
Poisson shot noise added post-sim: incident e/pattern = dose·step², streamed two-pass. The
**NL70 target and the NL70 PSF grids are NOISELESS** (no Poisson). The dose series and the
reviewer-2 PSFs are dosed at the stated e/Å².

### 6.4 Stage 3 — PtychoShelves multislice reconstruction
- Engine: Yu Lei MultiHollowPtycho, GPU **LSQ-ML/MLs**, hollow angle 0; two-engine presolve
  (Ndp 178→356) then full-resolution.
- `custom_data_flip = [0,0,1]` → the recon is **TRANSPOSED** vs the model (row↔Y, col↔X — the
  §5 coordinate map already accounts for this).
- **NL** depth layers over ~70 Å: NL70→dz 0.999, NL42→1.665, NL105→0.666 Å.
- Probe: NL70/NL42 target + their PSFs = **1 fixed mode, layer-reg OFF**; reviewer-2 PSFs =
  **released mixed-state probes** (so the incoherent TDS has somewhere to go instead of
  corrupting the object).
- Output = the reconstructed multislice **object**: a complex transmission function per depth
  layer, `(nL, Ny, Nx)` complex64 — the array you load.

### 6.5 What this means when you handle a volume
- It is a **multislice object, not a projected image**: `V[l]` is the complex transmission of
  depth slice `l` (entrance→exit). `np.angle(V[l])` ≈ that slice's projected potential; atoms
  are **positive-phase peaks after per-layer median subtraction**.
- **z** = beam axis at `dz` Å/layer; **in-plane** is transposed vs GT (use the §5 map).
- **Surface layers are artifact-prone**: entrance/exit planes accumulate a "dumping-ground"
  residual (hence `zmax_show_A=66` trims the exit, and why the z=10/64 depth PSFs broke).
  Trust the interior.
- **Noise**: NL70 = noiseless/coherent (cleanest, your best bet); dosed volumes carry Poisson
  noise → cap deconvolution iterations (your `effective_rl_iters()` already scales for this).

---

## 7. How the kernel is consumed (why the broken ones fail)

`psf.empirical_psf()` loads `cfg.single_atom_vol`, takes `np.angle`, subtracts the per-layer
median, finds the **global argmax**, and crops `±psf_half_z_A/dz` × `±psf_half_xy_A/dx` around
it. So the kernel *must* have its brightest voxel on the atom. The GOLD/clean kernels do; the
broken ones have argmax on an edge/surface artifact, which is why they're unusable as-is. If you
ever re-extract, `analysis/atomfind/extract_psf.py <recon_dir> <tag>` regenerates one (it
force-orients the atom to a positive peak and prints the same argmax/cleanliness check).

**Bottom line:** start with `NL70_coherent` + `psf_Pb_NL70_vol.npy`. That pair is the one we're
confident in.
