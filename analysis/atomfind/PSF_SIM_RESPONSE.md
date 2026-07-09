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
| `psf_O_NL70_vol.npy` | NL70, coherent | ⚠️ off-centre (24,47,46) | **marginal** — O is weakly resolved; blob sits in a noisy field. Prefer the Pb kernel as proxy. |
| `psf_Pb_rev2_d1e10_vol.npy` | NL105, 16-ph, dose 1e10 | ✅ (43,30,30) | clean — for the reviewer-2 data (§4) |
| `psf_Pb_rev2_d1e8_vol.npy` | NL105, 16-ph, dose 1e8 | ✅ (43,30,30) | clean — for the reviewer-2 data (§4) |
| `psf_O_rev2_d1e10/d1e8` | NL105, 16-ph, dosed | ❌ corner artifact | **broken — do not use.** O does not localise as an isolated atom under TDS+dose. |
| `psf_Ti_rev2_d1e10/d1e8` | NL105, 16-ph, dosed | ❌ corner artifact | **broken — do not use.** Same reason. |
| `psf_Pb_NL70_z10`, `psf_Pb_NL70_z64` | NL70 depth series | ❌ argmax at surface (layer 0) | **broken** — the z=10/64 atoms sit inside the surface-trim zone; re-extract with a smaller `--zdrop` if you want the depth dependence. |

Every kernel also has a `psf_<tag>_check.png` next to it (in-plane max-projection + axial
slice) — glance at these; the good ones look like a single compact blob, the broken ones are
obviously noise/edge junk.

**Physics takeaway (real, not a bug):** under 16-phonon thermal diffuse scattering + finite
dose, only the **heavy Pb** column reconstructs as an isolated atom. O (light) and Ti (medium)
are attenuated by the Debye–Waller factor and buried in the incoherent TDS background, so an
*isolated* O/Ti atom is too under-constrained to reconstruct there. Detect O/Ti with the Pb
system-PSF (the shape is what matters for deconvolution/matched-filtering).

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

## 6. How the kernel is consumed (why the broken ones fail)

`psf.empirical_psf()` loads `cfg.single_atom_vol`, takes `np.angle`, subtracts the per-layer
median, finds the **global argmax**, and crops `±psf_half_z_A/dz` × `±psf_half_xy_A/dx` around
it. So the kernel *must* have its brightest voxel on the atom. The GOLD/clean kernels do; the
broken ones have argmax on an edge/surface artifact, which is why they're unusable as-is. If you
ever re-extract, `analysis/atomfind/extract_psf.py <recon_dir> <tag>` regenerates one (it
force-orients the atom to a positive peak and prints the same argmax/cleanliness check).

**Bottom line:** start with `NL70_coherent` + `psf_Pb_NL70_vol.npy`. That pair is the one we're
confident in.
