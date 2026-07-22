# Data request — empirical single-atom PSF (for the atom-finding pipeline)

**Why:** every atom reconstructs as the system's 3-D point-spread function. If we
reconstruct a *single isolated atom* through the **identical** sim+recon pipeline, the
reconstructed blob **is** the PSF (anisotropy, missing cone, channeling and all) — a
*measured*, paper-defensible kernel, not an assumed one. The atom-finding pipeline
(`analysis/atomfind/`) then deconvolves / fits with it. It already runs with a data-derived
stand-in (averaged Pb blobs) whose axial tails are neighbour-contaminated and tapered; this
sim removes that limitation. Drop-in: `run_atomfind.py --single-atom-vol <the .npy>`.

This is small and cheap (one atom → fast sim, fast recon), unlike the labyrinth runs.

## What to simulate
Same geometry as the real sim so `Ndpx / d_alpha / dx` and the axial propagation match —
reuse `sim/simulate_4dstem.py`'s `build_phantom_atoms()` pattern but with **one atom**:

- **Structure:** a single atom at the **scan centre**, in the **same square in-plane box**
  as production: `cell = [70.008, 70.008, 74.0]` Å (2 Å vacuum + ~70 Å beam path), atom at
  `(x, y, z) = (40.0, 20.0, 37.0)` Å (mid-depth). PBC as in production.
  ```python
  from ase import Atoms
  atoms = Atoms("Pb", positions=[(40.0, 20.0, 37.0)], cell=[70.008, 70.008, 74.0], pbc=True)
  ```
- **Element (priority order):** ① **Pb** (Z=82) — bright reference, validates the
  data-derived PSF; ② **O** (Z=8) — the O-specific PSF (light atoms channel differently, so
  the O kernel may differ slightly and matters most for the O detector); ③ **Ti** (Z=22).
- **Coherent (no phonons)** — we want the instrument+reconstruction PSF (the mean-position
  kernel we deconvolve), not the atom's own thermal smear.
- **Everything else identical to the target dataset it will be used on:**
  - For the **current NL70** work: `SCAN_STEP=0.15`, `SLICE_THICKNESS=2`, scan centre
    (40, 20), 20 Å window, 300 keV, 100 mrad, defocus −20 Å, 4×4 binning.
  - For the **0.05 Å / 16-phonon reviewer-2** work: match *that* run — `SCAN_STEP=0.05`,
    `SLICE_THICKNESS=0.5`, and apply the **same dose** (`add_poisson_noise.py --dose <D>`)
    so the PSF carries the same noise/regularisation regime.
- **Reconstruct through the identical driver/params** as the matching production recon
  (reg off, `PROBE_MODES=1`, fixed probe `PROBE_START=inf`, `BETA_LSQ=0.1`, same `NL`
  = 70 for NL70 / matching for reviewer-2).

## What I need back
The reconstructed object exported to the **same format** as `NL70_new_vol.npy` — a
`complex64` array shape `(nL, Ny, Nx)` (use the same `analysis/figures/dose_fig_common.py::_read_mat`
→ `np.save` path that produced `NL70_new_vol.npy`). Name them:
- `psf_Pb_vol.npy` (priority 1), `psf_O_vol.npy` (priority 2), `psf_Ti_vol.npy` (priority 3),
- and, for the better data, `psf_Pb_step0.05_dose<D>_vol.npy` etc.

Put them on `~/Desktop/`. The pipeline just needs the path via `--single-atom-vol`.

## Nice-to-have (not blocking)
- **Depth dependence:** the same single Pb at `z = 10, 37, 64` Å (entrance / mid / exit) →
  three PSFs to check whether the axial kernel varies with depth (defocus/propagation). If
  it does, we can use a depth-varying PSF; if not, mid-depth alone is fine.
- A **phonon** single-atom (σ matching the run) to quantify the thermal contribution to the
  blur vs the instrument PSF.

## How the pipeline consumes it
`psf.empirical_psf()` loads the `.npy`, finds the atom's blob, and crops it onto the standard
kernel grid — it becomes the default PSF (over data/synthetic) automatically. No code change;
just the flag.

---

# REQUEST 2 (2026-07-16) — in-situ vacancy-difference kernels

**Status of request 1:** delivered (see `PSF_SIM_RESPONSE.md`) — the isolated Pb/Ti kernels are
in use. The isolated-**O** kernel failed for a physical reason the response itself identified:
*isolating the atom removes the crystalline support that makes O visible*. This request fixes
that by measuring each species' contribution **inside the lattice**.

## The idea
Reconstruct the full labyrinth (already exists: `NL70_new_vol.npy`) and the same labyrinth
**minus exactly one atom**, through the byte-identical pipeline. Then

```
PSF_insitu(X) = angle(recon_full) − angle(recon_minus_X)
```

is the reconstructed contribution of that ONE atom — channeling, neighbour support and all.
This is the correct matched filter for weak in-crystal species (O especially), which the
isolated-atom sim cannot give. (Recon is nonlinear so the difference isn't a strict linear
kernel, but it is the best empirical in-situ response and exactly what the guided O-detector
needs.)

## What to simulate (each = full labyrinth with ONE atom deleted; NL70 pipeline: 0.15 Å step,
coherent, slice 2 Å, same scan window, reg-off / 1 fixed probe mode / NL70 recon params)
| # | delete | where | priority |
|---|---|---|---|
| 1 | one **O** | on a B–O column, mid-depth (z ≈ 35 Å), mid-field | **1 (the O matched filter)** |
| 2 | one **O** | on a pure-O column, mid-depth | 2 |
| 3 | one **Ti** | mid-depth | 3 |
| 4 | one **Pb** | mid-depth | 4 (sanity: should ≈ the isolated-Pb kernel) |

Pick atoms well inside the scan window (X,Y ∈ [35,45] Å) so the difference blob is clean of
edge effects. Record WHICH atom was deleted (element, X, Y, Z in the prepared-cell frame) in
the filename or a sidecar note.

## Deliverables
- `vac_O_BO_vol.npy`, `vac_O_pure_vol.npy`, `vac_Ti_vol.npy`, `vac_Pb_vol.npy` on `~/Desktop`
  — full recon volumes, same complex64 (70,404,404) format (NOT pre-differenced; we difference
  locally against `NL70_new_vol.npy` so registration/scale is handled in one place).
- Plus (small, from request-1 leftovers): re-extract the broken `psf_Pb_NL70_z10/z64` depth
  kernels with a smaller `--zdrop` (per your response §2) → depth-dependent kernels for the
  entrance/exit bands where recall is weakest.

## Cost note
Each run = one full NL70-style sim+recon (same cost as the original NL70 run). Priorities 1–2
are the ones that matter for oxygen; 3–4 are cheap insurance/validation if the queue allows.
