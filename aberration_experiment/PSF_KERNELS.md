# Matched PSF kernels for the thin aberration campaign

What the per-α single-atom kernels are, how they broke at high α, and what has to be true for
them to come out clean. Written 2026-09-10 after the first full atomfind pass on
`atomfind_results_20260910_0944`.

## Why they exist

atomfind models the phase volume as a superposition of single-atom responses, so it needs the
**system PSF**: how ONE atom appears in *this* reconstruction, with *these* aberrations. An
unmatched (synthetic) kernel is not a small error — it costs most of the recall. Every α has a
different probe, so every α needs its own kernel.

## How they are made

`campaign/run_thin_atomfind.sh` runs, per α, a **Pb** and a **Ti** grid leg: `build_atom_grid()`
puts one element on a 2-D grid at a **single depth plane**, through the same aberrated probe and
the same full box as the lab recon. `analysis/atomfind/extract_psf.py` then crops one clean
interior blob → `psf_<el>_a<α>_vol.npy`, which atomfind loads via `--single-atom-vol` /
`--ti-kernel-vol`.

## What went wrong

The grid recons degrade with α. `extract_psf.py` prints a cleanliness diagnostic
(`|phase|>0.5` fraction — **small = clean single blob**):

| α | Pb | Ti | argmax z | verdict |
|---|----|----|----------|---------|
| 50 | 0.10% | 0.07% | mid-crop | clean |
| 70 | 1.67% | 10.4% | mid-crop | marginal (Ti) |
| 90 | 67.5% | 61.6% | **last layer** | junk |
| 100 | 81.3% | (NaN'd) | **last layer** | junk |

At a90/a100 the extractor's argmax lands on the **last z layer** instead of the box centre where
the grid atom actually is, and the rendered volumes are speckle across the illuminated disc
rather than isolated blobs.

## Diagnosis: conditioning, not geometry

The obvious suspect is **blob overlap**: the balanced probe grows with α (d90 = 4.0 / 3.9 / 6.6 /
11.0 Å at 50/70/90/100) while `GRID_SPACING` is fixed at 4 Å, so at a100 an 11 Å probe spans
~3 grid atoms. `build_atom_grid`'s own docstring assumes "spacing >> the ~1 Å blob".

**That suspect is wrong.** Measured on the clean a50 kernel: in-plane **FWHM 0.20 Å** against a
**4.01 Å** probe — ptychography resolves ~20× below the probe size, which is the entire point of
the technique. Isolated atoms stay isolated at 4 Å spacing even at a100, so widening the grid
would not have fixed anything.

The real cause is an **under-constrained depth solve**:

- the grid object is a **single atomic plane**;
- it is reconstructed into **NL = 23 (a90) / 28 (a100) free depth layers** — NL is chosen as
  Nyquist for the *lab* object, which is far more depth structure than a one-plane object has;
- with **`regularize_layers = 0`**.

`ptycho/run_synthetic_recon_ML.m` sets `REGLAYER=0` deliberately, and its comment says exactly
why: layer regularization "was only stabilising the under-constrained deep solve, a job the fixed
probe + high overlap should now do." That is true for the dense lab slab and **false for a sparse
single-plane grid**, which *is* an under-constrained deep solve. With nothing tying the layers
together, the solver is free to smear one plane's phase across 28 slots — the speckle we see.

Secondary factor: `SCAN_WINDOW = 14 Å` is barely wider than the a100 probe (11 Å), so there is
little positional diversity to condition the sparse solve.

## Fix, cheapest first

1. **Re-recon only — no re-sim.** The grid sim data already exists; only the reconstruction was
   wrong. Re-run the a90/a100 Pb+Ti grid recons with `REGLAYER=0.5` (plus `BETA_LSQ=0.05`, which
   the a100 Ti leg needs anyway after its NaN divergence).
2. If still speckled: **re-sim with `SCAN_WINDOW=30`** (≳3× the probe) for diversity. Keep
   `GRID_SPACING=4` — spacing was never the problem.
3. Optional: **thin the PSF box** (`--grid-box-z ≈ 8`) so a one-plane object is solved in ~8
   layers at the *same* dz rather than 28. Better conditioned; slightly different axial tail
   because there is less free-space propagation either side.

## How to tell a kernel is good

`extract_psf.py` prints what you need:

- `|phase|>0.5` fraction of a few % at most;
- argmax at the **centre** of the crop in z (the grid atom sits at box centre), not at a z edge;
- in-plane FWHM sub-Å; axial FWHM comparable to the depth resolution δz = λ/α²;
- **`--zdrop` must be the vacuum band in LAYERS**, `round(4.0/dz)` → 1/2/3/4 at a50/70/90/100.
  The default of 12 silently erases a NL=7 volume (`z0=12, z1=-5`) and reports "no blobs found".

## Current status — FLAG

- **a50, a70**: matched grid kernels used (a70's Ti is marginal, 10.4%).
- **a90, a100**: grid kernels unusable, so the reported numbers were produced with the
  **data-derived PSF** (`psf.data_psf()` — real Pb columns averaged from the same lab recon).
  It is matched by construction (same reconstruction, same aberrations) and it performs well,
  but it is measured from the stacked labyrinth, so its axial tail carries neighbouring-column
  contributions that a true isolated-atom kernel would not. **Any a90/a100 kernel-dependent
  number should be labelled data-derived until the `REGLAYER` re-runs land.**
