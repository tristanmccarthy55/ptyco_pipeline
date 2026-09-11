# Matched PSF kernels for the thin aberration campaign

What the per-α single-atom kernels are, the rule they have to obey, how they broke at high α, and
which fixes are legitimate. Written 2026-09-10 after the first full atomfind pass on
`atomfind_results_20260910_0944`.

## The rule: byte-identical to the lab reconstruction

atomfind models the phase volume as a superposition of single-atom responses, so the kernel must
be **the system PSF of the very reconstruction it will be matched against**. That means the grid
leg and the lab leg must agree on *everything except the object*:

| must match | value |
|---|---|
| depth slices `NL` | **α-dependent** — Nyquist over the 27.525 Å box |
| slice spacing `dz` | `27.525 / NL` |
| box | 5 cells + 2×4 Å z-vacuum = 27.525 Å, full-box recon |
| probe | the true aberrated probe for that α, fixed |
| `regularize_layers` | **0** |
| `beta_LSQ`, `GROUPING`, `NITER`, scan step, scan window | same as the lab leg |

Only the **object** may differ (one element on a grid instead of the labyrinth). Anything else and
you are matching data against a kernel from a different experiment.

Per α:

| α (mrad) | NL | dz (Å) | δz = λ/α² (Å) |
|---|---|---|---|
| 50 | 7 | 3.932 | 7.88 |
| 70 | 14 | 1.966 | 4.02 |
| 90 | 23 | 1.197 | 2.43 |
| 100 | 28 | 0.983 | 1.97 |

## How they are made

`campaign/run_thin_atomfind.sh` runs, per α, a **Pb** and a **Ti** grid leg: `build_atom_grid()`
puts one element on a 2-D grid at a **single depth plane**, through the same aberrated probe and
box as the lab recon. `analysis/atomfind/extract_psf.py` removes the phase-ramp gauge, finds
the grid atoms in the interior of the scan field and **averages** them →
`psf_<el>_a<α>_vol.npy`, loaded by atomfind via `--single-atom-vol` / `--ti-kernel-vol`.
(Before 2026-09-11 it cropped one blob; `--single` still does, for reproducing old kernels.)

## What went wrong

`extract_psf.py` prints a cleanliness diagnostic (`|phase|>0.5` fraction — **small = clean**):

| α | Pb | Ti | argmax z | verdict |
|---|----|----|----------|---------|
| 50 | 0.10% | 0.07% | mid-crop | clean |
| 70 | 1.67% | 10.4% | mid-crop | marginal (Ti) |
| 90 | 67.5% | 61.6% | **last layer** | junk |
| 100 | 81.3% | (NaN'd) | **last layer** | junk |

At a90/a100 the extractor locks onto the **last z layer** instead of the box centre where the grid
atom is, and the volumes are speckle across the illuminated disc rather than isolated blobs.

## Diagnosis: conditioning, not geometry

The obvious suspect is **blob overlap** — the balanced probe grows with α (d90 = 4.0 / 3.9 / 6.6 /
11.0 Å at 50/70/90/100) against a fixed 4 Å grid spacing, and `build_atom_grid`'s docstring assumes
"spacing >> the ~1 Å blob".

**That suspect is wrong.** Measured on the clean a50 kernel: in-plane **FWHM 0.20 Å** against a
**4.01 Å** probe. Ptychography resolves ~20× below the probe size, so atoms 4 Å apart stay isolated
even at a100, and widening the grid would have wasted a sim round.

The real problem is an **under-constrained solve**: the grid object is a *single atomic plane*
carrying very little total scattering, reconstructed into 23–28 depth layers. The lab object
reconstructs fine at the same NL with the same settings — so NL=28 with no regularization is not
intrinsically broken; the grid object is simply too thin a constraint once the probe is large.
Compounding it, the grid legs scanned a **14 Å** window while the lab legs scan the sim default
**20 Å**, so the sparse solve also had less positional diversity than the dense one.

## Rejected fix: layer regularization

The first fix proposed here was `REGLAYER=0.5` on the grid legs. **It was wrong and must not be
used.** Two reasons, both fatal:

1. **It destroys the kernel's z shape.** Odstrčil's own parameter documentation
   (`ptycho/ptychography_exp_30nm.m`): *"0<R<<1 -> apply regularization on the reconstructed object
   layers … 0.01 == weak regularization that will slowly **symmetrize information content between
   layers**"* — i.e. it low-passes the depth axis, which is exactly the axial response the kernel
   exists to carry. `run_fusion_hollow.m` already encodes this: it warns that REGLAYER *"will
   low-pass the depth axis, which is the quantity this experiment measures"*, notes that
   **anything above 0.05 "destroys the signal"**, and prints `regularize_layers <- MUST be 0` in
   its preflight. The proposed 0.5 was ten times that ceiling.
2. **It breaks byte-identity.** Lab at 0 and grid at 0.5 are two different reconstruction
   operators, so the resulting blob is not the lab recon's PSF.

`REGLAYER` stays **0 on every leg**. Fix the object and the scan, never the operator.

## The fix ladder (all at REGLAYER = 0)

Ordered by how little they disturb byte-identity.

1. **S1 — restore scan identity and densify.** `WIN=20` (the lab's own window; now the driver
   default) and `GRIDSP=3`. Only the object and the scan change, and the scan change makes the two
   legs *more* alike. Try this first.
2. **S2 — densify further.** `GRIDSP=2.5`. There is ~20× headroom before blobs interfere, but keep
   spacing above `extract_psf --half-xy` (30 px = 1.48 Å) so the crop stays neighbour-free.
3. **S3 — sparse multi-plane grid** (the agreed last resort). Adds the depth density the single
   plane lacks while keeping every recon setting identical, so it is *more* faithful than thinning
   the box. Planes must be far enough apart that each blob's axial response is neighbour-free —
   spacing ≳ 4·δz. Viable only at high α, which is conveniently where it is needed:

   | α | δz | planes in the 19.5 Å atomic region |
   |---|----|------------------------------------|
   | 50 | 7.88 Å | none (too coarse) |
   | 70 | 4.02 Å | none |
   | 90 | 2.43 Å | 2 @ 10 Å |
   | 100 | 1.97 Å | 2 @ 8 Å |

**Not recommended:** thinning the PSF box (`--grid-box-z`) to cut the layer count. It changes `NL`,
which is one of the things that must match, so it buys conditioning at the cost of the rule.

## Extraction: the second half of the fix

S1 made the reconstructions clean — every grid is a 7×7 array of isolated atoms and every axial
profile peaks on the true atom depth (13.3–14.7 Å vs 13.76 Å) with flat vacuum bands. But the
original one-blob extractor still mis-picked three of eight kernels (Ti_a50 off the crop centre
with a spurious conjugation; both a70 at the exit edge), because it searched the whole object
including the junk scan edge, took a single blob, and never removed the phase-ramp gauge (0.6 rad
across the field at a50). A clean recon is necessary but not sufficient.

The extractor now:

1. **removes the ramp** with the `analyze_thin_campaign` model — one plane fitted to the
   depth-summed phase, spread evenly over the layers — masking atom pixels out of the fit;
2. **locates** the grid atoms in the inner 13 Å (clear of the scan edge) on a *background-removed*
   band, because a plane leaves a bowl in weak kernels that hid 22 of Ti_a50's 25 atoms;
3. **averages** them, with the kernel *values* taken from the deplaned phase, not the
   background-removed band — high-passing the kernel would reshape it, the same class of
   distortion as REGLAYER. Averaging N measurements of one response does not change the operator,
   so byte-identity holds; it is what the grid is for.

## How to tell a kernel is good

`extract_psf.py` prints what you need:

- number of atoms averaged (25 for `GRIDSP=3`, `WIN=20`) and the fraction that peak in the same
  layer — the grid is one plane, so they should agree;
- argmax at the **crop centre** (y = x = `--half-xy`), and in z on the atom plane, not a z edge;
- peak/background of tens at least; `|phase|>0.5` fraction a few % at most;
- in-plane FWHM sub-Å; axial FWHM comparable to δz = λ/α² for that α;
- **`--zdrop` must be the vacuum band in LAYERS**, `round(4.0/dz)` → 1/2/3/4 at a50/70/90/100. It is
  now required: the old default of 12 silently erased a NL=7 volume and reported "no blobs found".

## Running it

```bash
# rebuild ONLY the PSF kernels (lab legs untouched), S1 settings, all four alphas for uniformity
ALPHAS="50 70 90 100" MODES="Pb Ti" GRIDSP=3 WIN=20 OVERWRITE=1 BETA_LSQ=0.05 \
  bash campaign/run_thin_atomfind.sh
# then per alpha and element (R = round(4/dz))
python analysis/atomfind/extract_psf.py recon_af_a<A>_<El>_NL<NL> <El>_a<A> --zdrop R --out psf
```

## Current status (2026-09-11) — S1 kernels are good, no fallback needed

`atomfind_results_20260910_2109`, all eight kernels:

| kernel | atoms | argmax (y,x) | peak/bg | same-layer | in-plane FWHM |
|---|---|---|---|---|---|
| Pb_a50 | 25 | (30,30) | 178 | 100% | 0.15 Å |
| Ti_a50 | 25 | (30,30) | 44 | 92% | 0.25 Å |
| Pb_a70 | 25 | (30,30) | 379 | 100% | 0.10 Å |
| Ti_a70 | 25 | (30,30) | 61 | 100% | 0.20 Å |
| Pb_a90 | 25 | (30,30) | 1694 | 100% | 0.15 Å |
| Ti_a90 | 25 | (30,30) | 513 | 100% | 0.10 Å |
| Pb_a100 | 25 | (30,30) | 270 | 76%* | 0.15 Å |
| Ti_a100 | 25 | (30,30) | 70 | 100% | 0.15 Å |

\* The a100 grid atom sits at z = 13.762 Å and 14 × 0.983 = 13.76 Å — exactly on a layer boundary,
so its peak splits between two layers. Physical, not a fault.

The **data-derived flag is lifted**: every kernel now comes out of the byte-identical pipeline.
Matched kernels vs the data-derived stand-in on the same lab recons:

| α | bulk recall Pb / Ti / O | z-RMS | species confusion |
|---|---|---|---|
| 90 | 92/90/88% → **98/98/95%** | 0.72 → **0.41 Å** | → 0.5% |
| 100 | 96/88/84% → **100/98/86%** | 0.84 → **0.46 Å** | → 2.9% |

The kernel is worth a factor of ~2 in depth error. S2 and the multi-plane fallback were not needed.

**One open deviation from the rule — `beta_LSQ`.** The S1 kernels ran at `BETA_LSQ=0.05`; the
a50/a90/a100 lab recons ran at the `.m` default 0.1. `beta_LSQ` is a step size, not a penalty, so
unlike REGLAYER it does not change what is minimised — but with a fixed 200 iterations it can still
change where the solve lands, so the numbers above are not *strictly* byte-identical. The driver now
pins one `BETA_LSQ` (default 0.05) for every leg, and a70 / a110 run with it throughout. To close the
gap for the other three:
`ALPHAS="50 90 100" MODES=lab RECON_ONLY=1 bash campaign/run_thin_atomfind.sh`.
