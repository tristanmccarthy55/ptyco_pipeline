# PtychoShelves data contract (for generating synthetic 4D-STEM with abTEM)

This documents **exactly** how diffraction data, scan positions, and the initial
probe must be laid out and paired so this pipeline (Yu's MultiHollowPtycho fork)
ingests them correctly. It was reverse-engineered from the working PrScO3 example
and verified against the loader code. Source of truth, file by file:

- data load:      [`+detector/+prep_data/+matlab_aps/load_data.m`](ptycho/+detector/+prep_data/+matlab_aps/load_data.m)
- data process:   [`+detector/+prep_data/+matlab_aps/process_raw_data.m`](ptycho/+detector/+prep_data/+matlab_aps/process_raw_data.m)
- positions load: [`+scans/+positions/hdf5_pos.m`](ptycho/+scans/+positions/hdf5_pos.m)
- position gen:   [`ptycho/utils_EM/position_generate.m`](ptycho/utils_EM/position_generate.m)
- detector cfg:   [`+detector/+empad/empad.m`](ptycho/+detector/+empad/empad.m)

---

## 1. Files the pipeline reads, per scan

Everything for scan number `N` lives in a folder named with `p.scan.format`
(`'%02d'` → `01`, `02`, …), under `p.base_path`:

```
exampleData/
  01/
    data_dp.hdf5   (or data_dp.mat)   <- diffraction intensities
    data_position.hdf5                <- scan positions
    probe_initial.mat                 <- initial probe guess
```

`get_filenames_cu.m` looks for **`data_dp.hdf5` first, then `data_dp.mat`**.
Writing `data_dp.hdf5` from Python (h5py) is the easiest route — no MATLAB v7.3
gymnastics.

---

## 2. `data_dp` — diffraction patterns

### What it holds
- A real, **non-negative INTENSITY** stack (detector counts), NOT amplitude.
  The engine forms the modulus itself: `fmag = sqrt(dp / ADU)` downstream.
- MATLAB shape: **`[Ndet, Ndet, Npos]`** = `[detector_y, detector_x, scan_point]`.
  For the PrScO3 example: `[256, 256, 4096]`.

### Center / orientation
- The zero-order (bright-field) disk must sit at the **array centre pixel**
  `ctr = [floor(Ndet/2)+1, floor(Ndet/2)+1]` (= `[129,129]` for 256).
  i.e. patterns are **`fftshift`ed** (DC in the middle), not corner-origin.
  `process_raw_data.m` crops `ctr ± Ndet/2`; with a centred 256 pattern this is
  the identity.
- `empad` detector orientation is `[0 0 0]` (no transpose, no flip) — so `dp` is
  used **exactly as stored**. Any mirror/rotation you bake in is what the engine
  sees. If a recon comes out flipped/rotated, fix it either in your generator or
  by setting `det.orientation = [<transpose> <fliplr> <flipud>]` in `empad.m`.

### ADU
- `ADU` = counts per electron (580 for the EMPAD example). The driver divides by
  it. For synthetic data, set `ADU = 1` and emit electron-intensity directly
  (then set `ADU = 1` in the driver), or pick a dose and a realistic ADU.

### Scan-point ordering — **THE critical pairing rule**
The 3rd axis (`Npos`) must be ordered **Y-fastest** (column-major over a
`[npy, npx]` grid, with the row/Y index varying fastest):

```
k = iy + (ix-1)*npy        % iy = 1..npy (fast),  ix = 1..npx (slow)
```

Equivalently: build `cbed[dy, dx, iy, ix]` and `reshape` to `[dy, dx, npy*npx]`
in MATLAB (column-major) — `iy` runs fastest. This **must** match the position
ordering in §3 (it does, because `position_generate` uses the same convention).

> Verified from the example: in `data_position.hdf5`, x stays constant while y
> steps by 0.41 Å for the first 64 entries → y is the fast axis.

---

## 3. `data_position.hdf5` — scan positions

- One dataset: **`/probe_positions_0`**.
- MATLAB shape **`[Npos, 2]`**, units **Ångström**, column 1 = **x**, column 2 = **y**.
  (In h5py this reads back as shape `(2, Npos)` — row 0 = x, row 1 = y — because
  HDF5/MATLAB reverse axis order.)
- Same **Y-fastest** ordering as `dp` (point `k` of positions ↔ frame `k` of `dp`).
- Centred raster, per `position_generate` → `scan_position_rot`:
  ```
  ppx = linspace(-floor(npx/2), ceil(npx/2)-1, npx) * step_x;   % Å
  ppy = linspace(-floor(npy/2), ceil(npy/2)-1, npy) * step_y;   % Å
  [ppX, ppY] = meshgrid(ppx, ppy);
  positions  = [ppX(:), ppY(:)];     % [x, y], y-fastest
  ```

### Sign/axis convention applied downstream (don't pre-apply it)
`hdf5_pos.m` converts to the engine's internal real coordinates as:
```
positions_real(:,1) = -y    % row
positions_real(:,2) = -x    % col
```
i.e. it **negates and swaps** x/y. You do **not** apply this yourself — just emit
`[x, y]` in Å in the same scan order as `dp`, and let the loader do it. (Because
the raster is symmetric about 0, the negation only sets handedness; what matters
is that `dp` and positions share one ordering.)

---

## 4. `probe_initial.mat` — initial probe

- Variable **`probe`**: complex `[Ndet, Ndet]` (`complex64` is fine) — the
  real-space probe wavefunction.
- Variable **`p`**: a small struct with `p.binning = false` and
  `p.detector.binning = false`.
- Normalisation: scale so `sum(|probe|^2) ≈ Itot` (mean total counts per pattern).
  Not critical — `initial_probe_rescaling = true` fixes the global scale on
  iteration 1. You can use the bundled `generateProbeFunction` or export abTEM's
  probe directly.

---

## 5. Geometry — numbers that must be self-consistent

These tie the detector sampling, probe, and object pixel together. Pick them in
abTEM and mirror them in the driver:

```
lambda  = 12.3986 / sqrt((2*511.0 + V)*V)     % electron wavelength [Å], V in keV
dk      = (alpha0*1e-3) / rbf / lambda          % recip-space pixel [1/Å per det px]
d_alpha = alpha0 / rbf                           % [mrad per detector pixel]
p.z     = 1 / (d_alpha*1e-3)                     % [rad^-1]  (driver sets this)
dx_obj  = 1 / (Ndet * dk)                        % real-space/object pixel [Å]
```
where `alpha0` = probe convergence semi-angle [mrad], `rbf` = bright-field disk
radius in detector **pixels**. So your abTEM detector angular calibration must be
`d_alpha = alpha0/rbf` mrad/pixel, and the BF disk must span `rbf` px.

PrScO3 example values: `V=300`, `alpha0=21.4`, `rbf=26` → `d_alpha=0.823 mrad/px`,
`dk=0.823 1/Å`-equiv, `scanstep=0.41 Å`, `Ndet=256`, `ADU=580`, `thickness=210 Å`.

---

## 6. End-to-end flow

```
abTEM multislice 4D-STEM
   │   (Npos diffraction patterns, intensity, centred; probe; scan grid)
   ▼
write per-scan folder NN/:
   data_dp.hdf5 (/dp)         [Ndet,Ndet,Npos], y-fastest, centred intensity
   data_position.hdf5 (/probe_positions_0)  [Npos,2] Å, [x,y], y-fastest
   probe_initial.mat          probe[Ndet,Ndet] complex, p struct
   ▼
matlab_aps loader: load dp -> sum bursts -> crop to ctr±asize/2
   ▼  hdf5_pos: positions Å -> m, negate+swap -> positions_real
fmag = sqrt(dp/ADU)  ->  LSQ-ML / MLs engine  ->  object + probe
```

---

## 7. Python recipe (abTEM → these files)

```python
import numpy as np, h5py
from scipy.io import savemat

Ndet = 256
npy, npx = 64, 64          # scan grid (rows=y, cols=x)
step = 0.41                # Å

# --- diffraction: I has shape (npy, npx, Ndet_y, Ndet_x), centred INTENSITY ---
# (fftshift each pattern so the BF disk is at the centre pixel)
I = abtem_intensity            # <-- your abTEM output, real >= 0

# flatten scan to y-fastest linear order to match positions:
#   MATLAB column-major over [npy,npx] => reshape with order='F' on (npy,npx)
I = I.reshape(npy*npx, Ndet, Ndet, order='F')      # (Npos, dy, dx), y-fastest
Npos = I.shape[0]

# write data_dp.hdf5 so MATLAB h5read('/dp') returns [Ndet_y, Ndet_x, Npos].
# h5py shape (s0,s1,s2) reads in MATLAB as [s2,s1,s0]; and M[a,b,c]=D[c,b,a],
# so to get M[dy,dx,k] we store D[k,dx,dy] = I[k,dy,dx]  => transpose last two axes.
D = np.ascontiguousarray(I.transpose(0, 2, 1)).astype('float32')   # (Npos, dx, dy)
with h5py.File('01/data_dp.hdf5', 'w') as f:
    f.create_dataset('dp', data=D)

# --- positions: [x, y] in Å, y-fastest, centred raster ---
ppx = np.linspace(-(npx//2), (npx - npx//2) - 1, npx) * step
ppy = np.linspace(-(npy//2), (npy - npy//2) - 1, npy) * step
ppX, ppY = np.meshgrid(ppx, ppy)                    # shape (npy, npx)
pos = np.stack([ppX.ravel(order='F'), ppY.ravel(order='F')], axis=1)  # (Npos,2) [x,y]
# store as (2, Npos) so MATLAB sees [Npos, 2]:
with h5py.File('01/data_position.hdf5', 'w') as f:
    f.create_dataset('probe_positions_0', data=pos.T.astype('float32'))

# --- initial probe ---
probe = abtem_probe_realspace.astype('complex64')   # (Ndet, Ndet)
savemat('01/probe_initial.mat',
        {'probe': probe, 'p': {'binning': False, 'detector': {'binning': False}}})
```

> **Verify orientation once.** Detector dy/dx transpose and scan x/y handedness
> are the easy things to get wrong. Put a single, deliberately **off-centre** test
> feature in one pattern, run a few iterations, and confirm the reconstructed
> object isn't mirrored/rotated/transposed vs. ground truth. If it is, flip the
> last-two-axis transpose above, or set `det.orientation` in `empad.m`. The safest
> reference is the working `sample_data_PrScO3.mat` layout (h5py `dp` shape
> `(Npos, Ndet, Ndet)`), which this recipe reproduces.

---

## 8. Driver knobs to change for your data

In [`ptycho/run_baseline_PrScO3.m`](ptycho/run_baseline_PrScO3.m): `Ndpx`,
`d_alpha (= alpha0/rbf*1e-3)`, `HT`, `ADU`, `thick`, `Nlayers`, `scanstep`,
and the scan grid (`npx`, `npy` if you regenerate positions). Keep
`global mask1 = ones(Ndpx)` for the hollow-angle-0 baseline.
