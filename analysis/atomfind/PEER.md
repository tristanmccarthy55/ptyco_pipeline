# Reproducing the nominated result

One page. Two commands. Everything needed is in this repository plus one data tarball.

The result offered for independent reproduction is **three-dimensional atom localisation with
calibrated uncertainty**, and the polarisation map derived from it. The simulation and
ptychographic reconstruction that produce the input volume need GPU hours, so the volume is
supplied pre-computed; **everything downstream of it is reproduced from scratch**.

---

## 1. Inputs

| File | What it is | Size |
|---|---|---|
| `NL70_new_vol.npy` | the reconstructed phase volume, complex64 `(70, 404, 404)`, dz 0.999 Å, dx 0.0495 Å/px | 91 MB |
| `psf_Pb_NL70_vol.npy` | the measured single-lead-atom response (the forward-model kernel) | 1.3 MB |
| `psf_Ti_NL70_vol.npy` | the measured single-titanium response | 1.3 MB |
| `gt_prepared.npz` | the reference structure, already rotated/orthogonalised into the beam frame | 0.4 MB |
| `PTO6_STO6_18_18_labyrinthPoscar.vasp` | the same structure, unprepared (only needed to rebuild the cache) | 1.4 MB |

Unpack the tarball anywhere; point the pipeline at it with `--data-dir`, or drop the files in
`atomfind/data/`, or set `$ATOMFIND_DATA`. There are **no absolute paths anywhere in the code**.

## 2. Install and run

Use **Python 3.10 or newer** if you have it. The pipeline runs correctly on 3.9, but pip then
resolves scipy to 1.13, whose NNLS raises tens of thousands of spurious floating-point warnings
on this problem; the results are identical either way (checked: largest relative difference
across `report.json` is 1e-11, and no non-finite value reaches any export).

```bash
python3 -m venv venv && . venv/bin/activate
pip install -r atomfind/requirements.txt

tar xzf atomfind_data_v1.tar.gz -C atomfind/data/

python atomfind/run_atomfind.py --preset NL70_coherent --out ./out
python atomfind/polarisation.py --out ./out
```

Run from the directory that contains `atomfind/`. A few minutes on one core; no GPU.

## 3. What you should get

`run_atomfind.py` writes `found_atoms.csv` (element, x/y/z in ångström in the reference frame,
per-axis 95% conformal half-widths, model σ, species confidence, and a flag marking
lattice-constrained detections), `found_atoms.extxyz`, `uq_conformal.json` and `report.json`.
`polarisation.py` writes `polarisation.json` / `.npz`.

On this input the run is deterministic and should reproduce

| Quantity | Expected |
|---|---|
| atoms found | 1834 |
| precision | 0.97 |
| bulk recall, Pb / Ti / O | 95.6 / 96.4 / 95.5 % |
| species confusion | 1.1 % |
| in-plane RMS accuracy | 0.032 Å |
| depth RMS accuracy | 0.37 Å |
| in-plane polarisation error (median) | 0.007 Å, 0.9° |
| propagated σ on the along-beam component | 0.237 Å |

**The uncertainty is the point, not the coordinates.** The number that validates it is the
coverage: 96% of true errors fall inside the nominal 95% interval, per stratum, held out on a
50/50 split. Typical half-widths are 0.02 Å in-plane and 0.5 Å in depth, rising to 1.5 Å for
weak oxygen near the exit surface. They are strongly heteroscedastic, so propagate the
**per-atom** interval, never a global figure.

## 4. Tolerances

The only stochastic element is the conformal calibration split, which is seeded. Re-seeding
(`polarisation.py --seed N`) moves the per-stratum quantiles by under 5% and the coverage
figures by under one percentage point. Coverage is finite-sample, so a stratum with fewer
than ~20 atoms — the constrained-oxygen entrance band here — shows several percentage points
of conformal noise and should not be read as a failure.

## 5. Failure modes you are expected to see

These are reported, not absorbed, and both appear in the output:

- **~28% of bulk titanium do not acquire a complete oxygen cage.** They are excluded from the
  polarisation map rather than completed from the lattice.
- **~5% have a cage of the wrong composition.** Their in-plane error rises to a median of
  0.45 Å, and the positional interval covers only 12% of them. This is a *detection* failure,
  not a localisation failure: conformal intervals are conditional on correct detection. It is
  visible in the exported cage completeness and species confidence, not in the error bars.

If your numbers differ by more than the tolerances above, the most likely causes are a
truncated download (check `sha256sum -c SHA256SUMS`) or a different reference structure.
