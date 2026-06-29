#!/usr/bin/env python
"""Recon-vs-dose comparison figure: depth cross-section + in-plane slice per dose.

Run AFTER pulling the dose recons to the Mac. Point DATA_ROOT at the folder holding
the per-dose recon dirs (each a recon_*dose<D>*/.../Niter*.mat). Produces a 2xD grid:
  top row    = depth cross-section down the SAME Pb column at each dose
  bottom row = the SAME in-plane slice at each dose
so you can read off how depth and in-plane fidelity degrade as the dose drops.

  python dose_compare.py
"""
import os, glob, numpy as np, h5py
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import maximum_filter

DATA_ROOT = os.path.expanduser("~/Desktop/dose_series")
OUT       = os.path.expanduser("~/Desktop/dose_compare.png")
DOSES     = ["1e10", "1e8", "1e6", "1e4"]        # high -> low
DX, W     = 0.0492, 22                            # object px (A); half cross-section strip


def load_vol(dose):
    mats = glob.glob(os.path.join(DATA_ROOT, f"*dose{dose}[!0-9]*", "**", "Niter*.mat"), recursive=True) \
        or glob.glob(os.path.join(DATA_ROOT, f"*dose{dose}", "**", "Niter*.mat"), recursive=True)
    if not mats:
        return None
    m = sorted(mats, key=lambda p: int("".join(filter(str.isdigit, os.path.basename(p)))))[-1]
    with h5py.File(m, "r") as f:
        g = f["outputs"]
        layers = []
        for r in g["object_roi"][:, 0]:
            a = f[r][:]
            a = (a["real"] + 1j * a["imag"]) if a.dtype.names else a
            layers.append(a.T)
        vol = np.angle(np.array(layers)).astype(float)
        vol -= np.median(vol, (1, 2), keepdims=True)
        z = g["z_distance"][:, 0]
        dz = float(np.median(z[np.isfinite(z)])) * 1e10
    return vol, dz


vols = {d: load_vol(d) for d in DOSES}
have = [d for d in DOSES if vols[d] is not None]
if not have:
    raise SystemExit(f"no dose recons found under {DATA_ROOT} (expected recon_*dose*/.../Niter*.mat)")
ref, refdz = vols[have[0]]                          # cleanest available, for column/slice choice
nL = ref.shape[0]

# pick a strong Pb column (away from edges) on the reference, reuse for every dose
dm = ref.mean(0); dmn = dm - dm.min()
mask = np.zeros_like(dmn, bool); mask[40:-40, 40:-40] = True
pk = (dmn == maximum_filter(dmn, 25)) & (dmn > np.percentile(dmn[mask], 97)) & mask
ys, xs = np.where(pk); j = np.argmax(dmn[ys, xs]); yc, xc = ys[j], xs[j]
ksl = nL // 2                                       # a mid-depth in-plane slice

fig, ax = plt.subplots(2, len(DOSES), figsize=(3.6 * len(DOSES), 7.4))
for i, d in enumerate(DOSES):
    a0, a1 = ax[0, i], ax[1, i]
    if vols[d] is None:
        for a in (a0, a1): a.text(.5, .5, f"{d}\nmissing", ha="center", va="center"); a.axis("off")
        continue
    vol, dz = vols[d]
    cs = vol[:, yc-1:yc+2, xc-W:xc+W].mean(1)
    a0.imshow(cs, extent=[-W*DX, W*DX, (nL-0.5)*dz, 0.5*dz], aspect="auto", cmap="inferno",
              vmin=np.percentile(cs, 5), vmax=np.percentile(cs, 99.3))
    a0.set_title(f"{d} e/Å²", fontsize=12)
    a0.set_xlabel("in-plane (Å)")
    sl = vol[ksl]
    a1.imshow(sl, cmap="inferno", vmin=np.percentile(sl, 5), vmax=np.percentile(sl, 99.5))
    a1.axis("off")
ax[0, 0].set_ylabel("depth cross-section\nz (Å)  [entrance→exit]")
ax[1, 0].set_ylabel(f"in-plane slice {ksl}/{nL}")
ax[1, 0].axis("on"); ax[1, 0].set_xticks([]); ax[1, 0].set_yticks([])
fig.suptitle("Reconstruction vs electron dose — depth cross-section (top) & in-plane slice (bottom)",
             fontsize=13)
fig.tight_layout()
fig.savefig(OUT, dpi=130)
print("wrote", OUT, "| doses present:", have, "| column", (yc, xc), "| slice", ksl)
