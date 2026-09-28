#!/usr/bin/env python
"""@file render_phase.py
@brief Phase images of a reconstruction: the depth sum, every slice, and an x-z section, one PNG per recon.

PtychoShelves' own image saving is off in this pipeline (save.store_images=0), so nothing else draws the object.
For each *_recons.h5 (or recon dir, whose newest h5 outside analysis.prev_* is used) this writes <out>/<name>.png:
  left   the phase summed over all slices, cropped to the illuminated field (illum_sum > 0.15 max)
  right  every slice on ONE colour scale (set by the atomic slices), labelled with its centre depth; slices whose
         centre lies in the z-vacuum bands are marked "vac" -- a vacuum slice that shows atoms is depth leakage
  bottom an x-z section: the phase averaged over a 1 A band in y through the field centre, true aspect
The in-plane phase ramp is removed as extract_psf does (one plane fitted to the depth sum, spread over the slices);
point-sampled sims carry almost none. Depths come from sim_meta beam_thickness_A when the recon dir links it, else
--box.

  ~/hyperspy-bundle/bin/python analysis/render_phase.py <h5 or recon dir> [...] --out DIR [--z-vacuum 4] [--box 27.525]
"""
import argparse, glob, os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def find_h5(path):
    if os.path.isfile(path):
        return path
    hs = [h for h in glob.glob(os.path.join(path, "**", "*_recons.h5"), recursive=True) if "analysis.prev" not in h]
    return max(hs, key=os.path.getmtime) if hs else None


def box_thickness(h5, default):
    """beam_thickness_A from the leg's sim_meta.mat (recon dir 01/ link, or a packed sim_out_af_<leg>/01/)."""
    from scipy.io import loadmat
    rdir = h5.split(os.sep + "analysis" + os.sep)[0]
    leg = os.path.basename(rdir).replace("recon_af_", "", 1).rsplit("_NL", 1)[0]
    for m in (os.path.join(rdir, "01", "sim_meta.mat"),
              os.path.join(os.path.dirname(rdir), f"sim_out_af_{leg}", "01", "sim_meta.mat")):
        if os.path.exists(m):
            return float(loadmat(m, simplify_cells=True)["meta"]["beam_thickness_A"]), os.path.basename(rdir)
    return default, os.path.basename(rdir)


def load(h5):
    with h5py.File(h5, "r") as f:
        obj = np.asarray(f["reconstruction/object"]).squeeze()
        dx = float(np.ravel(f["reconstruction/p/dx_spec"][...])[0]) * 1e10
        il = np.asarray(f["reconstruction/p/illum_sum/illum_sum_0"]).squeeze()
    if obj.ndim == 2:
        obj = obj[None]
    return obj, dx, il


def render(h5, out, z_vac, box_default):
    obj, dx, il = load(h5)
    box, name = box_thickness(h5, box_default)
    nl = obj.shape[0]; dz = box / nl
    m = il > 0.15 * il.max()
    ys, xs = np.where(m)
    y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
    ph = np.angle(obj[:, y0:y1, x0:x1]).astype(np.float64)
    ph -= ph.mean(axis=(1, 2), keepdims=True)
    s = ph.sum(0)                                            # deplane, as extract_psf.deplane
    yy, xx = np.mgrid[0:s.shape[0], 0:s.shape[1]].astype(float)
    dev = np.abs(s - np.median(s)); w = dev < np.percentile(dev, 97)
    c, *_ = np.linalg.lstsq(np.c_[xx[w], yy[w], np.ones(w.sum())], s[w], rcond=None)
    ph -= ((c[0] * xx + c[1] * yy) / nl)[None]
    s = ph.sum(0)
    zc = (np.arange(nl) + 0.5) * dz
    vac = (zc < z_vac) | (zc > box - z_vac)
    atomic = ph[~vac] if (~vac).any() else ph
    lo, hi = np.percentile(atomic, [0.5, 99.8])
    lim = max(abs(lo), abs(hi))
    ext = [0, s.shape[1] * dx, s.shape[0] * dx, 0]

    ncol = int(np.ceil(np.sqrt(nl * 1.6)))
    nrow = int(np.ceil(nl / ncol))
    fig = plt.figure(figsize=(4.2 + 1.9 * ncol, max(5.2, 1.9 * nrow + 2.4)))
    gs = fig.add_gridspec(nrow + 1, ncol + 2, height_ratios=[1] * nrow + [0.9], width_ratios=[1.6, 1.6] + [1] * ncol)
    ax = fig.add_subplot(gs[:nrow, :2])
    im = ax.imshow(s, cmap="inferno", extent=ext, vmin=np.percentile(s, 0.5), vmax=np.percentile(s, 99.8))
    ax.set_title(f"{name}\nphase summed over {nl} slices (rad)", fontsize=9)
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    for l in range(nl):
        a = fig.add_subplot(gs[l // ncol, 2 + l % ncol])
        a.imshow(ph[l], cmap="RdBu_r", vmin=-lim, vmax=lim, extent=ext)
        a.set_title(f"{l + 1}: z {zc[l]:.1f} Å" + ("  vac" if vac[l] else "") + f"\nstd {ph[l][m[y0:y1, x0:x1]].std():.3f}",
                    fontsize=7, color="tab:blue" if vac[l] else "k")
        a.set_xticks([]); a.set_yticks([])
    band = max(1, int(round(0.5 / dx)))
    cy = s.shape[0] // 2
    xz = ph[:, cy - band:cy + band + 1, :].mean(1)
    a = fig.add_subplot(gs[nrow, :])
    a.imshow(xz, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto", extent=[0, s.shape[1] * dx, box, 0],
             interpolation="nearest")
    for zb in (z_vac, box - z_vac):
        a.axhline(zb, color="k", lw=0.6, ls="--")
    a.set_xlabel("x (Å)"); a.set_ylabel("z (Å)")
    a.set_title(f"x-z section (y band ±0.5 Å through the centre); dashed = vacuum bands; slice colour scale ±{lim:.2f} rad",
                fontsize=8)
    fig.tight_layout()
    os.makedirs(out, exist_ok=True)
    p = os.path.join(out, f"{name}.png")
    fig.savefig(p, dpi=110); plt.close(fig)
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="*_recons.h5 files or recon dirs")
    ap.add_argument("--out", required=True)
    ap.add_argument("--z-vacuum", type=float, default=4.0)
    ap.add_argument("--box", type=float, default=27.525, help="beam thickness (A) if no sim_meta.mat is found")
    a = ap.parse_args()
    for p in a.paths:
        h5 = find_h5(p)
        if not h5:
            print(f"{p}: no *_recons.h5"); continue
        print("wrote", render(h5, a.out, a.z_vacuum, a.box), flush=True)


if __name__ == "__main__":
    main()
