#!/usr/bin/env python
"""@file extract_psf.py
@brief Extract an empirical PSF kernel from a grid-PSF reconstruction.

A lone atom is too under-constrained to reconstruct, but a sparse grid of one element on a single
depth plane converges like the real sample while each blob stays isolated -- and every blob IS the
system PSF. So this finds the grid atoms in the interior of the scan field and AVERAGES them: the
grid exists precisely so the same response is measured many times, and N blobs cut the noise by
~sqrt(N) without touching the reconstruction operator (the kernel stays byte-identical to the lab
recon; see aberration_experiment/PSF_KERNELS.md). Before averaging it removes the in-plane phase
ramp that fixed-probe recons carry -- a gauge, not part of the response -- with the same model as
analyze_thin_campaign.crop_phase (one plane fitted to the depth-summed phase, spread evenly over
the layers). The result is written as a complex volume in the NL70_new_vol.npy format, with the
atom as the single positive phase maximum, so psf.empirical_psf() (angle -> argmax -> crop)
consumes it directly via cfg.single_atom_vol / cfg.ti_kernel_vol.

  python analysis/atomfind/extract_psf.py <recon_dir> <name> --zdrop R [--out DIR]

--zdrop is the z-vacuum band in LAYERS dropped at each end: round(z_vacuum / dz), i.e.
1/2/3/4 at a50/70/90/100 for the thin campaign (4 A vacuum). --single restores the original
one-blob extraction, for reproducing kernels made before 2026-09-11.
<recon_dir> holds .../Niter*.mat (newest used) or *_recons.h5. Writes psf_<name>_vol.npy (+ _check.png).
"""
import argparse, glob, os
import numpy as np
import h5py
from scipy.ndimage import gaussian_filter, maximum_filter

DX_DEFAULT = 0.0492          # A/px, the campaign object pixel (p/dx_spec)


def load_vol(recon_dir):
    mats = glob.glob(os.path.join(recon_dir, "**", "Niter*.mat"), recursive=True)
    if not mats:
        # newer PtychoShelves writes *_recons.h5 (reconstruction/object) instead of Niter*.mat
        rec = glob.glob(os.path.join(recon_dir, "**", "*_recons.h5"), recursive=True)
        if not rec:
            raise SystemExit(f"no Niter*.mat or *_recons.h5 under {recon_dir}")
        m = sorted(rec, key=os.path.getmtime)[-1]
        with h5py.File(m, "r") as f:
            obj = f["reconstruction"]["object"][:]           # (NL,1,Ny,Nx) or (NL,Ny,Nx) complex
        obj = np.asarray(obj).squeeze()
        if obj.ndim == 2:
            obj = obj[None]
        return obj.astype(np.complex64), m
    m = sorted(mats, key=lambda p: int("".join(filter(str.isdigit, os.path.basename(p)))))[-1]
    with h5py.File(m, "r") as f:
        L = []
        for r in f["outputs"]["object_roi"][:, 0]:
            a = f[r][:]
            a = (a["real"] + 1j * a["imag"]) if a.dtype.names else a
            L.append(a.T)
    return np.array(L), m


def deplane(ph):
    """Remove the in-plane phase-ramp gauge: one plane fitted to the depth-summed phase, spread
    evenly over the layers (the analyze_thin_campaign model). Atom pixels are masked out of the
    fit so a few dozen bright blobs cannot tilt it. Returns (ph, ramp span across the field)."""
    s = ph.sum(0)
    yy, xx = np.mgrid[0:s.shape[0], 0:s.shape[1]].astype(float)
    dev = np.abs(s - np.median(s))
    m = dev < np.percentile(dev, 97)
    A = np.c_[xx[m], yy[m], np.ones(m.sum())]
    c, *_ = np.linalg.lstsq(A, s[m], rcond=None)
    plane = (c[0] * xx + c[1] * yy) / ph.shape[0]
    return ph - plane[None], float(plane.max() - plane.min())


def find_sites(band, inner_px, sep_px):
    """Grid-atom sites: local maxima of the depth-mean phase well above a robust noise floor,
    inside the inner box (away from the scan-field edge, where the recon is least constrained)."""
    med = np.median(band); mad = 1.4826 * np.median(np.abs(band - med))
    pk = (band == maximum_filter(band, sep_px)) & (band > med + 6 * mad)
    cy, cx = band.shape[0] / 2, band.shape[1] / 2
    ys, xs = np.where(pk)
    keep = (np.abs(ys - cy) <= inner_px) & (np.abs(xs - cx) <= inner_px)
    return ys[keep], xs[keep]


def fwhm_px(line):
    """Contiguous half-max width around the peak of a 1-D profile (px)."""
    k = int(np.argmax(line)); h = line[k] / 2; lo = hi = k
    while lo > 0 and line[lo - 1] >= h: lo -= 1
    while hi < len(line) - 1 and line[hi + 1] >= h: hi += 1
    return hi - lo + 1


def average_blobs(cvol, a):
    """The default: deplane, find the interior grid atoms, average their crops. -> (phase, info)."""
    nL, Ny, Nx = cvol.shape
    z0, z1 = a.zdrop, nL - a.zdrop
    h = int(round(a.fov / 2 / a.dx)); W = a.half_xy
    ph = np.angle(cvol[:, Ny // 2 - h:Ny // 2 + h, Nx // 2 - h:Nx // 2 + h]).astype(float)
    ph -= np.median(ph, (1, 2), keepdims=True)
    ph, span = deplane(ph)
    ph -= np.median(ph, (1, 2), keepdims=True)

    # Sites are LOCATED on a background-removed band, because a plane does not take out all of a
    # weak kernel's smooth background (Ti at a50 keeps a bowl that inflates the noise floor 12x
    # and hides 22 of its 25 atoms). The kernel VALUES below still come from the deplaned phase:
    # high-passing the kernel itself would reshape it, the same distortion REGLAYER causes.
    band = ph[z0:z1].mean(0)
    det = band - gaussian_filter(band, 1.0 / a.dx)
    # sign from the strongest local extrema: grid atoms are one sign, noise is both
    ext = (np.abs(det) == maximum_filter(np.abs(det), int(round(1.0 / a.dx))))
    top = np.sort(np.abs(det[ext]))[-9:]
    sgn = np.sign(np.median(det[ext][np.isin(np.abs(det[ext]), top)]))
    ph *= sgn; det *= sgn

    ys, xs = find_sites(det, int(round(a.inner / 2 / a.dx)), int(round(1.0 / a.dx)))
    ok = (ys - W >= 0) & (ys + W <= det.shape[0]) & (xs - W >= 0) & (xs + W <= det.shape[1])
    ys, xs = ys[ok], xs[ok]
    if len(ys) < 3:
        raise SystemExit(f"only {len(ys)} grid atoms found in the inner {a.inner} A -- "
                         f"is this a grid-PSF recon? (use --single for a lone atom)")
    crops = np.stack([ph[z0:z1, y - W:y + W, x - W:x + W] for y, x in zip(ys, xs)])
    zpk = np.argmax(crops[:, :, W, W], axis=1)
    mode = np.bincount(zpk).argmax()
    if len(ys) >= 2:
        d = np.hypot(ys[:, None] - ys[None], xs[:, None] - xs[None]); d[d == 0] = np.inf
        spacing = float(np.median(d.min(1))) * a.dx
    else:
        spacing = np.nan
    K = crops.mean(0)
    info = dict(n=len(ys), sign="+" if sgn > 0 else "-", ramp=span, spacing=spacing,
                zmode=int(mode), zagree=float(np.mean(zpk == mode)),
                noise=float(crops[:, :, W, W].std(0)[mode] / np.sqrt(len(ys))))
    return K, info


def single_blob(cvol, a):
    """The original extraction (pre 2026-09-11): the one blob nearest the object centre, no ramp
    removal, sign by |min| vs |max|. Kept only so earlier kernels can be reproduced."""
    nL, Ny, Nx = cvol.shape
    ph = np.angle(cvol); ph -= np.median(ph, (1, 2), keepdims=True)
    z0, z1 = a.zdrop, nL - a.zdrop
    band = np.abs(ph[z0:z1]).mean(0); s = band - band.min()
    pk = (s == maximum_filter(s, 9)) & (s > np.percentile(s, 99))
    ys, xs = np.where(pk)
    if len(ys) == 0:
        raise SystemExit("no blobs found in the mid-depth band")
    d = (ys - Ny / 2) ** 2 + (xs - Nx / 2) ** 2
    W = a.half_xy
    yy = int(np.clip(ys[np.argmin(d)], W, Ny - W)); xx = int(np.clip(xs[np.argmin(d)], W, Nx - W))
    K = np.angle(cvol[z0:z1, yy - W:yy + W, xx - W:xx + W])
    K -= np.median(K, (1, 2), keepdims=True)
    sgn = -1.0 if abs(K.min()) > abs(K.max()) else 1.0
    return sgn * K, dict(n=1, sign="+" if sgn > 0 else "-")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("recon_dir")
    ap.add_argument("name", help="element/geometry tag, e.g. Pb_a100 -> psf_Pb_a100_vol.npy")
    ap.add_argument("--out", default=os.path.expanduser("~/Desktop"))
    ap.add_argument("--zdrop", type=int, required=True,
                    help="z-vacuum layers to drop at each end = round(z_vacuum/dz)")
    ap.add_argument("--half-xy", type=int, default=30, help="in-plane crop half-width (px)")
    ap.add_argument("--dx", type=float, default=DX_DEFAULT, help="object pixel (A)")
    ap.add_argument("--fov", type=float, default=16.0, help="central field processed (A)")
    ap.add_argument("--inner", type=float, default=13.0, help="atoms accepted within this box (A)")
    ap.add_argument("--single", action="store_true", help="legacy one-blob extraction")
    a = ap.parse_args()

    cvol, m = load_vol(a.recon_dir)
    nL = cvol.shape[0]
    if 2 * a.zdrop >= nL:
        raise SystemExit(f"--zdrop {a.zdrop} leaves nothing of a {nL}-layer volume "
                         f"(it should be round(z_vacuum/dz))")
    K, info = single_blob(cvol, a) if a.single else average_blobs(cvol, a)

    os.makedirs(a.out, exist_ok=True)
    out = os.path.join(a.out, f"psf_{a.name}_vol.npy")
    np.save(out, np.exp(1j * K).astype(np.complex64))

    l, r, c = np.unravel_index(np.argmax(K), K.shape)
    W = a.half_xy
    print(f"src : {m}")
    print(f"wrote {out}  shape {K.shape}")
    if a.single:
        print(f"[single] 1 blob, sign {info['sign']}")
    else:
        print(f"[avg] {info['n']} grid atoms (spacing {info['spacing']:.2f} A), sign {info['sign']}, "
              f"ramp removed {info['ramp']:.3f} rad across the field; peak layer agrees for "
              f"{info['zagree']:.0%} of atoms (layer {info['zmode']})")
    print(f"agent argmax at (z={l}, y={r}, x={c}) -- centre is y=x={W}; "
          f"|phase|>0.5 fraction {(np.abs(K) > 0.5).mean():.2%} (small = clean single blob)")
    fx = fwhm_px(K[l, r, :]) * a.dx; fz = fwhm_px(K[:, r, c])
    bg = np.std(K[l][np.hypot(*np.meshgrid(np.arange(K.shape[2]) - c, np.arange(K.shape[1]) - r)) > 0.75 * W])
    print(f"peak {K[l, r, c]:.3f} rad, bg sd {bg:.4f} -> peak/bg {K[l, r, c] / bg:.0f}; "
          f"in-plane FWHM {fx:.2f} A, axial FWHM {fz} layer(s)")
    if not a.single and (r, c) != (W, W):
        print(f"WARNING: argmax is off the crop centre by ({r - W},{c - W}) px -- inspect _check.png")

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        p = np.abs(K)
        fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
        ax[0].imshow(p.max(0), cmap="inferno"); ax[0].set_title(f"psf_{a.name}: in-plane (max-proj)")
        ax[1].imshow(p[:, r, :], aspect="auto", cmap="inferno"); ax[1].set_title("axial (z vs x)")
        ax[1].set_ylabel("z layer")
        ax[2].plot(K[:, r, c], "-o"); ax[2].axhline(0, color="0.6", lw=0.8)
        ax[2].set_title("axial profile through the peak"); ax[2].set_xlabel("z layer")
        fig.tight_layout(); fig.savefig(os.path.join(a.out, f"psf_{a.name}_check.png"), dpi=110)
        print(f"wrote psf_{a.name}_check.png")
    except Exception as e:
        print("(no figure:", e, ")")


if __name__ == "__main__":
    main()
