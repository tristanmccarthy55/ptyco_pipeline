#!/usr/bin/env python
"""Extract an empirical PSF kernel from a grid-PSF reconstruction.

The sims thread reconstructs a sparse GRID of one element at one depth — a lone atom is
too under-constrained to reconstruct (it comes back ~84% noise), but a 4 A grid converges
like the real sample while each blob stays isolated in-plane and axially. Every blob IS
the system PSF. This pulls ONE clean interior blob onto a volume in the same complex
object format as NL70_new_vol.npy, so analysis/atomfind's empirical_psf() (np.angle ->
argmax -> crop) consumes it directly via cfg.single_atom_vol.

  python analysis/atomfind/extract_psf.py <recon_dir> <name> [--out ~/Desktop]

<recon_dir> holds .../Niter*.mat (newest used). Writes psf_<name>_vol.npy (+ _check.png):
one interior blob, entrance/exit surface layers dropped, atom forced to a POSITIVE phase
peak so empirical_psf()'s argmax locks onto it (not a negative-phase element or noise).
"""
import argparse, glob, os
import numpy as np
import h5py
from scipy.ndimage import maximum_filter


def load_vol(recon_dir):
    mats = glob.glob(os.path.join(recon_dir, "**", "Niter*.mat"), recursive=True)
    if not mats:
        raise SystemExit(f"no Niter*.mat under {recon_dir}")
    m = sorted(mats, key=lambda p: int("".join(filter(str.isdigit, os.path.basename(p)))))[-1]
    with h5py.File(m, "r") as f:
        L = []
        for r in f["outputs"]["object_roi"][:, 0]:
            a = f[r][:]
            a = (a["real"] + 1j * a["imag"]) if a.dtype.names else a
            L.append(a.T)
    return np.array(L), m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("recon_dir")
    ap.add_argument("name", help="element/geometry tag, e.g. Pb_NL70 -> psf_Pb_NL70_vol.npy")
    ap.add_argument("--out", default=os.path.expanduser("~/Desktop"))
    ap.add_argument("--zdrop", type=int, default=12, help="surface layers to drop each end")
    ap.add_argument("--half-xy", type=int, default=30, help="in-plane crop half-width (px)")
    a = ap.parse_args()

    cvol, m = load_vol(a.recon_dir)
    nL, Ny, Nx = cvol.shape
    ph = np.angle(cvol); ph -= np.median(ph, (1, 2), keepdims=True)

    # interior blob nearest the volume centre, from a mid-depth band (surfaces excluded)
    z0, z1 = a.zdrop, nL - a.zdrop
    band = np.abs(ph[z0:z1]).mean(0); s = band - band.min()
    pk = (s == maximum_filter(s, 9)) & (s > np.percentile(s, 99))
    ys, xs = np.where(pk)
    if len(ys) == 0:
        raise SystemExit("no blobs found in the mid-depth band")
    d = (ys - Ny / 2) ** 2 + (xs - Nx / 2) ** 2
    yy, xx = int(ys[np.argmin(d)]), int(xs[np.argmin(d)])
    W = a.half_xy
    yy = int(np.clip(yy, W, Ny - W)); xx = int(np.clip(xx, W, Nx - W))
    crop = cvol[z0:z1, yy - W:yy + W, xx - W:xx + W].astype(np.complex64)

    # force the atom to a POSITIVE phase peak so empirical_psf()'s argmax finds it
    pc = np.angle(crop); pc -= np.median(pc, (1, 2), keepdims=True)
    if abs(pc.min()) > abs(pc.max()):
        crop = np.conj(crop)
        print("[sign] atom was a negative phase feature -> conjugated")

    os.makedirs(a.out, exist_ok=True)
    out = os.path.join(a.out, f"psf_{a.name}_vol.npy")
    np.save(out, crop)

    V = np.angle(crop); V -= np.median(V, (1, 2), keepdims=True)
    l, r, c = np.unravel_index(np.argmax(V), V.shape)
    print(f"src : {m}")
    print(f"wrote {out}  shape {crop.shape} {crop.dtype}")
    print(f"agent argmax at (z={l}, y={r}, x={c}); |phase|>0.5 fraction {(np.abs(V) > 0.5).mean():.2%} (small = clean single blob)")

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        p = np.abs(V)
        fig, ax = plt.subplots(1, 2, figsize=(9, 4.2))
        ax[0].imshow(p.max(0), cmap="inferno"); ax[0].set_title(f"psf_{a.name}: in-plane (max-proj)")
        ax[1].imshow(p[:, r, :], aspect="auto", cmap="inferno"); ax[1].set_title("axial (z vs x)"); ax[1].set_ylabel("z layer")
        fig.tight_layout(); fig.savefig(os.path.join(a.out, f"psf_{a.name}_check.png"), dpi=120)
        print(f"wrote psf_{a.name}_check.png")
    except Exception as e:
        print("(no figure:", e, ")")


if __name__ == "__main__":
    main()
