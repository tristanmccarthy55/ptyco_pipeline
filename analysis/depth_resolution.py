#!/usr/bin/env python
"""Depth-resolution analysis for the multislice synthetic reconstructions.

Loads PtychoShelves Niter*.mat outputs (object_roi cell array of complex layers)
and asks the central question: do we resolve atomic planes ALONG THE BEAM?

The rigorous test is the kz power spectrum of the phase along each atom column:
if individual planes are resolved, the column phase oscillates at the atom-plane
frequency (Pb-Pb ~4 A -> 0.25 A^-1; AB sublattice ~2 A -> 0.50 A^-1). We compare
on-column pixels to vacuum pixels (baseline) for both the 42- and 70-layer runs,
and mark the true GT plane frequencies.

Reusable: re-point DATA_ROOT at any future pull of recon_NL*/.../Niter*.mat.
Usage:  python depth_resolution.py
Outputs: figures + volumes under OUT_DIR (~/Desktop).
"""
import os, glob, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import h5py

DATA_ROOT = os.path.expanduser("~/Desktop/recon_new")
OUT_DIR   = os.path.expanduser("~/Desktop")
VASP      = "/Users/u2109287/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PtycoShelves/PTO6_STO6_18_18_labyrinthPoscar.vasp"
SCAN_CX, SCAN_CY, SCAN_W = 40.0, 20.0, 20.0   # sim scan window (A)

# ---------------------------------------------------------------- loader
def load_recon(run):
    mats = glob.glob(os.path.join(DATA_ROOT, run, "**", "Niter*.mat"), recursive=True)
    m = sorted(mats, key=lambda p: int("".join(filter(str.isdigit, os.path.basename(p)))))[-1]
    with h5py.File(m, "r") as f:
        g = f["outputs"]
        refs = g["object_roi"][:, 0]
        layers = []
        for r in refs:
            a = f[r][:]
            if a.dtype.names: a = a["real"] + 1j * a["imag"]
            layers.append(a.T)                      # MATLAB col-major -> [Ny,Nx]
        vol = np.array(layers)                      # [nL,Ny,Nx] complex
        dx  = float(g["pixel_size"][0, 0]) * 1e10   # A
        z   = g["z_distance"][:, 0]; z = z[np.isfinite(z)]
        dz  = float(np.median(z)) * 1e10            # A per layer
        ferr = np.array(f["fourier_error_out"]).ravel()
    return vol, dx, dz, ferr, os.path.basename(m)

# ---------------------------------------------------------------- GT depth structure
def gt_plane_freqs():
    import ase.io, abtem
    a = ase.io.read(VASP); a.rotate(-90, "y", rotate_cell=True); a = abtem.orthogonalize_cell(a)
    Lx, Ly, Lz = a.cell.lengths(); s = max(Lx, Ly); a.cell[0, 0] = s; a.cell[1, 1] = s
    a.center(axis=0); a.center(axis=1); a.center(axis=2, vacuum=2.0)
    pos = a.get_positions(); Z = a.get_atomic_numbers()
    win = (np.abs(pos[:, 0] - SCAN_CX) < SCAN_W/2) & (np.abs(pos[:, 1] - SCAN_CY) < SCAN_W/2)
    out = {}
    for nm, zz in [("Pb", 82), ("Ti", 22)]:
        zc = pos[win & (Z == zz), 2]
        # dominant plane spacing from nearest-neighbour z gaps
        zs = np.sort(np.unique(np.round(zc, 1)))
        gaps = np.diff(zs); gap = np.median(gaps[gaps > 0.5]) if len(gaps) else np.nan
        out[nm] = gap
    return out, pos[win], Z[win]

# ---------------------------------------------------------------- kz spectrum
def kz_spectrum(vol, dz):
    V = np.angle(vol).astype(float)
    V -= np.median(V, (1, 2), keepdims=True)        # per-layer offset
    dm = V.mean(0); dmn = dm - dm.min()
    col = dmn > np.percentile(dmn, 95)              # on-column pixels
    vac = dmn < np.percentile(dmn, 35)              # vacuum pixels
    nL = V.shape[0]
    win = np.hanning(nL)[:, None, None]
    Vc = (V - V.mean(0, keepdims=True)) * win        # remove DC, window along z
    P = np.abs(np.fft.rfft(Vc, axis=0))**2
    kz = np.fft.rfftfreq(nL, d=dz)                   # A^-1
    Pcol = P[:, col].mean(1); Pvac = P[:, vac].mean(1)
    return kz, Pcol, Pvac

# ---------------------------------------------------------------- run
runs = ["NL42", "NL70"]
data = {r: load_recon(r) for r in runs}
gtgap, gtpos, gtZ = gt_plane_freqs()
f_pb = 1.0 / gtgap["Pb"]; f_ti = 1.0 / gtgap["Ti"]
print(f"GT plane spacing along beam:  Pb {gtgap['Pb']:.2f} A ({f_pb:.3f} A^-1)   Ti {gtgap['Ti']:.2f} A ({f_ti:.3f} A^-1)")
for r in runs:
    vol, dx, dz, ferr, fn = data[r]
    np.save(os.path.join(OUT_DIR, f"{r}_new_vol.npy"), vol.astype(np.complex64))
    print(f"{r}: {vol.shape} dz={dz:.3f} A  kz_Nyq={1/(2*dz):.3f} A^-1  ferr {ferr[0]:.2f}->{ferr[-1]:.2f}  [{fn}]")

# ---- Figure 1: convergence (error logged once per check; drop zero padding)
fig, ax = plt.subplots(figsize=(7, 4.5))
for r in runs:
    ferr = data[r][3]; ferr = ferr[np.isfinite(ferr) & (ferr > 0)]
    ax.plot(np.arange(1, len(ferr)+1), ferr, marker="o", ms=3, label=f"{r} (dz={data[r][2]:.2f} A, end={ferr[-1]:.1f})")
ax.set_xlabel("logged check (~ iteration)"); ax.set_ylabel("Fourier error"); ax.set_yscale("log")
ax.legend(); ax.set_title("Convergence — reg off, single mode, fixed probe")
fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, "depth_convergence.png"), dpi=130)

# ---- Figure 2: kz power spectrum (THE depth-resolution test)
fig, ax = plt.subplots(1, 2, figsize=(14, 5.2), sharey=True)
for i, r in enumerate(runs):
    vol, dx, dz, ferr, fn = data[r]
    kz, Pcol, Pvac = kz_spectrum(vol, dz)
    a = ax[i]
    a.semilogy(kz, Pcol, lw=2, label="on-column")
    a.semilogy(kz, Pvac, lw=1.2, color="gray", label="vacuum (baseline)")
    for f0, nm, c in [(f_pb, "Pb planes 4A", "tab:red"), (f_ti, "Ti planes", "tab:green")]:
        a.axvline(f0, color=c, ls="--", lw=1.2, label=f"{nm} {f0:.2f}/A")
    a.axvline(1/(2*dz), color="k", ls=":", lw=1, label=f"Nyquist {1/(2*dz):.2f}/A")
    a.set_xlabel("kz  (A$^{-1}$)"); a.set_title(f"{r}   dz={dz:.2f} A"); a.set_xlim(0, 0.6)
    if i == 0: a.set_ylabel("z-power (on column)")
    a.legend(fontsize=8)
fig.suptitle("Depth-resolution test: does the column phase oscillate at the atom-plane frequency?", fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, "depth_kz_spectrum.png"), dpi=130)

# ---- Figure 3: column line profiles (phase vs depth), brightest Pb columns
from scipy.ndimage import maximum_filter
fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
for i, r in enumerate(runs):
    vol, dx, dz, ferr, fn = data[r]
    V = np.angle(vol).astype(float); V -= np.median(V, (1, 2), keepdims=True)
    dm = V.mean(0); dmn = dm - dm.min()
    pk = (dmn == maximum_filter(dmn, 25)) & (dmn > np.percentile(dmn, 96))
    ys, xs = np.where(pk); zaxis = (np.arange(V.shape[0]) + 0.5) * dz
    order = np.argsort(-dmn[ys, xs])
    for j in order[:4]:
        prof = V[:, ys[j], xs[j]]; prof = prof - prof.mean()
        ax[i].plot(zaxis, prof, lw=1.4)
    for zp in np.arange(zaxis[0], zaxis[-1], gtgap["Pb"]):
        ax[i].axvline(zp, color="tab:red", ls=":", lw=0.6, alpha=0.5)
    ax[i].set_xlabel("depth z (A)"); ax[i].set_ylabel("phase (mean-sub)")
    ax[i].set_title(f"{r}: brightest Pb columns vs depth (red = GT 4A planes)")
fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, "depth_line_profiles.png"), dpi=130)

# ---- verdict numbers: peak-at-0.25 excess over vacuum
print("\n=== depth-resolution verdict ===")
for r in runs:
    vol, dx, dz, ferr, fn = data[r]
    kz, Pcol, Pvac = kz_spectrum(vol, dz)
    def band(f0):
        m = np.abs(kz - f0) < 0.03
        return Pcol[m].max() / (Pvac[m].max() + 1e-30)
    print(f"{r}:  on-column/vacuum power  @Pb(0.25)={band(f_pb):.1f}x   @Ti({f_ti:.2f})={band(f_ti):.1f}x"
          f"   (>~2x = real plane signal)")
print("\nwrote depth_convergence.png, depth_kz_spectrum.png, depth_line_profiles.png + NL*_new_vol.npy")
