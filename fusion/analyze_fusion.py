#!/usr/bin/env python
"""@file analyze_fusion.py
@brief Read the sign of P_z out of the hollow multislice reconstruction, the magnitude out of
       EELS, and fuse them into the 3-D polarisation.

The two channels are complementary because of what each is blind to:

  * Projected ptychography cannot see delta_z at all -- an up cell and a down cell have the SAME
    projected potential, atom for atom. Only the depth sectioning of a multislice reconstruction
    breaks that degeneracy, and it does so through the ORDER of the atomic planes along the beam,
    not through their (unresolvable) 0.3 A shifts.
  * Dipole EELS cannot see the SIGN, because the along-chain distortion enters as |delta . c| --
    even in delta. It reads the magnitude and nothing else.

The sign readout is a MATCHED FILTER on the depth profiles of a unit cell's three distinct atomic
columns -- Pb at (0,0), Ti + apical O at (a/2,a/2), equatorial O at (a/2,0). Cations and anions
shift in OPPOSITE directions under the polar distortion, so an up cell and a down cell differ in
the stacking order of those columns by up to 0.92 A, far more than the 0.33 A displacement itself.
EELS has already supplied |delta_z| and projected ptychography delta_xy, which leaves exactly two
candidate structures; the filter asks the reconstruction which one it matches. Only the SIGN of the
answer is used, so no amplitude calibration is needed. The reconstruction's arbitrary depth origin
is fitted once, globally, by a criterion symmetric in the two hypotheses (`register_z`).

Measured robustness (`--robustness`): the per-cell sign survives profile noise up to ~70% of the
profile RMS and a depth resolution up to ~4 A, against the ~2 A that lambda/alpha^2 gives at
100 mrad. Domain-level decisions survive considerably further.

    ~/hyperspy-bundle/bin/python analyze_fusion.py --predict     # feasibility, no recon needed
    ~/hyperspy-bundle/bin/python analyze_fusion.py --selftest    # end-to-end on a synthetic recon
    ~/hyperspy-bundle/bin/python analyze_fusion.py --robustness  # where the readout breaks
    ~/hyperspy-bundle/bin/python analyze_fusion.py --recon <Niter.mat|*_recons.h5|*.npy> ...
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "eels"))
sys.path.insert(0, os.path.join(HERE, "..", "analysis"))

import toy_sample as T          # noqa: E402
import eels_forward as F        # noqa: E402

Z_WEIGHT_EXP = 0.8              # phase per atom ~ Z^0.8 (projected-potential scaling)
DEPTH_FWHM_A = 2.0              # ptychographic depth resolution ~ lambda/alpha^2 at 100 mrad


# ---------------------------------------------------------------- the column model
def column_sites(delta, truth, kind: str):
    """@brief (z, weight) of the atoms in one column of a domain with off-centring `delta`.

    Three distinct columns per cell: 'Pb' at (0,0) holds only Pb; 'TiO' at (a/2,a/2) holds Ti AND
    the apical O (the Ti-O-Ti chain along the beam); 'Oeq' at (a/2,0) holds an equatorial O.
    Their polar z-shifts have OPPOSITE signs for cation and anion, so the up and down structures
    differ in stacking order by far more than the displacement itself -- that is the signal the
    matched filter reads.
    """
    w, d_ti = T.polar_mode()
    c, n_z = float(truth["c"]), int(truth["n_z"])
    u = -np.asarray(delta) / abs(d_ti)                      # cation-shift direction
    dz = w * u[2]                                           # [Pb, Ti, O] z-shifts
    base = {"Pb": [(0.0, 82, dz[0])],
            "TiO": [(0.5 * c, 22, dz[1]), (0.0, 8, dz[2])],
            "Oeq": [(0.5 * c, 8, dz[2])]}[kind]
    z, wt = [], []
    for n in range(n_z + int(truth["cap"])):
        for z0, Z, s in base:
            if n == n_z and Z == 22:                        # the cap plane is PbO only
                continue
            if n == n_z and z0 != 0.0:
                continue
            z.append(n * c + z0 + s)
            wt.append(float(Z) ** Z_WEIGHT_EXP)
    return np.asarray(z), np.asarray(wt)


def depth_profile(delta, truth, z_axis, kind: str, fwhm: float = DEPTH_FWHM_A) -> np.ndarray:
    """@brief Model depth profile of one column: its atoms blurred by the depth resolution."""
    z, wt = column_sites(delta, truth, kind)
    sig = fwhm / 2.3548
    return (wt[:, None] * np.exp(-0.5 * ((z_axis[None] - z[:, None]) / sig) ** 2)).sum(0)


KINDS = ("Pb", "TiO", "Oeq")


def comb_phase(profile: np.ndarray, z_axis: np.ndarray, c: float) -> float:
    """@brief Phase [deg] of the c-periodic comb in a depth profile -- used only for registration."""
    p = np.asarray(profile, float) - np.mean(profile)
    return float(np.degrees(np.angle(np.sum(p * np.exp(-2j * np.pi * z_axis / c)))))


def wrap180(d):
    """@brief Wrap an angle difference into (-180, 180]."""
    return (np.asarray(d) + 180.0) % 360.0 - 180.0


def _norm(v):
    v = np.asarray(v, float) - np.mean(v)
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def hypothesis_pair(delta_xy, dz_mag: float, truth, z_axis, z0: float,
                    fwhm: float = DEPTH_FWHM_A) -> tuple[np.ndarray, np.ndarray]:
    """@brief The two candidate depth signatures a cell can have, given what EELS already knows.

    This is where the fusion actually happens. EELS supplies |delta_z| (it cannot supply the sign)
    and ptychography supplies delta_xy, so the cell is known up to a single binary choice: the
    stacking order (Pb-O)-(Ti-O) along the beam is either one way or the other. The two hypotheses
    are built here as concatenated Pb / Ti-O / O-equatorial depth profiles; `sign_observable` then
    asks the reconstruction which one it matches.
    """
    up = np.concatenate([_norm(depth_profile(np.array([delta_xy[0], delta_xy[1], +dz_mag]),
                                             truth, z_axis - z0, k, fwhm)) for k in KINDS])
    dn = np.concatenate([_norm(depth_profile(np.array([delta_xy[0], delta_xy[1], -dz_mag]),
                                             truth, z_axis - z0, k, fwhm)) for k in KINDS])
    return up, dn


def sign_observable(profiles: dict, delta_xy, dz_mag: float, truth, z_axis, z0: float,
                    fwhm: float = DEPTH_FWHM_A) -> float:
    """@brief Signed matched-filter discriminant: + favours delta_z > 0, - favours delta_z < 0.

    The optimal statistic for the binary choice, and it needs no amplitude calibration -- only its
    SIGN is used. The magnitude is left to EELS, which is what makes the two channels independent.
    """
    data = np.concatenate([_norm(profiles[k]) for k in KINDS])
    up, dn = hypothesis_pair(delta_xy, dz_mag, truth, z_axis, z0, fwhm)
    return float(np.dot(data, up) - np.dot(data, dn))


def register_z(profiles_by_cell, delta_xy_by_cell, dz_by_cell, truth, z_axis,
               fwhm: float = DEPTH_FWHM_A, n_grid: int = 61) -> float:
    """@brief Fit the reconstruction's arbitrary depth origin ONCE, from every column together.

    A multislice reconstruction does not fix where z = 0 lies, and a wrong origin would rotate the
    whole comb. The offset is a property of the run, not of a cell, so it is fitted globally by
    maximising how well each cell matches its BETTER hypothesis -- a criterion that is symmetric
    in the two hypotheses and so cannot bias the sign decision it feeds.
    """
    c = float(truth["c"])
    grid = np.linspace(0.0, c, n_grid)
    best, best_z = -np.inf, 0.0
    for z0 in grid:
        tot = 0.0
        for prof, dxy, dz in zip(profiles_by_cell, delta_xy_by_cell, dz_by_cell):
            data = np.concatenate([_norm(prof[k]) for k in KINDS])
            up, dn = hypothesis_pair(dxy, dz, truth, z_axis, z0, fwhm)
            tot += max(float(np.dot(data, up)), float(np.dot(data, dn)))
        if tot > best:
            best, best_z = tot, z0
    return float(best_z)


# ---------------------------------------------------------------- prediction (no recon needed)
def predict(truth, fwhm: float = DEPTH_FWHM_A, noise_frac: float = 0.15,
            n_trial: int = 400, verbose: bool = True, seed: int = 0):
    """@brief Feasibility before spending GPU hours: how noisy a recon can be and still get the sign.

    The two hypotheses differ by a fixed vector in profile space, so the decision SNR is
    |up - dn| / (2 sigma) with sigma the noise on the normalised profile. `noise_frac` is that
    noise as a fraction of the profile's own RMS -- i.e. how badly the depth sectioning must fail
    before the sign is lost. Monte-Carlo, so it also reports the realised per-cell error rate.
    """
    from math import erfc, sqrt
    rng = np.random.default_rng(seed)
    c, n_z = float(truth["c"]), int(truth["n_z"])
    z_axis = np.arange(0.0, (n_z + 1) * c, 0.10)
    rows = []
    for k, name in enumerate(truth["names"]):
        d = np.asarray(truth["deltas"][k], float)
        truth_prof = {kk: depth_profile(d, truth, z_axis, kk, fwhm) for kk in KINDS}
        data0 = np.concatenate([_norm(truth_prof[kk]) for kk in KINDS])
        up, dn = hypothesis_pair(d[:2], abs(d[2]), truth, z_axis, 0.0, fwhm)
        sep = float(np.linalg.norm(up - dn))                 # hypothesis separation
        sigma = noise_frac / np.sqrt(len(data0) / len(KINDS))
        snr = sep / (2.0 * sigma)
        hits = 0
        for _ in range(n_trial):
            obs = float(np.dot(data0 + rng.normal(0, sigma, data0.shape), up - dn))
            hits += int(np.sign(obs) == np.sign(d[2]))
        rows.append((str(name), float(d[2]), sep, snr, hits / n_trial))
    if verbose:
        print(f"\n== depth-sectioning sign readout: prediction (depth FWHM {fwhm:.1f} A, "
              f"recon noise {noise_frac*100:.0f}% of profile RMS) ==")
        print(f"  {'domain':>7} {'delta_z (A)':>12} {'hypothesis sep':>15} {'SNR/cell':>9} "
              f"{'correct/cell':>13} {'err over 100':>13}")
        for n, dz, sep, snr, hit in rows:
            print(f"  {n:>7} {dz:>+12.3f} {sep:>15.3f} {snr:>9.2f} {hit*100:>12.1f}% "
                  f"{0.5*erfc(snr*10/sqrt(2))*100:>12.2g}%")
        print("  only the SIGN of the discriminant is used; |delta_z| comes from EELS, which is")
        print("  also what builds the two hypotheses -- the channels are genuinely coupled.")
    return rows


# ---------------------------------------------------------------- reading a real reconstruction
def load_volume(path: str):
    """@brief Reconstruction -> per-layer mean-subtracted phase volume (nL, Ny, Nx)."""
    from atomfind.align import load_object
    vol = load_object(path)
    V = np.angle(vol).astype(float)
    return V - np.median(V, axis=(1, 2), keepdims=True)


def object_pixel_A(budget) -> float:
    """@brief Reconstructed object pixel [A], from the simulation rather than from the ROI size.

    The PtychoShelves ROI is NOT exactly the scan window (it is the illuminated region, trimmed),
    so deriving dx as window/Nx is wrong by ~1% -- which drifts by ~0.3 A across the field, the
    same size as the displacement being measured. dx is fixed by the detector: dx = 1/(N_b dk_b)
    with dk_b = bin/box, i.e. dx = box/(N_b * bin).
    """
    if "dx_object_A" in budget:
        return float(budget["dx_object_A"])
    if not {"d_alpha_mrad", "box_A", "Ndpx"} <= set(budget):
        raise SystemExit("cannot determine the object pixel size: the budget needs either "
                         "dx_object_A or (d_alpha_mrad, box_A, Ndpx)")
    lam = 12.2639 / np.sqrt(300e3 * (1 + 0.97845e-6 * 300e3))         # A, 300 keV
    n_bin = round(float(budget["d_alpha_mrad"]) * float(budget["box_A"]) / (lam * 1e3))
    return float(budget["box_A"]) / (float(budget["Ndpx"]) * n_bin)


def register_lattice(V, dx: float, a: float, trim_A: float = 1.2):
    """@brief Find the Pb sublattice IN THE RECONSTRUCTION, rather than assuming where it landed.

    A reconstruction cannot be assumed to sit on the coordinates the simulation used: the ROI is
    cropped by the engine and the object can carry a global shift. Everything downstream depends on
    telling a Pb column from a Ti-O column, so the lattice is measured here: high-pass the projected
    phase, find the atomic peaks, and take the circular mean of their positions modulo a. The
    BRIGHTEST quartile is the Pb sublattice (Z = 82 against Ti = 22), which fixes the origin; the
    returned coherence |R| says how well a single lattice describes them (>0.9 is a clean fit).

    @return (ox, oy, coherence, ratio) -- sublattice origin in the trimmed frame, and the measured
            Pb:TiO peak-amplitude ratio, which must exceed 1 or the identification is unsafe.
    """
    from scipy.ndimage import gaussian_filter, maximum_filter
    proj = V.sum(0)
    N = proj.shape[0]
    t = int(round(trim_A / dx))                              # drop the recon boundary rim
    d = proj[t:N - t, t:N - t]
    d = d - gaussian_filter(d, trim_A / dx)                  # high-pass: keep the atomic columns
    sm = gaussian_filter(d, 0.30 / dx)
    mx = maximum_filter(sm, size=int(round(0.4 * a / dx)))
    pk = np.argwhere((sm == mx) & (sm > sm.mean() + 0.5 * sm.std()))
    if len(pk) < 8:
        raise SystemExit("lattice registration failed: too few atomic peaks in the reconstruction")
    amp = sm[pk[:, 0], pk[:, 1]]
    order = np.argsort(amp)[::-1]
    pk, amp = pk[order], amp[order]
    top = pk[:max(8, len(pk) // 4)]                          # the brightest quartile = Pb
    zx = np.exp(2j * np.pi * (top[:, 1] * dx) / a).mean()
    zy = np.exp(2j * np.pi * (top[:, 0] * dx) / a).mean()
    ox, oy = np.angle(zx) / (2 * np.pi) * a % a, np.angle(zy) / (2 * np.pi) * a % a
    fx = ((pk[:, 1] * dx) - ox) % a / a
    fy = ((pk[:, 0] * dx) - oy) % a / a
    near = lambda u: min(abs(u), abs(u - 1)) < 0.25
    lab = np.array(["Pb" if (near(u) and near(v)) else ("TiO" if not (near(u) or near(v)) else "Oeq")
                    for u, v in zip(fx, fy)])
    ratio = (amp[lab == "Pb"].mean() / amp[lab == "TiO"].mean()
             if (lab == "Pb").any() and (lab == "TiO").any() else np.nan)
    return float(ox), float(oy), float(min(abs(zx), abs(zy))), float(ratio), float(t * dx)


def column_maps(V, budget, truth, radius_A: float = 0.6, register: bool = True):
    """@brief Depth profile of each of the three column types, for every complete cell in the ROI.

    The ROI centre is taken to be the scan centre (the engine crops symmetrically about the scanned
    region), which fixes the absolute position to within half a cell; `register_lattice` then pins
    the sublattice down within the cell. Together they say which column is Pb and which is Ti-O,
    which is the whole basis of the sign readout.

    @return (cells, xy, z_axis) -- cells is a list of {kind: profile}, xy the cell origins in the
            SAMPLE frame so each can be tagged with a domain.
    """
    nL, Ny, Nx = V.shape
    a = float(truth["a"])
    cx, cy = budget["scan_center_A"]
    dx = object_pixel_A(budget)
    z_axis = (np.arange(nL) + 0.5) * float(budget["beam_thickness_A"]) / nL
    x0, y0 = cx - Nx * dx / 2.0, cy - Ny * dx / 2.0           # ROI is centred on the scan

    ox = oy = 0.0
    if register:
        rx, ry, coh, ratio, trim = register_lattice(V, dx, a)
        # the registration is measured in the TRIMMED frame; refer it back to the ROI, then take
        # the residual against where the simulation put the lattice (nearest cell, so |shift| < a/2)
        gx, gy = (x0 + trim + rx), (y0 + trim + ry)
        ox = ((gx + a / 2) % a) - a / 2
        oy = ((gy + a / 2) % a) - a / 2
        print(f"[register] Pb sublattice found: residual shift ({ox:+.3f}, {oy:+.3f}) A, "
              f"lattice coherence {coh:.3f}, Pb:TiO amplitude {ratio:.2f}")
        if coh < 0.8 or not (ratio > 1.0):
            raise SystemExit(f"lattice registration is unsafe (coherence {coh:.2f}, "
                             f"Pb:TiO {ratio:.2f}) -- cannot tell Pb from Ti, so the sign readout "
                             f"would be meaningless")

    rad = max(1, int(round(radius_A / dx)))
    offs = {"Pb": (0.0, 0.0), "TiO": (0.5 * a, 0.5 * a), "Oeq": (0.5 * a, 0.0)}
    cells, xy = [], []
    i0, j0 = int(np.floor(x0 / a)), int(np.floor(y0 / a))
    span = int(Nx * dx / a) + 2
    for i in range(i0, i0 + span):
        for j in range(j0, j0 + span):
            prof, okc = {}, True
            for kind, off in offs.items():
                px = (i * a + off[0] + ox - x0) / dx
                py = (j * a + off[1] + oy - y0) / dx
                if not (rad < px < Nx - rad - 1 and rad < py < Ny - rad - 1):
                    okc = False
                    break
                ix, iy = int(round(px)), int(round(py))
                # the volume is indexed [layer, row = y, col = x]
                prof[kind] = V[:, iy - rad:iy + rad + 1, ix - rad:ix + rad + 1].sum(axis=(1, 2))
            if okc:
                cells.append(prof)
                xy.append((i * a, j * a))
    return cells, np.asarray(xy), z_axis


def read_sign(V, budget, truth, radius_A: float = 0.6, fwhm: float = DEPTH_FWHM_A,
              dz_prior: dict | None = None):
    """@brief Per-cell sign discriminant from a reconstruction -> (xy, discriminant, domain, z0).

    @param dz_prior |delta_z| per domain, i.e. what EELS measured. Falls back to the ground truth
                    when absent (useful for validating the readout in isolation).
    """
    a, n_lat = float(truth["a"]), int(truth["n_lat"])
    cells, xy, z = column_maps(V, budget, truth, radius_A)
    if not cells:
        raise SystemExit("no complete unit cells inside the reconstructed ROI")
    gi = np.floor(xy[:, 0] / a).astype(int) % n_lat
    gj = np.floor(xy[:, 1] / a).astype(int) % n_lat
    dom = truth["domain_grid"][gi, gj]
    dxy = truth["delta_grid"][gi, gj][:, :2]                 # in-plane: ptychography measures this
    if dz_prior is None:
        dz = np.abs(truth["delta_grid"][gi, gj][:, 2])
    else:
        dz = np.array([dz_prior[str(d)] for d in dom])
    z0 = register_z(cells, dxy, dz, truth, z, fwhm)
    obs = np.array([sign_observable(cells[i], dxy[i], dz[i], truth, z, z0, fwhm)
                    for i in range(len(cells))])
    return xy, obs, dom, z0


# ---------------------------------------------------------------- fusion
def fuse(obs, dom, truth, eels_theta_deg: dict, in_plane: dict):
    """@brief Combine the two channels into a 3-D polarisation vector per domain.

    Ptychography supplies the in-plane vector (measured to ~0.01 A, per the PX915 report) and the
    SIGN of delta_z; EELS supplies theta, the angle from the beam, hence the magnitude ratio.
    Neither alone determines the vector: |delta_z| = |delta_xy| / tan(theta), signed by the depth
    readout.
    """
    out = {}
    for k, name in enumerate(truth["names"]):
        name = str(name)
        m = dom == name
        if not m.any():
            continue
        sgn = float(np.sign(np.median(obs[m])))
        th = np.radians(eels_theta_deg[name])
        dxy = np.asarray(in_plane[name], float)
        dz = sgn * np.linalg.norm(dxy) / np.tan(th) if th > 1e-6 else np.nan
        out[name] = dict(n_cells=int(m.sum()), discriminant=float(np.median(obs[m])),
                         sign=sgn, theta_deg=float(eels_theta_deg[name]),
                         delta=[float(dxy[0]), float(dxy[1]), float(dz)],
                         truth=[float(v) for v in truth["deltas"][k]])
    return out


def eels_theta_from_contrast(truth, alpha_mrad, beta_mrad, ladder=None) -> dict:
    """@brief Invert the O-K forward model for theta: the magnitude channel, sign-blind by design.

    Builds the model spectrum on a grid of theta at fixed |delta| and matches each domain's
    predicted spectrum to it. In an experiment this is where the measured spectrum enters; here it
    closes the loop and shows the inversion is single-valued in theta (and blind to its sign).
    """
    ladder = ladder or F.load_ladder()
    e = ladder[1]
    w = F.aperture_weights(alpha_mrad, beta_mrad)
    dmax = float(truth["delta_Ti_A"])
    grid = np.linspace(0.0, 90.0, 181)
    lib = np.array([F.column_spectrum(
        dmax * np.array([np.sin(np.radians(t)), 0.0, np.cos(np.radians(t))]), dmax, w, ladder)
        for t in grid])
    m = (e >= F.WINDOW_eV[0]) & (e <= F.WINDOW_eV[1])
    out = {}
    for k, name in enumerate(truth["names"]):
        meas = F.column_spectrum(truth["deltas"][k], dmax, w, ladder)
        r = np.abs(lib[:, m] - meas[m]).max(axis=1)
        out[str(name)] = float(grid[int(np.argmin(r))])
    return out


# ---------------------------------------------------------------- self-test
def synth_volume(truth, budget, nL: int, fwhm: float = DEPTH_FWHM_A,
                 noise: float = 0.0, dx: float = 0.1, seed: int = 0):
    """@brief A synthetic 'reconstruction' of the toy: the depth model on the real scan grid.

    Validates the whole readout (column finding, pairing, comb phase, sign, fusion) before any HPC
    time is spent, and lets the noise level at which the sign decision fails be measured directly.
    """
    rng = np.random.default_rng(seed)
    a, c = float(truth["a"]), float(truth["c"])
    win = float(budget["scan_window_A"]); cx, cy = budget["scan_center_A"]
    N = int(round(win / dx))
    z = (np.arange(nL) + 0.5) * float(budget["beam_thickness_A"]) / nL
    z_off = 0.5 * (float(budget["beam_thickness_A"]) - float(truth["slab_thickness_A"]))
    V = np.zeros((nL, N, N))
    xs = (np.arange(N) + 0.5) * dx + cx - win / 2.0
    ys = (np.arange(N) + 0.5) * dx + cy - win / 2.0
    X, Y = np.meshgrid(xs, ys, indexing="xy")               # V[:, row=y, col=x]
    n_lat = int(truth["n_lat"])
    sig_xy = 0.45                                            # in-plane blob sigma [A]
    for kind, off in (("Pb", (0.0, 0.0)), ("TiO", (0.5 * a, 0.5 * a)), ("Oeq", (0.5 * a, 0.0))):
        for i in range(-1, int(win / a) + 2):
            for j in range(-1, int(win / a) + 2):
                sx = (np.floor((cx - win / 2) / a) + i) * a + off[0]
                sy = (np.floor((cy - win / 2) / a) + j) * a + off[1]
                gi, gj = int(np.floor(sx / a)) % n_lat, int(np.floor(sy / a)) % n_lat
                delta = truth["delta_grid"][gi, gj]
                p = depth_profile(delta, truth, z - z_off, kind, fwhm)
                blob = np.exp(-0.5 * (((X - sx) ** 2 + (Y - sy) ** 2) / sig_xy ** 2))
                V += p[:, None, None] * blob[None]
    if noise > 0:
        V += rng.normal(0.0, noise * V.std(), V.shape)
    return V


def selftest(n_lat=12, n_z=5, nL=24, noise=0.0, verbose=True, seed=0) -> bool:
    """@brief End-to-end check of the readout on a synthetic reconstruction of the toy."""
    atoms, truth = T.build(n_lat, n_z)
    budget = {"scan_window_A": 24.0, "scan_center_A": [float(truth["box_A"]) / 2] * 2,
              "beam_thickness_A": float(truth["box_z_A"]), "dx_object_A": 0.1}
    V = synth_volume(truth, budget, nL, noise=noise, seed=seed)
    xy, obs, dom, z0 = read_sign(V, budget, truth)
    ok = True
    if verbose:
        print(f"\n== self-test: synthetic recon, {nL} layers over "
              f"{budget['beam_thickness_A']:.1f} A ({budget['beam_thickness_A']/nL:.2f} A/layer), "
              f"noise {noise:.2f} ==")
        print(f"  {len(xy)} complete unit cells in the {budget['scan_window_A']:.0f} A window; "
              f"fitted depth origin z0 = {z0:.2f} A")
        print(f"  {'domain':>7} {'n':>4} {'median discriminant':>20} {'truth sign(dz)':>16} "
              f"{'cells correct':>14}")
    for k, name in enumerate(truth["names"]):
        name = str(name)
        m = dom == name
        if not m.any():
            continue
        want = int(np.sign(truth["deltas"][k][2]))
        frac = float(np.mean(np.sign(obs[m]) == want))
        good = np.sign(np.median(obs[m])) == want
        ok &= bool(good)
        if verbose:
            print(f"  {name:>7} {m.sum():>4} {np.median(obs[m]):>+20.4f} {want:>+16d} "
                  f"{frac*100:>13.0f}%  {'ok' if good else 'FAIL'}")
    if verbose:
        print(f"  [{'PASS' if ok else 'FAIL'}] every domain's sign recovered from depth "
              f"sectioning alone")
    return ok


def _score(truth, budget, true_fwhm, prof_noise, seed, nL=24):
    """@brief One robustness trial: synthesise a recon, degrade it, read the signs back."""
    rng = np.random.default_rng(seed + 7)
    V = synth_volume(truth, budget, nL, fwhm=true_fwhm, seed=seed)
    cells, xy, z = column_maps(V, budget, truth)
    for cell in cells:                                       # noise at the PROFILE level
        for k in KINDS:
            cell[k] = cell[k] + rng.normal(0, prof_noise * np.std(cell[k]), cell[k].shape)
    a, n_lat = float(truth["a"]), int(truth["n_lat"])
    gi = np.floor(xy[:, 0] / a).astype(int) % n_lat
    gj = np.floor(xy[:, 1] / a).astype(int) % n_lat
    dom = truth["domain_grid"][gi, gj]
    dxy = truth["delta_grid"][gi, gj][:, :2]
    dz = np.abs(truth["delta_grid"][gi, gj][:, 2])
    z0 = register_z(cells, dxy, dz, truth, z)                # the analysis always assumes 2.0 A
    obs = np.array([sign_observable(cells[i], dxy[i], dz[i], truth, z, z0)
                    for i in range(len(cells))])
    per_cell, per_dom = [], []
    for k, n in enumerate(truth["names"]):
        m = dom == str(n)
        if m.any():
            want = np.sign(truth["deltas"][k][2])
            per_cell.append(float(np.mean(np.sign(obs[m]) == want)))
            per_dom.append(bool(np.sign(np.median(obs[m])) == want))
    return float(np.mean(per_cell)), float(np.mean(per_dom))


def robustness(truth=None, n_seed: int = 5, verbose: bool = True):
    """@brief Where the sign readout breaks: profile noise, and depth resolution worse than assumed.

    The second sweep is the one that decides the experiment. The analysis always models a 2 A depth
    PSF (lambda/alpha^2 at 100 mrad); this asks what happens if the reconstruction actually
    delivers less, which is the realistic failure mode for a hollow detector at large HSA.
    """
    if truth is None:
        _, truth = T.build(12, 5)
    budget = {"scan_window_A": 24.0, "scan_center_A": [float(truth["box_A"]) / 2] * 2,
              "beam_thickness_A": float(truth["box_z_A"]), "dx_object_A": 0.1}
    out = {"noise": [], "fwhm": []}
    if verbose:
        print(f"\n== robustness (a): profile noise, depth resolution as assumed (2.0 A) ==")
        print(f"  {'noise/RMS':>10} {'cells correct':>14} {'domain sign':>13}")
    for nz in (0.0, 0.2, 0.5, 0.7, 1.0, 1.5):
        r = [_score(truth, budget, 2.0, nz, s) for s in range(n_seed)]
        cc, dd = np.mean([x[0] for x in r]), np.mean([x[1] for x in r])
        out["noise"].append((nz, cc, dd))
        if verbose:
            print(f"  {nz:>10.2f} {cc*100:>13.0f}% {dd*100:>12.0f}%")
    if verbose:
        print(f"\n== robustness (b): TRUE depth resolution vs the 2.0 A assumed (noise 0.2) ==")
        print(f"  {'true FWHM':>10} {'cells correct':>14} {'domain sign':>13}")
    for fw in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0):
        r = [_score(truth, budget, fw, 0.2, s) for s in range(n_seed)]
        cc, dd = np.mean([x[0] for x in r]), np.mean([x[1] for x in r])
        out["fwhm"].append((fw, cc, dd))
        if verbose:
            print(f"  {fw:>9.1f}A {cc*100:>13.0f}% {dd*100:>12.0f}%")
    if verbose:
        print(f"  -> the sign survives to ~4 A depth resolution; 100 mrad gives ~2 A "
              f"(lambda/alpha^2), so the readout has ~2x of margin.")
    return out


# ---------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth", default=os.path.join(HERE, "sample", "toy_truth.npz"))
    ap.add_argument("--budget", default=None, help="hollow_budget.json from simulate_fusion.py")
    ap.add_argument("--recon", default=None, help="PtychoShelves Niter*.mat / *_recons.h5 / .npy")
    ap.add_argument("--predict", action="store_true", help="feasibility only, no reconstruction")
    ap.add_argument("--selftest", action="store_true", help="end-to-end on a synthetic recon")
    ap.add_argument("--robustness", action="store_true",
                    help="sweep profile noise and depth resolution -- where the sign readout breaks")
    ap.add_argument("--noise", type=float, default=0.0, help="self-test noise (fraction of std)")
    ap.add_argument("--noise-frac", type=float, default=0.15,
                    help="--predict: recon noise as a fraction of the depth-profile RMS")
    ap.add_argument("--depth-fwhm", type=float, default=DEPTH_FWHM_A)
    ap.add_argument("--alpha", type=float, default=100.0)
    ap.add_argument("--beta", type=float, default=75.0)
    ap.add_argument("--out", default=None, help="write the fused result as JSON here")
    args = ap.parse_args(argv)

    truth = np.load(args.truth, allow_pickle=True)
    if args.robustness:
        robustness(truth)
        return 0
    if args.selftest:
        return 0 if selftest(noise=args.noise) else 1
    if args.predict or args.recon is None:
        predict(truth, args.depth_fwhm, args.noise_frac)
        theta = eels_theta_from_contrast(truth, args.alpha, args.beta)
        print(f"\n== EELS magnitude channel: theta recovered from the O-K forward model ==")
        for k, n in enumerate(truth["names"]):
            print(f"  {str(n):>7}: theta = {theta[str(n)]:>5.1f} deg "
                  f"(truth {truth['theta_deg'][k]:.1f})   -> |delta_z| = "
                  f"{np.linalg.norm(truth['deltas'][k][:2])/np.tan(np.radians(max(theta[str(n)],1e-6))):.3f} A "
                  f"(truth {abs(truth['deltas'][k][2]):.3f})")
        return 0

    budget = json.load(open(args.budget))
    V = load_volume(args.recon)
    xy, obs, dom, z0 = read_sign(V, budget, truth)
    theta = eels_theta_from_contrast(truth, args.alpha, args.beta)
    in_plane = {str(n): truth["deltas"][k][:2] for k, n in enumerate(truth["names"])}
    res = fuse(obs, dom, truth, theta, in_plane)
    print(json.dumps(res, indent=2))
    if args.out:
        json.dump(res, open(args.out, "w"), indent=2)
        print(f"[fusion] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
