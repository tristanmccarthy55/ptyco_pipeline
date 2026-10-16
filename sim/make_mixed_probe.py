#!/usr/bin/env python
"""
@file make_mixed_probe.py
@brief Write a MULTI-MODE start probe: the top eigenmodes of the partially coherent state the sim's own settings describe.

A probe-mode reconstruction seeded the default way starts every extra mode at mode_start_pow (2 %), and no mode run of
the campaign has moved mode 1 far from its seed (mono300/PHASE1.md). This writes instead the expected mixed state of a
real gun with the data's focal spread and source size (read from sim_meta.mat): a continuous Gaussian spread of defocus
(rms focal spread) and of probe position (source FWHM / 2.355 per axis), sampled on symmetric stratified quantiles (a
product grid for the source), applied to the TRUE probe (probe_initial_true.mat) as k-space phases. The mixture's
eigenmodes follow from the K x K Gram matrix (only the aperture's |F|^2 and the offsets enter it: the aberration phase is
common to every configuration and cancels), each mode = F(k) x sum_i c_i exp(-i theta_i(k)). Because the distributions are
symmetric and isotropic, the state does not depend on the array's orientation or on the sign conventions of defocus and
shift, so the engine's data flip cannot put it in the wrong frame.

Output: probe_initial.mat with probe of shape (N, N, 1, M) -- the live loader (core/ptycho_model_probe.m) keeps a 4-D
probe as modes (its 3-D reshape is commented out), and prepare_initial_probes then uses the file's modes instead of
seeding Hermite modes -- scaled so the M modes hold the true probe's total power with the state's relative occupancies.
probe_initial.json: occupancies (the M kept and the whole spectrum's head), the share of the state they hold, the share
of the sim's OWN 16-configuration state their span holds (both orientations of its offsets), mode 1's overlap with the
true probe (ov_start for analysis/kick_probe.py), and self-checks (powers from the modes vs the Gram eigenvalues, mutual
orthogonality). The writer exits non-zero if a self-check fails.

    python sim/make_mixed_probe.py --sim-dir <recon>/01 --modes 8 --out <recon>/01/probe_initial.mat

SEARCHED START (2026-10-16): --base <probe .mat> builds the modes on THAT coherent probe (e.g. an outer search's start,
written by sim/make_probe.py) instead of the true one -- what an experimenter would do with their measured probe and known
energy spread / source size. The truth is still read, for the record only: mode 1's overlap with it, and how much of the
sim's own 16-configuration state (built on the true probe) the modes' span holds.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat
from scipy.stats import norm

E0_EV = 510998.95
POWER_TOL = 2e-3          # |power from the modes / Gram eigenvalue - 1|, every kept mode
ORTHO_TOL = 2e-3          # max |<m|n>| / sqrt(p_m p_n), m != n
GRAM_CHUNK = 40000        # k-pixels per block of the Gram sum (memory: configurations x block complex numbers)


def wavelength_A(e_ev):
    return 12.2643 / np.sqrt(e_ev * (1 + 0.97845e-6 * e_ev))


def quantiles(n):
    q = norm.ppf((np.arange(n) + 0.5) / n)
    return q / q.std()                                   # symmetric, rms exactly 1 (as the sim's coherence_offsets)


def offsets(fs_A, src_fwhm_A, n_fs, n_src):
    """(defocus offsets [A], (K, 2) position offsets [A]) of the continuous state's equal-weight configurations."""
    ddf = quantiles(n_fs) * fs_A if fs_A > 0 else np.zeros(1)
    if src_fwhm_A > 0:
        g = quantiles(n_src) * src_fwhm_A / (2 * np.sqrt(2 * np.log(2)))
        X, Y = np.meshgrid(g, g, indexing="ij"); d = np.c_[X.ravel(), Y.ravel()]
    else:
        d = np.zeros((1, 2))
    return np.repeat(ddf, len(d)), np.tile(d, (len(ddf), 1))


def phases(ddf, dr, k0, k1, lam):
    """theta_i(k) for every configuration (rows) at the given k-pixels (columns)."""
    return (np.pi * lam * ddf[:, None] * (k0 ** 2 + k1 ** 2)[None, :]
            + 2 * np.pi * (dr[:, :1] * k0[None, :] + dr[:, 1:] * k1[None, :]))


def sim16_offsets(fs_A, src_fwhm_A, n, seed=7):
    """The sim's own coherence_offsets (sim/simulate_4dstem.py), reproduced: same draws, same rescaling."""
    rng = np.random.default_rng(seed)
    ddf = np.zeros(n); dr = np.zeros((n, 2))
    if fs_A > 0 and n > 1:
        ddf = rng.permutation(quantiles(n)) * fs_A
    if src_fwhm_A > 0 and n > 1:
        d = rng.standard_normal((n, 2)); d -= d.mean(0); d /= d.std(0); dr = d * src_fwhm_A / (2 * np.sqrt(2 * np.log(2)))
    return ddf, dr


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sim-dir", type=Path, required=True, help="01/ holding probe_initial_true.mat and sim_meta.mat")
    ap.add_argument("--modes", type=int, required=True, help="modes to write (the reconstruction's PROBE_MODES)")
    ap.add_argument("--out", type=Path, required=True, help="probe_initial.mat to write")
    ap.add_argument("--base", type=Path, default=None,
                    help="coherent probe to build the modes on (default: <sim-dir>/probe_initial_true.mat)")
    ap.add_argument("--focal-spread", type=float, default=None, help="override the sim's focal spread [A rms]")
    ap.add_argument("--source-fwhm", type=float, default=None, help="override the sim's source FWHM [A]")
    ap.add_argument("--n-fs", type=int, default=None, help="defocus quantiles (default 200 alone, 24 with a source)")
    ap.add_argument("--n-src", type=int, default=None, help="source quantiles per axis (default 20 alone, 8 with a focal spread)")
    a = ap.parse_args(argv)

    meta = loadmat(str(a.sim_dir / "sim_meta.mat"), simplify_cells=True)["meta"]
    fs = float(meta.get("focal_spread_A", 0.0)) if a.focal_spread is None else a.focal_spread
    src = float(meta.get("source_size_fwhm_A", 0.0)) if a.source_fwhm is None else a.source_fwhm
    lam = wavelength_A(float(meta["energy_kev"]) * 1e3)
    dx = float(meta["dx_object_A"])
    T = loadmat(str(a.sim_dir / "probe_initial_true.mat"))["probe"].astype(np.complex128)
    n = T.shape[0]
    assert T.shape == (n, n), f"true probe {T.shape}"
    B = T if a.base is None else np.squeeze(loadmat(str(a.base))["probe"]).astype(np.complex128)   # read before --out may overwrite it
    if B.shape != T.shape:
        raise SystemExit(f"[mixed] base probe {B.shape} != true probe {T.shape} (a multi-mode file is not a base)")
    if a.modes < 1:
        raise SystemExit("[mixed] --modes must be >= 1")
    if fs <= 0 and src <= 0:
        raise SystemExit("[mixed] the sim is fully coherent (no focal spread, no source size): nothing to mix")
    n_fs = a.n_fs or (200 if src <= 0 else 24)
    n_src = a.n_src or (20 if fs <= 0 else 8)
    ddf, dr = offsets(fs, src, n_fs, n_src)
    K = len(ddf); w = 1.0 / K

    F = np.fft.fft2(B)                                   # the modes are built on the base (the true probe by default)
    Ft = F if a.base is None else np.fft.fft2(T)         # the sim's own configurations are built on the truth
    kk = np.fft.fftfreq(n, d=dx)
    K0, K1 = np.meshgrid(kk, kk, indexing="ij")
    a2 = np.abs(F) ** 2
    sup = a2 > 1e-10 * a2.max()                           # the aperture (the probe's whole k-space support)
    # Gram matrix over the whole support, in blocks: G_ij = w sum_k |F|^2 exp(i(theta_i - theta_j)). (A strided subset
    # left the full-grid modes 0.5 % from orthogonal; the exact sum makes them orthogonal to rounding.)
    k0, k1, amp = K0[sup], K1[sup], np.sqrt(a2[sup] * w)
    G = np.zeros((K, K), np.complex128)
    for s0 in range(0, len(k0), GRAM_CHUNK):
        sl = slice(s0, s0 + GRAM_CHUNK)
        E = np.exp(-1j * phases(ddf, dr, k0[sl], k1[sl], lam)) * amp[sl][None, :]           # (K, block)
        G += E.conj() @ E.T
    lamb, V = np.linalg.eigh(G)
    order = np.argsort(lamb)[::-1]; lamb = np.clip(lamb[order].real, 0, None); V = V[:, order]
    occ = lamb / lamb.sum()
    M = min(a.modes, K)
    # the kept modes on the full support: F_m(k) = F(k) sum_i sqrt(w) V_im exp(-i theta_i(k))
    acc = np.zeros((M, int(sup.sum())), np.complex128)
    for i in range(K):
        ph = np.exp(-1j * phases(ddf[i:i + 1], dr[i:i + 1], k0, k1, lam))[0]
        acc += (np.sqrt(w) * V[i, :M])[:, None] * ph[None, :]
    modes = np.zeros((M, n, n), np.complex128)
    for m in range(M):
        Fm = np.zeros((n, n), np.complex128); Fm[sup] = F[sup] * acc[m]
        modes[m] = np.fft.ifft2(Fm)
    p = (np.abs(modes) ** 2).sum(axis=(1, 2))
    # self-checks: the modes' powers follow the Gram eigenvalues; the modes are mutually orthogonal
    p_rel = p / (np.linalg.norm(B) ** 2)                                # share of the base probe's power per mode
    power_err = float(np.max(np.abs(p_rel / occ[:M] - 1))) if M else 0.0
    U = modes.reshape(M, -1) / np.sqrt(p)[:, None]
    O = np.abs(U.conj() @ U.T); np.fill_diagonal(O, 0)
    ortho_err = float(O.max()) if M > 1 else 0.0
    # how much of the sim's OWN 16-configuration state the kept span holds (both orientations of its position offsets)
    cover = {}
    n16 = int(meta.get("coherence_samples", 16)) or 16
    for name, swap in (("as_written", False), ("transposed", True)):
        d16, r16 = sim16_offsets(fs, src, n16)
        r16 = r16[:, ::-1] if swap else r16
        tot = held = 0.0
        for i in range(n16):
            psi = np.zeros((n, n), np.complex128)
            psi[sup] = Ft[sup] * np.exp(-1j * phases(d16[i:i + 1], r16[i:i + 1], k0, k1, lam))[0]
            psi = np.fft.ifft2(psi)
            tot += np.linalg.norm(psi) ** 2
            held += float(np.sum(np.abs(U.conj() @ psi.ravel()) ** 2))
        cover[name] = held / tot
    # normalise: the M modes carry the base probe's total power, with the state's relative occupancies
    modes *= np.sqrt(np.linalg.norm(B) ** 2 / p.sum())
    P = np.moveaxis(modes, 0, -1)[:, :, None, :].astype(np.complex64)   # (N, N, 1, M)
    if a.out.is_symlink():
        raise SystemExit(f"[mixed] {a.out} is a symlink (a known-probe link?) -- refusing to write through it")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    savemat(str(a.out), {"probe": P, "p": {"binning": False, "detector": {"binning": False}}})

    m1 = modes[0]
    ov1 = float(np.abs(np.vdot(T.ravel(), m1.ravel())) / (np.linalg.norm(T) * np.linalg.norm(m1)))
    ok = power_err <= POWER_TOL and ortho_err <= ORTHO_TOL
    info = dict(kind="mixed-state start (top eigenmodes of a continuous Gaussian focal spread / source)",
                focal_spread_A=fs, source_fwhm_A=src, n_fs=n_fs if fs > 0 else 0, n_src=n_src if src > 0 else 0,
                configurations=K, modes=M, occupancy_kept=[float(x) for x in occ[:M]],
                occupancy_head=[float(x) for x in occ[:16]], state_held=float(occ[:M].sum()),
                sim16_state_held=cover, overlap_with_sim_probe=ov1, c1_A=None, c1_sim_A=None, scaled_terms={},
                self_check=dict(power_err=power_err, ortho_err=ortho_err, result="PASS" if ok else "FAIL"),
                n_b=n, dx_A=dx, wavelength_A=lam, sim_dir=str(a.sim_dir.resolve()),
                base="probe_initial_true.mat (the truth)" if a.base is None else str(a.base.resolve()),
                base_overlap_with_sim_probe=float(np.abs(np.vdot(T.ravel(), B.ravel())) / (np.linalg.norm(T) * np.linalg.norm(B))))
    js = a.out.with_suffix(".json")
    if js.is_symlink():
        js.unlink()
    js.write_text(json.dumps(info, indent=2))
    print(f"[mixed] focal spread {fs:g} A rms, source {src:g} A FWHM: {K} configurations, {M} modes written to {a.out}")
    print(f"[mixed] occupancies {' '.join(f'{x:.3f}' for x in occ[:M])}  (state held {occ[:M].sum():.3f}; "
          f"sim's own 16-config state held {cover['as_written']:.3f} / {cover['transposed']:.3f})")
    print(f"[mixed] mode 1 overlap with the true probe {ov1:.4f}; SELF-CHECK {info['self_check']['result']} "
          f"(power err {power_err:.1e}, orthogonality {ortho_err:.1e})")
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
