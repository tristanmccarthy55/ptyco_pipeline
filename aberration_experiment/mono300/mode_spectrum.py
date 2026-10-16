"""Exact mixed-state spectrum of the phase-1 probes: how many probe modes each focal spread / source size needs.

The sim's partially coherent probe is an equal-weight mixture of COHERENCE_SAMPLES coherent probes, each with its own
defocus offset and source (position) offset (sim/simulate_4dstem.coherence_offsets, seed 7). The common aberration
phase is a unitary diagonal in Fourier space and does not change the spectrum, so only the aperture, the wavelength and
the offsets matter. Occupancies = eigenvalues of the mixture, normalised. 'sim16' uses the sim's own 16 offsets;
'cont' approximates the continuous Gaussian state with 1000 samples (what a real gun gives).
"""
import numpy as np
from scipy.stats import norm

E = 300e3; m0c2 = 510998.95
LAM = 12.2643 / np.sqrt(E * (1 + 0.97845e-6 * E))      # A
ALPHA = 0.080
DK = 0.04                                               # 1/A
kmax = ALPHA / LAM
n = int(np.ceil(kmax / DK)) + 1
k1 = np.arange(-n, n + 1) * DK
KX, KY = np.meshgrid(k1, k1, indexing="ij")
K2 = KX ** 2 + KY ** 2
ap = K2 <= kmax ** 2
kx, ky, k2 = KX[ap], KY[ap], K2[ap]

def offsets_sim(fs, src, nconf=16, seed=7):
    rng = np.random.default_rng(seed)
    ddf = np.zeros(nconf); dr = np.zeros((nconf, 2))
    if fs > 0:
        q = norm.ppf((np.arange(nconf) + 0.5) / nconf); q = q / q.std(); ddf = rng.permutation(q) * fs
    if src > 0:
        s = src / (2 * np.sqrt(2 * np.log(2))); d = rng.standard_normal((nconf, 2)); d -= d.mean(0); d /= d.std(0); dr = d * s
    return ddf, dr

def offsets_cont(fs, src, nconf=1000, seed=1):
    rng = np.random.default_rng(seed)
    ddf = rng.standard_normal(nconf) * fs
    dr = rng.standard_normal((nconf, 2)) * src / (2 * np.sqrt(2 * np.log(2)))
    return ddf, dr

def spectrum(ddf, dr):
    ph = np.pi * LAM * ddf[:, None] * k2[None, :] + 2 * np.pi * (dr[:, :1] * kx[None, :] + dr[:, 1:] * ky[None, :])
    P = np.exp(-1j * ph) / np.sqrt(ap.sum())
    G = P @ P.conj().T / len(ddf)                          # Gram matrix: same nonzero spectrum as the mixture
    w = np.sort(np.linalg.eigvalsh(G).real)[::-1]; w = np.clip(w, 0, None)
    return w / w.sum()

def need(w, f):
    return int(np.searchsorted(np.cumsum(w), f) + 1)

print(f"lambda {LAM:.5f} A, kmax {kmax:.3f} 1/A, aperture px {ap.sum()}")
print(f"{'point':<16}{'state':<6}{'top 6 occupancies':<44}{'n90':>5}{'n95':>5}{'n99':>5}")
pts = [("fs 7.29", 7.29, 0), ("fs 5", 5, 0), ("fs 3.5", 3.5, 0), ("fs 2.5", 2.5, 0), ("fs 1.5", 1.5, 0),
       ("fs 0.75", 0.75, 0), ("src 0.05", 0, 0.05), ("src 0.1", 0, 0.1), ("src 0.15", 0, 0.15), ("src 0.2", 0, 0.2), ("src 0.4", 0, 0.4), ("src 0.8", 0, 0.8), ("_coh 7.29+0.4", 7.29, 0.4),
       ("fs 2.5+src 0.1", 2.5, 0.1), ("fs 2.5+src 0.2", 2.5, 0.2), ("fs 3.5+src 0.1", 3.5, 0.1), ("fs 3.5+src 0.2", 3.5, 0.2)]
for name, fs, src in pts:
    for tag, fn in (("sim16", offsets_sim), ("cont", offsets_cont)):
        w = spectrum(*fn(fs, src))
        print(f"{name:<16}{tag:<6}{' '.join(f'{x:.3f}' for x in w[:6]):<44}{need(w, .90):>5}{need(w, .95):>5}{need(w, .99):>5}")
