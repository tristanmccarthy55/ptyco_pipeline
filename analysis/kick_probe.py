#!/usr/bin/env python
"""@file kick_probe.py
@brief Did a probe-update reconstruction, started from a kicked probe, find its way back to the true probe?

A "kick" leg (campaign/run_thin_atomfind.sh KICK_LABELS) reconstructs the lab simulation with the probe released,
several modes, and a start probe built from the true tableau with C1, C3 and A1 scaled by KICK (sim/make_probe.py
--scale). This compares, on the reconstruction's own grid:
  ov_start    |<start, true>|, invariant to a real-space shift and a global phase (a shift is a gauge the object
              absorbs): how far the start probe was from the truth
  ov_final    the same for the dominant recovered mode: how close the reconstruction got back
  mode_power  the share of the probe power in each recovered mode (one coherent probe should leave ~all in mode 1)
Frames: probe_initial_true.mat and probe_initial.mat are in the abTEM frame the simulator wrote; the h5's
reconstruction/probes reads back in that frame untransposed (probe_refit_check.py, checked 2026-09-23).

    python analysis/kick_probe.py <recon_af_<label>_kick_NL<n>> --out <dir>
"""
from __future__ import annotations

import argparse
import glob
import json
import os

import numpy as np


def shift_invariant_overlap(a, b):
    """max over real-space shifts of |<a, b>| / (||a|| ||b||)."""
    x = np.fft.ifft2(np.fft.fft2(a) * np.conj(np.fft.fft2(b)))
    return float(np.abs(x).max() / (np.linalg.norm(a) * np.linalg.norm(b)))


def modes_of(P):
    """(n_modes, N, N) from the h5's probes array, whichever axis carries the modes."""
    P = np.asarray(P).squeeze()
    if P.ndim == 2:
        return P[None]
    return np.moveaxis(P, -1, 0) if P.shape[-1] < P.shape[0] else P


def compare(recon_dir, out=None):
    from scipy.io import loadmat
    import h5py
    h5 = sorted(glob.glob(os.path.join(recon_dir, "analysis", "**", "*_recons.h5"), recursive=True),
                key=os.path.getmtime)
    T = loadmat(os.path.join(recon_dir, "01", "probe_initial_true.mat"))["probe"].astype(np.complex128)
    S = loadmat(os.path.join(recon_dir, "01", "probe_initial.mat"))["probe"].astype(np.complex128)
    res = dict(recon_dir=os.path.basename(os.path.normpath(recon_dir)), ov_start=shift_invariant_overlap(S, T))
    js = os.path.join(recon_dir, "01", "probe_initial.json")
    if os.path.exists(js):
        info = json.load(open(js))
        res.update(kick_c1_A=info.get("c1_A"), sim_c1_A=info.get("c1_sim_A"), scaled_terms=info.get("scaled_terms"))
    if not h5:
        res["status"] = "NO H5"
        return res
    with h5py.File(h5[-1], "r") as f:
        M = modes_of(f["reconstruction/probes"][...]).astype(np.complex128)
    if not np.isfinite(M).all():
        res["status"] = "PROBE NOT FINITE"
        return res
    pw = (np.abs(M) ** 2).sum(axis=(1, 2)); order = np.argsort(pw)[::-1]; M, pw = M[order], pw[order]
    res.update(status="ok", n_modes=int(len(M)), mode_power=[float(x) for x in pw / pw.sum()],
               ov_final=shift_invariant_overlap(M[0], T),
               ov_final_modes=[shift_invariant_overlap(m, T) for m in M])
    if out:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        os.makedirs(out, exist_ok=True)
        n = T.shape[0]; c = n // 2; h = n // 4
        fig, ax = plt.subplots(2, 3, figsize=(9.6, 6.6))
        for j, (P, name) in enumerate(((T, "true probe"), (S, f"kicked start ({res['ov_start']:.3f})"),
                                       (M[0], f"recovered, mode 1 ({res['ov_final']:.3f})"))):
            cy, cx = np.unravel_index((np.abs(P) ** 2).argmax(), P.shape)
            P = np.roll(np.roll(P, c - cy, 0), c - cx, 1)
            ax[0, j].imshow(np.abs(P[c - h:c + h, c - h:c + h]) ** 0.8, cmap="magma"); ax[0, j].set_title(name, fontsize=10)
            ax[1, j].imshow(np.angle(P[c - h:c + h, c - h:c + h]), cmap="twilight", vmin=-np.pi, vmax=np.pi)
            for a_ in ax[:, j]:
                a_.set_xticks([]); a_.set_yticks([])
        ax[0, 0].set_ylabel("amplitude"); ax[1, 0].set_ylabel("phase")
        fig.suptitle(f"{res['recon_dir']}: overlap with the true probe, start {res['ov_start']:.3f} -> recovered "
                     f"{res['ov_final']:.3f}; mode powers " + " / ".join(f"{x:.2f}" for x in res["mode_power"]), fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(out, res["recon_dir"] + ".png"), dpi=110)
        plt.close(fig)
        json.dump(res, open(os.path.join(out, res["recon_dir"] + ".json"), "w"), indent=1)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("recon_dirs", nargs="+")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    for d in a.recon_dirs:
        r = compare(d, a.out)
        print(json.dumps({k: v for k, v in r.items() if k != "ov_final_modes"}))


if __name__ == "__main__":
    main()
