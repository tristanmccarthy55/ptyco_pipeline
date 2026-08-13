#!/usr/bin/env python
"""@file aberration_compare.py
@brief JEOL-ARM aberration-retrieval result: depth resolution + probe recovery.

Compares the three legs of run_aberration_experiment.sh (all NL70, 70 mrad, reg off):
  perfect70      -- aberration-free 70-mrad reference
  ab_fitprobe    -- aberrated data, probe FITTED from a nominal start (the experiment)
  ab_knownprobe  -- aberrated data, TRUE probe fixed (control)

Depth resolution = on-column Pb-plane kz prominence (same metric as figures/dose_fig_common).
Probe recovery = the fitted probe (outputs.probe of ab_fitprobe) vs the true aberrated probe
(sim_out_aberrated70/01/probe_initial_true.mat, pulled to ab_experiment/).

Read: fit-probe prominence ~ perfect  => ptycho retrieved the ARM residual (70-mrad depth
resolution from a 28-mrad corrector). fit << perfect but known ~ perfect => the info is
there, the fit failed (raise PMODES / lower PSTART).

  python analysis/aberration_compare.py
Reads ~/Desktop/ab_experiment/{perfect70,ab_fitprobe,ab_knownprobe}/Niter*.mat (+
probe_initial_true.mat). Writes ~/Desktop/ab_experiment/aberration_compare.png.
"""
import os, glob
import numpy as np
import h5py
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import maximum_filter

ROOT   = os.path.expanduser("~/Desktop/ab_experiment")
LEGS   = ["perfect70", "ab_fitprobe", "ab_knownprobe"]
TITLE  = {"perfect70": "PERFECT-70 (reference)",
          "ab_fitprobe": "AB + FIT-PROBE (experiment)",
          "ab_knownprobe": "AB + KNOWN-PROBE (control)"}
DX     = 0.04916          # in-plane Å/px (BIN=1: 70.008/1426)
PERIOD = 3.9; F_PB = 1.0 / PERIOD          # Pb-plane spacing / frequency (Å⁻¹)
W      = 22               # half in-plane strip (px) for the cross-section


def _newest_mat(d):
    m = glob.glob(os.path.join(d, "**", "Niter*.mat"), recursive=True)
    if not m:
        raise SystemExit(f"no Niter*.mat under {d} — pull the recon first (see run_aberration_experiment.sh)")
    return sorted(m, key=lambda p: int("".join(filter(str.isdigit, os.path.basename(p)))))[-1]


def _resolve(f, node):
    """Return a complex array from an h5 dataset that may be an object-ref or a compound."""
    a = f[node[0, 0]][:] if node.dtype == object or node.dtype.kind == "O" else node[:]
    return (a["real"] + 1j * a["imag"]) if a.dtype.names else a


def load_vol(leg):
    m = _newest_mat(os.path.join(ROOT, leg))
    with h5py.File(m, "r") as f:
        g = f["outputs"]
        layers = []
        for r in g["object_roi"][:, 0]:
            a = f[r][:]; a = (a["real"] + 1j * a["imag"]) if a.dtype.names else a
            layers.append(a.T)
        V = np.angle(np.array(layers)).astype(float); V -= np.median(V, (1, 2), keepdims=True)
        z = g["z_distance"][:, 0]; dz = float(np.median(z[np.isfinite(z)])) * 1e10
    return V, dz


def load_recon_probe(leg):
    m = _newest_mat(os.path.join(ROOT, leg))
    with h5py.File(m, "r") as f:
        pr = f["outputs"]["probe"]
        try:
            a = f[pr[0, 0]][:]                 # array of mode refs -> first mode
        except Exception:
            a = pr[:]
        a = (a["real"] + 1j * a["imag"]) if a.dtype.names else a
    return np.asarray(a).T


def kz_prominence(V, dz):
    dm = V.mean(0); dmn = dm - dm.min(); col = dmn > np.percentile(dmn, 95)
    nL = V.shape[0]; win = np.hanning(nL)[:, None, None]
    P = np.abs(np.fft.rfft((V - V.mean(0, keepdims=True)) * win, axis=0)) ** 2
    kz = np.fft.rfftfreq(nL, d=dz); Pcol = P[:, col].mean(1)
    m = np.abs(kz - F_PB) < 0.03
    base = np.median(Pcol[(kz > F_PB - 0.12) & (kz < F_PB + 0.12)]) + 1e-30
    prom = float(Pcol[m].max() / base) if m.any() else np.nan
    return kz, Pcol, prom


def main():
    vols = {leg: load_vol(leg) for leg in LEGS}
    ref = vols["perfect70"][0]

    # strong interior Pb column on the reference, reused for all legs
    dm = ref.mean(0); dmn = dm - dm.min()
    mask = np.zeros_like(dmn, bool); mask[40:-40, 40:-40] = True
    pk = (dmn == maximum_filter(dmn, 25)) & (dmn > np.percentile(dmn[mask], 97)) & mask
    ys, xs = np.where(pk); j = int(np.argmax(dmn[ys, xs])); yc, xc = int(ys[j]), int(xs[j])

    fig = plt.figure(figsize=(14, 8)); gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.28)
    prom = {}
    for i, leg in enumerate(LEGS):
        V, dz = vols[leg]; nL = V.shape[0]
        cs = V[:, yc-1:yc+2, xc-W:xc+W].mean(1)
        ax = fig.add_subplot(gs[0, i])
        ax.imshow(cs, extent=[-W*DX, W*DX, (nL-0.5)*dz, 0.5*dz], aspect="auto", cmap="inferno",
                  vmin=np.percentile(cs, 5), vmax=np.percentile(cs, 99.3))
        _, _, prom[leg] = kz_prominence(V, dz)
        ax.set_title(f"{TITLE[leg]}\nPb-plane kz prominence {prom[leg]:.1f}×")
        ax.set_xlabel("in-plane (Å)")
        if i == 0:
            ax.set_ylabel("depth z (Å)  [entrance→exit]")

    # kz spectra
    axk = fig.add_subplot(gs[1, 0])
    for leg, c in zip(LEGS, ["k", "tab:red", "tab:blue"]):
        V, dz = vols[leg]; kz, Pcol, _ = kz_prominence(V, dz)
        axk.plot(kz, Pcol / (np.median(Pcol) + 1e-30), color=c, lw=1.6, label=TITLE[leg].split(" (")[0])
    axk.axvline(F_PB, ls=":", color="gray"); axk.set_xlim(0, 0.5)
    axk.set_xlabel("kz (Å⁻¹)"); axk.set_ylabel("on-column power (norm.)")
    axk.set_title(f"depth spectrum — Pb planes at {F_PB:.3f} Å⁻¹"); axk.legend(fontsize=8)

    # probe recovery: true vs fitted
    ov = np.nan
    tp = os.path.join(ROOT, "probe_initial_true.mat")
    try:
        import scipy.io as sio
        tru = np.asarray(sio.loadmat(tp)["probe"])
        rec = load_recon_probe("ab_fitprobe")
        if rec.shape == tru.shape:
            ov = float(np.abs(np.vdot(rec, tru)) / (np.linalg.norm(rec) * np.linalg.norm(tru)))
        for k, (Pw, ttl) in enumerate([(tru, "TRUE aberrated probe"),
                                       (rec, f"RECOVERED probe  (overlap {ov:.2f})")]):
            ax = fig.add_subplot(gs[1, 1+k]); I = np.abs(Pw) ** 2
            cy, cx = np.unravel_index(I.argmax(), I.shape); w = min(220, cy, cx, I.shape[0]-cy, I.shape[1]-cx)
            ax.imshow(I[cy-w:cy+w, cx-w:cx+w] ** 0.35, cmap="inferno"); ax.set_title(ttl); ax.axis("off")
    except Exception as e:
        fig.add_subplot(gs[1, 1]).annotate(f"probe compare skipped:\n{e}", (.5, .5), ha="center");

    fig.suptitle("JEOL ARM (Cs 1 µm + C5/C56 1 mm) opened to 70 mrad — can ptychography fit it?",
                 fontsize=13.5, y=0.98)
    out = os.path.join(ROOT, "aberration_compare.png"); fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)
    print(f"\nPb-plane kz prominence (depth resolution):")
    for leg in LEGS:
        print(f"  {leg:14s} {prom[leg]:6.1f}×")
    if not np.isnan(prom['perfect70']) and prom['perfect70'] > 0:
        print(f"\n  fit-probe / perfect  = {prom['ab_fitprobe']/prom['perfect70']:.2f}   "
              f"(→1.0 = ptycho fully recovered the ARM aberration)")
        print(f"  known-probe / perfect = {prom['ab_knownprobe']/prom['perfect70']:.2f}   "
              f"(the info ceiling if the probe were known)")
    print(f"  recovered-vs-true probe overlap = {ov:.2f}")


if __name__ == "__main__":
    main()
