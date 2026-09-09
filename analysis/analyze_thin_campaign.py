#!/usr/bin/env python
"""@file analyze_thin_campaign.py
@brief Sweep analysis for the thin-slab aberration campaign — depth x-z, in-plane, probe recovery.

Reads the PtychoShelves *_recons.h5 for recon_round_a<alpha>_{perfect,known,fitprobe}_NL* under a
given directory (the extracted results tarball). Per alpha: object depth (x-z) section, in-plane
depth-summed phase, and shift-invariant probe overlap (blind-recovered vs the known/true probe).
The in-plane phase RAMP (a gauge artifact on fixed-probe legs) is removed by a plane fit.

    ~/hyperspy-bundle/bin/python analysis/analyze_thin_campaign.py <dir-with-recon_round_*> [--out DIR]

Figures default into aberration_experiment/figs/<ISO-week>/ (committed, dated) — not the Desktop.
"""
from __future__ import annotations
import sys, os, re, argparse, datetime
import numpy as np, h5py, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from glob import glob

THICK = 11.715                                        # atomic slab beam path [Å] (3-cell base)
PLANES = np.array([2.07, 4.00, 5.95, 7.89, 9.84])     # AO/BO2 cation-plane depths (build_thin_sample)

def week_dir(sub="figs"):
    repo = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    d = os.path.join(repo, "aberration_experiment", sub, datetime.date.today().strftime("%G-W%V"))
    os.makedirs(d, exist_ok=True); return d

def load(h5):
    f = h5py.File(h5, "r"); r = f["reconstruction"]
    obj = np.asarray(r["object"]).squeeze()
    if obj.ndim == 2: obj = obj[None]
    prb = np.asarray(r["probes"]).squeeze()
    prb0 = np.asarray(r["p"]["probe_initial"]).squeeze()
    dx = float(np.asarray(r["p"]["dx_spec"]).ravel()[0]) * 1e10
    il = np.asarray(r["p"]["illum_sum"]["illum_sum_0"]); f.close()
    return dict(obj=obj, prb=prb, prb0=prb0, dx=dx, nl=obj.shape[0], il=il)

def crop_phase(d, half_A=9.0, deplane=True):
    """per-layer DC-removed phase, cropped to +-half_A about the illuminated centre, ramp removed."""
    il = d["il"]; m = il > 0.15*il.max(); ys, xs = np.where(m)
    cy, cx = int((ys.min()+ys.max())/2), int((xs.min()+xs.max())/2); h = int(half_A/d["dx"])
    ph = np.angle(d["obj"]); ph = ph - np.median(ph, axis=(1, 2), keepdims=True)
    ph = ph[:, max(0,cy-h):cy+h, max(0,cx-h):cx+h]
    if deplane and ph.shape[1] > 4 and ph.shape[2] > 4:
        s = ph.sum(0); yy, xx = np.mgrid[0:s.shape[0], 0:s.shape[1]].astype(float)
        A = np.c_[xx.ravel(), yy.ravel(), np.ones(xx.size)]
        c, *_ = np.linalg.lstsq(A, s.ravel(), rcond=None)
        ph = ph - ((c[0]*xx + c[1]*yy)/ph.shape[0])[None]
    return ph

def overlap(a, b):
    """shift- and global-phase-invariant probe coherence (max normalised complex xcorr)."""
    from numpy.fft import fft2, ifft2
    a = a/(np.linalg.norm(a)+1e-30); b = b/(np.linalg.norm(b)+1e-30)
    return float(np.abs(ifft2(fft2(a)*np.conj(fft2(b)))).max())

def find_round(base):
    out = {}
    for h in glob(f"{base}/**/recon_round_a*_*_NL*/analysis/*/*/*_recons.h5", recursive=True):
        m = re.search(r"recon_round_a(\d+)_(perfect|known|fitprobe)_NL(\d+)", h)
        if m: out.setdefault(int(m.group(1)), {})[m.group(2)] = h
    return dict(sorted(out.items()))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base", help="directory holding recon_round_* (the extracted tarball)")
    ap.add_argument("--out", default=None, help="figure dir (default aberration_experiment/figs/<week>)")
    a = ap.parse_args()
    out = a.out or week_dir("figs"); os.makedirs(out, exist_ok=True)
    R = find_round(a.base); alphas = list(R.keys())
    if not alphas: sys.exit(f"no recon_round_* found under {a.base}")
    print("alphas:", alphas, "-> figs to", out)
    D = {al: {leg: load(h) for leg, h in R[al].items()} for al in alphas}
    legs = ["perfect", "known", "fitprobe"]; titles = {"perfect":"perfect","known":"known (true)","fitprobe":"fit (blind)"}

    # probe recovery
    ov = {}
    for al in alphas:
        if "known" in D[al] and "fitprobe" in D[al]:
            ov[al] = (overlap(D[al]["known"]["prb"], D[al]["fitprobe"]["prb"]),
                      overlap(D[al]["known"]["prb"], D[al]["fitprobe"]["prb0"]))
            print(f"  a{al}: NL={D[al]['known']['nl']:2d}  probe overlap fit={ov[al][0]:.3f} (start {ov[al][1]:.3f})")

    # FIG depth x-z
    fig, ax = plt.subplots(len(alphas), 3, figsize=(13, 2.0*len(alphas)), squeeze=False, constrained_layout=True)
    for i, al in enumerate(alphas):
        for j, leg in enumerate(legs):
            aa = ax[i][j]
            if leg not in D[al]: aa.axis("off"); continue
            ph = crop_phase(D[al][leg]); yb = ph.shape[1]//2
            band = ph[:, max(0,yb-8):yb+8, :].mean(1); A = ph.shape[2]*D[al][leg]["dx"]
            v = np.percentile(np.abs(band), 98)+1e-9
            aa.imshow(band, cmap="magma", vmin=-v, vmax=v, extent=[-A/2, A/2, 0, THICK], origin="lower", aspect="auto", interpolation="bilinear")
            for zp in PLANES: aa.axhline(zp, color="cyan", lw=0.4, alpha=0.35)
            aa.set_xticks([])
            if j == 0: aa.set_ylabel(f"a{al} NL{D[al][leg]['nl']}\nz(Å)", fontsize=8)
            else: aa.set_yticks([])
            if i == 0: aa.set_title(titles[leg], fontsize=11)
    fig.suptitle("ROUND sweep — depth (x–z) vs α  (cyan = AO/BO₂ planes)", fontsize=13)
    fig.savefig(os.path.join(out, "round_depth.png"), dpi=140); plt.close(fig)

    # FIG in-plane
    fig, ax = plt.subplots(len(alphas), 3, figsize=(11, 2.0*len(alphas)), squeeze=False, constrained_layout=True)
    for i, al in enumerate(alphas):
        for j, leg in enumerate(legs):
            aa = ax[i][j]
            if leg not in D[al]: aa.axis("off"); continue
            ph = crop_phase(D[al][leg]).sum(0); A = ph.shape[1]*D[al][leg]["dx"]
            v0, v1 = np.percentile(ph, [1, 99])
            aa.imshow(ph, cmap="gray", vmin=v0, vmax=v1, extent=[-A/2, A/2, -A/2, A/2], origin="lower")
            aa.set_xticks([]); aa.set_yticks([])
            if j == 0: aa.set_ylabel(f"a{al}", fontsize=9)
            if i == 0: aa.set_title(titles[leg], fontsize=11)
    fig.suptitle("ROUND sweep — in-plane (depth-summed phase) vs α", fontsize=13)
    fig.savefig(os.path.join(out, "round_inplane.png"), dpi=140); plt.close(fig)

    # FIG probe recovery
    if ov:
        al = sorted(ov); fig, axc = plt.subplots(figsize=(7, 4.4), constrained_layout=True)
        axc.plot(al, [ov[a][0] for a in al], "o-", label="blind-recovered vs true")
        axc.plot(al, [ov[a][1] for a in al], "s--", color="gray", label="perfect start vs true")
        axc.axvline(70, color="r", ls=":", lw=1); axc.set_ylim(0, 1.02)
        axc.set_xlabel("α (mrad)"); axc.set_ylabel("shift-invariant |⟨true, probe⟩|"); axc.legend(fontsize=9)
        axc.set_title("blind probe retrieval vs α", fontsize=12)
        fig.savefig(os.path.join(out, "round_probe_overlap.png"), dpi=140); plt.close(fig)
    print("wrote round_depth.png, round_inplane.png" + (", round_probe_overlap.png" if ov else ""))

if __name__ == "__main__":
    main()
