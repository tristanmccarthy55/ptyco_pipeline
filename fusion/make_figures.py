#!/usr/bin/env python
"""@file make_figures.py
@brief Figures from the reconstructed volumes: the phase volume, atoms in depth, and the EELS axis.

Three figures, all from real reconstructions and real CASTEP spectra:

  1  the phase volume      -- projected lattice, depth profiles against the true atom planes, and
                              the k_z spectra that show what scan sampling does to the depth axis
  2  atoms in depth        -- z-centroids fitted per column against ground truth, and the
                              Pb-to-Ti-O offset that carries sign(delta_z)
  3  the EELS axis         -- the CASTEP displacement ladder, the four domain spectra, and the
                              inversion that turns O-K contrast into |delta_z|

    ~/hyperspy-bundle/bin/python make_figures.py --out-dir figs
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "analysis"))
import analyze_fusion as AF     # noqa: E402
import eels_forward as F        # noqa: E402
import toy_sample as T          # noqa: E402

DESK = "/Users/u2109287/Desktop/fusion_recons"
COL = {"A": "#c0392b", "B": "#2980b9", "C": "#e67e22", "D": "#16a085"}
SPEC = {"Pb": "#6D4AA6", "TiO": "#16786A", "Oeq": "#B3701A"}
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "white"})


def vol(path):
    return AF.load_volume(glob.glob(path, recursive=True)[0])


def budget_for(nl, dx=None):
    """@brief Geometry for a fusion recon; dx must come from the file for a presolve object."""
    b = json.load(open(os.path.join(DESK, "round2", "hollow_budget_production.json")))
    # derive dx from the detector geometry unless the caller knows better (a presolve object has
    # its own, coarser pixel size and must pass it in)
    b["dx_object_A"] = dx if dx is not None else AF.object_pixel_A(b)
    return b


def pixel_size(path):
    import h5py
    with h5py.File(glob.glob(path, recursive=True)[0], "r") as f:
        return float(np.asarray(f["outputs"]["pixel_size"]).ravel()[0]) * 1e10


def kz_spectrum(P, dz):
    R = P.reshape(P.shape[0], -1)
    R = R - R.mean(0)
    t = np.polynomial.polynomial.polyvander(np.arange(P.shape[0]) * dz, 3)
    R = R - t @ np.linalg.lstsq(t, R, rcond=None)[0]
    S = (np.abs(np.fft.rfft(R, axis=0)) ** 2).mean(1)
    f = np.fft.rfftfreq(P.shape[0], d=dz)
    return f, S / np.median(S[f > 0.12])


# ---------------------------------------------------------------- figure 1
def figure_volume(truth, out):
    """@brief The reconstructed phase volume, and what scan sampling does to its depth axis."""
    fig = plt.figure(figsize=(13.2, 7.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.05, 1.25, 1.25])
    a_lat, c_lat = float(truth["a"]), float(truth["c"])

    good = os.path.join(DESK, "round3", "**", "recon_hsa0_NL16_known", "01", "*", "Niter200.mat")
    fine = os.path.join(DESK, "round4", "**", "recon_hsa0_NL26", "01", "*step01*", "Niter200.mat")
    Vg, Vf = vol(good), vol(fine)
    bg = budget_for(16)
    bf = budget_for(26, dx=pixel_size(fine))

    # (a) projected phase, zoomed, with the three column types marked
    ax = fig.add_subplot(gs[0, 0])
    P = Vg.sum(0)
    dx = bg["dx_object_A"]
    n = int(round(3.2 * a_lat / dx))
    o = P.shape[0] // 2
    sub = P[o:o + n, o:o + n]
    ax.imshow(sub, cmap="bone", extent=[0, n * dx, n * dx, 0])
    # anchor the sublattices on the brightest Pb column IN THIS CROP, so the markers cannot
    # drift away from the image the way a coordinate round-trip can
    from scipy.ndimage import gaussian_filter
    sm = gaussian_filter(sub, 0.25 / dx)
    py0, px0 = np.unravel_index(np.argmax(sm), sm.shape)
    step = a_lat / dx
    for lab, (sx, sy), key, mk in (("Pb", (0, 0), "Pb", "o"),
                                   ("Ti + O$_{ap}$", (.5, .5), "TiO", "s"),
                                   ("O$_{eq}$", (.5, 0), "Oeq", "^")):
        gx = (px0 + (np.arange(-2, 4) + sx) * step) * dx
        gy = (py0 + (np.arange(-2, 4) + sy) * step) * dx
        X, Y = np.meshgrid(gx, gy)
        k = (X > .5) & (X < n * dx - .5) & (Y > .5) & (Y < n * dx - .5)
        ax.scatter(X[k], Y[k], s=30, marker=mk, facecolors="none", edgecolors=SPEC[key],
                   linewidths=1.4, label=lab)
    ax.set_xlim(0, n * dx); ax.set_ylim(n * dx, 0)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax.set_title("(a) projected phase — the lattice", loc="left")
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)")

    # (b, c) depth profiles against the true atom planes
    z_off = 0.5 * (float(bg["beam_thickness_A"]) - float(truth["slab_thickness_A"]))
    for col, (lbl, V, b) in enumerate([("0.3 Å step — blind", vol(os.path.join(DESK, "hsa0_Niter200.mat")),
                                        budget_for(24)),
                                       ("depth-capable reference", Vg, bg)]):
        cells, xy, z = AF.column_maps(V, b, truth)
        gi = np.floor(xy[:, 0] / a_lat).astype(int) % int(truth["n_lat"])
        gj = np.floor(xy[:, 1] / a_lat).astype(int) % int(truth["n_lat"])
        dom = truth["domain_grid"][gi, gj]
        ax = fig.add_subplot(gs[col, 1])
        for name in ("A", "B"):
            m = dom == name
            if not m.any():
                continue
            prof = np.mean([AF._norm(cells[i]["TiO"]) for i in np.where(m)[0]], axis=0)
            ax.plot(z, prof, color=COL[name], lw=1.9, label=f"domain {name}")
        d = truth["deltas"][0]
        zs, _ = AF.column_sites(d, truth, "TiO")
        for zz in zs + z_off:
            ax.axvline(zz, color="0.75", lw=0.8, zorder=0)
        ax.set_title(f"({'bc'[col]}) Ti–O column depth profile · {lbl}", loc="left")
        ax.set_xlabel("z (Å)"); ax.set_ylabel("phase (norm.)")
        ax.legend(frameon=False, fontsize=8, loc="upper left", ncol=2)
        ax.text(0.5, -0.19, "grey = true atom planes (domain A)", transform=ax.transAxes,
                ha="center", fontsize=7.5, color="0.45")

    # (d, e) kz spectra
    ax = fig.add_subplot(gs[:, 2])
    runs = [("0.3 Å step, NL24", vol(os.path.join(DESK, "hsa0_Niter200.mat")), 25.368 / 24, "#B3701A"),
            ("0.15 Å step, NL26", Vf, 25.368 / 26, "#16786A"),
            ("70 Å reference (NL70)", None, 0.999, "#6D4AA6")]
    for lbl, V, dz, c in runs:
        if V is None:
            W = np.load("/Users/u2109287/Desktop/NL70_new_vol.npy")
            P = np.angle(W).astype(float)
        else:
            P = V
        P = P - np.median(P, axis=(1, 2), keepdims=True)
        f, S = kz_spectrum(P, dz)
        m = f > 0.10
        ax.semilogy(f[m], S[m], color=c, lw=1.8, label=lbl)
        ax.axvline(1 / (2 * dz), color=c, lw=0.9, ls=":", alpha=0.8)
    ax.axvline(1 / c_lat, color="0.35", lw=1.2, ls="--")
    ax.annotate(f"PbTiO₃ {c_lat:.2f} Å", (1 / c_lat, ax.get_ylim()[1]), xytext=(4, -12),
                textcoords="offset points", fontsize=8, color="0.35")
    ax.set_xlabel("depth spatial frequency $k_z$ (Å$^{-1}$)")
    ax.set_ylabel("power / median band power")
    ax.set_title("(d) the depth axis: sampling decides", loc="left")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.text(0.02, 0.02, "dotted = each run's layer-grid Nyquist", transform=ax.transAxes,
            fontsize=7.5, color="0.45")

    fig.savefig(out, dpi=170)
    print(f"wrote {out}")


# ---------------------------------------------------------------- figure 2
def fit_z(profile, z, z0, half):
    """@brief Sub-layer z centroid of the peak nearest z0, by parabolic fit on the maximum."""
    m = np.abs(z - z0) <= half
    if not m.any():
        return np.nan
    idx = np.where(m)[0]
    k = idx[np.argmax(profile[idx])]
    if k == 0 or k == len(z) - 1:
        return float(z[k])
    y0, y1, y2 = profile[k - 1], profile[k], profile[k + 1]
    den = y0 - 2 * y1 + y2
    d = 0.5 * (y0 - y2) / den if abs(den) > 1e-12 else 0.0
    return float(z[k] + np.clip(d, -1, 1) * (z[1] - z[0]))


def figure_atoms(truth, out):
    """@brief Fitted atom depths against ground truth, and the offset that carries the sign."""
    good = os.path.join(DESK, "round3", "**", "recon_hsa0_NL16_known", "01", "*", "Niter200.mat")
    V, b = vol(good), budget_for(16)
    a_lat, c_lat = float(truth["a"]), float(truth["c"])
    z_off = 0.5 * (float(b["beam_thickness_A"]) - float(truth["slab_thickness_A"]))
    cells, xy, z = AF.column_maps(V, b, truth)
    gi = np.floor(xy[:, 0] / a_lat).astype(int) % int(truth["n_lat"])
    gj = np.floor(xy[:, 1] / a_lat).astype(int) % int(truth["n_lat"])
    dom = truth["domain_grid"][gi, gj]

    rows = []
    for i in range(len(cells)):
        d = truth["delta_grid"][gi[i], gj[i]]
        for kind in ("Pb", "TiO"):
            zs, w = AF.column_sites(d, truth, kind)
            if kind == "TiO":                       # the Ti sub-peak only (drop the apical O)
                zs = zs[::2]
            for zt in zs + z_off:
                if zt < z.min() + 1.2 or zt > z.max() - 1.2:
                    continue
                zf = fit_z(cells[i][kind], z, zt, 0.45 * c_lat)
                rows.append((dom[i], kind, zt, zf))
    R = np.array([(d, k, a, b_) for d, k, a, b_ in rows if np.isfinite(b_)],
                 dtype=object)

    fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.0), constrained_layout=True)

    # (a) found vs true
    for kind in ("Pb", "TiO"):
        m = R[:, 1] == kind
        ax[0].scatter(R[m, 2].astype(float), R[m, 3].astype(float), s=16, alpha=0.65,
                      color=SPEC[kind], label={"Pb": "Pb column", "TiO": "Ti (Ti–O column)"}[kind],
                      edgecolors="none")
    lim = [z.min(), z.max()]
    ax[0].plot(lim, lim, color="0.4", lw=1, ls="--", zorder=0)
    ax[0].set_xlabel("true atom depth (Å)"); ax[0].set_ylabel("fitted depth (Å)")
    ax[0].set_title("(a) atom depths recovered", loc="left")
    ax[0].legend(frameon=False, fontsize=8)

    # (b) residuals
    res = R[:, 3].astype(float) - R[:, 2].astype(float)
    for kind in ("Pb", "TiO"):
        m = R[:, 1] == kind
        ax[1].hist(res[m], bins=np.linspace(-1.2, 1.2, 33), alpha=0.65, color=SPEC[kind],
                   label=f"{kind}  RMS {np.sqrt((res[m]**2).mean()):.2f} Å")
    ax[1].axvline(0, color="0.4", lw=1, ls="--")
    ax[1].set_xlabel("fitted − true depth (Å)"); ax[1].set_ylabel("atoms")
    ax[1].set_title("(b) depth residual", loc="left")
    ax[1].legend(frameon=False, fontsize=8)

    # (c) the readout that actually works: matched filter over the whole depth profile
    theta = AF.eels_theta_from_contrast(truth, 100.0, 75.0)
    dzp = {str(n): float(np.linalg.norm(truth["deltas"][k][:2])
                         / np.tan(np.radians(theta[str(n)])))
           for k, n in enumerate([str(x) for x in truth["names"]])}
    xy2, obs, dom2, z0 = AF.read_sign(V, b, truth, dz_prior=dzp)
    names = [str(x) for x in truth["names"]]
    for k, name in enumerate(names):
        m = dom2 == name
        if not m.any():
            continue
        rng = np.random.default_rng(k)
        ax[2].scatter(np.full(m.sum(), k) + rng.normal(0, .07, m.sum()), obs[m], s=26,
                      alpha=0.75, color=COL[name], edgecolors="none")
        ax[2].plot([k - .3, k + .3], [np.median(obs[m])] * 2, color=COL[name], lw=2.6)
    ax[2].axhline(0, color="0.35", lw=1.1)
    ax[2].set_xticks(range(4))
    ax[2].set_xticklabels([f"{n}\n{'δz > 0' if truth['deltas'][i][2] > 0 else 'δz < 0'}"
                           for i, n in enumerate(names)])
    ax[2].set_xlabel("domain"); ax[2].set_ylabel("matched-filter discriminant")
    ax[2].set_title("(c) sign(δz) from the whole depth profile", loc="left")
    ax[2].text(0.02, 0.97, "above 0 → up", transform=ax[2].transAxes, fontsize=8,
               va="top", color="0.4")
    ax[2].text(0.02, 0.03, "below 0 → down", transform=ax[2].transAxes, fontsize=8,
               va="bottom", color="0.4")

    fig.suptitle("Atoms placed in depth — depth-capable reconstruction, "
                 f"{len(R)} atoms over {len(cells)} unit cells. A single Ti–Pb offset is too noisy "
                 "for the sign;\nthe matched filter over the whole profile gets every cell right.",
                 fontsize=10)
    fig.savefig(out, dpi=170)
    print(f"wrote {out}")


# ---------------------------------------------------------------- figure 3
def figure_eels(truth, out, alpha=100.0, beta=75.0):
    """@brief The CASTEP ladder, the domain spectra, and the inversion to |delta_z|."""
    lad = F.load_ladder()
    s_grid, e, S_par, S_perp = lad
    w = F.aperture_weights(alpha, beta)
    dmax = float(truth["delta_Ti_A"])
    names = [str(x) for x in truth["names"]]
    m = (e >= -3) & (e <= 30)

    fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.1), constrained_layout=True)

    # (a) the displacement ladder, q perp c -- the channel a 100 mrad probe lands on
    ramp = plt.get_cmap("viridis")(np.linspace(0.05, 0.78, len(s_grid)))   # skip the pale yellow
    for i, sv in enumerate(s_grid):
        ax[0].plot(e[m], S_perp[i][m] / S_perp.max(), color=ramp[i], lw=1.7,
                   label=f"s = {sv:.2f}")
    ax[0].set_xlabel("energy − O-K onset (eV)"); ax[0].set_ylabel("intensity (norm.)")
    ax[0].set_title("(a) CASTEP ladder, q ⊥ c", loc="left")
    ax[0].legend(frameon=False, fontsize=7.5, title="along-chain distortion", title_fontsize=7.5,
                 loc="center right")
    ax[0].annotate("π* grows with the\npolar displacement", (8.9, 1.0), xytext=(-4, -30),
                   textcoords="offset points", fontsize=8, color="0.4", ha="left",
                   arrowprops=dict(arrowstyle="-", color="0.55", lw=0.8))

    # (b) the four domains, and the difference
    spec = {n: F.column_spectrum(truth["deltas"][k], dmax, w, lad) for k, n in enumerate(names)}
    ref = spec[names[0]]
    for n in names[1:]:
        ax[1].plot(e[m], (spec[n] - ref)[m] / ref[m].max() * 100, color=COL[n],
                   lw=2.6 if n == "C" else 1.5,
                   label=f"{n} − A   |δz| = {abs(truth['deltas'][names.index(n)][2]):.3f} Å")
    ax[1].annotate("C − A and D − A coincide", (0.97, 0.92), xycoords="axes fraction",
                   ha="right", fontsize=8, color="0.4")
    ax[1].axhline(0, color="0.4", lw=0.8)
    ax[1].set_xlabel("energy − O-K onset (eV)"); ax[1].set_ylabel("difference (% of edge max)")
    ax[1].set_title("(b) domain spectra, differenced", loc="left")
    ax[1].legend(frameon=False, fontsize=7.5)
    ax[1].annotate("B − A ≡ 0 : sign-blind", (0.04, 0.16), xycoords="axes fraction",
                   fontsize=8.5, color=COL["B"], weight="bold")

    # (c) the inversion: contrast against |delta_z|
    th = np.linspace(0, 90, 46)
    lib = [F.column_spectrum(dmax * np.array([np.sin(np.radians(t)), 0, np.cos(np.radians(t))]),
                             dmax, w, lad) for t in th]
    base = lib[0]
    c = np.array([F.contrast(x, base, e) * 100 for x in lib])
    dz = dmax * np.cos(np.radians(th))
    ax[2].plot(dz, c, color="#6D4AA6", lw=2)
    for k, n in enumerate(names):
        dzk = abs(truth["deltas"][k][2])
        ck = np.interp(-dzk, -dz, c)
        ax[2].scatter([dzk], [ck], s=48, color=COL[n], zorder=3, edgecolors="white", linewidths=1.2)
        ax[2].annotate(n, (dzk, ck), xytext=(7 if n in ("A", "C") else -13, 5),
                       textcoords="offset points", fontsize=9, color=COL[n], weight="bold")
    ax[2].set_xlabel("|δz| (Å)"); ax[2].set_ylabel("O-K contrast vs fully along-beam (%)")
    ax[2].set_title("(c) EELS is monotonic in |δz|", loc="left")
    for lbl, k in (("A, B coincide — sign-blind", 0), ("C, D coincide", 2)):
        dzk = abs(truth["deltas"][k][2])
        ax[2].annotate(lbl, (dzk, np.interp(-dzk, -dz, c)), xytext=(-14, 30 if k == 0 else -26),
                       textcoords="offset points", fontsize=8, color="0.4", ha="right",
                       arrowprops=dict(arrowstyle="-", color="0.55", lw=0.8))

    fig.suptitle("The EELS axis: O-K near-edge structure carries |δz|, and nothing about its sign",
                 fontsize=10)
    fig.savefig(out, dpi=170)
    print(f"wrote {out}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default=os.path.join(HERE, "figs"))
    ap.add_argument("--truth", default=os.path.join(HERE, "sample", "toy_truth.npz"))
    ap.add_argument("--fig", nargs="*", default=["1", "2", "3"])
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    truth = np.load(args.truth, allow_pickle=True)
    if "1" in args.fig:
        figure_volume(truth, os.path.join(args.out_dir, "fig1_phase_volume.png"))
    if "2" in args.fig:
        figure_atoms(truth, os.path.join(args.out_dir, "fig2_atoms_depth.png"))
    if "3" in args.fig:
        figure_eels(truth, os.path.join(args.out_dir, "fig3_eels_axis.png"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
