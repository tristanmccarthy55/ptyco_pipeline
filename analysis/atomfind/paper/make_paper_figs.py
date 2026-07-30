#!/usr/bin/env python
"""@file make_paper_figs.py
@brief Publication figures for the atom-localisation methods section.

Produces three vector PDFs (+ PNG previews) in paper/figs/, each backing a specific claim in
atomfind_methods.tex:

  fig_method.pdf      -- forward model: measured anisotropic PSF; a superposition of that PSF
                         at the GT depths reproduces a real column profile (corr 0.94); blind
                         atoms + calibrated depth intervals on a raw cross-section.
  fig_uncertainty.pdf -- (a) variance inflation factor vs interatomic separation (why the
                         per-peak CRLB is optimistic on a dense lattice); (b) HELD-OUT
                         reliability of the raw model sigma vs the conformal interval;
                         (c) per-stratum coverage, pooled vs Mondrian.
  fig_comparison.pdf  -- oxygen / titanium recall by method (read from report.json).

Run:  ~/hyperspy-bundle/bin/python atomfind/paper/make_paper_figs.py
(one finder run ~1-2 min; everything else is cheap).
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))          # -> analysis/
from atomfind import config, align, psf as psfmod, find, validate, uncertainty
from atomfind.fit import _render
from atomfind.deconv import crop_kernel_inplane

FIGDIR = os.path.join(HERE, "figs")
os.makedirs(FIGDIR, exist_ok=True)
COL = {82: "#2166ac", 22: "#8a8f98", 8: "#d6604d"}                  # Pb / Ti / O
LAB = {82: "Pb", 22: "Ti", 8: "O"}

plt.rcParams.update({
    "font.size": 8.5, "font.family": "serif", "mathtext.fontset": "cm",
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "axes.labelsize": 9, "axes.titlesize": 9, "legend.fontsize": 7.5,
    "figure.dpi": 150, "savefig.dpi": 300,
})


def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"{stem}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", stem + ".pdf/.png")


# ---------------------------------------------------------------- shared compute
def compute():
    cfg = config.preset("NL70_coherent")
    V, dx = align.load_phase(cfg)
    pos, Z = align.load_gt(cfg)
    al = align.register(V, dx, pos, Z, cfg)
    kernels = psfmod.species_kernels(cfg, dx)
    found, seeds = find.find_atoms_v3(V, cfg, dx, kernels)
    al = align.refine_with_atoms(al, found, pos, Z, cfg)
    rep, m = validate.finder_report(found, pos, Z, al, cfg)
    return dict(cfg=cfg, V=V, dx=dx, pos=pos, Z=Z, al=al, kernels=kernels,
                found=found, rep=rep, m=m)


# ---------------------------------------------------------------- VIF vs separation
def vif_curve(cfg, dx, K):
    K = find._unit_norm(crop_kernel_inplane(K))
    hz = (K.shape[0]-1)//2; hxy = (K.shape[1]-1)//2
    Kz, Ky, Kx = np.gradient(K)
    box = (K.shape[0]+30, K.shape[1]+8, K.shape[2]+8)
    cz, cr, cc = box[0]/2, box[1]/2, box[2]/2

    def cols(l):
        Kp = _render(box, (l, cr, cc), K, hz, hxy)
        Gz = _render(box, (l, cr, cc), Kz, hz, hxy)
        Gy = _render(box, (l, cr, cc), Ky, hz, hxy)
        Gx = _render(box, (l, cr, cc), Kx, hz, hxy)
        return np.column_stack([Kp, -Gz, -Gy, -Gx])
    seps = np.arange(0.7, 4.6, 0.1)
    vz, vb = [], []
    for s in seps:
        J1 = cols(cz); J2 = cols(cz + s/cfg.dz)
        J = np.hstack([J1, J2]); F = J.T @ J
        joint = np.linalg.inv(F); cond = np.linalg.inv(F[:4, :4])
        vb.append(joint[0, 0]/cond[0, 0]); vz.append(joint[1, 1]/cond[1, 1])
    return seps, np.array(vz), np.array(vb)


# ---------------------------------------------------------------- reliability (held-out)
def _q_finite(scores, alpha):
    s = np.sort(scores[np.isfinite(scores)])
    n = s.size
    if n < 10:
        return np.nan
    k = int(np.ceil((n+1)*(1-alpha)))
    return s[min(k, n)-1]


def reliability(found, m, cfg, axis="z", seed=0):
    """Split matched atoms 50/50; calibrate q on one half, measure coverage on the other,
    for model-sigma (Gaussian) and for conformal, across nominal levels."""
    from scipy.stats import norm
    matched = np.where(m["match_gi"] >= 0)[0]
    err = np.abs(m[f"match_d{axis}"]); sig = found[f"s{axis}_A"]
    ok = matched[np.isfinite(err[matched]) & (sig[matched] > 0)]
    rng = np.random.default_rng(seed); rng.shuffle(ok)
    cal, test = ok[:len(ok)//2], ok[len(ok)//2:]
    strat = uncertainty.strata_array(found, cfg)
    nominal = np.linspace(0.50, 0.99, 22)
    cov_model, cov_conf, cov_mond = [], [], []
    sc_cal = err[cal]/sig[cal]
    for nom in nominal:
        a = 1-nom
        # raw model sigma, treated as a Gaussian 1-sigma -> two-sided interval
        k = norm.ppf(1-a/2)
        cov_model.append(np.mean(err[test] <= k*sig[test]))
        # pooled conformal
        q = _q_finite(sc_cal, a)
        cov_conf.append(np.mean(err[test] <= q*sig[test]))
        # Mondrian conformal (per-stratum q)
        hit = []
        for i in test:
            cs = err[cal[strat[cal] == strat[i]]] / sig[cal[strat[cal] == strat[i]]]
            qi = _q_finite(cs, a)
            if not np.isfinite(qi):
                qi = q
            hit.append(err[i] <= qi*sig[i])
        cov_mond.append(np.mean(hit))
    return nominal, np.array(cov_model), np.array(cov_conf), np.array(cov_mond)


def per_stratum_cov(found, m, cfg, alpha=0.05, seed=0, axis="z"):
    matched = np.where(m["match_gi"] >= 0)[0]
    err = np.abs(m[f"match_d{axis}"]); sig = found[f"s{axis}_A"]
    ok = matched[np.isfinite(err[matched]) & (sig[matched] > 0)]
    rng = np.random.default_rng(seed); rng.shuffle(ok)
    cal, test = ok[:len(ok)//2], ok[len(ok)//2:]
    strat = uncertainty.strata_array(found, cfg)
    qp = _q_finite(err[cal]/sig[cal], alpha)
    rows = []
    for st in sorted(set(strat[test])):
        ti = test[strat[test] == st]
        if ti.size < 8:
            continue
        ci = cal[strat[cal] == st]
        qm = _q_finite(err[ci]/sig[ci], alpha)
        if not np.isfinite(qm):
            qm = qp
        cov_p = np.mean(err[ti] <= qp*sig[ti])
        cov_m = np.mean(err[ti] <= qm*sig[ti])
        rows.append((st, ti.size, cov_p, cov_m))
    return rows


# ================================================================ FIGURE 1
def fig_method(D):
    cfg, V, dx, al = D["cfg"], D["V"], D["dx"], D["al"]
    pos, Z, found = D["pos"], D["Z"], D["found"]
    K3 = psfmod.empirical_psf(cfg, dx)
    fig = plt.figure(figsize=(7.2, 2.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.85, 1.1, 1.0], wspace=0.42)

    # (a) measured PSF, x-z slice
    axa = fig.add_subplot(gs[0])
    hz = (K3.shape[0]-1)//2; hxy = (K3.shape[1]-1)//2
    fz, fxy = psfmod.measure_fwhm(K3, cfg, dx)
    ext = [-hxy*dx, hxy*dx, hz*cfg.dz, -hz*cfg.dz]
    axa.imshow(K3[:, hxy, :], extent=ext, aspect="auto", cmap="magma")
    axa.set_xlim(-1.2, 1.2); axa.set_ylim(3.5, -3.5)
    axa.set_xlabel(r"$x$ (\AA)"); axa.set_ylabel(r"depth $z$ (\AA)")
    axa.set_title(r"(a) measured PSF $K$")
    axa.text(0.04, 0.03, f"FWHM\n$z$={fz:.1f}, $xy$={fxy:.2f} " + r"\AA",
             transform=axa.transAxes, color="w", fontsize=7, va="bottom")

    # (b) forward-model validation: profile vs superposition of K at GT depths
    axb = fig.add_subplot(gs[1])
    from atomfind.run_atomfind import _pick_BO_column
    r, c, on = _pick_BO_column(pos, Z, cfg, al, V)
    pr, pc, prof = find.track_column(V, r, c, cfg)
    prof = np.clip(prof - np.percentile(prof, 10), 0, None)
    k1d, zoff = psfmod.axial_kernel(K3, navg=cfg.find_profile_navg)
    nL = V.shape[0]
    gr, gc, gl = al.site_to_index(pos[on, 0], pos[on, 1], pos[on, 2])
    A = np.zeros((nL, len(on)))
    for j, lf in enumerate(gl):
        A[:, j] = np.interp(np.arange(nL) - lf, zoff, k1d, left=0, right=0)
    from scipy.optimize import nnls
    amp, _ = nnls(A, prof)
    model = A @ amp
    zc = (np.arange(nL)+0.5)*cfg.dz
    keep = (zc > 2) & (zc < 66)
    cc_corr = np.corrcoef(prof[keep], model[keep])[0, 1]
    axb.plot(zc[keep], prof[keep]/prof[keep].max(), "o", ms=2.4, color="0.25",
             label="recon profile", zorder=3)
    axb.plot(zc[keep], model[keep]/prof[keep].max(), "-", lw=1.3, color="#d6604d",
             label=r"$\sum_i\beta_i K$ at GT depths")
    for j in range(len(on)):
        zz = gl[j]*cfg.dz
        if 2 < zz < 66:
            axb.plot([zz], [-0.05], marker="^", ms=3,
                     color=COL.get(int(Z[on[j]]), "k"), clip_on=False)
    axb.set_xlim(2, 66); axb.set_ylim(-0.08, 1.15)
    axb.set_xlabel(r"depth $z$ (\AA)"); axb.set_ylabel("norm. phase")
    axb.set_title("(b) forward model on a B–O column")
    axb.text(0.97, 0.93, rf"$r={cc_corr:.2f}$", transform=axb.transAxes,
             ha="right", va="top", fontsize=9)
    axb.legend(loc="upper left", handlelength=1.4)

    # (c) cross-section with blind atoms + calibrated 95% depth bars
    axc = fig.add_subplot(gs[2])
    _, m = validate.finder_report(found, pos, Z, al, cfg)
    qtab = uncertainty.calibrate(found, m, cfg, alphas=cfg.uq_alphas, min_n=cfg.uq_min_stratum)
    hw = uncertainty.apply(found, qtab, cfg, 0.05)
    W = int(round(2.6/dx))
    cs = V[:, r-1:r+2, c-W:c+W].mean(1)
    l0 = int(round(2/cfg.dz)); l1 = int(round(66/cfg.dz))
    ext = [(-W-0.5)*dx, (W-0.5)*dx, (cs.shape[0]-0.5)*cfg.dz, -0.5*cfg.dz]
    axc.imshow(cs, extent=ext, aspect="auto", cmap="inferno",
               vmin=max(np.percentile(cs[l0:l1], 5), 0), vmax=np.percentile(cs[l0:l1], 99.5))
    strip = np.where((np.abs(found["row"]-r)*dx < 0.35) & (np.abs(found["col"]-c)*dx < 2.6))[0]
    for i in strip:
        zz = found["layer"][i]*cfg.dz
        if zz > 66:
            continue
        sp = int(found["species"][i])
        axc.errorbar([(found["col"][i]-c)*dx], [zz], yerr=hw["z"][i], xerr=hw["x"][i],
                     fmt="o", ms=3, color=COL.get(sp, "w"), ecolor=COL.get(sp, "w"),
                     elinewidth=0.8, capsize=1.5, markeredgecolor="k", markeredgewidth=0.3)
    axc.set_xlim((-W-0.5)*dx, (W-0.5)*dx); axc.set_ylim(66, -0.5)
    axc.set_xlabel(r"$x$ (\AA)"); axc.set_ylabel(r"depth $z$ (\AA)")
    axc.set_title("(c) blind atoms + 95\\% interval")
    h = [plt.Line2D([], [], marker="o", ls="", mfc=COL[z], mec="k", ms=4, label=LAB[z])
         for z in (22, 8)]
    axc.legend(handles=h, loc="lower right", framealpha=0.85)
    save(fig, "fig_method")


# ================================================================ FIGURE 2
def fig_uncertainty(D):
    cfg, dx, found, m = D["cfg"], D["dx"], D["found"], D["m"]
    fig = plt.figure(figsize=(7.2, 2.4))
    gs = fig.add_gridspec(1, 3, wspace=0.42)

    # (a) VIF vs separation
    axa = fig.add_subplot(gs[0])
    seps, vz, vb = vif_curve(cfg, dx, D["kernels"][82])
    axa.plot(seps, vz, "-", lw=1.5, color="#d6604d", label=r"$\sigma_z$ (depth)")
    axa.plot(seps, vb, "-", lw=1.2, color="#2166ac", label=r"$\beta$ (amplitude)")
    axa.axhline(1, color="0.6", lw=0.7, ls=":")
    for sx, lab in [(1.95, "Ti–O"), (3.90, "Pb–Pb")]:
        axa.axvline(sx, color="0.6", lw=0.7, ls="--")
        axa.text(sx, axa.get_ylim()[1], f" {lab}", fontsize=6.5, va="top", rotation=90)
    axa.set_xlabel(r"neighbour separation (\AA)"); axa.set_ylabel("variance inflation factor")
    axa.set_title("(a) VIF of the joint fit"); axa.legend(handlelength=1.4)
    axa.set_ylim(0.8, None)

    # (b) held-out reliability
    axb = fig.add_subplot(gs[1])
    nom, cmod, cconf, cmond = reliability(found, m, cfg, axis="z")
    axb.plot([0.5, 1], [0.5, 1], "-", color="0.6", lw=0.8, label="ideal")
    axb.plot(nom, cmod, "s-", ms=3, color="#b2182b", label=r"raw model $\sigma$ (Gaussian)")
    axb.plot(nom, cmond, "o-", ms=3, color="#1b7837", label="conformal (Mondrian)")
    axb.set_xlabel("nominal coverage"); axb.set_ylabel("empirical coverage (held-out)")
    axb.set_title("(b) depth-interval reliability"); axb.legend(loc="upper left")
    axb.set_xlim(0.5, 1); axb.set_ylim(0, 1.02)

    # (c) per-stratum coverage, pooled vs Mondrian
    axc = fig.add_subplot(gs[2])
    rows = per_stratum_cov(found, m, cfg, alpha=0.05)
    rows = sorted(rows, key=lambda x: x[2])[:8]                     # worst 8 for pooled
    labels = [r[0].replace("|", "\n") for r in rows]
    xx = np.arange(len(rows))
    axc.bar(xx-0.2, [r[2] for r in rows], 0.38, color="#b2182b", label="pooled q")
    axc.bar(xx+0.2, [r[3] for r in rows], 0.38, color="#1b7837", label="Mondrian q")
    axc.axhline(0.95, color="0.3", lw=0.9, ls="--")
    axc.set_xticks(xx); axc.set_xticklabels(labels, fontsize=5.2)
    axc.set_ylabel("held-out coverage @95\\%"); axc.set_ylim(0.6, 1.02)
    axc.set_title("(c) coverage per stratum"); axc.legend(loc="lower left")
    save(fig, "fig_uncertainty")


# ================================================================ FIGURE 3
def fig_comparison(D):
    path = os.path.join(D["cfg"].out_dir, "report.json")
    if not os.path.exists(path):
        print("skip fig_comparison: no report.json (run run_atomfind.py first)")
        return
    fnd = json.load(open(path))["finder"]
    methods = [("3-D peaks (raw)", "peaks3d_raw"), ("RL + peaks", "peaks3d_rl"),
               ("MEM + peaks", "peaks3d_mem"), ("1-D spike", "v1_spike"),
               ("model fit (this work)", "v3")]
    methods = [(lbl, k) for lbl, k in methods if k in fnd]
    fig, ax = plt.subplots(figsize=(3.6, 2.5))
    xx = np.arange(len(methods))
    o = [fnd[k]["O"]["recall"]*100 for _, k in methods]
    ti = [fnd[k]["Ti"]["recall"]*100 for _, k in methods]
    ax.bar(xx-0.2, ti, 0.38, color="#8a8f98", label="Ti")
    ax.bar(xx+0.2, o, 0.38, color="#d6604d", label="O")
    for i, v in enumerate(o):
        ax.text(i+0.2, v+1, f"{v:.0f}", ha="center", fontsize=6.5)
    ax.set_xticks(xx); ax.set_xticklabels([lbl for lbl, _ in methods], rotation=25,
                                          ha="right", fontsize=6.8)
    ax.set_ylabel("recall (\\%)"); ax.set_ylim(0, 100)
    ax.set_title("Oxygen recovery by method"); ax.legend(loc="upper left")
    save(fig, "fig_comparison")


if __name__ == "__main__":
    D = compute()
    print(f"[compute] {len(D['found'])} atoms; precision {D['rep']['precision']:.2f}")
    fig_method(D)
    fig_uncertainty(D)
    fig_comparison(D)
