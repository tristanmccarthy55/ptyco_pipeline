#!/usr/bin/env python
"""End-to-end atom-finding pipeline: PSF -> deconvolution -> per-site fit -> validation.

Runs UNCHANGED on the current NL70 volume and on the better data (switch with --preset,
or point config.preset('reviewer2') at the new .npy). Requires the abtem/skimage stack:

    ~/hyperspy-bundle/bin/python run_atomfind.py                 # NL70, data PSF (default)
    ~/hyperspy-bundle/bin/python run_atomfind.py --psf all       # compare data vs synthetic
    ~/hyperspy-bundle/bin/python run_atomfind.py --preset reviewer2 --dose 1e8 \
        --single-atom-vol ~/Desktop/psf_atom_vol.npy             # when the sim PSF lands

Outputs (to cfg.out_dir, default ~/Desktop/atomfind_out): figures + report.json + a
printed CAN / CAN'T-SHOW summary.
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from atomfind import config, align, psf as psfmod, deconv, fit, validate, find


# ---------------------------------------------------------------- null sites
def make_null(pos, cfg, n=400, min_sep_A=0.9, seed=1):
    rng = np.random.default_rng(seed)
    win = align.in_window(pos, cfg)
    A = pos[win]
    cx, cy = cfg.scan_center_xy
    h = cfg.scan_window_A / 2
    out = []
    tries = 0
    while len(out) < n and tries < 50 * n:
        tries += 1
        p = np.array([rng.uniform(cx - h, cx + h), rng.uniform(cy - h, cy + h),
                      rng.uniform(3, cfg.zmax_show_A - 3)])
        if np.hypot(A[:, 0] - p[0], A[:, 1] - p[1]).min() > min_sep_A:
            out.append(p)
    return np.array(out)


# ---------------------------------------------------------------- figures
def fig_psf(psfs, cfg, dx, path):
    names = [n for n in ("empirical", "data", "synthetic") if n in psfs]
    fig, axes = plt.subplots(2, len(names), figsize=(4.2*len(names), 7), squeeze=False)
    for j, nm in enumerate(names):
        k = psfs[nm]; hz = (k.shape[0]-1)//2; hxy = (k.shape[1]-1)//2
        fz, fxy = psfmod.measure_fwhm(k, cfg, dx)
        ax = axes[0][j]
        ext = [-hxy*dx, hxy*dx, hz*cfg.dz, -hz*cfg.dz]
        ax.imshow(k[:, hxy, :], extent=ext, aspect="auto", cmap="inferno")
        ax.set_title(f"{nm} PSF (x-z)\nFWHM z={fz:.1f} A  xy={fxy:.2f} A")
        ax.set_xlabel("in-plane x (A)"); ax.set_ylabel("depth z (A)")
        ax2 = axes[1][j]
        zc = k[:, hxy, hxy]; zc = zc/zc.max()
        xc = k[hz, hxy, :]; xc = xc/xc.max()
        ax2.plot(np.arange(-hz, hz+1)*cfg.dz, zc, "-o", ms=3, label="along z")
        ax2.plot(np.arange(-hxy, hxy+1)*dx, xc, "-", label="in-plane")
        ax2.axhline(0.5, color="gray", ls=":", lw=0.8)
        ax2.set_xlim(-4, 4); ax2.set_xlabel("offset (A)"); ax2.set_ylabel("norm."); ax2.legend(fontsize=8)
    fig.suptitle("Point-spread functions (measured vs synthetic) — anisotropic, missing-cone along z")
    fig.tight_layout(); fig.savefig(path, dpi=130, bbox_inches="tight"); plt.close(fig)


def fig_amplitude(rec, olabel, cfg, path):
    groups = [("Pb", rec["Z"] == 82, "#00f5ff"),
              ("Ti", rec["Z"] == 22, "#39ff14"),
              ("O (in-plane\nisolated)", (rec["Z"] == 8) & (olabel == 2), "#ff21ff"),
              ("O (axial\noverlap)", (rec["Z"] == 8) & (olabel == 1), "#ff8c00"),
              ("null", rec["Z"] == -1, "#888888")]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    data = [np.clip(rec["beta"][m], 1e-3, None) for _, m, _ in groups]
    parts = ax.violinplot(data, showmedians=True, showextrema=False)
    for pc, (_, _, col) in zip(parts["bodies"], groups):
        pc.set_facecolor(col); pc.set_alpha(0.6)
    nullb = rec["beta"][rec["Z"] == -1]
    tau = np.quantile(nullb, 0.95)
    ax.axhline(tau, color="red", ls="--", lw=1.2, label=f"null p95 = {tau:.2f} (detect >)")
    ax.set_yscale("log"); ax.set_xticks(range(1, len(groups)+1))
    ax.set_xticklabels([g[0] for g in groups], fontsize=9)
    ax.set_ylabel("fitted amplitude beta (log)"); ax.legend(fontsize=9)
    ax.set_title("Per-site fitted amplitude by species")
    # amplitude vs Z
    pts, slope = validate.amplitude_vs_Z(rec)
    ax2.plot(pts[:, 0], pts[:, 1], "ko-")
    for z, b in pts:
        ax2.annotate({0: "null", 8: "O", 22: "Ti", 82: "Pb"}.get(int(z), ""),
                     (z, b), textcoords="offset points", xytext=(6, 4))
    if np.isfinite(slope):
        zz = np.linspace(0, 85, 10); ax2.plot(zz, pts[pts[:,0]==0,1][0] + slope*zz, "r:",
                                              label=f"light-atom Z-scaling (slope {slope:.3f})")
    ax2.set_xlabel("atomic number Z"); ax2.set_ylabel("median beta")
    ax2.set_title("Amplitude vs Z (does O follow the lattice?)"); ax2.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=130, bbox_inches="tight"); plt.close(fig)


def fig_roc(rec, olabel, path):
    nullb = rec["beta"][rec["Z"] == -1]
    fig, ax = plt.subplots(figsize=(6.2, 6))
    for nm, m, col in [("O all", rec["Z"] == 8, "#c000c0"),
                       ("O in-plane isolated", (rec["Z"] == 8) & (olabel == 2), "#ff21ff"),
                       ("O axial overlap", (rec["Z"] == 8) & (olabel == 1), "#ff8c00"),
                       ("Ti (reference)", rec["Z"] == 22, "#39aa14")]:
        b = rec["beta"][m]
        fpr, tpr, auc = validate.roc(b, nullb)
        ax.plot(fpr, tpr, color=col, lw=2, label=f"{nm}  AUC={auc:.2f} (n={len(b)})")
    ax.plot([0, 1], [0, 1], "k:", lw=1, label="chance")
    ax.set_xlabel("false-positive rate (null)"); ax.set_ylabel("true-positive rate")
    ax.set_title("Detection ROC vs off-lattice null"); ax.legend(fontsize=9, loc="lower right")
    fig.tight_layout(); fig.savefig(path, dpi=130, bbox_inches="tight"); plt.close(fig)


def _pick_BO_column(pos, Z, cfg, al, V):
    """The BRIGHTEST Ti (B-O) column in the recon -> (row, col, Ti&O atoms on it).

    (Fix: the old version picked the column nearest the field centre, which was 2.4x
    dimmer than the validated overlay's -> a blurry, unrepresentative cross-section.)"""
    win = align.in_window(pos, cfg)
    dm = np.clip(V, 0, None).mean(0)
    ti = np.where(win & (Z == 22))[0]
    best = None
    for j in ti:
        rf, cf, _ = al.site_to_index(pos[j, 0], pos[j, 1], pos[j, 2])
        r, c = int(round(rf)), int(round(cf))
        if 3 <= r < V.shape[1]-3 and 3 <= c < V.shape[2]-3:
            b = dm[r-1:r+2, c-1:c+2].mean()
            if best is None or b > best[0]:
                best = (b, pos[j])
    seed = best[1]
    same = win & (np.hypot(pos[:, 0]-seed[0], pos[:, 1]-seed[1]) < 0.8)
    on = np.where(same & ((Z == 22) | (Z == 8)))[0]
    r, c, _ = al.site_to_index(seed[0], seed[1], seed[2])
    return int(round(r)), int(round(c)), on


SPCOL = {82: "#00f5ff", 22: "#39ff14", 8: "#ff21ff"}       # Pb / Ti / O marker colours


def fig_overlay(V, dec, found, al, pos, Z, cfg, dx, path, W=22):
    """Brightest B-O column: raw + RL-deconvolved cross-sections (GT markers), plus the
    money panel -- the blind v2 atoms (species-coloured, with z error bars) vs the GT comb."""
    r, c, on = _pick_BO_column(pos, Z, cfg, al, V)
    nL = V.shape[0]; zrec = (np.arange(nL)+0.5)*cfg.dz
    ext = [-W*dx, W*dx, zrec[-1], zrec[0]]
    l0 = int(round(cfg.trim_z_A[0]/cfg.dz)); l1 = int(round(cfg.trim_z_A[1]/cfg.dz))
    fig, axes = plt.subplots(1, 3, figsize=(15, 8), gridspec_kw={"width_ratios": [1, 1, 1.2]})
    for ax, img, ttl in [(axes[0], V, "raw phase"), (axes[1], dec, "RL-deconvolved (interior)")]:
        cs = img[:, r-1:r+2, c-W:c+W].mean(1)
        shown = cs[l0:l1]
        ax.imshow(cs, extent=ext, aspect="auto", cmap="inferno",
                  vmin=max(np.percentile(shown, 5), 0), vmax=np.percentile(shown, 99.5))
        for j in on:
            zz = al.SGN*pos[j, 2] + al.OFF
            if zz > cfg.zmax_show_A: continue
            _, gcj, _ = al.site_to_index(pos[j, 0], pos[j, 1], pos[j, 2])
            xin = (gcj - c) * dx
            if Z[j] == 22:
                mk = ax.scatter([xin], [zz], s=60, marker="x", c="#39ff14", linewidths=1.8)
            else:
                mk = ax.scatter([xin], [zz], s=60, marker="D",
                                facecolors="none", edgecolors="#ff21ff", linewidths=1.6)
            mk.set_path_effects([pe.withStroke(linewidth=2.2, foreground="black")])
        ax.set_xlim(-W*dx, W*dx); ax.set_ylim(cfg.zmax_show_A, zrec[0])
        ax.set_xlabel("in-plane x (A)"); ax.set_ylabel("depth z (A)"); ax.set_title(ttl)
    # panel 3: blind v2 atoms on this column (species-coloured, z error bars) vs GT ticks
    ax = axes[2]
    fx = found["col"]*dx; fy = found["row"]*dx
    cx0, cy0 = (cfg.X0 + c*dx - cfg.X0), (cfg.Y0 + r*dx - cfg.Y0)   # column centre in A offsets
    oncol = np.where(np.hypot(fx - c*dx, fy - r*dx) < 0.8)[0]
    for i in oncol:
        zz = found["z_A"][i]
        sp = int(found["species"][i])
        ax.errorbar([0.5], [zz], yerr=[found["sz_A"][i]], fmt="o", ms=6,
                    color=SPCOL.get(sp, "w"), ecolor=SPCOL.get(sp, "w"),
                    capsize=2, elinewidth=1.2)
    for zz_lab, cc in [("Pb", 82), ("Ti", 22), ("O", 8)]:
        ax.plot([], [], "o", color=SPCOL[cc], label=f"found {zz_lab}")
    for j in on:
        zz = al.SGN*pos[j, 2] + al.OFF
        if zz > cfg.zmax_show_A: continue
        col = "#39aa14" if Z[j] == 22 else "#c000c0"
        ax.plot([-0.1], [zz], marker="_", color=col, ms=16, mew=2.5)
    ax.plot([], [], "_", color="#39aa14", label="GT Ti"); ax.plot([], [], "_", color="#c000c0", label="GT O")
    ax.set_ylim(cfg.zmax_show_A, zrec[0]); ax.set_xlim(-0.25, 1.0); ax.set_xticks([])
    ax.set_title("blind v2 atoms (colour=species,\nbar=z sigma) vs GT ticks")
    ax.legend(fontsize=8, loc="lower right")
    fig.suptitle("Brightest B-O column -- 3-D tube-CLEAN atoms + error bars vs ground truth")
    fig.tight_layout(); fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)


def fig_zaccuracy(methods, cfg, path):
    """z-error histograms across finder methods (dict name->match-dict): the payoff of
    3-D fitting over 1-D spike over raw peak-picking, in both count matched and z-RMS."""
    colours = {"v3 lattice+guided": "tab:green", "v1 spike": "tab:blue", "raw peak-pick": "gray"}
    fig, ax = plt.subplots(figsize=(8, 5))
    for lab, m in methods.items():
        dz = m["match_dz"][m["match_gi"] >= 0]
        dz = dz[np.isfinite(dz)]
        rms = np.sqrt(np.mean(dz**2)) if len(dz) else np.nan
        ax.hist(dz, bins=np.linspace(-cfg.match_tol_z_A, cfg.match_tol_z_A, 41),
                histtype="stepfilled", alpha=0.45, color=colours.get(lab, None),
                label=f"{lab}  (n={len(dz)}, z-RMS={rms:.2f} A)")
    ax.axvline(0, color="k", lw=0.8, ls=":")
    ax.set_xlabel("z error  (found - GT)  [A]"); ax.set_ylabel("matched atoms")
    ax.set_title(f"Depth-localisation error vs ground truth (dz={cfg.dz:.2f} A/layer)")
    ax.legend()
    fig.tight_layout(); fig.savefig(path, dpi=130, bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", default="NL70_coherent")
    ap.add_argument("--psf", default="auto",
                    choices=["auto", "data", "synthetic", "empirical", "all"])
    ap.add_argument("--single-atom-vol", default=None)
    ap.add_argument("--dose", type=float, default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--n-null", type=int, default=400)
    args = ap.parse_args()

    cfg = config.preset(args.preset)
    if args.single_atom_vol: cfg.single_atom_vol = args.single_atom_vol
    if args.dose is not None: cfg.dose_e_per_A2 = args.dose
    if args.out: cfg.out_dir = args.out
    os.makedirs(cfg.out_dir, exist_ok=True)

    print(f"[atomfind] preset={cfg.name}  vol={cfg.recon_vol}  dose={cfg.dose_e_per_A2}")
    V, dx = align.load_phase(cfg)
    pos, Z = align.load_gt(cfg)
    al = align.register(V, dx, pos, Z, cfg)
    print("[align] " + align.summarize(al))

    psfs = psfmod.build_psfs(V, dx, al, pos, Z, cfg)
    print(f"[psf] built {[n for n in psfs if not n.startswith('_')]}  default={psfs['_default']}"
          f"  (data averaged {psfs.get('_data_n','-')} Pb blobs)")
    fig_psf(psfs, cfg, dx, os.path.join(cfg.out_dir, "psf_compare.png"))

    default_psf = psfs[psfs["_default"]]
    dec, dinfo = deconv.richardson_lucy_3d(V, default_psf, cfg)
    np.save(os.path.join(cfg.out_dir, "deconvolved_vol.npy"), dec)
    print(f"[deconv] RL {dinfo}")

    olab = validate.classify_oxygen(pos, Z, cfg)

    # ---- BLIND atom finder (no GT input) ---------------------------------
    # v3 (default): preprocess -> 3-D tube CLEAN -> lattice-aware species (column typing +
    # comb parity) -> guided re-detection at empty comb slots (prior = the column's OWN
    # atoms; guided atoms tagged) -> xyz + 3-axis error bars + species -> ASE object.
    kernels = psfmod.species_kernels(cfg, dx)
    found, seeds = find.find_atoms_v3(V, cfg, dx, kernels)
    # refine the recon<->GT map's affine on matched heavy atoms (fiducial refinement --
    # calibration only, the finder never sees the map). Kills the ~1.6 px peak-vs-fit
    # convention bias + the 0.6% dx scale residual.
    al = align.refine_with_atoms(al, found, pos, Z, cfg)
    print("[align] refined: " + align.summarize(al))
    frep, m_v3 = validate.finder_report(found, pos, Z, al, cfg, olabel=olab)
    blind = found[found["guided"] == 0]
    frep_blind, _ = validate.finder_report(blind, pos, Z, al, cfg, olabel=olab)
    np.save(os.path.join(cfg.out_dir, "found_atoms.npy"), found)
    csv_p, xyz_p = find.export_atoms(found, al, cfg, os.path.join(cfg.out_dir, "found_atoms"))
    print(f"[find] v3: {len(seeds)} columns -> {frep['n_found']} atoms "
          f"({(found['guided']==1).sum()} guided), prec={frep['precision']:.2f}, "
          f"xy-RMS={frep.get('xy_rms_A', float('nan')):.2f}A, z-RMS={frep['z_rms_A']:.2f}A  "
          f"(Pb {frep['Pb']['recall']:.0%} / Ti {frep['Ti']['recall']:.0%} / O {frep['O']['recall']:.0%})")
    print(f"[export] {os.path.basename(csv_p)} + {os.path.basename(xyz_p)} (ASE object)")

    found_spike, _, _ = find.find_atoms(V, cfg, dx, default_psf, method="spike")
    found_raw, _, _ = find.find_atoms(V, cfg, dx, default_psf, method="raw")
    frep_spike, m_spike = validate.finder_report(found_spike, pos, Z, al, cfg, olabel=olab)
    frep_raw, m_raw = validate.finder_report(found_raw, pos, Z, al, cfg, olabel=olab)
    fig_zaccuracy({"v3 lattice+guided": m_v3, "v1 spike": m_spike, "raw peak-pick": m_raw},
                  cfg, os.path.join(cfg.out_dir, "z_accuracy.png"))

    null = make_null(pos, cfg, n=args.n_null)

    if args.psf == "all":
        which = ["empirical", "data", "synthetic"]
    elif args.psf == "auto":
        which = [psfs["_default"]]
    else:
        which = [args.psf]
    which = [w for w in which if w in psfs] or [psfs["_default"]]
    reports = {}
    rec_default = None
    for w in which:
        rec = fit.fit_amplitudes(V, dx, al, pos, Z, psfs[w], cfg, extra_sites=null)
        olabel = np.zeros(len(rec), int)
        m = rec["idx"] >= 0; olabel[m] = olab[rec["idx"][m]]
        rep = validate.detection_report(rec, olabel)
        rep["amp_summary"] = validate.amplitude_summary(rec)
        reports[w] = rep
        if w == psfs["_default"] or rec_default is None:
            rec_default, olabel_default = rec, olabel
        print(f"[fit:{w}] tau={rep['tau']:.3f}  " +
              "  ".join(f"{k}:TPR={rep[k]['tpr']:.2f}/AUC={rep[k]['auc']:.2f}"
                        for k in ("Pb", "Ti", "O_all", "O_axial_overlap", "O_inplane_isolated") if k in rep))

    fig_amplitude(rec_default, olabel_default, cfg, os.path.join(cfg.out_dir, "amplitude_vs_Z.png"))
    fig_roc(rec_default, olabel_default, os.path.join(cfg.out_dir, "roc_oxygen.png"))
    fig_overlay(V, dec, found, al, pos, Z, cfg, dx, os.path.join(cfg.out_dir, "detection_overlay.png"))

    report = dict(preset=cfg.name, vol=cfg.recon_vol, dose=cfg.dose_e_per_A2,
                  dx=dx, dz=cfg.dz,
                  align=dict(SGN=al.SGN, OFF=al.OFF, CAL_X=al.CAL_X, CAL_Y=al.CAL_Y,
                             mX=al.mX, bX=al.bX, mY=al.mY, bY=al.bY,
                             corr_depth=al.corr_depth),
                  psf=dict(default=psfs["_default"], data_n=psfs.get("_data_n"),
                           fwhm={n: list(map(float, psfmod.measure_fwhm(psfs[n], cfg, dx)))
                                 for n in psfs if not n.startswith("_")}),
                  deconv=dinfo,
                  finder=dict(v3=frep, v3_blind_only=frep_blind, v1_spike=frep_spike,
                              raw=frep_raw),
                  reports=reports)
    with open(os.path.join(cfg.out_dir, "report.json"), "w") as f:
        json.dump(report, f, indent=2, default=float)

    _print_verdict(reports, frep, frep_blind, frep_spike, frep_raw, psfs, cfg)
    print(f"\n[atomfind] wrote figures + report.json + found_atoms.{{csv,extxyz,npy}} to {cfg.out_dir}")


def _print_verdict(reports, frep, frep_blind, frep_spike, frep_raw, psfs, cfg):
    d = reports.get(psfs["_default"], next(iter(reports.values())))
    print("\n" + "=" * 78)
    print(f"VERDICT  ({cfg.name}, dose={cfg.dose_e_per_A2},  PSF={psfs['_default']})")
    print("=" * 78)
    print("BLIND FINDER v3 (tube CLEAN + lattice species + guided; no GT anywhere):")
    print(f"  {'':16s}   raw   v1-spike   v3-blind   v3-full   v3-BULK   z-RMS")
    for k, lab in [("Pb", "Pb (Z=82)"), ("Ti", "Ti (Z=22)"), ("O", "O (Z=8)")]:
        rr = f"{frep_raw[k]['recall']:.0%}"; rs = f"{frep_spike[k]['recall']:.0%}"
        rb0 = f"{frep_blind[k]['recall']:.0%}"; rv = f"{frep[k]['recall']:.0%}"
        rb = f"{frep[k]['recall_bulk']:.0%}"
        print(f"  {lab:16s} {rr:>5} {rs:>8} {rb0:>9} {rv:>8} {rb:>8}   {frep[k]['z_rms_A']:.2f}A")
    print(f"  overall precision  raw {frep_raw['precision']:.0%} / v1 {frep_spike['precision']:.0%}"
          f" / v3 {frep['precision']:.0%}    xy-RMS {frep.get('xy_rms_A', float('nan')):.2f}A"
          f"    z-RMS raw {frep_raw['z_rms_A']:.2f} -> v3 {frep['z_rms_A']:.2f}A")
    if "sigma_coverage_1s" in frep:
        cv = frep["sigma_coverage_1s"]
        pz = "  ".join(f"{nm} {frep[nm]['z_cov_1s']:.0%}" for nm in ("Pb", "Ti", "O")
                       if "z_cov_1s" in frep.get(nm, {}))
        print(f"  error-bar 1-sigma COVERAGE (want ~68%):  x={cv['x']:.0%} y={cv['y']:.0%} "
              f"z={cv['z']:.0%}   per-species z: {pz}")
    if "confusion" in frep:
        cf = frep["confusion"]
        offd = sum(cf[f"{a}->{b}"] for a in (82, 22, 8) for b in (82, 22, 8) if a != b)
        diag = sum(cf[f"{a}->{a}"] for a in (82, 22, 8))
        print(f"  species confusion (found->matched GT), off-diagonal {offd}/{offd+diag} = {offd/(offd+diag):.1%}:")
        for pz, nm in [(82, "found Pb"), (22, "found Ti"), (8, "found O ")]:
            print(f"    {nm}: ->Pb {cf[f'{pz}->82']:4d}  ->Ti {cf[f'{pz}->22']:4d}  "
                  f"->O {cf[f'{pz}->8']:4d}  ->none {cf[f'{pz}->none']:4d}")
    print("\nGT-SEEDED O AMPLITUDE DETECTOR (calibrated contrast, vs off-lattice null):")
    for k, lab in [("O_all", "O all"), ("O_inplane_isolated", "O in-plane isolated"),
                   ("O_axial_overlap", "O axial-overlap")]:
        if k in d:
            r = d[k]
            print(f"  {lab:22s}  n={r['n']:4d}  detect(FPR5%)={r['tpr']:.0%}  AUC={r['auc']:.2f}")


if __name__ == "__main__":
    main()
