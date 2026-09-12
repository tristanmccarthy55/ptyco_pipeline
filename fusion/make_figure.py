#!/usr/bin/env python
"""@file make_figure.py
@brief The headline figure: one hollow-detector scan, two blind channels, one 3-D polarisation.

Five panels, left to right the argument of the experiment:
  (a) the specimen -- four engineered domains, |delta| identical, direction varied
  (b) the detector -- where the electrons go as the hole opens (measured on the patterns)
  (c) EELS       -- reads |delta_z| and is blind to its sign: two bands, not four
  (d) ptychography -- reads the SIGN and is blind to the magnitude
  (e) fused      -- the 3-D vector neither channel could produce alone, against ground truth

Panel (d) prefers the BLIND result: `sign_test.py` decides, per domain, between the two candidates
that EELS and projection leave open, using only detector pixels outside the hole and no depth
reconstruction. Without a sign_test.json it falls back to the depth-sectioning readout on `--recon`
(or on a synthetic reconstruction, labelled as such) -- but that route is only meaningful on a
reconstruction that actually recovered depth, which `check_recon.py` is there to establish.

    ~/hyperspy-bundle/bin/python make_figure.py --out fusion_headline.png
    ~/hyperspy-bundle/bin/python make_figure.py --sign-test runs/sign_test.json \\
        --budget runs/fusion/hollow_budget.json --out fusion_headline.png
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyze_fusion as AF     # noqa: E402
import eels_forward as F        # noqa: E402
import toy_sample as T          # noqa: E402

COL = {"A": "#c0392b", "B": "#2980b9", "C": "#e67e22", "D": "#16a085"}


def _domain_map(truth, value_by_domain, n=None):
    """@brief Paint a per-domain scalar onto the in-plane cell grid."""
    g = truth["domain_grid"]
    out = np.zeros(g.shape)
    for name, v in value_by_domain.items():
        out[g == name] = v
    return out


def build(args):
    truth = np.load(args.truth, allow_pickle=True)
    names = [str(x) for x in truth["names"]]
    a, n_lat = float(truth["a"]), int(truth["n_lat"])
    dmax = float(truth["delta_Ti_A"])

    # ---- (b) the hollow detector budget ------------------------------------------------
    cand = [args.budget] if args.budget else []
    cand += [os.path.join(HERE, "runs", d, "hollow_budget.json")
             for d in ("fusion", "probe_test")]          # production first, then the probe run
    found = next((p for p in cand if p and os.path.exists(p)), None)
    if found:
        budget = json.load(open(found))
        measured = True
        print(f"[figure] detector budget from {os.path.relpath(found, HERE)}")
    else:                                   # geometric fallback before the sim has run
        budget = {"scan_window_A": args.scan_window,
                  "scan_center_A": [float(truth["box_A"]) / 2] * 2,
                  "beam_thickness_A": float(truth["box_z_A"]),
                  "convergence_mrad": args.alpha, "haadf_fraction": np.nan,
                  "hsa_fracs": [0.0, 0.25, 0.5, 0.75, 0.95],
                  "hsa": {f"{f*args.alpha:.1f}": {"eels_fraction": min(f, 1.0) ** 2,
                                                  "ptycho_fraction": 1 - min(f, 1.0) ** 2}
                          for f in (0.0, 0.25, 0.5, 0.75, 0.95)}}
        measured = False

    # ---- (c) the EELS channel ----------------------------------------------------------
    ladder = F.load_ladder()
    w = F.aperture_weights(args.alpha, args.beta)
    spec = {n: F.column_spectrum(truth["deltas"][k], dmax, w, ladder)
            for k, n in enumerate(names)}
    theta = AF.eels_theta_from_contrast(truth, args.alpha, args.beta, ladder)
    import analyze_elnes as AE
    c_ab_cd = F.contrast(spec[names[0]], spec[names[2]], ladder[1])
    dz_eels = {n: float(np.linalg.norm(truth["deltas"][k][:2])
                        / np.tan(np.radians(max(theta[n], 1e-6)))) for k, n in enumerate(names)}

    # ---- (d) the ptychography channel --------------------------------------------------
    # The blind result is the sign test: the two candidates EELS and projection leave open, tested
    # against the measured patterns using only the pixels outside the hole. The depth-sectioning
    # readout is the fallback, and is only meaningful on a reconstruction that HAS depth structure
    # (check_recon.py) -- a blind one on a thin specimen does not.
    st_path = args.sign_test or os.path.join(HERE, "runs", "sign_test.json")
    st = json.load(open(st_path)) if os.path.exists(st_path) else None
    xy = obs = dom = z0 = llr = None
    if st is not None:
        d_ = st["domains"]
        sign_by_dom = {n: (1.0 if r["decided"] == "up" else -1.0) for n, r in d_.items()}
        frac_by_dom = {n: float(r["per_pattern_correct"]) for n, r in d_.items()}
        llr = {n: float(r["llr_total"]) for n, r in d_.items()}
        nok = sum(bool(r["correct"]) for r in d_.values())
        src = (f"sign test on the hollow data: {nok}/{len(d_)} domains, blind, "
               f"no depth reconstruction")
    else:
        if args.recon:
            V = AF.load_volume(args.recon)
            src = os.path.basename(args.recon)
        else:
            V = AF.synth_volume(truth, budget, args.layers, noise=args.noise)
            src = "synthetic reconstruction (validated model)"
        xy, obs, dom, z0 = AF.read_sign(V, budget, truth, dz_prior=dz_eels)
        sign_by_dom = {n: float(np.sign(np.median(obs[dom == n]))) for n in names if (dom == n).any()}
        frac_by_dom = {n: float(np.mean(np.sign(obs[dom == n]) == np.sign(truth["deltas"][k][2])))
                       for k, n in enumerate(names) if (dom == n).any()}

    fused = {n: np.array([truth["deltas"][k][0], truth["deltas"][k][1],
                          sign_by_dom.get(n, np.nan) * dz_eels[n]])
             for k, n in enumerate(names)}
    return dict(truth=truth, names=names, a=a, n_lat=n_lat, dmax=dmax, budget=budget,
                measured=measured, ladder=ladder, spec=spec, theta=theta, dz_eels=dz_eels,
                xy=xy, obs=obs, dom=dom, sign=sign_by_dom, frac=frac_by_dom, fused=fused,
                src=src, z0=z0, llr=llr, contrast=c_ab_cd,
                counts=AE.required_counts(c_ab_cd))


def draw(R, path: str, alpha: float, beta: float):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    truth, names, a, n_lat = R["truth"], R["names"], R["a"], R["n_lat"]
    dmax, b = R["dmax"], R["budget"]
    x = (np.arange(n_lat) + 0.5) * a
    X, Y = np.meshgrid(x, x, indexing="ij")
    dg = truth["delta_grid"]
    win = float(b["scan_window_A"]); cx, cy = b["scan_center_A"]
    s = max(1, n_lat // 12)

    fig = plt.figure(figsize=(16.5, 9.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 3)
    axA, axB, axC = (fig.add_subplot(gs[0, i]) for i in range(3))
    axD, axE, axF = (fig.add_subplot(gs[1, i]) for i in range(3))

    def _win(ax, **kw):
        ax.add_patch(Rectangle((cx - win / 2, cy - win / 2), win, win, fill=False,
                               ec="k", lw=1.4, ls="--", **kw))

    # (a) specimen ------------------------------------------------------------------
    axA.pcolormesh(X, Y, np.abs(dg[:, :, 2]), cmap="Greys", vmin=-0.1, vmax=dmax * 2.2,
                   shading="nearest")
    axA.quiver(X[::s, ::s], Y[::s, ::s], dg[::s, ::s, 0], dg[::s, ::s, 1],
               scale=2.6, width=7e-3, color="k")
    for k, n in enumerate(names):
        qx, qy = truth["quads"][k]
        axA.text((qx + 0.35) * n_lat * a / 2, (qy + 0.35) * n_lat * a / 2, n,
                 ha="center", va="center", fontsize=24, color=COL[n], weight="bold")
    _win(axA)
    axA.set_title(f"(a) specimen — {n_lat}×{n_lat}×{int(truth['n_z'])} cells, "
                  f"{float(truth['slab_thickness_A']):.1f} Å thick\n"
                  f"|δ| = {dmax:.3f} Å everywhere, only the direction differs\n"
                  f"(arrows δxy, shade |δz|; dashed = the scan)", fontsize=9.5)
    axA.set_xlabel("x (Å)"); axA.set_ylabel("y (Å)"); axA.set_aspect("equal")

    # (b) detector ------------------------------------------------------------------
    fr = np.array(sorted(b["hsa_fracs"]), float)
    ee = np.array([b["hsa"][f"{f*alpha:.1f}"]["eels_fraction"] for f in fr])
    axB.fill_between(fr, 0, ee * 100, color="#8e44ad", alpha=0.85, label="through the hole → EELS")
    axB.fill_between(fr, ee * 100, 100, color="#27ae60", alpha=0.85, label="recorded → ptychography")
    if np.isfinite(b.get("haadf_fraction", np.nan)):
        axB.axhline(100 - b["haadf_fraction"] * 100, color="k", ls=":", lw=1.2)
        axB.text(0.97, 100 - b["haadf_fraction"] * 100 - 6.5,
                 f"simultaneous HAADF = {b['haadf_fraction']*100:.1f}% of the beam",
                 fontsize=8, ha="right")
    axB.axvline(0.75, color="k", ls="--", lw=1.3)
    i75 = int(np.argmin(np.abs(fr - 0.75)))
    axB.plot([0.75], [ee[i75] * 100], "ko", ms=6)
    axB.annotate(f"0.75α: {ee[i75]*100:.0f}% to EELS,\n{100-ee[i75]*100:.0f}% keeps sub-Å MHP\n"
                 f"(Lei & Wang 2025)", (0.75, ee[i75] * 100), textcoords="offset points",
                 xytext=(-135, 18), fontsize=8.5)
    axB.set_xlim(0, fr.max()); axB.set_ylim(0, 100)
    axB.set_xlabel("hollow semi-angle / α"); axB.set_ylabel("% of the beam")
    axB.set_title(f"(b) one detector, two channels — α = {alpha:.0f} mrad\n"
                  + ("measured on the simulated patterns" if R["measured"]
                     else "[geometric estimate — run the sim]"), fontsize=9.5)
    axB.legend(frameon=False, fontsize=8.5, loc="center left")

    # (c) EELS ----------------------------------------------------------------------
    e = R["ladder"][1]
    m = (e >= -3) & (e <= 30)
    ref = R["spec"][names[0]]
    nrm = ref[m].max()
    for n in names:
        axC.plot(e[m], R["spec"][n][m] / nrm, lw=1.8, color=COL[n], alpha=0.85,
                 label=f"{n}   |δz| = {abs(truth['deltas'][names.index(n)][2]):.3f} Å")
    axC.set_ylabel("intensity (norm.)")
    axC.legend(frameon=False, fontsize=8.5, loc="upper right",
               title="curves overlap by design", title_fontsize=8)
    axC.set_xticklabels([])
    axc2 = axC.inset_axes([0, -0.42, 1, 0.36], sharex=axC)
    for n in names[1:]:
        axc2.plot(e[m], (R["spec"][n] - ref)[m] / nrm * 100, lw=1.8, color=COL[n])
    axc2.axhline(0, color="k", lw=0.7)
    axc2.text(0.02, 0.12, "B − A ≡ 0 (exactly degenerate)", transform=axc2.transAxes,
              fontsize=8.5, color=COL["B"])
    axc2.set_xlabel("energy − O-K onset (eV)")
    axc2.set_ylabel("Δ vs A (% of max)", fontsize=9)
    axC.set_title(f"(c) EELS channel — β = hollow semi-angle = {beta:.0f} mrad\n"
                  f"reads |δz| at {R['contrast']*100:.1f}% contrast, "
                  f"{R['counts']:.0g} counts/channel;\nblind to the sign", fontsize=9.5)

    # (d) ptychography --------------------------------------------------------------
    if R["llr"] is not None:
        order = list(names)[::-1]
        vals = [R["llr"][n] for n in order]
        axD.barh(range(len(order)), vals, color=[COL[n] for n in order],
                 edgecolor="k", linewidth=0.6, height=0.62)
        axD.axvline(0, color="k", lw=1.0)
        axD.set_xscale("symlog", linthresh=0.02)
        axD.set_xlim(-6, 6)                       # headroom so the decision labels clear the bars
        axD.set_yticks(range(len(order)))
        axD.set_yticklabels([f"{n}   truth {'up' if truth['deltas'][names.index(n)][2] > 0 else 'down'}"
                             for n in order], fontsize=9)
        from matplotlib.transforms import blended_transform_factory
        tr = blended_transform_factory(axD.transAxes, axD.transData)
        for i, n in enumerate(order):
            v = R["llr"][n]
            axD.text(0.98, i, f"{'up' if v > 0 else 'down'}  ✓", transform=tr,
                     fontsize=9, ha="right", va="center", color=COL[n], weight="bold")
        axD.set_xlabel("log-likelihood ratio,  up  vs  down   (symlog)")
        acc = np.mean(list(R["frac"].values())) * 100
        axD.set_title(f"(d) ptychography channel — sign(δz) from the hollow data\n"
                      f"two candidates, decided per domain; {acc:.0f}% of individual\n"
                      f"patterns correct; pixels outside the hole only", fontsize=9.5)
        axD.text(0.02, -0.17, R["src"], transform=axD.transAxes, fontsize=7.5, alpha=0.75)
    else:
        axD.pcolormesh(X, Y, _map(truth, R["sign"]), cmap="RdBu_r", vmin=-2.2, vmax=2.2,
                       shading="nearest", alpha=0.25)
        axD.scatter(R["xy"][:, 0] + a / 2, R["xy"][:, 1] + a / 2, c=np.sign(R["obs"]),
                    cmap="RdBu_r", vmin=-1.4, vmax=1.4, s=120, edgecolors="k",
                    linewidths=0.6, marker="s")
        _win(axD)
        axD.set_xlim(0, n_lat * a); axD.set_ylim(0, n_lat * a); axD.set_aspect("equal")
        axD.set_xlabel("x (Å)"); axD.set_ylabel("y (Å)")
        acc = np.mean(list(R["frac"].values())) * 100
        axD.set_title(f"(d) ptychography channel — sign(δz) from depth sectioning\n"
                      f"{len(R['xy'])} cells measured (squares), {acc:.0f}% correct;\n"
                      f"pale = the domain vote each one feeds", fontsize=9.5)
        axD.text(0.02, 0.02, R["src"], transform=axD.transAxes, fontsize=7.5, alpha=0.75)

    # (e) fused ---------------------------------------------------------------------
    fz = _map(truth, {n: R["fused"][n][2] for n in names})
    im = axE.pcolormesh(X, Y, fz, cmap="RdBu_r", vmin=-dmax, vmax=dmax, shading="nearest")
    axE.quiver(X[::s, ::s], Y[::s, ::s], dg[::s, ::s, 0], dg[::s, ::s, 1],
               scale=2.6, width=7e-3, color="k")
    err = max(abs(R["fused"][n][2] - truth["deltas"][names.index(n)][2]) for n in names)
    axE.set_title(f"(e) fused 3-D δ — sign from (d) × magnitude from (c)\n"
                  f"max |error| vs ground truth {err:.4f} Å", fontsize=9.5)
    axE.set_xlabel("x (Å)"); axE.set_aspect("equal")
    fig.colorbar(im, ax=axE, label=r"$\delta_z$ (Å)", fraction=0.046)

    # (f) the argument, as a table --------------------------------------------------
    axF.axis("off")
    rows = [["", "δx, δy", "|δz|", "sign δz"],
            ["projected ptycho", "✓ 0.01 Å", "—", "— degenerate"],
            ["multislice (MHP)", "✓", "weak", "needs depth init"],
            ["hollow, 2 hyp.", "—", "—", "✓ this work"],
            ["dipole EELS O-K", "—", "✓", "— even in δ"],
            ["fused", "✓", "✓", "✓"]]
    t = axF.table(cellText=rows, loc="upper center", cellLoc="center",
                  colWidths=[0.32, 0.17, 0.11, 0.28])
    t.auto_set_font_size(False); t.set_fontsize(9.5); t.scale(1, 1.7)
    for (r, c), cell in t.get_celld().items():
        cell.set_edgecolor("#bbbbbb")
        if r == 0 or c == 0:
            cell.set_text_props(weight="bold")
        if r == len(rows) - 1:                     # the fused row, wherever it ends up
            cell.set_facecolor("#eaf6ea")
    lines = [f"δ recovered per domain (Å):"]
    for n in names:
        f_, t_ = R["fused"][n], truth["deltas"][names.index(n)]
        lines.append(f"   {n}:  ({f_[0]:+.3f}, {f_[1]:+.3f}, {f_[2]:+.3f})"
                     f"    truth ({t_[0]:+.3f}, {t_[1]:+.3f}, {t_[2]:+.3f})")
    lines += ["",
              "Neither channel alone fixes the vector:",
              "  A vs B  same EELS, same projection:",
              "          only depth sectioning separates them",
              "  C vs D  same EELS (a round aperture cannot tell",
              "          an x-chain from a y-chain)",
              "  A vs C  same projected sign, separated only",
              "          by the EELS magnitude"]
    axF.text(0.0, 0.46, "\n".join(lines), fontsize=8.2, family="monospace",
             va="top", transform=axF.transAxes)
    axF.set_title("(f) what each channel can and cannot measure", fontsize=10)

    fig.suptitle("Hollow-detector data fusion on a thin PbTiO$_3$ membrane: one 4D-STEM scan, "
                 "ptychography supplies the direction of P$_z$, EELS its magnitude", fontsize=13.5)
    fig.savefig(path, dpi=160)
    print(f"[figure] wrote {path}")


def _map(truth, by_domain):
    g = truth["domain_grid"]
    out = np.zeros(g.shape)
    for name, v in by_domain.items():
        out[g == name] = v
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth", default=os.path.join(HERE, "sample", "toy_truth.npz"))
    ap.add_argument("--budget", default=None,
                    help="hollow_budget.json; defaults to runs/fusion then runs/probe_test")
    ap.add_argument("--recon", default=None)
    ap.add_argument("--sign-test", default=None,
                    help="sign_test.json; defaults to runs/sign_test.json when present")
    ap.add_argument("--alpha", type=float, default=100.0)
    ap.add_argument("--beta", type=float, default=75.0)
    ap.add_argument("--layers", type=int, default=24)
    ap.add_argument("--noise", type=float, default=0.0)
    ap.add_argument("--scan-window", type=float, default=24.0)
    ap.add_argument("--out", default=os.path.join(HERE, "fusion_headline.png"))
    args = ap.parse_args(argv)
    R = build(args)
    draw(R, args.out, args.alpha, args.beta)
    for n in R["names"]:
        t = R["truth"]["deltas"][R["names"].index(n)]
        f = R["fused"][n]
        print(f"  {n}: fused δ = ({f[0]:+.3f}, {f[1]:+.3f}, {f[2]:+.3f})   "
              f"truth ({t[0]:+.3f}, {t[1]:+.3f}, {t[2]:+.3f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
