#!/usr/bin/env python
"""@file eels_forward.py
@brief Per-column O-K spectra for the toy membrane: which polarisation the EELS channel can read.

The spectroscopy half of the fusion. Each O in a perovskite sits on a linear Ti-O-Ti chain, and
its K edge is dichroic about THAT chain axis -- a lattice direction, not the polarisation. What
the polar displacement does is change the two Ti-O bond lengths along the chain, which changes the
size of the sigma*/pi* dichroism. That dependence is exactly the M5 displacement ladder
(`eels/runs/scan_{0.00..0.75}_Oap` + `tet_Pz_Oap` at s = 1), so one existing set of CASTEP+OptaDOS
runs determines the whole model with no new DFT:

    sigma_site(E) = S_par(s) <cos^2 (q, c_chain)>  +  S_perp(s) (1 - <cos^2 (q, c_chain)>)
    s = |delta . c_chain| / |delta|_max     (the along-chain component of the polar displacement)

The aperture average <cos^2> is the M6c/M6d convergent-probe model generalised to an arbitrary
axis: q_t = k(theta_f - theta_i) is a VECTOR difference, so the convergence cone alpha spreads q
on its own and the collection cone beta (= the hollow semi-angle, since the hole is what feeds the
spectrometer) adds to it. With the beam along z, <cos^2(q,z)> is small and <cos^2(q,x)> =
<cos^2(q,y)> = (1 - <cos^2(q,z)>)/2.

Summing the three O per formula unit gives the measurable column spectrum. The contrast between
domains is the fusion's magnitude channel; it is SIGN-BLIND by construction (s depends on |delta.c|),
which is what ptychography has to supply.

    ~/hyperspy-bundle/bin/python eels_forward.py --cache    # extract the M5 ladder once (slow)
    ~/hyperspy-bundle/bin/python eels_forward.py            # domain spectra + contrast + dose
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "eels"))
import analyze_elnes as A     # noqa: E402
import config as C            # noqa: E402

LADDER = [0.00, 0.25, 0.50, 0.75, 1.00]
LADDER_SEED = {0.00: "scan_0.00_Oap", 0.25: "scan_0.25_Oap", 0.50: "scan_0.50_Oap",
               0.75: "scan_0.75_Oap", 1.00: "tet_Pz_Oap"}
CACHE = os.path.join(HERE, "sample", "elnes_ladder.npz")
WINDOW_eV = A.ANALYSIS_WINDOW_eV                     # (-2, 30) relative to the edge onset


# ---------------------------------------------------------------- the CASTEP ladder (cached)
def build_cache(runs_dir: str | None = None, out: str = CACHE) -> str:
    """@brief Extract the `:exc` block of every ladder spectrum into one small npz.

    The raw OptaDOS `_core_edge.dat` files are ~208 MB each (one block per atom); only the
    excited atom's block is the physical ELNES. Parsing all ten every run is minutes of I/O, so
    they are reduced once to (energy, S_par, S_perp) per ladder step on a common energy grid.
    """
    runs_dir = runs_dir or os.path.join(HERE, "..", "eels", "runs")
    e_ref, par, perp = None, [], []
    for s in LADDER:
        seed = LADDER_SEED[s]
        qc = os.path.join(runs_dir, f"{seed}.qc_core_edge.dat")
        qp = os.path.join(runs_dir, f"{seed}.qperp_core_edge.dat")
        for p in (qc, qp):
            if not os.path.exists(p):
                raise SystemExit(f"missing {p} -- the M5 ladder must be present (see eels/RESULTS.md)")
        e1, y1 = A.load_optados_core(qc)
        e2, y2 = A.load_optados_core(qp)
        if e_ref is None:
            e_ref = e1
        par.append(np.interp(e_ref, e1, y1))
        perp.append(np.interp(e_ref, e2, y2))
        print(f"[ladder] s={s:.2f}  {seed}: {len(e1)} pts, "
              f"max S_par {y1.max():.3g}  max S_perp {y2.max():.3g}")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out, s=np.array(LADDER), energy=e_ref,
             S_par=np.array(par), S_perp=np.array(perp))
    print(f"[ladder] wrote {out}")
    return out


def load_ladder(path: str = CACHE):
    """@brief (s, energy, S_par(s,E), S_perp(s,E)) for the O-K displacement ladder."""
    if not os.path.exists(path):
        raise SystemExit(f"{path} not found -- run `eels_forward.py --cache` first")
    d = np.load(path)
    return d["s"], d["energy"], d["S_par"], d["S_perp"]


def ladder_at(s_query, s, S_par, S_perp):
    """@brief Linearly interpolate the ladder to an arbitrary along-chain distortion s in [0,1]."""
    q = np.clip(np.atleast_1d(s_query), s.min(), s.max())
    par = np.stack([np.interp(q, s, S_par[:, i]) for i in range(S_par.shape[1])], axis=-1)
    perp = np.stack([np.interp(q, s, S_perp[:, i]) for i in range(S_perp.shape[1])], axis=-1)
    return par, perp


# ---------------------------------------------------------------- aperture average of cos^2(q,n)
def q_cos2(alpha_rad: float, beta_rad: float, theta_E: float, axis: str = "z",
           n: int = 400_000, seed: int = 0) -> float:
    """@brief Dipole-weighted <cos^2 angle(q, axis)> over the convergence AND collection cones.

    Generalises `analyze_elnes.parallel_weight` (axis = beam) to an arbitrary crystal axis, which
    is what a per-site model needs: the z-chain O is probed by q_z, the x- and y-chain O by the
    transverse q the convergence cone supplies. q = (q_t, q_z) with q_z = k theta_E and
    q_t = k(theta_f - theta_i); the dipole weight is |q.r|^2/q^4 -> component^2 / q^4.

    @param axis 'z' (along the beam) or 'x'/'y' (in-plane; equal by the round-cone symmetry)
    @return <cos^2> in [0,1]; the three axes sum to 1 by construction.
    """
    rng = np.random.default_rng(seed)
    ri = alpha_rad * np.sqrt(rng.random(n)); pi_ = 2 * np.pi * rng.random(n)
    rf = beta_rad * np.sqrt(rng.random(n)); pf = 2 * np.pi * rng.random(n)
    qx = rf * np.cos(pf) - ri * np.cos(pi_)
    qy = rf * np.sin(pf) - ri * np.sin(pi_)
    qz = np.full_like(qx, theta_E)                    # the minimum (on-axis) momentum transfer
    q2 = qx**2 + qy**2 + qz**2
    num = {"x": qx**2, "y": qy**2, "z": qz**2}[axis]
    return float(np.sum(num / q2**2) / np.sum(1.0 / q2))


def aperture_weights(alpha_mrad: float, beta_mrad: float, edge_eV: float = 532.0,
                     energy_keV: float = 300.0, n: int = 400_000) -> dict:
    """@brief <cos^2(q, axis)> for the three chain axes at this (alpha, beta), plus theta_E."""
    tE = A.characteristic_angle_rad(edge_eV, energy_keV)
    wz = q_cos2(alpha_mrad * 1e-3, beta_mrad * 1e-3, tE, "z", n=n)
    return {"z": wz, "x": 0.5 * (1 - wz), "y": 0.5 * (1 - wz), "theta_E_mrad": tE * 1e3}


# ---------------------------------------------------------------- per-column spectrum
CHAIN_VEC = {"z": np.array([0., 0., 1.]), "x": np.array([1., 0., 0.]), "y": np.array([0., 1., 0.])}


def column_spectrum(delta, delta_max: float, w: dict, ladder) -> np.ndarray:
    """@brief O-K spectrum of one column: the three O per formula unit, each on its own chain.

    @param delta      the column's Ti off-centring vector (A); only its along-chain components
                      matter, and only through |.| -- hence the sign blindness.
    @param delta_max  |delta| of the fully polarised reference (the s = 1 ladder step)
    @param w          aperture_weights() for this (alpha, beta)
    """
    s_grid, e, S_par, S_perp = ladder
    total = np.zeros_like(e)
    for ax in ("z", "x", "y"):
        s = abs(float(np.dot(delta, CHAIN_VEC[ax]))) / delta_max
        par, perp = ladder_at(s, s_grid, S_par, S_perp)
        total += w[ax] * par[0] + (1.0 - w[ax]) * perp[0]
    return total


def contrast(sa: np.ndarray, sb: np.ndarray, e: np.ndarray, window=WINDOW_eV) -> float:
    """@brief max|Sa - Sb| / max(Sa, Sb) over the near-edge window -- RESULTS.md's convention."""
    m = (e >= window[0]) & (e <= window[1])
    return float(np.abs(sa[m] - sb[m]).max() / max(sa[m].max(), sb[m].max()))


# ---------------------------------------------------------------- report
def report(truth, alpha_mrad: float, beta_mrad: float, ladder=None, quiet: bool = False) -> dict:
    """@brief Predicted domain spectra, the pairwise contrasts and the dose each one needs."""
    ladder = ladder or load_ladder()
    s_grid, e, S_par, S_perp = ladder
    w = aperture_weights(alpha_mrad, beta_mrad)
    names = [str(n) for n in truth["names"]]
    dmax = float(truth["delta_Ti_A"])
    spec = {n: column_spectrum(truth["deltas"][k], dmax, w, ladder) for k, n in enumerate(names)}

    if not quiet:
        print(f"\n== O-K forward model  (alpha = {alpha_mrad:.0f} mrad, "
              f"beta = hollow semi-angle = {beta_mrad:.0f} mrad, theta_E = "
              f"{w['theta_E_mrad']:.2f} mrad) ==")
        print(f"  aperture weights  <cos^2(q,z)> = {w['z']:.4f}   "
              f"<cos^2(q,x)> = <cos^2(q,y)> = {w['x']:.4f}")
        print(f"  (beam || z: the along-beam chain is probed by q_perp, so its sigma*/pi* "
              f"dichroism enters INVERTED and at ~1/3 strength -- the M6c/M6d result)")
        print(f"\n  {'domain':>7} {'s_z':>6} {'s_x':>6} {'s_y':>6}   along-chain distortion "
              f"of each Ti-O-Ti chain")
        for k, n in enumerate(names):
            d = truth["deltas"][k]
            print(f"  {n:>7} {abs(d[2])/dmax:>6.3f} {abs(d[0])/dmax:>6.3f} "
                  f"{abs(d[1])/dmax:>6.3f}")
        print(f"\n  {'pair':>9} {'contrast':>10} {'counts/chan (SNR 3)':>22}  verdict")
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                cij = contrast(spec[names[i]], spec[names[j]], e)
                nc = A.required_counts(cij) if cij > 1e-9 else float("inf")
                verd = ("DEGENERATE -> ptychography must supply it" if cij < 1e-6
                        else "measurable")
                print(f"  {names[i]}-{names[j]:>7} {cij*100:>9.3f}% "
                      f"{nc:>22.3g}  {verd}")
    return {"energy": e, "spectra": spec, "weights": w}


def sweep(truth, alpha_mrad: float = 100.0,
          fracs=(0.10, 0.25, 0.50, 0.60, 0.75, 0.90, 0.95, 1.00), ladder=None) -> dict:
    """@brief Does the hollow hole cost the spectroscopy anything? Contrast vs hollow semi-angle.

    The design question from Lei & Wang: how large a hole can we get away with. Their answer is a
    PTYCHOGRAPHY answer (lateral resolution holds to ~0.75 alpha, degrades past ~0.90). This is the
    other half: the hole is also the EELS collection aperture, so does opening it wash the
    polarisation signal out? It does not -- beyond alpha ~ 5 mrad the convergence cone alone
    already sets the momentum-transfer distribution (M6c), so <cos^2(q,z)> and every contrast
    derived from it are FLAT in beta. The hole radius is a pure dose splitter.
    """
    ladder = ladder or load_ladder()
    e = ladder[1]
    dmax = float(truth["delta_Ti_A"])
    names = [str(n) for n in truth["names"]]
    zhat, xhat = np.array([0., 0., 1.]), np.array([1., 0., 0.])
    print(f"\n== how large a hollow hole?  (alpha = {alpha_mrad:.0f} mrad) ==")
    print(f"  {'HSA/alpha':>10} {'beta (mrad)':>12} {'<cos2(q,z)>':>12} "
          f"{'AB-CD contrast':>15} {'counts/chan':>13} {'||/perp limit':>14}")
    out = []
    for f in fracs:
        b = alpha_mrad * f
        w = aperture_weights(alpha_mrad, b)
        cAB = contrast(column_spectrum(truth["deltas"][0], dmax, w, ladder),
                       column_spectrum(truth["deltas"][2], dmax, w, ladder), e)
        clim = contrast(column_spectrum(dmax * zhat, dmax, w, ladder),
                        column_spectrum(dmax * xhat, dmax, w, ladder), e)
        print(f"  {f:>10.2f} {b:>12.1f} {w['z']:>12.4f} {cAB*100:>14.2f}% "
              f"{A.required_counts(cAB):>13.2e} {clim*100:>13.2f}%")
        out.append((f, b, w["z"], cAB, clim))
    print(f"  -> contrast varies by <{(max(o[3] for o in out)/min(o[3] for o in out)-1)*100:.0f}% "
          f"across the whole range: the hole costs the spectroscopy nothing.")
    print(f"     Size it from the PTYCHOGRAPHY side alone (Lei & Wang: safe to 0.75 alpha, "
          f"0.90 at >=1e5 e/A^2).")
    print(f"  -> near-parallel reference (alpha = beta = 2 mrad): ||/perp limit "
          f"{contrast(column_spectrum(dmax*zhat, dmax, aperture_weights(2,2), ladder), column_spectrum(dmax*xhat, dmax, aperture_weights(2,2), ladder), e)*100:.2f}%"
          f"  -- WORSE than the 100 mrad probe, so the convergent probe is not the problem.")
    return {"sweep": np.array(out), "names": names}


def plot(res, truth, path: str, alpha_mrad: float, beta_mrad: float):
    """@brief The EELS panel of the headline figure: domain spectra + what they can/cannot separate."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    e, spec = res["energy"], res["spectra"]
    names = list(spec)
    m = (e >= WINDOW_eV[0] - 3) & (e <= WINDOW_eV[1])
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    col = {"A": "#c0392b", "B": "#2980b9", "C": "#e67e22", "D": "#16a085"}
    for n in names:
        ax[0].plot(e[m], spec[n][m], lw=1.6, color=col.get(n), label=f"domain {n}")
    ax[0].set_xlabel("energy relative to O-K onset (eV)")
    ax[0].set_ylabel("intensity (arb.)")
    ax[0].set_title(f"predicted O-K, α={alpha_mrad:.0f} / β={beta_mrad:.0f} mrad")
    ax[0].legend(frameon=False, fontsize=9)

    ref = spec[names[0]]
    for n in names[1:]:
        ax[1].plot(e[m], (spec[n] - ref)[m] / ref[m].max() * 100, lw=1.6, color=col.get(n),
                   label=f"{n} − {names[0]}")
    ax[1].axhline(0, color="k", lw=0.6)
    ax[1].set_xlabel("energy relative to O-K onset (eV)")
    ax[1].set_ylabel("difference (% of edge max)")
    ax[1].set_title("A−B is identically zero: EELS is sign-blind")
    ax[1].legend(frameon=False, fontsize=9)
    fig.savefig(path, dpi=150)
    print(f"[eels] wrote {path}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", action="store_true", help="extract the M5 ladder into sample/ (slow, once)")
    ap.add_argument("--runs-dir", default=None, help="where the OptaDOS *_core_edge.dat live")
    ap.add_argument("--truth", default=os.path.join(HERE, "sample", "toy_truth.npz"))
    ap.add_argument("--alpha", type=float, default=100.0, help="probe convergence semi-angle [mrad]")
    ap.add_argument("--beta", type=float, default=75.0,
                    help="EELS collection = the HOLLOW semi-angle [mrad] (default 0.75*alpha)")
    ap.add_argument("--plot", default=None, help="write the spectra figure here")
    ap.add_argument("--sweep", action="store_true",
                    help="contrast vs hollow semi-angle -- the 'how big a hole' answer")
    args = ap.parse_args(argv)

    if args.cache:
        build_cache(args.runs_dir)
        return 0
    truth = np.load(args.truth, allow_pickle=True)
    res = report(truth, args.alpha, args.beta)
    if args.sweep:
        sweep(truth, args.alpha)
    if args.plot:
        plot(res, truth, args.plot, args.alpha, args.beta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
