#!/usr/bin/env python
"""@file atomfind_sign.py
@brief sign(delta_z) from located atoms: atomfind on a reconstruction, then Ti against its oxygens.

The atom-level route to the sign. Every atom in the volume is located with atomfind's blind v3
finder, and each located Ti is referred to the oxygens around it along the beam. Up and down
domains differ only in the sign of that offset. Three references, from strictest to most lenient:

  O6 cage       z_Ti - mean z of the six nearest located O (`atomfind.polarisation.offcentring`)
  apical pair   z_Ti - mean z of the nearest located O above and below it in its own column
  equatorial O  z_Ti - mean z of the located O in the neighbouring columns at its height (>= 2)

All O move rigidly in the polar mode, so the three agree on the true structure (checked here on
the true atoms). Each located Ti is scored against ITS OWN true offset (nearest true Ti).

Error bars are atomfind's own: model sigma from the finder (joint CRLB + kernel mismatch), turned
into 95 % intervals by `atomfind.uncertainty` split-conformal calibration, with the matched true
atoms of the same volume as the calibration set -- exactly how atomfind calibrates on NL70.

Two kinds of volume go through the identical finder, kernels and analysis:

  ideal  the TRUE atoms rendered with the measured single-atom kernels (the forward model atomfind
         inverts) on the reconstruction's own grid and lattice registration, scaled to its phase
         level, plus white noise at a chosen fraction of its background level: a reconstruction
         with every atom at its true depth. The positive control -- NOT a reconstruction.
  recon  the real reconstruction.

Kernels are the NL70 measured Pb and Ti responses: the same 300 keV / 100 mrad / -20 A probe,
bin 4 and engine settings as the fusion reconstructions, at dz 0.999 against 0.976 here.

    ~/hyperspy-bundle/bin/python atomfind_sign.py --recon <Niter200.mat> --budget <hollow_budget.json> \\
        --ideal-noise 0 0.1 1 --cache <prefix> --reuse
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "analysis"))
import analyze_fusion as AF                                                  # noqa: E402
from atomfind import config as AFC, find as AFF, psf as AFP, polarisation as AFPOL  # noqa: E402
from atomfind import uncertainty as AFU                                      # noqa: E402

KERNEL_PB = os.path.expanduser("~/Desktop/psf_Pb_NL70_vol.npy")
KERNEL_TI = os.path.expanduser("~/Desktop/psf_Ti_NL70_vol.npy")
Z_OF = {"Pb": 82, "Ti": 22, "O": 8}
EDGE_A = 1.2                    # stay off the reconstruction rim when counting and scoring
ESTIMATORS = ("O6 cage", "apical pair", "equatorial O")


# ---------------------------------------------------------------- set-up
def finder_config(dx, dz, thickness_A, c_A):
    """@brief atomfind settings for this slab: whole depth analysed, no exit-band inflation.

    `bulk_z_A` only sets atomfind's UQ depth strata (entrance / bulk / exit); the finder ignores it.
    """
    return AFC.Config(name="fusion", recon_vol="", dx=dx, dz=dz, X0=None, Y0=None,
                      trim_z_A=(0.0, thickness_A), zmax_show_A=thickness_A,
                      exit_band_z_A=thickness_A, clean_max_atoms=16, comb_period_A=c_A,
                      bulk_z_A=(4.0, thickness_A - 4.0), convergence_mrad=100.0,
                      single_atom_vol=KERNEL_PB, single_atom_species=82, ti_kernel_vol=KERNEL_TI)


def context(recon, budget_path, truth_path, vasp):
    """@brief Everything both the analysis and the figures need about one reconstruction."""
    import ase.io
    import h5py
    truth = np.load(truth_path, allow_pickle=True)
    budget = json.load(open(budget_path))
    with h5py.File(recon, "r") as f:
        budget["dx_object_A"] = float(np.asarray(f["outputs"]["pixel_size"]).ravel()[0]) * 1e10
    V = AF.load_volume(recon)
    T = float(budget["beam_thickness_A"])
    dz = T / V.shape[0]
    frame = AF.lattice_frame(V, budget, truth)
    cfg = finder_config(frame[4], dz, T, float(truth["c"]))
    atoms = ase.io.read(vasp)
    return dict(recon=recon, truth=truth, budget=budget, V=V, T=T, dz=dz, frame=frame,
                dx=frame[4], cfg=cfg, kernels=AFP.species_kernels(cfg, frame[4]),
                pos_t=atoms.get_positions(), sym_t=np.array(atoms.get_chemical_symbols()))


def render_ideal(shape, frame, pos, sym, kernels, dz):
    """@brief True atoms x measured kernels on the reconstruction grid (O shares the Ti shape)."""
    from scipy.ndimage import shift as nd_shift
    x0, y0, ox, oy, dx = frame
    nL, Ny, Nx = shape
    K = {82: kernels[82], 22: kernels[22], 8: kernels[22]}
    hz, hxy = (K[82].shape[0] - 1) // 2, (K[82].shape[1] - 1) // 2
    V = np.zeros(shape)
    for (X, Y, Zc), s in zip(pos, sym):
        z = Z_OF[s]
        l, r, c = Zc / dz - 0.5, (Y + oy - y0) / dx, (X + ox - x0) / dx
        if not (-hz <= l < nL + hz and -hxy <= r < Ny + hxy and -hxy <= c < Nx + hxy):
            continue
        li, ri, ci = int(np.floor(l)), int(np.floor(r)), int(np.floor(c))
        k = nd_shift(K[z], (l - li, r - ri, c - ci), order=1, mode="constant") * z ** AF.Z_WEIGHT_EXP
        o = (li - hz, ri - hxy, ci - hxy)
        dst = tuple(slice(max(0, o[i]), min(shape[i], o[i] + k.shape[i])) for i in range(3))
        src = tuple(slice(dst[i].start - o[i], dst[i].stop - o[i]) for i in range(3))
        V[dst] += k[src]
    return V


def level_and_background(V_ideal, V_recon, cfg, dx):
    """@brief Scale that puts the ideal volume at the reconstruction's phase level, and the
    reconstruction's background sigma (robust, per-layer high-passed)."""
    from scipy.ndimage import gaussian_filter
    scale = (np.percentile(AFF.preprocess(V_recon, cfg, dx), 99.9)
             / np.percentile(AFF.preprocess(V_ideal, cfg, dx), 99.9))
    hp = np.stack([L - gaussian_filter(L, cfg.bg_smooth_A / dx) for L in V_recon])
    return float(scale), 1.4826 * float(np.median(np.abs(hp - np.median(hp))))


def ideal_volume(ctx, noise, seed=0):
    """@brief The ideal control at `noise` x the reconstruction background (one draw per level)."""
    V0 = render_ideal(ctx["V"].shape, ctx["frame"], ctx["pos_t"], ctx["sym_t"], ctx["kernels"],
                      ctx["dz"])
    scale, sigma = level_and_background(V0, ctx["V"], ctx["cfg"], ctx["dx"])
    return V0 * scale + np.random.default_rng(seed).normal(0.0, noise * sigma, V0.shape), scale, sigma


def find_or_load(ctx, V, cache_path=None, reuse=False):
    """@brief atomfind v3 on V, cached to `cache_path` (.npy) when given."""
    if reuse and cache_path and os.path.exists(cache_path):
        return np.load(cache_path)
    found, _ = AFF.find_atoms_v3(V, ctx["cfg"], ctx["dx"], ctx["kernels"])
    if cache_path:
        np.save(cache_path, found)
    return found


def located(found, frame, dz):
    """@brief Finder records -> sample-frame xyz (A) and species."""
    x0, y0, ox, oy, dx = frame
    xyz = np.c_[found["col"] * dx + x0 - ox, found["row"] * dx + y0 - oy,
                (found["layer"] + 0.5) * dz]
    return xyz, np.asarray(found["species"])


def in_field(xyz, frame, shape):
    x0, y0, ox, oy, dx = frame
    X0, Y0 = x0 - ox + EDGE_A, y0 - oy + EDGE_A
    X1, Y1 = x0 - ox + shape[2] * dx - EDGE_A, y0 - oy + shape[1] * dx - EDGE_A
    return (xyz[:, 0] > X0) & (xyz[:, 0] < X1) & (xyz[:, 1] > Y0) & (xyz[:, 1] < Y1)


# ---------------------------------------------------------------- matching + atomfind error bars
def match_to_truth(xyz, pos_t, tol_xy=0.6, tol_z=2.0):
    """@brief One-to-one match of located atoms to true atoms (any species), nearest first.

    Returns the dict `atomfind.uncertainty.calibrate` expects: match_gi (-1 = unmatched) and the
    signed x/y/z errors (located - true). atomfind's own tolerances (cfg.match_tol_*)."""
    from scipy.spatial import cKDTree
    n = len(xyz)
    out = dict(match_gi=np.full(n, -1), match_dx=np.full(n, np.nan),
               match_dy=np.full(n, np.nan), match_dz=np.full(n, np.nan))
    if not n:
        return out
    tree = cKDTree(pos_t[:, :2])
    pairs = []
    for i, p in enumerate(xyz):
        for g in tree.query_ball_point(p[:2], tol_xy):
            dzv = abs(p[2] - pos_t[g, 2])
            if dzv <= tol_z:
                pairs.append((np.hypot(p[0] - pos_t[g, 0], p[1] - pos_t[g, 1]) + dzv, i, g))
    used = set()
    for _, i, g in sorted(pairs):
        if out["match_gi"][i] >= 0 or g in used:
            continue
        used.add(g)
        out["match_gi"][i] = g
        out["match_dx"][i], out["match_dy"][i], out["match_dz"][i] = xyz[i] - pos_t[g]
    return out


def atomfind_intervals(found, xyz, ctx, alpha=0.05):
    """@brief atomfind's calibrated half-widths (x, y, z) at level 1-alpha, and the coverage check."""
    m = match_to_truth(xyz, ctx["pos_t"])
    qtab = AFU.calibrate(found, m, ctx["cfg"], alphas=(alpha,), min_n=ctx["cfg"].uq_min_stratum)
    hw = AFU.apply(found, qtab, ctx["cfg"], alpha)
    _, cov = AFU.coverage_table(found, m, ctx["cfg"], qtab, alpha)
    return hw, m, cov


# ---------------------------------------------------------------- Ti against its oxygens
def ti_offsets(ti, ox):
    """@brief Per-Ti z offset against its oxygens, by each estimator (NaN where undefined)."""
    from scipy.spatial import cKDTree
    out = {k: np.full(len(ti), np.nan) for k in ESTIMATORS}
    if not len(ti) or not len(ox):
        return out
    if len(ox) >= 6:
        d6, ok, _ = AFPOL.offcentring(ti, ox)
        out["O6 cage"][ok] = d6[ok, 2]
    tree = cKDTree(ox[:, :2])
    for i, p in enumerate(ti):
        near = tree.query_ball_point(p[:2], 2.4)
        if not near:
            continue
        q = ox[near]
        rxy = np.hypot(q[:, 0] - p[0], q[:, 1] - p[1])
        dzq = q[:, 2] - p[2]
        col = rxy < 0.6
        up, dn = col & (dzq > 0.8) & (dzq < 2.8), col & (dzq < -0.8) & (dzq > -2.8)
        if up.any() and dn.any():
            out["apical pair"][i] = p[2] - 0.5 * (q[up, 2][np.argmin(dzq[up])]
                                                  + q[dn, 2][np.argmax(dzq[dn])])
        eq = (rxy > 1.5) & (rxy < 2.4) & (np.abs(dzq) < 1.0)
        if eq.sum() >= 2:
            out["equatorial O"][i] = p[2] - q[eq, 2].mean()
    return out


def equatorial_offsets(found, ctx, alpha=0.05):
    """@brief Per located Ti: z_Ti - mean z(equatorial O) with atomfind's propagated interval.

    @return dict of per-Ti arrays: xyz, dz, hw (half-width at 1-alpha), n_o, truth_dz, domain;
            plus 'coverage' (atomfind's per-axis coverage on this volume) and 'z_rms'.
    """
    from scipy.spatial import cKDTree
    truth, frame, V = ctx["truth"], ctx["frame"], ctx["V"]
    a, n = float(truth["a"]), int(truth["n_lat"])
    xyz, spec = located(found, frame, ctx["dz"])
    hw, m, cov = atomfind_intervals(found, xyz, ctx, alpha)
    sig = hw["z"] / 1.959964                    # the half-width back to a sigma, as polarisation.py
    io = np.where(spec == 8)[0]
    tree = cKDTree(xyz[io, :2]) if len(io) else None

    pos_t, sym_t = ctx["pos_t"], ctx["sym_t"]
    tti = np.where(sym_t == "Ti")[0]
    true_eq = ti_offsets(pos_t[tti], pos_t[sym_t == "O"])["equatorial O"]
    ttree = cKDTree(pos_t[tti, :2])

    rows = []
    for i in np.where((spec == 22) & in_field(xyz, frame, V.shape))[0]:
        p = xyz[i]
        near = io[tree.query_ball_point(p[:2], 2.4)] if tree is not None else []
        if not len(near):
            continue
        rxy = np.hypot(xyz[near, 0] - p[0], xyz[near, 1] - p[1])
        eq = near[(rxy > 1.5) & (rxy < 2.4) & (np.abs(xyz[near, 2] - p[2]) < 1.0)]
        if len(eq) < 2:
            continue
        dz = p[2] - xyz[eq, 2].mean()
        s = np.sqrt(sig[i] ** 2 + np.sum(sig[eq] ** 2) / len(eq) ** 2)
        tdz = np.nan
        cand = ttree.query_ball_point(p[:2], 0.6)
        if cand:
            g = cand[int(np.argmin(np.abs(pos_t[tti[cand], 2] - p[2])))]
            if abs(pos_t[tti[g], 2] - p[2]) < 2.0:
                tdz = true_eq[g]
        rows.append((p[0], p[1], p[2], dz, 1.959964 * s, len(eq), tdz,
                     str(truth["domain_grid"][int(np.floor(p[0] / a)) % n,
                                              int(np.floor(p[1] / a)) % n])))
    R = {k: np.array([r[j] for r in rows]) for j, k in enumerate(
        ("x", "y", "z", "dz", "hw", "n_o", "truth_dz", "domain"))}
    R["xyz_all"], R["spec_all"], R["hw_all"] = xyz, spec, hw
    R["coverage"] = cov
    ok = m["match_gi"] >= 0
    R["z_rms"] = {s: float(np.sqrt(np.nanmean(m["match_dz"][ok & (spec == Z_OF[s])] ** 2)))
                  for s in Z_OF if (ok & (spec == Z_OF[s])).any()}
    return R


def domain_summary(R, truth):
    """@brief Per domain: n, inverse-variance mean with its interval, and how the Ti bars fall."""
    names = [str(x) for x in truth["names"]]
    out = {}
    for k, nm in enumerate(names):
        mm = R["domain"] == nm if len(R["dz"]) else np.array([], bool)
        if not mm.any():
            continue
        v, h = R["dz"][mm], R["hw"][mm]
        w = 1.0 / h ** 2
        s = float(np.sign(truth["deltas"][k][2]))
        out[nm] = dict(n=int(mm.sum()), truth=s, mean=float(np.sum(w * v) / np.sum(w)),
                       hw_mean=float(1.0 / np.sqrt(np.sum(w))),
                       right_side=float(np.mean(np.sign(v) == s)),
                       excl_zero_right=int(np.sum((np.abs(v) > h) & (np.sign(v) == s))),
                       excl_zero_wrong=int(np.sum((np.abs(v) > h) & (np.sign(v) != s))),
                       median_hw=float(np.median(h)))
    return out


# ---------------------------------------------------------------- the three-estimator score
def score(xyz, spec, truth, pos_t, sym_t, frame, shape):
    """@brief Located atoms -> counts, z accuracy, and the sign by each estimator, per domain."""
    from scipy.spatial import cKDTree
    a, n = float(truth["a"]), int(truth["n_lat"])
    grid, names = truth["domain_grid"], [str(x) for x in truth["names"]]
    sign_of = {nm: float(np.sign(truth["deltas"][k][2])) for k, nm in enumerate(names)}

    ff, ft = in_field(xyz, frame, shape), in_field(pos_t, frame, shape)
    m = match_to_truth(xyz[ff], pos_t[ft])
    ok = m["match_gi"] >= 0
    sp_f = spec[ff]
    R = dict(counts={s: [int(((spec == Z_OF[s]) & ff).sum()), int(((sym_t == s) & ft).sum())]
                     for s in Z_OF},
             z_rms_A={s: (float(np.sqrt(np.mean(m["match_dz"][ok & (sp_f == Z_OF[s])] ** 2)))
                          if (ok & (sp_f == Z_OF[s])).any() else None) for s in Z_OF},
             estimators={})

    ti = xyz[(spec == 22) & ff]
    found_off = ti_offsets(ti, xyz[spec == 8])
    tti = pos_t[sym_t == "Ti"]
    true_off = ti_offsets(tti, pos_t[sym_t == "O"])
    dist, j = cKDTree(tti).query(ti, distance_upper_bound=1.0) if len(ti) else (np.array([]),) * 2
    gi, gj = np.floor(ti[:, 0] / a).astype(int), np.floor(ti[:, 1] / a).astype(int)
    dom = grid[gi % n, gj % n] if len(ti) else np.array([])

    for est in ESTIMATORS:
        v = found_off[est]
        okv = np.isfinite(v)
        truth_v = np.full(len(ti), np.nan)
        mm = np.isfinite(dist)
        truth_v[mm] = true_off[est][j[mm]]
        scored = okv & np.isfinite(truth_v)
        E = dict(n_ti=int(okv.sum()), of_ti=int(len(ti)),
                 true_abs_median_A=float(np.nanmedian(np.abs(true_off[est]))),
                 ti_correct=(float(np.mean(np.sign(v[scored]) == np.sign(truth_v[scored])))
                             if scored.any() else None),
                 error_rms_A=(float(np.sqrt(np.mean((v[scored] - truth_v[scored]) ** 2)))
                              if scored.any() else None),
                 domains={})
        for nm in names:
            md = okv & (dom == nm)
            if not md.any():
                continue
            w = v[md]
            se = float(w.std(ddof=1) / np.sqrt(len(w))) if len(w) > 1 else float("nan")
            E["domains"][nm] = dict(truth=int(sign_of[nm]), n=int(len(w)), mean_A=float(w.mean()),
                                    se_A=se, decided=int(np.sign(w.mean())),
                                    correct=bool(np.sign(w.mean()) == sign_of[nm]))
        R["estimators"][est] = E
    return R


def report(label, R):
    c = R["counts"]
    zr = {k: (v if v is not None else float("nan")) for k, v in R["z_rms_A"].items()}
    L = [f"[{label}] located/true in field: Pb {c['Pb'][0]}/{c['Pb'][1]}  Ti {c['Ti'][0]}/{c['Ti'][1]}  "
         f"O {c['O'][0]}/{c['O'][1]};  z error RMS  Pb {zr['Pb']:.2f}  Ti {zr['Ti']:.2f}  O {zr['O']:.2f} A"]
    for est, E in R["estimators"].items():
        tc = "  n/a" if E["ti_correct"] is None else f"{100 * E['ti_correct']:4.0f}%"
        er = "  n/a" if E["error_rms_A"] is None else f"{E['error_rms_A']:.2f}"
        cells = [f"{nm} {d['mean_A']:+.2f}±{d['se_A']:.2f} (n={d['n']}) "
                 f"{'OK' if d['correct'] else 'WRONG'}" for nm, d in E["domains"].items()]
        L.append(f"  {est:13s} Ti used {E['n_ti']:3d}/{E['of_ti']:<3d} | per-Ti sign right {tc} | "
                 f"offset error RMS {er} A (true |dz| {E['true_abs_median_A']:.2f}) | "
                 + ("  ".join(cells) if cells else "no Ti usable"))
    return "\n".join(L)


def run(recon, budget_path, truth_path, vasp, noises=(1.0,), only="both", cache=None,
        reuse=False, out=None, seed=0):
    ctx = context(recon, budget_path, truth_path, vasp)
    volumes = []
    if only in ("both", "ideal"):
        for fr in noises:
            V, scale, sigma = ideal_volume(ctx, fr, seed)
            volumes.append((f"ideal, noise {fr:g}x background", f"ideal_noise{fr:g}", V))
        print(f"[atomfind-sign] ideal volume x{scale:.3g}; reconstruction background sigma "
              f"{sigma:.3g} rad")
    if only in ("both", "recon"):
        volumes.append(("blind reconstruction", "recon", ctx["V"]))

    results = dict(recon=recon, dx_A=ctx["dx"], dz_A=ctx["dz"])
    for label, tag, V in volumes:
        found = find_or_load(ctx, V, f"{cache}_{tag}_found.npy" if cache else None, reuse)
        xyz, spec = located(found, ctx["frame"], ctx["dz"])
        results[tag] = score(xyz, spec, ctx["truth"], ctx["pos_t"], ctx["sym_t"], ctx["frame"],
                             V.shape)
        E = equatorial_offsets(found, ctx)
        results[tag]["equatorial_with_atomfind_95"] = dict(
            coverage95=E["coverage"], domains=domain_summary(E, ctx["truth"]))
        print(report(label, results[tag]))
        for nm, d in results[tag]["equatorial_with_atomfind_95"]["domains"].items():
            print(f"    {nm}: weighted mean {d['mean']:+.3f} ± {d['hw_mean']:.3f} A (95 %), "
                  f"n={d['n']}, median Ti bar ±{d['median_hw']:.2f} A, bars clear of zero: "
                  f"{d['excl_zero_right']} right / {d['excl_zero_wrong']} wrong")
        print(f"    atomfind 95 % coverage on this volume: " +
              "  ".join(f"{k} {100 * v:.0f}%" for k, v in E["coverage"].items()), flush=True)
    if out:
        json.dump(results, open(out, "w"), indent=1, default=str)
        print(f"[atomfind-sign] wrote {out}")
    return results


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--recon", required=True)
    ap.add_argument("--budget", required=True)
    ap.add_argument("--truth", default=os.path.join(HERE, "sample", "toy_truth.npz"))
    ap.add_argument("--vasp", default=os.path.join(HERE, "sample", "toy_membrane.vasp"))
    ap.add_argument("--ideal-noise", type=float, nargs="+", default=[1.0],
                    help="noise on the ideal volume, as multiples of the reconstruction background")
    ap.add_argument("--only", choices=("both", "ideal", "recon"), default="both")
    ap.add_argument("--cache", default=None, help="prefix for cached finder output")
    ap.add_argument("--reuse", action="store_true", help="reuse cached finder output")
    ap.add_argument("--out", default=None)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)
    run(args.recon, args.budget, args.truth, args.vasp, tuple(args.ideal_noise), args.only,
        args.cache, args.reuse, args.out, args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
