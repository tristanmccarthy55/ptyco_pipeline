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

Two kinds of volume go through the identical finder, kernels and analysis:

  ideal  the TRUE atoms rendered with the measured single-atom kernels (the forward model atomfind
         inverts) on the reconstruction's own grid and lattice registration, scaled to its phase
         level, plus white noise at a chosen fraction of its background level: a reconstruction
         with every atom at its true depth. The positive control.
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

KERNEL_PB = os.path.expanduser("~/Desktop/psf_Pb_NL70_vol.npy")
KERNEL_TI = os.path.expanduser("~/Desktop/psf_Ti_NL70_vol.npy")
Z_OF = {"Pb": 82, "Ti": 22, "O": 8}
EDGE_A = 1.2                    # stay off the reconstruction rim when counting and scoring
ESTIMATORS = ("O6 cage", "apical pair", "equatorial O")


def finder_config(dx, dz, thickness_A, c_A):
    """@brief atomfind settings for this slab: whole depth analysed, no exit-band inflation."""
    return AFC.Config(name="fusion", recon_vol="", dx=dx, dz=dz, X0=None, Y0=None,
                      trim_z_A=(0.0, thickness_A), zmax_show_A=thickness_A,
                      exit_band_z_A=thickness_A, clean_max_atoms=16, comb_period_A=c_A,
                      convergence_mrad=100.0, single_atom_vol=KERNEL_PB, single_atom_species=82,
                      ti_kernel_vol=KERNEL_TI)


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


def located(found, frame, dz):
    """@brief Finder records -> sample-frame xyz (A) and species."""
    x0, y0, ox, oy, dx = frame
    xyz = np.c_[found["col"] * dx + x0 - ox, found["row"] * dx + y0 - oy,
                (found["layer"] + 0.5) * dz]
    return xyz, np.asarray(found["species"])


def _in_field(xyz, frame, shape):
    x0, y0, ox, oy, dx = frame
    X0, Y0 = x0 - ox + EDGE_A, y0 - oy + EDGE_A
    X1, Y1 = x0 - ox + shape[2] * dx - EDGE_A, y0 - oy + shape[1] * dx - EDGE_A
    return (xyz[:, 0] > X0) & (xyz[:, 0] < X1) & (xyz[:, 1] > Y0) & (xyz[:, 1] < Y1)


def _z_error(xyz_f, xyz_t, tol_xy=0.6):
    """@brief Signed z error of each located atom against the nearest true atom in its column."""
    from scipy.spatial import cKDTree
    if not len(xyz_f) or not len(xyz_t):
        return np.array([])
    tree = cKDTree(xyz_t[:, :2])
    err = []
    for p in xyz_f:
        idx = tree.query_ball_point(p[:2], tol_xy)
        if idx:
            dzs = p[2] - xyz_t[idx, 2]
            err.append(float(dzs[np.argmin(np.abs(dzs))]))
    return np.asarray(err)


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


def score(xyz, spec, truth, pos_t, sym_t, frame, shape):
    """@brief Located atoms -> counts, z accuracy, and the sign by each estimator, per domain."""
    from scipy.spatial import cKDTree
    a, n = float(truth["a"]), int(truth["n_lat"])
    grid, names = truth["domain_grid"], [str(x) for x in truth["names"]]
    sign_of = {nm: float(np.sign(truth["deltas"][k][2])) for k, nm in enumerate(names)}

    ff, ft = _in_field(xyz, frame, shape), _in_field(pos_t, frame, shape)
    R = dict(counts={s: [int(((spec == Z_OF[s]) & ff).sum()), int(((sym_t == s) & ft).sum())]
                     for s in Z_OF},
             z_rms_A={s: (lambda e: float(np.sqrt(np.mean(e ** 2))) if len(e) else None)(
                 _z_error(xyz[(spec == Z_OF[s]) & ff], pos_t[(sym_t == s) & ft]))
                 for s in Z_OF}, estimators={})

    ti = xyz[(spec == 22) & ff]
    found_off = ti_offsets(ti, xyz[spec == 8])
    tti = pos_t[sym_t == "Ti"]
    true_off = ti_offsets(tti, pos_t[sym_t == "O"])
    dist, j = cKDTree(tti).query(ti, distance_upper_bound=1.0) if len(ti) else (np.array([]),) * 2
    gi, gj = np.floor(ti[:, 0] / a).astype(int), np.floor(ti[:, 1] / a).astype(int)
    dom = grid[gi % n, gj % n] if len(ti) else np.array([])

    for est in ESTIMATORS:
        v = found_off[est]
        ok = np.isfinite(v)
        truth_v = np.full(len(ti), np.nan)
        m = np.isfinite(dist)
        truth_v[m] = true_off[est][j[m]]
        scored = ok & np.isfinite(truth_v)
        E = dict(n_ti=int(ok.sum()), of_ti=int(len(ti)),
                 true_abs_median_A=float(np.nanmedian(np.abs(true_off[est]))),
                 ti_correct=(float(np.mean(np.sign(v[scored]) == np.sign(truth_v[scored])))
                             if scored.any() else None),
                 error_rms_A=(float(np.sqrt(np.mean((v[scored] - truth_v[scored]) ** 2)))
                              if scored.any() else None),
                 domains={})
        for nm in names:
            mm = ok & (dom == nm)
            if not mm.any():
                continue
            w = v[mm]
            se = float(w.std(ddof=1) / np.sqrt(len(w))) if len(w) > 1 else float("nan")
            E["domains"][nm] = dict(truth=int(sign_of[nm]), n=int(len(w)), mean_A=float(w.mean()),
                                    se_A=se, decided=int(np.sign(w.mean())),
                                    correct=bool(np.sign(w.mean()) == sign_of[nm]),
                                    significance=float(abs(w.mean()) / se) if se > 0 else None)
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
        cells = []
        for nm, d in E["domains"].items():
            cells.append(f"{nm} {d['mean_A']:+.2f}±{d['se_A']:.2f} (n={d['n']}) "
                         f"{'OK' if d['correct'] else 'WRONG'}")
        L.append(f"  {est:13s} Ti used {E['n_ti']:3d}/{E['of_ti']:<3d} | per-Ti sign right {tc} | "
                 f"offset error RMS {er} A (true |dz| {E['true_abs_median_A']:.2f}) | "
                 + ("  ".join(cells) if cells else "no Ti usable"))
    return "\n".join(L)


def run(recon, budget_path, truth_path, vasp, noises=(1.0,), only="both", cache=None,
        reuse=False, out=None, seed=0):
    import ase.io
    import h5py
    truth = np.load(truth_path, allow_pickle=True)
    budget = json.load(open(budget_path))
    with h5py.File(recon, "r") as f:
        budget["dx_object_A"] = float(np.asarray(f["outputs"]["pixel_size"]).ravel()[0]) * 1e10
    V_recon = AF.load_volume(recon)
    T = float(budget["beam_thickness_A"])
    dz = T / V_recon.shape[0]
    frame = AF.lattice_frame(V_recon, budget, truth)
    dx = frame[4]
    cfg = finder_config(dx, dz, T, float(truth["c"]))
    kernels = AFP.species_kernels(cfg, dx)
    atoms = ase.io.read(vasp)
    pos_t, sym_t = atoms.get_positions(), np.array(atoms.get_chemical_symbols())

    volumes = []
    if only in ("both", "ideal"):
        V0 = render_ideal(V_recon.shape, frame, pos_t, sym_t, kernels, dz)
        scale, sigma = level_and_background(V0, V_recon, cfg, dx)
        print(f"[atomfind-sign] ideal volume x{scale:.3g}; reconstruction background sigma "
              f"{sigma:.3g} rad")
        rng = np.random.default_rng(seed)
        for fr in noises:
            volumes.append((f"ideal, noise {fr:g}x background",
                            f"ideal_noise{fr:g}", V0 * scale + rng.normal(0.0, fr * sigma, V0.shape)))
    if only in ("both", "recon"):
        volumes.append(("blind reconstruction", "recon", V_recon))

    results = dict(recon=recon, dx_A=dx, dz_A=dz)
    for label, tag, V in volumes:
        cpath = f"{cache}_{tag}_found.npy" if cache else None
        if reuse and cpath and os.path.exists(cpath):
            found = np.load(cpath)
        else:
            found, _ = AFF.find_atoms_v3(V, cfg, dx, kernels)
            if cpath:
                np.save(cpath, found)
        xyz, spec = located(found, frame, dz)
        results[tag] = score(xyz, spec, truth, pos_t, sym_t, frame, V.shape)
        print(report(label, results[tag]), flush=True)
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
