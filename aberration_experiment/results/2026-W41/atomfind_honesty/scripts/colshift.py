"""Blind per-column depth evidence: does the image prefer the found depths over the same atoms moved half a spacing?
Per column (found col_id): tube as the finder cut it (preprocessed, clipped, interior depth); model = sum of the Pb kernel at
the found positions with NNLS amplitudes; the same with every atom moved by +-s (s = half the column's own median atom
spacing). evidence e = (min chi2(+-s) - chi2(found)) / ||tube - depth-mean(tube)||^2: the share of the column's depth
structure that the found depths explain better than the shifted ones. ~0 on a featureless rod."""
import os, sys, json, copy, numpy as np
sys.argv_ = sys.argv
sys.argv = ['x']; exec(open(os.path.join(os.environ['S'], 'refind.py')).read().split('if __name__')[0])
from scipy.optimize import nnls
from atomfind.fit import _render
from atomfind.deconv import crop_kernel_inplane

def evid(tube, at, s, dz, K, hz, hxy, S2=None):
    """(e, chi2 at the found depths, design, amplitudes) for atoms at tube-frame (l, r, c)."""
    def fit(shift_A):
        M = np.column_stack([_render(tube.shape, (l + shift_A / dz, r, cc), K, hz, hxy) for (l, r, cc) in at])
        b, rn = nnls(M, tube.ravel())
        return rn ** 2, M, b
    c0, M0, b0 = fit(0.0)
    cs = min(fit(s)[0], fit(-s)[0])
    if S2 is None:
        S2 = float(((tube - tube.mean(0, keepdims=True)) ** 2).sum())
    return ((cs - c0) / S2 if S2 > 0 else np.nan), c0, M0, b0

def column_evidence(V, c, dx, found, K):
    Vp = find.preprocess(V, c, dx) if c.preprocess_bg else V
    l0 = int(round(c.trim_z_A[0] / c.dz)); l1 = min(int(round(c.trim_z_A[1] / c.dz)), V.shape[0])
    HW = c.tube_halfwidth_px
    hz = (K.shape[0] - 1) // 2; hxy = (K.shape[1] - 1) // 2
    out = {}
    for cid in np.unique(found["col_id"]):
        A = found[found["col_id"] == cid]
        ri, ci = int(round(np.median(A["row"]))), int(round(np.median(A["col"])))
        if not (HW <= ri < V.shape[1] - HW and HW <= ci < V.shape[2] - HW):
            continue
        tube = np.clip(Vp[l0:l1, ri - HW:ri + HW + 1, ci - HW:ci + HW + 1], 0, None)
        zl = np.sort(A["layer"]) * c.dz
        s = 0.5 * float(np.median(np.diff(zl))) if len(zl) >= 3 else 0.5 * c.comb_period_A
        s = float(np.clip(s, 0.5, 2.5))
        at = [(a["layer"] - l0, a["row"] - (ri - HW), a["col"] - (ci - HW)) for a in A]
        e, c0, M0, b0 = evid(tube, at, s, c.dz, K, hz, hxy)
        # NO-DEPTH NULL: the tube with the depth structure the found atoms explain removed (depth mean + the fit residual),
        # put through the finder's own tube core (CLEAN + refine + quality cut); e of what it finds there
        model = (M0 @ b0).reshape(tube.shape)
        S2r = float(((tube - tube.mean(0, keepdims=True)) ** 2).sum())
        null = tube - model + model.mean(0, keepdims=True)      # the found atoms' depth modulation replaced by its depth mean
        null = np.clip(null, 0, None)
        na = find.clean_tube(null, K, c, floor=find.clean_floor_for(null, K, c, dx))
        en = np.nan
        if na:
            rf = [d for d in find.refine_tube(null, na, K, c) if d["quality"] >= c.quality_min_corr and d["amp"] > 0]
            if rf:
                zn = np.sort([d["l"] for d in rf]) * c.dz
                sn = float(np.clip(0.5 * np.median(np.diff(zn)), 0.5, 2.5)) if len(zn) >= 3 else 0.5 * c.comb_period_A
                en = evid(null, [(d["l"], d["r"], d["c"]) for d in rf], sn, c.dz, K, hz, hxy, S2=S2r)[0]
        out[int(cid)] = dict(e=e, e_null=en, n=len(A), s=s, r2=1 - c0 / float((tube ** 2).sum()))
    return out

if __name__ == "__main__":
    for key in sys.argv_[1:]:
        info = R.find(key)
        V, cfg, pos, Z, found, al = dh.frame(info)
        pb, ti = kernels_for(info)
        c = copy.copy(cfg); c.single_atom_vol, c.ti_kernel_vol = pb, ti
        kern = psfmod.species_kernels(c, al.dx)
        K = find._unit_norm(crop_kernel_inplane(kern[82]))
        ev = column_evidence(V.astype(np.float64), c, al.dx, found, K)
        e = np.array([v["e"] for v in ev.values()]); n = np.array([v["n"] for v in ev.values()])
        r2 = np.array([v["r2"] for v in ev.values()])
        en = np.array([v["e_null"] for v in ev.values()]); thr = np.nanpercentile(en, 95); atoms_w = np.repeat(e, n)
        print(f"{'':11s} null e: median {np.nanmedian(en):+.3f} p95 {thr:+.3f} (finder found atoms in {np.isfinite(en).sum()}/{len(en)} null columns) | "
              f"atoms in columns above the null p95: {100*(atoms_w > thr).mean():3.0f}%")
        atoms_w = np.repeat(e, n)
        print(f"{key:11s} cols {len(e):3d} | evidence e: median {np.median(e):+.3f}  atom-weighted median {np.median(atoms_w):+.3f}  "
              f"share of atoms in columns with e>0.05: {100*(atoms_w>0.05).mean():3.0f}%  e<0.01: {100*(atoms_w<0.01).mean():3.0f}% | fit r2 median {np.median(r2):.2f}", flush=True)
        json.dump({str(k): v for k, v in ev.items()}, open(os.path.join(os.environ['S'], f"colshift_{key}.json"), "w"))
