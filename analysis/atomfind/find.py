#!/usr/bin/env python
"""Blind atom finder -- no ground truth anywhere in this module.

The physics: in-plane is easy (~0.1 A, 2 px -> columns crisp, sub-pixel xy trivial);
DEPTH is the crutch (dz ~ 1 A, atoms ~4 layers apart overlap through the measured kernel).

v2 (default, `find_atoms_v2`) -- 3-D TUBE CLEAN. The v1 1-D column-profile reduction was
measured to fail on dim B-O columns: lean-tracking wanders 0.5-0.65 A off-column in the
weak-signal stretches, distorting the profile (Ti recall capped at 65%). v2 removes the
reduction entirely: per detected column, take a fixed tube (half-width ~0.5 A, contains
the lean) and run MATCHING PURSUIT with the measured 3-D system kernel -- correlate, take
the max, subtract amp*K, repeat (CLEAN, as in radio astronomy). Every atom gets its OWN
(x, y, z). Then a joint Gauss-Newton refinement per tube gives sub-voxel positions and
1-sigma ERROR BARS from the nonlinear-LS covariance, and each atom is species-classified
(amplitude bands + Pb-vs-Ti kernel shape, both measured from the single-atom sims) with
junk rejected by FIT QUALITY, not amplitude (an amplitude floor alone was measured to
trade Ti recall against precision).

v1 (`find_atoms`, kept as the comparison baseline): per-column z-profile + 1-D NNLS spike
deconvolution; `raw_peak_z` is the naive baseline below that.

Outputs are in the RECON frame -- the validator maps GT into this frame to score, so no GT
leaks here. `export_atoms` converts to the prepared-cell physical frame (calibration
constants, not GT positions) and writes CSV + an ASE extxyz object.
"""
from __future__ import annotations
import numpy as np
from scipy.ndimage import maximum_filter
from scipy.signal import fftconvolve

from . import psf as _psf
from .fit import _nnls_gram, _render
from .deconv import crop_kernel_inplane


# ---------------------------------------------------------------- in-plane columns
def detect_columns(V, cfg, dx):
    """Local maxima of the interior depth-mean -> sub-pixel (row, col) column seeds."""
    l0, l1 = int(round(cfg.trim_z_A[0]/cfg.dz)), int(round(cfg.trim_z_A[1]/cfg.dz))
    dm = np.clip(V[max(l0, 0):l1], 0, None).mean(0)
    sep = max(3, int(round(cfg.find_min_sep_A / dx)))
    thr = np.percentile(dm, cfg.find_col_pct)
    peaks = (dm == maximum_filter(dm, size=sep)) & (dm > thr)
    ys, xs = np.where(peaks)
    seeds = []
    for r, c in zip(ys, xs):
        if not (1 <= r < dm.shape[0]-1 and 1 <= c < dm.shape[1]-1):
            continue
        # parabolic sub-pixel refine on the 3x3 patch
        rr = r + _parab(dm[r-1, c], dm[r, c], dm[r+1, c])
        cc = c + _parab(dm[r, c-1], dm[r, c], dm[r, c+1])
        seeds.append((rr, cc, float(dm[r, c])))
    return seeds, dm


def _parab(ym, y0, yp):
    d = ym - 2*y0 + yp
    return 0.0 if abs(d) < 1e-12 else np.clip(0.5*(ym - yp)/d, -0.5, 0.5)


def track_column(V, r0, c0, cfg):
    """Follow a column down z by per-layer centroid in a small window -> path r(l), c(l)."""
    nL = V.shape[0]
    w = cfg.find_track_win_px
    navg = cfg.find_profile_navg
    r, c = r0, c0
    path_r, path_c, prof = np.zeros(nL), np.zeros(nL), np.zeros(nL)
    for l in range(nL):
        ri, ci = int(round(r)), int(round(c))
        if not (w <= ri < V.shape[1]-w and w <= ci < V.shape[2]-w):
            path_r[l], path_c[l] = r, c
            continue
        win = np.clip(V[l, ri-w:ri+w+1, ci-w:ci+w+1], 0, None)
        if win.sum() > 1e-6:                                    # update centroid if there is signal
            gy, gx = np.indices(win.shape)
            r = ri - w + (gy*win).sum()/win.sum()
            c = ci - w + (gx*win).sum()/win.sum()
        path_r[l], path_c[l] = r, c
        rp, cp = int(round(r)), int(round(c))
        prof[l] = np.clip(V[l, rp-navg:rp+navg+1, cp-navg:cp+navg+1], 0, None).mean()
    return path_r, path_c, prof


# ---------------------------------------------------------------- axial spike deconvolution
def _spike_design(nL, k1d, zoff, cand, dz):
    """Design matrix A[l, j] = axial kernel placed at fractional-layer candidate cand[j]."""
    A = np.zeros((nL, len(cand)))
    L = np.arange(nL)
    for j, p in enumerate(cand):
        A[:, j] = np.interp(L - p, zoff, k1d, left=0.0, right=0.0)
    return A


def spike_deconv(prof, k1d, zoff, cfg):
    """NNLS spike deconvolution of a 1-D column profile. Returns atoms [(z_layer, amp)]."""
    nL = len(prof)
    l0 = int(round(cfg.trim_z_A[0]/cfg.dz)); l1 = int(round(cfg.trim_z_A[1]/cfg.dz))
    step = cfg.spike_grid_A / cfg.dz
    cand = np.arange(l0, l1 + 1e-9, step)
    A = _spike_design(nL, k1d, zoff, cand, cfg.dz)
    Aa = np.column_stack([A, np.ones(nL)])                     # + constant background
    beta = _nnls_gram(Aa, np.clip(prof, 0, None), ridge=cfg.fit_ridge)
    amp = beta[:-1]
    if amp.max() <= 0:
        return []
    # merge adjacent nonzero candidates within spike_merge_A -> one atom (amp-weighted z)
    keep = amp > cfg.spike_min_frac * amp.max()
    atoms, i = [], 0
    idx = np.where(keep)[0]
    merge_lay = cfg.spike_merge_A / cfg.dz
    groups = []
    for k in idx:
        if groups and (cand[k] - cand[groups[-1][-1]]) <= merge_lay:
            groups[-1].append(k)
        else:
            groups.append([k])
    for g in groups:
        w = amp[g]
        zc = np.sum(cand[g] * w) / w.sum()
        atoms.append((float(zc), float(w.sum())))
    return atoms


def raw_peak_z(prof, cfg, k1d=None):
    """Baseline: local maxima of the raw (optionally lightly matched-filtered) profile."""
    nL = len(prof)
    p = np.clip(prof, 0, None)
    l0 = int(round(cfg.trim_z_A[0]/cfg.dz)); l1 = int(round(cfg.trim_z_A[1]/cfg.dz))
    sep = max(1, int(round(3.5 / cfg.dz)))                     # ~one plane spacing
    mx = maximum_filter(p, size=2*sep+1)
    atoms = []
    for l in range(max(l0, 1), min(l1, nL-1)):
        if p[l] == mx[l] and p[l] > 0.15 * p.max():
            zc = l + _parab(p[l-1], p[l], p[l+1])              # sub-layer refine
            atoms.append((float(zc), float(p[l])))
    return atoms


# ================================================================= v2: 3-D tube CLEAN
def _unit_norm(K):
    K = np.clip(K, 0, None)
    return K / np.sqrt((K ** 2).sum())


def clean_tube(tube, K, cfg):
    """Matching pursuit (CLEAN) with the unit-L2-norm kernel K in one column tube.

    Correlate residual with K, take the max (= that atom's amplitude, since ||K||=1),
    subtract amp*K, non-max-suppress around accepted atoms, stop at cfg.clean_floor.
    The floor is deliberately LOW (completeness); junk is culled later by fit quality."""
    hz = (K.shape[0]-1)//2
    hxy = (K.shape[1]-1)//2
    Km = K[::-1, ::-1, ::-1]
    R = tube.copy()
    atoms = []
    nz = cfg.clean_nms_z_layers
    nxy = cfg.clean_nms_xy_px
    for _ in range(cfg.clean_max_atoms):
        C = fftconvolve(R, Km, mode="same")
        for (l, r, c, _a) in atoms:                      # suppress re-detections
            z0, z1 = max(int(l-nz), 0), min(int(l+nz)+1, C.shape[0])
            y0, y1 = max(int(r)-nxy, 0), min(int(r)+nxy+1, C.shape[1])
            x0, x1 = max(int(c)-nxy, 0), min(int(c)+nxy+1, C.shape[2])
            C[z0:z1, y0:y1, x0:x1] = -1
        i = np.unravel_index(np.argmax(C), C.shape)
        pk = float(C[i])
        if pk < cfg.clean_floor:
            break
        l, r, c = i
        z0, z1 = max(l-hz, 0), min(l+hz+1, R.shape[0])
        y0, y1 = max(r-hxy, 0), min(r+hxy+1, R.shape[1])
        x0, x1 = max(c-hxy, 0), min(c+hxy+1, R.shape[2])
        kz0, ky0, kx0 = z0-(l-hz), y0-(r-hxy), x0-(c-hxy)
        R[z0:z1, y0:y1, x0:x1] -= pk * K[kz0:kz0+z1-z0, ky0:ky0+y1-y0, kx0:kx0+x1-x0]
        atoms.append((float(l), float(r), float(c), pk))
    return atoms


def _model_tube(shape, atoms, K, hz, hxy):
    """Render the full tube model: sum of amp * K at each (fractional) atom position."""
    M = np.zeros(shape)
    for (l, r, c, a) in atoms:
        M += a * _render(shape, (l, r, c), K, hz, hxy).reshape(shape)
    return M


def _contrib(shape, l, r, c, a, K, hz, hxy):
    return a * _render(shape, (l, r, c), K, hz, hxy).reshape(shape)


def refine_tube(tube, atoms, K, cfg):
    """Joint Gauss-Newton refinement of every atom in a tube -> sub-voxel positions,
    amplitudes, 1-sigma errors, and a fit-quality score.

    Per sweep, per atom: subtract all OTHER atoms' model (kept incrementally as
    total - own contribution), then fit (amp, dz, dy, dx) of amp*K((z,y,x)-p) on the
    kernel-support patch. Errors from the standard nonlinear-LS covariance
    sigma^2 = resid_var * diag((J^T J)^-1); on noiseless data resid_var is model
    mismatch, i.e. an honest systematic error bar. Quality = normalised correlation of
    the (others-subtracted) patch with the fitted kernel."""
    hz = (K.shape[0]-1)//2
    hxy = (K.shape[1]-1)//2
    Kz, Ky, Kx = np.gradient(K)
    cur = [list(a) for a in atoms]                        # [l, r, c, amp]
    contribs = [_contrib(tube.shape, *a, K, hz, hxy) for a in atoms]
    total = np.sum(contribs, axis=0)

    def patch_and_bases(l, r, c):
        li, ri, ci = int(round(l)), int(round(r)), int(round(c))
        z0, z1 = max(li-hz, 0), min(li+hz+1, tube.shape[0])
        y0, y1 = max(ri-hxy, 0), min(ri+hxy+1, tube.shape[1])
        x0, x1 = max(ci-hxy, 0), min(ci+hxy+1, tube.shape[2])
        sl_ = (slice(z0, z1), slice(y0, y1), slice(x0, x1))
        pshape = (z1-z0, y1-y0, x1-x0)
        ctr = (l-z0, r-y0, c-x0)
        Kp = _render(pshape, ctr, K,  hz, hxy).reshape(pshape)
        Gz = _render(pshape, ctr, Kz, hz, hxy).reshape(pshape)
        Gy = _render(pshape, ctr, Ky, hz, hxy).reshape(pshape)
        Gx = _render(pshape, ctr, Kx, hz, hxy).reshape(pshape)
        return sl_, Kp, Gz, Gy, Gx

    for _sweep in range(cfg.refine_sweeps):
        for i in range(len(cur)):
            l, r, c, a = cur[i]
            if a <= 0:
                continue
            sl_, Kp, Gz, Gy, Gx = patch_and_bases(l, r, c)
            patch = (tube - total + contribs[i])[sl_]     # others removed, self kept
            # model m = a*K(u), u = x - p  =>  dm/dp = -a * K'(u)
            J = np.column_stack([Kp.ravel(), -a*Gz.ravel(), -a*Gy.ravel(), -a*Gx.ravel()])
            rvec = (patch - a*Kp).ravel()
            try:
                dp = np.linalg.solve(J.T @ J + 1e-9*np.eye(4), J.T @ rvec)
            except np.linalg.LinAlgError:
                continue
            da, dl, dr, dc = dp
            cur[i][3] = max(a + da, 0.0)
            cur[i][0] = float(np.clip(l + np.clip(dl, -1, 1), 0, tube.shape[0]-1))
            cur[i][1] = float(np.clip(r + np.clip(dr, -1, 1), 0, tube.shape[1]-1))
            cur[i][2] = float(np.clip(c + np.clip(dc, -1, 1), 0, tube.shape[2]-1))
            new = _contrib(tube.shape, *cur[i], K, hz, hxy)
            total += new - contribs[i]
            contribs[i] = new

    # tube-level noise floor for the error bars: on noiseless data the per-atom patch
    # residual after CLEAN+GN is near-zero (overfit) -> the honest positional uncertainty is
    # set by the UNMODELLED signal (kernel-vs-real mismatch), estimated over the column core.
    hw = tube.shape[1] // 2
    core = (slice(None), slice(max(hw-hxy, 0), hw+hxy+1), slice(max(hw-hxy, 0), hw+hxy+1))
    core_res = (tube - total)[core]
    rvar = float(np.mean(core_res**2)) + 1e-12

    # final pass: covariance + quality per surviving atom
    out = []
    for i in range(len(cur)):
        l, r, c, a = cur[i]
        if a <= 1e-6:
            continue
        sl_, Kp, Gz, Gy, Gx = patch_and_bases(l, r, c)
        patch = (tube - total + contribs[i])[sl_]
        J = np.column_stack([Kp.ravel(), -a*Gz.ravel(), -a*Gy.ravel(), -a*Gx.ravel()])
        try:
            cov = rvar * np.linalg.inv(J.T @ J + 1e-9*np.eye(4))
            sl, sr, sc = np.sqrt(np.clip(np.diag(cov)[1:], 0, None))
        except np.linalg.LinAlgError:
            sl = sr = sc = np.nan
        den = float(np.linalg.norm(patch) * np.linalg.norm(a*Kp))
        quality = float((patch * (a*Kp)).sum())/den if den > 0 else 0.0
        out.append(dict(l=l, r=r, c=c, amp=a, sl=sl, sr=sr, sc=sc, quality=quality))
    return out


def _corr_with_kernel(tube, others_model, l, r, c, a, K, hz, hxy):
    """Normalised correlation of the (others-subtracted) patch with kernel K at (l,r,c)."""
    li, ri, ci = int(round(l)), int(round(r)), int(round(c))
    z0, z1 = max(li-hz, 0), min(li+hz+1, tube.shape[0])
    y0, y1 = max(ri-hxy, 0), min(ri+hxy+1, tube.shape[1])
    x0, x1 = max(ci-hxy, 0), min(ci+hxy+1, tube.shape[2])
    patch = (tube - others_model)[z0:z1, y0:y1, x0:x1]
    Kp = _render(patch.shape, (l-z0, r-y0, c-x0), K, hz, hxy).reshape(patch.shape)
    den = float(np.linalg.norm(patch) * np.linalg.norm(Kp))
    return float((patch*Kp).sum())/den if den > 0 else 0.0


def _kmeans1d(vals, k=3, iters=50):
    """Tiny 1-D k-means (no sklearn). Returns sorted cluster centres + boundaries."""
    cen = np.percentile(vals, np.linspace(10, 90, k))
    for _ in range(iters):
        lab = np.argmin(np.abs(vals[:, None] - cen[None, :]), axis=1)
        new = np.array([vals[lab == j].mean() if (lab == j).any() else cen[j] for j in range(k)])
        if np.allclose(new, cen):
            break
        cen = new
    cen = np.sort(cen)
    bounds = 0.5*(cen[:-1] + cen[1:])
    return cen, bounds


def find_atoms_v2(V, cfg, dx, kernels):
    """The v2 blind finder. kernels = {82: K_Pb, [22: K_Ti]} (measured, from psf.species_kernels).

    Returns a record array: row, col, layer (recon frame, sub-voxel), z_A, amp,
    sx_A, sy_A, sz_A (1-sigma), quality, species (best guess Z), col_id."""
    KPb = _unit_norm(crop_kernel_inplane(kernels[82]))
    KTi = _unit_norm(crop_kernel_inplane(kernels[22])) if 22 in kernels else None
    hz = (KPb.shape[0]-1)//2
    hxy = (KPb.shape[1]-1)//2
    l0 = int(round(cfg.trim_z_A[0]/cfg.dz))
    l1 = min(int(round(cfg.trim_z_A[1]/cfg.dz)), V.shape[0])
    HW = cfg.tube_halfwidth_px

    seeds, _ = detect_columns(V, cfg, dx)
    raw = []                                             # per-atom dicts + tube refs
    for cid, (r0, c0, _b) in enumerate(seeds):
        ri, ci = int(round(r0)), int(round(c0))
        if not (HW <= ri < V.shape[1]-HW and HW <= ci < V.shape[2]-HW):
            continue
        tube = np.clip(V[l0:l1, ri-HW:ri+HW+1, ci-HW:ci+HW+1], 0, None)
        atoms = clean_tube(tube, KPb, cfg)
        if not atoms:
            continue
        refined = refine_tube(tube, atoms, KPb, cfg)
        # Ti-kernel correlation (species shape evidence), others-subtracted via the
        # incremental total model (own contribution added back per atom)
        contribs = [_contrib(tube.shape, d["l"], d["r"], d["c"], d["amp"], KPb, hz, hxy)
                    for d in refined]
        total = np.sum(contribs, axis=0) if contribs else 0.0
        for i, d in enumerate(refined):
            om = total - contribs[i]
            d["corr_pb"] = d["quality"]
            d["corr_ti"] = (_corr_with_kernel(tube, om, d["l"], d["r"], d["c"], d["amp"],
                                              KTi, (KTi.shape[0]-1)//2, (KTi.shape[1]-1)//2)
                            if KTi is not None else np.nan)
            d["cid"] = cid
            d["ri"], d["ci"], d["l0"] = ri, ci, l0
            raw.append(d)

    # ---- junk cut by FIT QUALITY (not amplitude) ----
    # Use the Pb-kernel fit correlation (`quality`): it is the honest "does this look like an
    # atom" score. (corr_ti is only for species SHAPE classification -- folding it into the
    # cut admitted shape-plausible junk and cost ~5% precision.)
    keep = [d for d in raw if d["quality"] >= cfg.quality_min_corr
            and d["amp"] > 0 and np.isfinite(d["sl"])]

    # ---- species by blind amplitude bands + kernel-shape veto ----
    if keep:
        la = np.log(np.array([d["amp"] for d in keep]))
        _, bounds = _kmeans1d(la, k=3)
        b1, b2 = bounds                                  # O/Ti and Ti/Pb boundaries (log-amp)
        for d in keep:
            v = np.log(d["amp"])
            sp = 82 if v >= b2 else (22 if v >= b1 else 8)
            # shape veto: light atoms have the broader (Ti) response -- demote a boundary
            # "Pb" whose patch clearly prefers the Ti kernel shape
            if sp == 82 and np.isfinite(d["corr_ti"]) and d["corr_ti"] > d["corr_pb"] + 0.05 \
               and v < b2 + 0.5*(la.max()-b2):
                sp = 22
            d["species"] = sp

    # error bars: formal CRB (from refine_tube) in quadrature with the resolution floors
    fxy, fz = cfg.sigma_floor_xy_A, cfg.sigma_floor_z_A
    dt = np.dtype([("row", float), ("col", float), ("layer", float), ("z_A", float),
                   ("amp", float), ("sx_A", float), ("sy_A", float), ("sz_A", float),
                   ("quality", float), ("species", int), ("col_id", int)])
    recs = [(d["ri"]-HW+d["r"], d["ci"]-HW+d["c"], d["l0"]+d["l"], (d["l0"]+d["l"])*cfg.dz,
             d["amp"],
             np.hypot(d["sc"]*dx, fxy), np.hypot(d["sr"]*dx, fxy), np.hypot(d["sl"]*cfg.dz, fz),
             d["quality"], d["species"], d["cid"]) for d in keep]
    return np.array(recs, dtype=dt), seeds


# ================================================================= v3: lattice-aware + guided
def preprocess(V, cfg, dx):
    """Stage C: remove the smooth per-layer background (depth haze).

    Per layer, subtract a heavy Gaussian blur (bg_smooth_A ~ the column pitch, >> atom
    width) and clip >= 0. The blur barely samples the razor-thin atomic peaks, so their
    amplitude survives while the broad haze that biases weak-O amplitudes and feeds junk
    detections is removed."""
    from scipy.ndimage import gaussian_filter
    sig = cfg.bg_smooth_A / dx
    W = np.empty_like(V)
    for l in range(V.shape[0]):
        W[l] = V[l] - gaussian_filter(V[l], sig)
    return np.clip(W, 0, None)


def _comb_fit(zs, amps, cfg):
    """Fit the BRIGHT comb of a column: period p (~3.9 A) and phase phi (A).

    Period from the median spacing of the brighter half; phase from the amp^2-weighted
    circular mean of ALL atoms (Ti ~1.5 dominates O ~0.65 quadratically, so the phase
    locks onto the bright comb even when more O than Ti were detected)."""
    zs = np.asarray(zs, float); amps = np.asarray(amps, float)
    p = cfg.comb_period_A
    if len(zs) >= 4:
        zb = np.sort(zs[amps >= np.median(amps)])
        d = np.diff(zb)
        d = d[(d > 0.6*p) & (d < 1.4*p)]
        if len(d) >= 2:
            p = float(np.median(d))
    w = amps**2
    ang = 2*np.pi*zs/p
    phi = (np.angle(np.sum(w*np.exp(1j*ang))) % (2*np.pi)) * p/(2*np.pi)
    return p, float(phi)


def _comb_dist(z, phi, p):
    r = (z - phi) % p
    return min(r, p - r)


def _patch_at(arr, l, r, c, hz, hxy):
    """Clipped patch slices around a (possibly fractional) centre; returns (view, ctr)."""
    li, ri, ci = int(round(l)), int(round(r)), int(round(c))
    z0, z1 = max(li-hz, 0), min(li+hz+1, arr.shape[0])
    y0, y1 = max(ri-hxy, 0), min(ri+hxy+1, arr.shape[1])
    x0, x1 = max(ci-hxy, 0), min(ci+hxy+1, arr.shape[2])
    return (slice(z0, z1), slice(y0, y1), slice(x0, x1)), (l-z0, r-y0, c-x0)


def _fit_single_guided(resid, K, Kgrads, l_s, r_s, c_s, cfg, dx, min_corr=None):
    """Fit ONE atom on the residual at a predicted comb slot: integer matched-filter search
    inside the position gates, then 2 Gauss-Newton refinements. Returns a dict or None.
    min_corr overrides cfg.guided_min_corr (species-dependent bar)."""
    hz = (K.shape[0]-1)//2; hxy = (K.shape[1]-1)//2
    Kz, Ky, Kx = Kgrads
    gz = max(1, int(np.ceil(cfg.guided_gate_z_A/cfg.dz)))
    gxy = max(1, int(round(cfg.guided_gate_xy_A/dx)))
    best = None
    for dl in range(-gz, gz+1):
        for dr in range(-gxy, gxy+1, 2):
            for dc in range(-gxy, gxy+1, 2):
                l, r, c = int(round(l_s))+dl, int(round(r_s))+dr, int(round(c_s))+dc
                if not (0 <= l < resid.shape[0] and 0 <= r < resid.shape[1] and 0 <= c < resid.shape[2]):
                    continue
                sl_, ctr = _patch_at(resid, l, r, c, hz, hxy)
                patch = resid[sl_]
                Kp = _render(patch.shape, ctr, K, hz, hxy).reshape(patch.shape)
                den = float((Kp*Kp).sum())
                if den <= 0:
                    continue
                a = float((patch*Kp).sum())/den
                if best is None or a > best[0]:
                    best = (a, float(l), float(r), float(c))
    if best is None or best[0] <= 0:
        return None
    a, l, r, c = best
    for _ in range(2):                                    # GN sub-voxel refine
        sl_, ctr = _patch_at(resid, l, r, c, hz, hxy)
        patch = resid[sl_]
        Kp = _render(patch.shape, ctr, K,  hz, hxy).reshape(patch.shape)
        Gz = _render(patch.shape, ctr, Kz, hz, hxy).reshape(patch.shape)
        Gy = _render(patch.shape, ctr, Ky, hz, hxy).reshape(patch.shape)
        Gx = _render(patch.shape, ctr, Kx, hz, hxy).reshape(patch.shape)
        J = np.column_stack([Kp.ravel(), -a*Gz.ravel(), -a*Gy.ravel(), -a*Gx.ravel()])
        try:
            dp = np.linalg.solve(J.T @ J + 1e-9*np.eye(4), J.T @ (patch - a*Kp).ravel())
        except np.linalg.LinAlgError:
            break
        a = max(a + dp[0], 0.0)
        l = float(np.clip(l + np.clip(dp[1], -1, 1), 0, resid.shape[0]-1))
        r = float(np.clip(r + np.clip(dp[2], -1, 1), 0, resid.shape[1]-1))
        c = float(np.clip(c + np.clip(dp[3], -1, 1), 0, resid.shape[2]-1))
    # gates + quality on the final fit
    if a <= 0 or abs(l-l_s)*cfg.dz > cfg.guided_gate_z_A or \
       np.hypot(r-r_s, c-c_s)*dx > 2.0*cfg.guided_gate_xy_A:
        return None
    sl_, ctr = _patch_at(resid, l, r, c, hz, hxy)
    patch = resid[sl_]
    Kp = _render(patch.shape, ctr, K, hz, hxy).reshape(patch.shape)
    den = float(np.linalg.norm(patch)*np.linalg.norm(a*Kp))
    quality = float((patch*(a*Kp)).sum())/den if den > 0 else 0.0
    if quality < (min_corr if min_corr is not None else cfg.guided_min_corr):
        return None
    # covariance for the error bars (same machinery as refine_tube's final pass)
    Gz = _render(patch.shape, ctr, Kz, hz, hxy).reshape(patch.shape)
    Gy = _render(patch.shape, ctr, Ky, hz, hxy).reshape(patch.shape)
    Gx = _render(patch.shape, ctr, Kx, hz, hxy).reshape(patch.shape)
    J = np.column_stack([Kp.ravel(), -a*Gz.ravel(), -a*Gy.ravel(), -a*Gx.ravel()])
    rvec = (patch - a*Kp).ravel()
    rvar = float(rvec @ rvec) / max(rvec.size - 4, 1)
    try:
        cov = rvar * np.linalg.inv(J.T @ J + 1e-9*np.eye(4))
        sl, sr, sc = np.sqrt(np.clip(np.diag(cov)[1:], 0, None))
    except np.linalg.LinAlgError:
        sl = sr = sc = np.nan
    return dict(l=l, r=r, c=c, amp=a, sl=sl, sr=sr, sc=sc, quality=quality)


def _lean_predict(kept, l):
    """Predict (r, c) at layer l from the column's own atoms (linear fit vs layer)."""
    ls = np.array([d["l"] for d in kept]); rs = np.array([d["r"] for d in kept])
    cs = np.array([d["c"] for d in kept])
    if len(kept) >= 3:
        pr = np.polyval(np.polyfit(ls, rs, 1), l)
        pc = np.polyval(np.polyfit(ls, cs, 1), l)
    else:
        pr, pc = rs.mean(), cs.mean()
    return float(pr), float(pc)


def find_atoms_v3(V, cfg, dx, kernels):
    """v3 blind finder: preprocess -> tube CLEAN + GN refine (v2 core) -> junk cut ->
    LATTICE-AWARE species (column typing + B-O comb parity) -> GUIDED re-detection at
    empty comb slots (position prior from the column's OWN atoms -- still no GT).

    Returns records with v2 fields + guided (0=blind, 1=guided)."""
    Vp = preprocess(V, cfg, dx) if cfg.preprocess_bg else V
    KPb = _unit_norm(crop_kernel_inplane(kernels[82]))
    KTi = _unit_norm(crop_kernel_inplane(kernels[22])) if 22 in kernels else KPb
    KPb_g = np.gradient(KPb); KTi_g = np.gradient(KTi)
    hz = (KPb.shape[0]-1)//2; hxy = (KPb.shape[1]-1)//2
    l0 = int(round(cfg.trim_z_A[0]/cfg.dz))
    l1 = min(int(round(cfg.trim_z_A[1]/cfg.dz)), V.shape[0])
    HW = cfg.tube_halfwidth_px

    # ---- v2 core per tube, keeping the tube for the guided stage ----
    seeds, _ = detect_columns(Vp, cfg, dx)
    tubes = []
    for cid, (r0, c0, _b) in enumerate(seeds):
        ri, ci = int(round(r0)), int(round(c0))
        if not (HW <= ri < V.shape[1]-HW and HW <= ci < V.shape[2]-HW):
            continue
        tube = np.clip(Vp[l0:l1, ri-HW:ri+HW+1, ci-HW:ci+HW+1], 0, None)
        atoms = clean_tube(tube, KPb, cfg)
        if not atoms:
            continue
        refined = refine_tube(tube, atoms, KPb, cfg)
        kept = [d for d in refined if d["quality"] >= cfg.quality_min_corr
                and d["amp"] > 0 and np.isfinite(d["sl"])]
        if kept:
            tubes.append(dict(cid=cid, ri=ri, ci=ci, tube=tube, kept=kept))

    # ---- Stage A: column typing (k-means on per-column amp p75; measured zero-overlap) ----
    p75 = np.array([np.percentile([d["amp"] for d in t["kept"]], 75) for t in tubes])
    if len(p75) >= 3:
        _, cb = _kmeans1d(np.log(p75), k=3)
        cb1, cb2 = cb                                    # O/BO and BO/A boundaries
    else:
        cb1, cb2 = np.log(1.0), np.log(2.5)
    for t, v in zip(tubes, np.log(p75)):
        t["ctype"] = "A" if v >= cb2 else ("BO" if v >= cb1 else "O")

    # ---- Stage A: species per atom + comb parameters per column ----
    for t in tubes:
        zs = np.array([d["l"] for d in t["kept"]]) * cfg.dz     # tube-frame z (A)
        amps = np.array([d["amp"] for d in t["kept"]])
        p, phi = _comb_fit(zs, amps, cfg)
        t["p"], t["phi"] = p, phi
        if t["ctype"] == "A":
            # an A column contains ONLY Pb: quality-passing blobs at a fraction of the
            # column amplitude are haze/duplicates, not atoms -- drop them (measured: ~60
            # such junk detections otherwise inherit the Pb label).
            med = float(np.median(amps))
            t["kept"] = [d for d in t["kept"] if d["amp"] >= 0.25*med]
            for d in t["kept"]:
                d["species"] = 82
        elif t["ctype"] == "O":
            for d in t["kept"]:
                d["species"] = 8
        else:
            # B-O: amplitude split with a dead zone, then LOCAL parity for the ambiguous.
            # A global comb phase accumulates period error toward the column ends and flips
            # parity there (measured: nearly all residual Ti<->O confusion sat at z>54 A).
            # Local parity is translation-invariant: distance of (z_i - z_j) from multiples
            # of p, against the nearest CONFIDENT atoms -- no phase to drift.
            la = np.log(amps)
            if len(la) >= 4:
                cen2, bnd2 = _kmeans1d(la, k=2)
                b, spread = bnd2[0], max(cen2[1]-cen2[0], 0.2)
            else:
                b, spread = np.log(1.0), 0.8
            dead = 0.25 * spread
            conf = np.where(np.abs(la - b) >= dead)[0]
            lab = np.where(la >= b, 22, 8)                    # amplitude vote
            for i in range(len(zs)):
                if abs(la[i] - b) >= dead or len(conf) == 0:
                    continue                                   # confident by amplitude
                d_ti = [(_comb_dist(zs[i]-zs[j], 0, p)) for j in conf if lab[j] == 22]
                d_o = [(_comb_dist(zs[i]-zs[j], 0, p)) for j in conf if lab[j] == 8]
                dt_ = min(d_ti) if d_ti else np.inf
                do_ = min(d_o) if d_o else np.inf
                if np.isfinite(min(dt_, do_)):
                    lab[i] = 22 if dt_ <= do_ else 8
            for d, s in zip(t["kept"], lab):
                d["species"] = int(s)

    # ---- Stage B: guided re-detection at empty comb slots ----
    Lz = (l1 - l0 - 1) * cfg.dz
    for t in tubes:
        total = np.sum([_contrib(t["tube"].shape, d["l"], d["r"], d["c"], d["amp"],
                                 KPb, hz, hxy) for d in t["kept"]], axis=0)
        resid = t["tube"] - total
        p = t["p"]
        species_list = ([82] if t["ctype"] == "A" else [8] if t["ctype"] == "O" else [22, 8])
        species_list = [s for s in species_list if s in cfg.guided_species]
        t["guided"] = []
        for sp in species_list:
            K, Kg = (KPb, KPb_g) if sp == 82 else (KTi, KTi_g)
            # species-dependent evidence bar: the lower guided bar is only justified for the
            # contrast-limited O (position prior compensates weak signal); heavy atoms are
            # already found blind at ~97%, so guided Pb/Ti must clear the BLIND bar or they
            # mostly add junk (measured: +36 spurious Pb at the low bar).
            bar = cfg.guided_min_corr if sp == 8 else cfg.quality_min_corr
            same_amp = [d["amp"] for d in t["kept"] if d["species"] == sp]
            amp_med = float(np.median(same_amp)) if same_amp else None
            same_z = [d["l"]*cfg.dz for d in t["kept"] if d["species"] == sp]
            # LOCAL slot prediction (translation-invariant, no global phase to drift):
            # extend +-k*p from each same-species atom; O slots on B-O columns also from
            # the Ti atoms' positions +- p/2 when O anchors are sparse near the ends.
            base = list(same_z)
            if t["ctype"] == "BO" and sp == 8:
                ti_z = [d["l"]*cfg.dz for d in t["kept"] if d["species"] == 22]
                base += [z + p/2 for z in ti_z] + [z - p/2 for z in ti_z]
            if not base:
                continue
            nmax = int(np.ceil(Lz / p)) + 1
            cand = {}                                     # 0.5-A-binned, keep nearest base
            for zb in base:
                for k in range(-nmax, nmax + 1):
                    zc = zb + k*p
                    if -0.4 <= zc <= Lz + 0.4:
                        key = int(round(zc / 0.5))
                        if key not in cand or abs(k) < cand[key][1]:
                            cand[key] = (zc, abs(k))
            for z_s, _k in sorted(cand.values()):
                if same_z and min(abs(np.array(same_z) - z_s)) < cfg.slot_empty_A:
                    continue                              # slot already filled
                l_s = z_s / cfg.dz
                r_s, c_s = _lean_predict(t["kept"], l_s)
                fit = _fit_single_guided(resid, K, Kg, l_s, r_s, c_s, cfg, dx, min_corr=bar)
                if fit is None:
                    continue
                # lattice-consistency amplitude gate: a guided atom must look like its
                # column-mates. Kills low-amp junk AND cross-species catches (a guided-Ti
                # slot grabbing an O has ~0.4x the Ti median -> rejected, and vice versa).
                # Pb gets a lower bound of 0.25 (A-columns hold no confusable species, and
                # domain-wall/edge Pb can be dim -- this restores real Pb that the A-column
                # junk cut removed, at proper lattice slots only). Wide upper bound absorbs
                # the KTi-vs-KPb matched-filter scale difference.
                lo = 0.25 if sp == 82 else 0.45
                hi = 1.6 if sp == 8 else 2.2      # tight O upper bound: a dim Ti (~1.2)
                #                                   must NOT pass as O (measured leak)
                if amp_med is not None and not (lo*amp_med <= fit["amp"] <= hi*amp_med):
                    continue
                fit["species"] = sp
                t["guided"].append(fit)
                same_z.append((fit["l"])*cfg.dz)
                resid -= _contrib(t["tube"].shape, fit["l"], fit["r"], fit["c"],
                                  fit["amp"], K, hz, hxy)

    # ---- assemble records (per-species 68%-coverage floors; guided get the wider scale) ----
    fxy = cfg.sigma_floor_xy_A
    dt = np.dtype([("row", float), ("col", float), ("layer", float), ("z_A", float),
                   ("amp", float), ("sx_A", float), ("sy_A", float), ("sz_A", float),
                   ("quality", float), ("species", int), ("col_id", int), ("guided", int)])
    recs = []
    for t in tubes:
        for d, g in [(d, 0) for d in t["kept"]] + [(d, 1) for d in t["guided"]]:
            s = cfg.guided_sigma_scale if g else 1.0
            z_A = (l0 + d["l"]) * cfg.dz
            if z_A > cfg.exit_band_z_A:                  # exit-artifact band: wider bars
                s *= cfg.exit_sigma_scale
            fz = cfg.sigma_floor_z_species.get(int(d["species"]), cfg.sigma_floor_z_A)
            recs.append((t["ri"]-HW+d["r"], t["ci"]-HW+d["c"], l0+d["l"], z_A,
                         d["amp"],
                         np.hypot(d["sc"]*dx, fxy)*s, np.hypot(d["sr"]*dx, fxy)*s,
                         np.hypot(d["sl"]*cfg.dz, fz)*s,
                         d["quality"], d["species"], t["cid"], g))
    return np.array(recs, dtype=dt), seeds


# ---------------------------------------------------------------- export
def export_atoms(found, al, cfg, out_prefix):
    """Write the found atoms as CSV + an ASE extxyz object in the prepared-cell frame.

    Uses ONLY the calibration constants of the recon<->prepared-cell map (the affine
    in-plane map + depth registration) -- no GT atom positions. The frame matches the GT
    VASP after the sec.11 prep, so the exported object overlays the model directly."""
    X, Y = al.index_to_site(found["row"], found["col"])
    Zc = al.layer_to_z(found["layer"])
    sym = {82: "Pb", 22: "Ti", 8: "O"}
    guided = found["guided"] if "guided" in found.dtype.names else np.zeros(len(found), int)
    # CSV
    csv_path = out_prefix + ".csv"
    with open(csv_path, "w") as f:
        f.write("element,X_A,Y_A,Z_A,sigma_x_A,sigma_y_A,sigma_z_A,amplitude,quality,col_id,guided\n")
        for i in range(len(found)):
            f.write(f"{sym.get(int(found['species'][i]),'X')},{X[i]:.4f},{Y[i]:.4f},{Zc[i]:.4f},"
                    f"{found['sx_A'][i]:.4f},{found['sy_A'][i]:.4f},{found['sz_A'][i]:.4f},"
                    f"{found['amp'][i]:.4f},{found['quality'][i]:.3f},{found['col_id'][i]},"
                    f"{int(guided[i])}\n")
    # ASE extxyz (per-atom arrays for the uncertainties)
    import ase, ase.io
    side, height = 70.008, 74.0                          # prepared-box constants (sec.3.1)
    atoms = ase.Atoms([sym.get(int(s), "X") for s in found["species"]],
                      positions=np.column_stack([X, Y, Zc]),
                      cell=[side, side, height], pbc=False)
    atoms.set_array("sigma_A", np.column_stack([found["sx_A"], found["sy_A"], found["sz_A"]]))
    atoms.set_array("amplitude", found["amp"].astype(float))
    atoms.set_array("quality", found["quality"].astype(float))
    atoms.set_array("guided", guided.astype(int))
    xyz_path = out_prefix + ".extxyz"
    ase.io.write(xyz_path, atoms, format="extxyz")
    return csv_path, xyz_path


# ---------------------------------------------------------------- orchestration
def find_atoms(V, cfg, dx, psf3d, method="spike"):
    """Blind finder. method in {'spike','raw'}. Returns a record array of found atoms in
    the recon frame: row, col, layer (z in layer units), z_A, amp, col_id."""
    k1d, zoff = _psf.axial_kernel(psf3d, navg=cfg.find_profile_navg)
    seeds, dm = detect_columns(V, cfg, dx)
    recs = []
    for cid, (r0, c0, _) in enumerate(seeds):
        pr, pc, prof = track_column(V, r0, c0, cfg)
        if method == "spike":
            atoms = spike_deconv(prof, k1d, zoff, cfg)
        else:
            atoms = raw_peak_z(prof, cfg, k1d)
        for (zl, amp) in atoms:
            li = int(np.clip(round(zl), 0, V.shape[0]-1))
            recs.append((pr[li], pc[li], zl, zl*cfg.dz, amp, cid))
    dt = np.dtype([("row", float), ("col", float), ("layer", float),
                   ("z_A", float), ("amp", float), ("col_id", int)])
    return np.array(recs, dtype=dt), seeds, dm
