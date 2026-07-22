#!/usr/bin/env python
"""@file validate.py
@brief Validation harness -- score found atoms and fitted amplitudes against ground truth.

Every atom is scored AT its GT-mapped position (align.py); the detection floor is a NULL
built by running the same PSF fit at off-lattice sites. Oxygen is split into its two
physically distinct failure modes:
  * AXIALLY-OVERLAPPED O -- a heavy atom at nearly the same (x,y) within a short z-gap
    (apical O ~1.9 A under Ti); unresolvable-from-Ti, the HARD case.
  * IN-PLANE-ISOLATED O -- nearest heavy displaced in-plane (equatorial / pure-O columns);
    resolvable in-plane, a contrast/SNR case.
health_warnings flags the canaries (confusion, sigma-coverage, collapsed recall).
"""
from __future__ import annotations
import numpy as np
from . import align as _align

OVERLAP_XY_A = 0.7      # "same column" in-plane tolerance
OVERLAP_Z_A = 2.6      # z-gap under which a heavy neighbour axially buries the O


# ---------------------------------------------------------------- blind-finder scoring
def match_found_to_gt(found, pos, Z, al, cfg):
    """Greedy nearest-neighbour match of blind-found atoms to GT, in the RECON frame.

    found: record array from find.find_atoms (row, col, layer, ...). GT atoms are mapped
    into the recon frame with al.site_to_index. A match requires in-plane <= match_tol_xy
    AND |dz| <= match_tol_z. Returns dicts of arrays: per found atom -> matched GT idx
    (or -1) and dz (A); plus GT-side recall bookkeeping. Brightest-first greedy so strong
    atoms claim their GT partner before weak/spurious ones."""
    win = _align.in_window(pos, cfg)
    gi = np.where(win)[0]
    gr, gc, gl = al.site_to_index(pos[gi, 0], pos[gi, 1], pos[gi, 2])
    gx, gy, gz = gc * al.dx, gr * al.dx, gl * cfg.dz          # A units in recon frame
    claimed = np.zeros(len(gi), bool)
    fx = found["col"] * al.dx; fy = found["row"] * al.dx; fz = found["layer"] * cfg.dz
    order = np.argsort(-found["amp"])
    match_gi = np.full(len(found), -1)
    match_dz = np.full(len(found), np.nan)
    match_dx = np.full(len(found), np.nan)
    match_dy = np.full(len(found), np.nan)
    txy, tz = cfg.match_tol_xy_A, cfg.match_tol_z_A
    for f in order:
        d_xy = np.hypot(gx - fx[f], gy - fy[f])
        d_z = np.abs(gz - fz[f])
        ok = (~claimed) & (d_xy <= txy) & (d_z <= tz)
        if ok.any():
            cand = np.where(ok)[0]
            j = cand[np.argmin(d_xy[cand] + d_z[cand])]
            claimed[j] = True
            match_gi[f] = gi[j]
            match_dz[f] = (gz[j] - fz[f])
            match_dx[f] = (gx[j] - fx[f])
            match_dy[f] = (gy[j] - fy[f])
    return dict(match_gi=match_gi, match_dz=match_dz, match_dx=match_dx, match_dy=match_dy,
                gt_idx=gi, gt_claimed=claimed)


def finder_report(found, pos, Z, al, cfg, olabel=None):
    """Per-species precision / recall / z-error RMS for a blind-found atom set.

    If `found` carries the v2 fields (species, sx_A/sy_A/sz_A), the report also includes
    the BULK recall split (cfg.bulk_z_A band), the species CONFUSION matrix, and the
    ERROR-BAR CALIBRATION (median |GT error| / sigma per axis; ~1 = honest sigmas)."""
    m = match_found_to_gt(found, pos, Z, al, cfg)
    gi, claimed = m["gt_idx"], m["gt_claimed"]
    matched = m["match_gi"] >= 0
    rep = {"n_found": int(len(found)), "n_matched": int(matched.sum()),
           "precision": float(matched.mean()) if len(found) else 0.0,
           "z_rms_A": float(np.sqrt(np.nanmean(m["match_dz"][matched]**2))) if matched.any() else np.nan,
           "z_bias_A": float(np.nanmean(m["match_dz"][matched])) if matched.any() else np.nan}
    if matched.any():
        rep["xy_rms_A"] = float(np.sqrt(np.nanmean(m["match_dx"][matched]**2
                                                   + m["match_dy"][matched]**2)))
    # bulk mask on GT (depth band away from entrance/exit artifacts)
    _, _, gl = al.site_to_index(pos[gi, 0], pos[gi, 1], pos[gi, 2])
    gz_A = gl * cfg.dz
    bulk = (gz_A >= cfg.bulk_z_A[0]) & (gz_A <= cfg.bulk_z_A[1])
    # per-species recall + species-resolved z-error (using the matched found atoms' GT Z)
    fz_species = np.array([Z[g] if g >= 0 else -1 for g in m["match_gi"]])
    for zz, nm in [(82, "Pb"), (22, "Ti"), (8, "O")]:
        sel = Z[gi] == zz
        gt_n = int(sel.sum())
        rec_n = int((sel & claimed).sum())
        sm = (fz_species == zz) & matched
        rep[nm] = dict(gt=gt_n, recall=(rec_n/gt_n if gt_n else np.nan),
                       n_matched=int(sm.sum()),
                       z_rms_A=float(np.sqrt(np.nanmean(m["match_dz"][sm]**2))) if sm.any() else np.nan,
                       recall_bulk=(float((sel & bulk & claimed).sum() / max((sel & bulk).sum(), 1))
                                    if (sel & bulk).any() else np.nan),
                       recall_edge=(float((sel & ~bulk & claimed).sum() / max((sel & ~bulk).sum(), 1))
                                    if (sel & ~bulk).any() else np.nan))
    if olabel is not None:
        for lab, nm in [(1, "O_axial_overlap"), (2, "O_inplane_isolated")]:
            sel = (Z[gi] == 8) & (olabel[gi] == lab)
            gt_n = int(sel.sum()); rec_n = int((sel & claimed).sum())
            rep[nm] = dict(gt=gt_n, recall=(rec_n/gt_n if gt_n else np.nan))
    # ---- v2 extras: confusion matrix + error-bar calibration ----
    if "species" in found.dtype.names:
        conf = {}
        for pz in (82, 22, 8):
            for tz_ in (82, 22, 8):
                conf[f"{pz}->{tz_}"] = int(((found["species"] == pz) & (fz_species == tz_)).sum())
            conf[f"{pz}->none"] = int(((found["species"] == pz) & ~matched).sum())
        rep["confusion"] = conf
    if "sz_A" in found.dtype.names and matched.any():
        cal, cov = {}, {}
        for err, sig, ax in [(m["match_dx"], found["sx_A"], "x"),
                             (m["match_dy"], found["sy_A"], "y"),
                             (m["match_dz"], found["sz_A"], "z")]:
            ok = matched & np.isfinite(err) & (sig > 1e-6)
            cal[ax] = float(np.median(np.abs(err[ok]) / sig[ok])) if ok.any() else np.nan
            # COVERAGE is the honest error-bar metric: fraction of atoms whose true error
            # is within +-1 sigma (target ~0.68). Median-ratio calibration alone hides
            # heavy tails (measured: median ratio 1.0 but coverage only ~50%).
            cov[ax] = float((np.abs(err[ok]) <= sig[ok]).mean()) if ok.any() else np.nan
        rep["sigma_calibration"] = cal
        rep["sigma_coverage_1s"] = cov
        for zz, nm in [(82, "Pb"), (22, "Ti"), (8, "O")]:
            s = matched & (found["species"] == zz) & np.isfinite(m["match_dz"])
            if s.any():
                rep[nm]["z_cov_1s"] = float((np.abs(m["match_dz"][s]) <= found["sz_A"][s]).mean())
    return rep, m


# ---------------------------------------------------------------- health metrics
def confusion_rate(rep):
    """Off-diagonal fraction of the species confusion matrix (0 if unavailable)."""
    cf = rep.get("confusion")
    if not cf:
        return float("nan")
    off = sum(cf[f"{a}->{b}"] for a in (82, 22, 8) for b in (82, 22, 8) if a != b)
    diag = sum(cf[f"{a}->{a}"] for a in (82, 22, 8))
    return off / (off + diag) if (off + diag) else float("nan")


def confusion_line(rep):
    """One-line confusion summary for run-time health printing."""
    cf = rep.get("confusion")
    if not cf:
        return "confusion: n/a"
    r = confusion_rate(rep)
    return (f"confusion {r:.1%}  "
            f"[Ti->O {cf['22->8']}, O->Ti {cf['8->22']}, Ti->Pb {cf['22->82']}, "
            f"Pb->Ti {cf['82->22']}]")


def health_warnings(rep, conf_max=0.05, cov_lo=0.55, cov_hi=0.95, recall_min=0.30):
    """Run-time sanity checks that PRECISION DOES NOT GIVE YOU.

    Measured failure mode (noise sweep): precision read 0.95 while Ti/O recall was 0% and
    81% of species labels were wrong -- the surviving bright Pb were placed correctly, so
    precision looked fine. Confusion rate and sigma-coverage are the canaries."""
    warns = []
    r = confusion_rate(rep)
    if r == r and r > conf_max:
        warns.append(f"species confusion {r:.1%} > {conf_max:.0%} — labels unreliable "
                     f"(check depth registration: a ~1 A OFF error swaps Ti/O)")
    # NOTE: coverage of the MODEL sigma is intentionally NOT a health check any more --
    # the model sigma is a pre-calibration quantity (joint CRLB (+) kernel mismatch) and is
    # expected to under-cover before the conformal step corrects it. Coverage of the
    # CALIBRATED interval is audited in uncertainty.coverage_table instead. What CAN warn
    # here is a degenerate conformal q (a stratum whose calibrated interval fails to hit
    # its own target on held-in data implies too few points / a broken stratum).
    for s in ("Pb", "Ti", "O"):
        if s in rep and rep[s].get("recall", 1.0) < recall_min:
            warns.append(f"{s} recall {rep[s]['recall']:.0%} < {recall_min:.0%} — species "
                         f"effectively not detected")
    return warns


# ---------------------------------------------------------------- O classification
def classify_oxygen(pos, Z, cfg):
    """Return an int label per atom: 1 = axially-overlapped O, 2 = in-plane-isolated O,
    0 = not an in-field O. Uses 3-D geometry of the true structure."""
    lab = np.zeros(len(pos), int)
    win = _align.in_window(pos, cfg)
    heavy = pos[(Z == 82) | (Z == 22)]
    oi = np.where(win & (Z == 8))[0]
    for j in oi:
        p = pos[j]
        dxy = np.hypot(heavy[:, 0] - p[0], heavy[:, 1] - p[1])
        dz = np.abs(heavy[:, 2] - p[2])
        buried = (dxy < OVERLAP_XY_A) & (dz < OVERLAP_Z_A) & (dz > 1e-3)
        lab[j] = 1 if buried.any() else 2
    return lab


# ---------------------------------------------------------------- amplitude stats
def amplitude_summary(rec):
    """Per-species beta stats. Returns dict keyed by label -> (n, median, mean, std, q)."""
    out = {}
    groups = {"Pb": rec["Z"] == 82, "Ti": rec["Z"] == 22, "O": rec["Z"] == 8,
              "null": rec["Z"] == -1}
    for nm, m in groups.items():
        b = rec["beta"][m]
        if len(b) == 0:
            continue
        out[nm] = dict(n=int(len(b)), median=float(np.median(b)), mean=float(b.mean()),
                       std=float(b.std()), p90=float(np.percentile(b, 90)),
                       frac_pos=float(np.mean(b > 1e-6)))
    return out


def amplitude_vs_Z(rec):
    """Median beta vs Z (with null at Z=0). Fit beta = a*Z through the light species to
    show whether O follows the heavier atoms' Z-scaling."""
    pts = []
    for zz in (0, 8, 22, 82):
        m = rec["Z"] == (-1 if zz == 0 else zz)
        if m.any():
            pts.append((zz, float(np.median(rec["beta"][m]))))
    pts = np.array(pts)
    # slope through Ti & O (light, unsaturated) referenced to null
    zt = {int(z): b for z, b in pts}
    slope = ((zt.get(22, np.nan) - zt.get(0, 0)) / 22.0) if 22 in zt else np.nan
    return pts, slope


# ---------------------------------------------------------------- detection
def _auc(pos_b, neg_b):
    """AUC = P(pos beta > neg beta) via rank statistic (Mann-Whitney)."""
    if len(pos_b) == 0 or len(neg_b) == 0:
        return np.nan
    allv = np.concatenate([pos_b, neg_b])
    order = allv.argsort(kind="mergesort")
    ranks = np.empty_like(order, float)
    ranks[order] = np.arange(1, len(allv) + 1)
    # average ranks for ties
    _, inv, counts = np.unique(allv, return_inverse=True, return_counts=True)
    csum = np.cumsum(counts)
    avg = {i: (csum[i] - counts[i] + 1 + csum[i]) / 2.0 for i in range(len(counts))}
    ranks = np.array([avg[i] for i in inv])
    r_pos = ranks[:len(pos_b)].sum()
    n1, n2 = len(pos_b), len(neg_b)
    return (r_pos - n1 * (n1 + 1) / 2.0) / (n1 * n2)


def roc(pos_b, neg_b, n=200):
    """ROC curve sweeping threshold; returns fpr, tpr, auc."""
    if len(pos_b) == 0 or len(neg_b) == 0:
        return np.array([0, 1]), np.array([0, 1]), np.nan
    thr = np.linspace(min(pos_b.min(), neg_b.min()), max(pos_b.max(), neg_b.max()), n)
    tpr = np.array([(pos_b >= t).mean() for t in thr])
    fpr = np.array([(neg_b >= t).mean() for t in thr])
    return fpr, tpr, _auc(pos_b, neg_b)


def detection_report(rec, olabel, fpr_target=0.05):
    """Threshold at the null's (1-fpr_target) quantile; report per-species TPR + AUC.

    olabel is classify_oxygen()[rec['idx']] aligned to rec (pass via caller)."""
    nullb = rec["beta"][rec["Z"] == -1]
    if len(nullb) == 0:
        raise ValueError("no null sites in records; pass extra_sites to fit_amplitudes")
    tau = float(np.quantile(nullb, 1 - fpr_target))
    rep = {"tau": tau, "fpr_target": fpr_target, "null_n": int(len(nullb))}
    species = {"Pb": rec["Z"] == 82, "Ti": rec["Z"] == 22, "O_all": rec["Z"] == 8,
               "O_axial_overlap": (rec["Z"] == 8) & (olabel == 1),
               "O_inplane_isolated": (rec["Z"] == 8) & (olabel == 2)}
    for nm, m in species.items():
        b = rec["beta"][m]
        if len(b) == 0:
            continue
        rep[nm] = dict(n=int(len(b)), tpr=float((b >= tau).mean()),
                       auc=float(_auc(b, nullb)), median=float(np.median(b)))
    return rep
