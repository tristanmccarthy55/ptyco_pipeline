#!/usr/bin/env python
"""@file uncertainty.py
@brief Uncertainty quantification: model sigma -> calibrated prediction intervals.

Two stages, deliberately kept separate:

  MODEL sigma (from find.py, no ground truth): sigma^2 = sigma_fit^2 + sigma_kernel^2.
    sigma_fit is the JOINT Cramer-Rao bound (full-tube Fisher inverse, not the per-atom
    block, which understates overlapping atoms); sigma_kernel is the kernel-mismatch
    systematic (find.kernel_mismatch_sigma), the component that transfers without GT.

  CALIBRATION (this module, GT required): stratified split-conformal (Mondrian). Score
    s = |Delta|/sigma_model, empirical quantile per stratum, interval = q * sigma_model, so
    coverage holds per stratum by construction -- replacing the old per-species tuned floors.
    Strata = species x mode (blind/guided) x depth band.

The q-table is written to JSON; applying one from another dataset assumes exchangeability,
which is REPORTED (uncertainty_report), not assumed. Full derivation and measured VIF /
coverage tables: RESULTS.md §6.
"""
from __future__ import annotations
import json
import numpy as np

AXES = ("x", "y", "z")


# ---------------------------------------------------------------- strata
def stratum_of(rec, cfg):
    """(species, mode, depth-band) label for one found atom."""
    sp = {82: "Pb", 22: "Ti", 8: "O"}.get(int(rec["species"]), "other")
    mode = "guided" if rec["guided"] else "blind"
    lo, hi = cfg.bulk_z_A
    z = rec["z_A"]
    band = "bulk" if lo <= z <= hi else ("entrance" if z < lo else "exit")
    return f"{sp}|{mode}|{band}"


def strata_array(found, cfg):
    return np.array([stratum_of(r, cfg) for r in found])


# ---------------------------------------------------------------- calibration
def _conformal_q(scores, alpha, min_n=20):
    """Split-conformal quantile with the finite-sample correction ceil((n+1)(1-alpha))/n.

    Returns (q, n). If the stratum is too small to support the requested level, returns
    NaN so the caller can fall back to a pooled q rather than fabricate a tight interval
    from three points."""
    s = np.asarray([v for v in scores if np.isfinite(v)])
    n = s.size
    if n < min_n:
        return np.nan, n
    k = int(np.ceil((n + 1) * (1 - alpha)))
    if k > n:                      # n too small for this alpha even if >= min_n
        return np.nan, n
    return float(np.sort(s)[k - 1]), n


def calibrate(found, match, cfg, alphas=(0.32, 0.05), min_n=20):
    """Fit the Mondrian conformal q-table from matched atoms.

    match: the dict returned by validate.finder_report / match_found_to_gt, containing
    match_gi and the per-axis signed errors. Only matched atoms carry a residual, so only
    they can calibrate -- that is a property of any GT-based calibration, and unmatched
    (spurious) detections are excluded from coverage claims by construction.

    Returns {"q": {alpha: {stratum: {axis: q}}}, "pooled": ..., "n": {...}}."""
    matched = match["match_gi"] >= 0
    strat = strata_array(found, cfg)
    err = {"x": match["match_dx"], "y": match["match_dy"], "z": match["match_dz"]}
    sig = {"x": found["sx_A"], "y": found["sy_A"], "z": found["sz_A"]}
    out = {"q": {}, "pooled": {}, "n": {}, "alphas": list(alphas)}
    for a in alphas:
        ak = f"{a:g}"
        out["q"][ak] = {}
        out["pooled"][ak] = {}
        for ax in AXES:
            ok = matched & np.isfinite(err[ax]) & (sig[ax] > 0)
            sc_all = np.abs(err[ax][ok]) / sig[ax][ok]
            qp, _ = _conformal_q(sc_all, a, min_n=min_n)
            out["pooled"][ak][ax] = qp
            for st in np.unique(strat):
                m = ok & (strat == st)
                q, n = _conformal_q(np.abs(err[ax][m]) / sig[ax][m], a, min_n=min_n)
                out["q"][ak].setdefault(st, {})[ax] = q
                out["n"][st] = int(m.sum())
    return out


def apply(found, qtab, cfg, alpha):
    """Half-widths for every atom at the requested level. Falls back to the pooled q where
    a stratum is under-populated (and never silently to 1.0)."""
    ak = f"{alpha:g}"
    if ak not in qtab["q"]:
        raise KeyError(f"no calibration at alpha={alpha}; have {list(qtab['q'])}")
    strat = strata_array(found, cfg)
    sig = {"x": found["sx_A"], "y": found["sy_A"], "z": found["sz_A"]}
    out = {}
    for ax in AXES:
        q = np.array([qtab["q"][ak].get(s, {}).get(ax, np.nan) for s in strat], float)
        pooled = qtab["pooled"][ak][ax]
        q = np.where(np.isfinite(q), q, pooled)
        out[ax] = sig[ax] * q
    return out


# ---------------------------------------------------------------- reporting
def coverage_table(found, match, cfg, qtab, alpha):
    """Empirical coverage per stratum at the given level -- the audit of the calibration."""
    matched = match["match_gi"] >= 0
    strat = strata_array(found, cfg)
    hw = apply(found, qtab, cfg, alpha)
    err = {"x": match["match_dx"], "y": match["match_dy"], "z": match["match_dz"]}
    rows = []
    for st in sorted(np.unique(strat)):
        m = matched & (strat == st)
        if m.sum() < 5:
            continue
        cov = {ax: float(np.mean(np.abs(err[ax][m]) <= hw[ax][m])) for ax in AXES}
        rows.append((st, int(m.sum()), cov))
    allm = matched
    cov_all = {ax: float(np.mean(np.abs(err[ax][allm]) <= hw[ax][allm])) for ax in AXES}
    return rows, cov_all


def uncertainty_report(found, match, cfg, qtab, alphas=(0.32, 0.05)):
    """Human-readable UQ block for the run verdict."""
    lines = ["UNCERTAINTY (model sigma = joint CRLB (+) kernel mismatch; "
             "intervals = split-conformal, Mondrian strata)"]
    for a in alphas:
        rows, cov_all = coverage_table(found, match, cfg, qtab, a)
        tgt = 1 - a
        lines.append(f"  target {tgt:.0%}:  overall coverage "
                     f"x {cov_all['x']:.0%}  y {cov_all['y']:.0%}  z {cov_all['z']:.0%}")
        worst = sorted(rows, key=lambda r: min(r[2].values()))[:3]
        for st, n, cov in worst:
            lines.append(f"      worst stratum {st:<22s} n={n:<5d} "
                         f"x {cov['x']:.0%}  y {cov['y']:.0%}  z {cov['z']:.0%}")
    ak = f"{alphas[0]:g}"
    qz = {s: v.get("z") for s, v in qtab["q"][ak].items() if np.isfinite(v.get("z", np.nan))}
    if qz:
        lo = min(qz, key=qz.get); hi = max(qz, key=qz.get)
        lines.append(f"  q(z) spread across strata: {qz[lo]:.2f} ({lo}) .. {qz[hi]:.2f} ({hi})"
                     f"   <- a single floor cannot represent this")
    return lines


def species_probability(found):
    """Posterior probability of the ASSIGNED species from the fitted amplitude.

    The species call is amplitude-driven. The discriminating uncertainty is NOT the per-atom
    statistical CRLB on beta (~0.016; on noiseless data it saturates the posterior to a hard
    0/1 and is meaningless as a probability) -- it is the WITHIN-CLASS AMPLITUDE SPREAD: how
    much genuine Ti (or O) amplitudes vary within a column, which is what actually makes the
    two classes confusable. We use the robust per-class spread on the column (MAD, floored by
    the CRLB beta so a single-member class still has a scale), equal priors, Gaussian
    likelihood against the class centres present in that column.

    This gives a graded, honest posterior: a faint Ti sitting at the O amplitude centre
    correctly gets LOW p for its Ti-or-O call. (It cannot catch the residual confusion where
    an atom's amplitude is genuinely anomalous -- a real Ti reconstructed at an O amplitude --
    because by amplitude that atom truly looks like O; that is a reconstruction-quality error,
    not an uncertainty the amplitude carries.)

    On single-species columns (A-site Pb, pure-O) there is no within-column alternative, so
    p = 1 by construction -- a statement about the within-column decision only; it does NOT
    cover column-type misassignment, the dominant error mode there."""
    p = np.ones(len(found))
    amp = found["amp"].astype(float)
    samp = (found["samp"].astype(float) if "samp" in found.dtype.names
            else np.full(len(found), 0.0))
    for cid in np.unique(found["col_id"]):
        m = np.where(found["col_id"] == cid)[0]
        sp = found["species"][m]
        classes = list(np.unique(sp))
        if len(classes) < 2:
            continue
        centres, spreads = {}, {}
        for c in classes:
            a_c = amp[m[sp == c]]
            centres[c] = float(np.median(a_c))
            mad = 1.4826 * np.median(np.abs(a_c - centres[c])) if a_c.size > 2 else 0.0
            spreads[c] = mad
        # common discrimination scale: pooled within-class spread, floored by the CRLB
        base = float(np.median([s for s in spreads.values() if s > 0]) or 0.0)
        for i in m:
            sp_i = int(found["species"][i])
            s = max(base, samp[i], 0.05*abs(centres.get(sp_i, amp[i])) or 1e-6)
            ll = np.array([-0.5*((amp[i]-centres[c])/s)**2 for c in classes])
            ll -= ll.max()
            w = np.exp(ll); w /= w.sum()
            p[i] = float(w[classes.index(sp_i)])
    return p


def save(qtab, path):
    with open(path, "w") as f:
        json.dump(qtab, f, indent=2, default=lambda o: None if o != o else float(o))


def load(path):
    with open(path) as f:
        return json.load(f)
