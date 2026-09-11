#!/usr/bin/env python
"""@file test_fusion.py
@brief Gate the whole fusion proof of concept without touching a GPU.

Runs in a few seconds on the laptop or on Blythe, and covers every step that a reconstruction
cannot check for itself: the specimen really carries the designed polarisation, the two channels
really are degenerate in the ways the argument depends on, the hollow mask really lands on the
pixel PtychoShelves calls the centre, and the sign readout really recovers all four domains.

It also guards the working constraint: `sim/simulate_4dstem.py` and `ptycho/run_synthetic_recon_ML.m`
carry the validated PTO/STO labyrinth geometry and must stay byte-identical to HEAD.

    ~/hyperspy-bundle/bin/python test_fusion.py
"""
from __future__ import annotations

import os
import subprocess
import sys
import traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "eels"))

import analyze_elnes as AE      # noqa: E402
import analyze_fusion as AF     # noqa: E402
import eels_forward as F        # noqa: E402
import simulate_fusion as SF    # noqa: E402
import toy_sample as T          # noqa: E402

RESULTS = []


def test(fn):
    RESULTS.append(fn)
    return fn


# ---------------------------------------------------------------- the specimen
@test
def test_polar_mode_matches_pbtio3():
    """The rotated distortion reproduces the literature off-centrings and shifts nothing rigidly."""
    w, d_ti = T.polar_mode()
    assert abs(abs(d_ti) - 0.331) < 0.002, f"|delta_Ti| = {abs(d_ti):.4f}, expected 0.331"
    assert abs(abs(w[0] - w[2]) - 0.487) < 0.002, f"|delta_Pb| = {abs(w[0]-w[2]):.4f}, expected 0.487"
    net = w[0] + w[1] + 3 * w[2]
    assert abs(net) < 1e-12, f"net cell translation {net:.2e} -- domains would step at the wall"


@test
def test_build_gates():
    """Every structural gate in toy_sample.validate passes."""
    atoms, truth = T.build(12, 5)
    assert T.validate(atoms, truth, verbose=False), "toy_sample.validate failed (run it for detail)"


@test
def test_surfaces_symmetric():
    """Both surfaces are PbO: no composition up/down asymmetry to fake a polarisation signal."""
    atoms, truth = T.build(8, 4)
    Z, z = atoms.get_atomic_numbers(), atoms.get_positions()[:, 2]
    bot, top = set(Z[z < z.min() + 1.0]), set(Z[z > z.max() - 1.0])
    assert bot == top == {8, 82}, f"terminations differ: bottom {sorted(bot)}, top {sorted(top)}"


@test
def test_columns_homogeneous_along_beam():
    """Every bulk column carries ONE polarisation over its full depth -- no deconvolution needed."""
    import build_cells as B
    atoms, truth = T.build(8, 4)
    d = B.offcentering(atoms, 22, 6)
    dom, ix, iy, bulk = T.domain_index(atoms, truth)
    key = ix * int(truth["n_lat"]) + iy
    for k in np.unique(key[bulk]):
        m = (key == k) & bulk
        assert np.abs(d[m] - d[m].mean(0)).max() < 1e-9, "a bulk column is not homogeneous"


# ---------------------------------------------------------------- the EELS channel
@test
def test_q_weights_sum_to_one():
    """<cos^2(q,x)> + <cos^2(q,y)> + <cos^2(q,z)> = 1 -- the aperture average is a proper projection."""
    tE = AE.characteristic_angle_rad(532.0, 300.0)
    tot = sum(F.q_cos2(0.1, 0.075, tE, ax, n=200_000) for ax in "xyz")
    assert abs(tot - 1.0) < 5e-3, f"weights sum to {tot:.4f}"


@test
def test_parallel_limit_reproduces_analytic_model():
    """As alpha -> 0 the convergent aperture average must reduce to analyze_elnes.parallel_weight."""
    tE = AE.characteristic_angle_rad(532.0, 300.0)
    for beta in (2e-3, 5e-3):
        mc = F.q_cos2(1e-9, beta, tE, "z", n=400_000)
        an = AE.parallel_weight(beta, tE)
        assert abs(mc - an) < 5e-3, f"beta {beta*1e3:.0f} mrad: MC {mc:.4f} vs analytic {an:.4f}"


@test
def test_eels_is_blind_to_the_sign():
    """The central claim: A and B (opposite delta_z, all else equal) give IDENTICAL spectra."""
    _, truth = T.build(8, 4)
    lad = F.load_ladder()
    w = F.aperture_weights(100.0, 75.0)
    dmax = float(truth["delta_Ti_A"])
    sA = F.column_spectrum(truth["deltas"][0], dmax, w, lad)
    sB = F.column_spectrum(truth["deltas"][1], dmax, w, lad)
    assert F.contrast(sA, sB, lad[1]) < 1e-12, "A and B are not degenerate -- the premise is wrong"
    sC = F.column_spectrum(truth["deltas"][2], dmax, w, lad)
    sD = F.column_spectrum(truth["deltas"][3], dmax, w, lad)
    assert F.contrast(sC, sD, lad[1]) < 1e-12, "C and D are not degenerate under a round aperture"


@test
def test_eels_magnitude_channel_is_measurable():
    """AB vs CD (different |delta_z|) must be separable at a dose that is actually achievable."""
    _, truth = T.build(8, 4)
    lad = F.load_ladder()
    w = F.aperture_weights(100.0, 75.0)
    dmax = float(truth["delta_Ti_A"])
    c = F.contrast(F.column_spectrum(truth["deltas"][0], dmax, w, lad),
                   F.column_spectrum(truth["deltas"][2], dmax, w, lad), lad[1])
    assert c > 0.03, f"contrast {c*100:.2f}% is too small to be the magnitude channel"
    assert AE.required_counts(c) < 1e5, f"needs {AE.required_counts(c):.1e} counts/channel"


@test
def test_hole_size_does_not_cost_the_spectroscopy():
    """Beyond alpha ~ 5 mrad the convergence sets the q spread, so the contrast is flat in beta."""
    _, truth = T.build(8, 4)
    lad = F.load_ladder()
    dmax = float(truth["delta_Ti_A"])
    cs = []
    for f in (0.1, 0.5, 0.75, 0.95):
        w = F.aperture_weights(100.0, 100.0 * f)
        cs.append(F.contrast(F.column_spectrum(truth["deltas"][0], dmax, w, lad),
                             F.column_spectrum(truth["deltas"][2], dmax, w, lad), lad[1]))
    spread = max(cs) / min(cs) - 1.0
    assert spread < 0.20, f"contrast varies by {spread*100:.0f}% across hole size -- not flat"


@test
def test_eels_inverts_to_the_right_angle():
    """The magnitude channel inverts single-valuedly for theta, hence for |delta_z|."""
    _, truth = T.build(8, 4)
    th = AF.eels_theta_from_contrast(truth, 100.0, 75.0)
    for k, n in enumerate(truth["names"]):
        assert abs(th[str(n)] - float(truth["theta_deg"][k])) < 1.0, \
            f"domain {n}: recovered {th[str(n)]:.1f} deg, truth {truth['theta_deg'][k]:.1f}"


# ---------------------------------------------------------------- the hollow detector
@test
def test_hollow_mask_geometry():
    """The hole is centred on the pixel PtychoShelves is told about and has the right radius."""
    n_b, d_alpha = 238, 1.6810
    theta = SF.detector_axes(n_b, d_alpha)
    m = SF.hollow_mask(theta, 75.0)
    ys, xs = np.nonzero(m == 0)
    assert abs(xs.mean() - n_b / 2) < 1e-9 and abs(ys.mean() - n_b / 2) < 1e-9, \
        "hole is not centred on index N/2 (MATLAB p.ctr = fix(N/2)+1)"
    r = np.hypot(xs - n_b / 2, ys - n_b / 2).max() * d_alpha
    assert 74.0 < r <= 75.0, f"hole radius {r:.2f} mrad, expected 75"
    assert set(np.unique(m)) <= {0.0, 1.0}, "mask must be binary (1 = use, 0 = ignore)"


@test
def test_dose_budget_is_monotonic_and_physical():
    """More hole -> more EELS, less ptychography, and the two always sum to the whole beam."""
    n_b, d_alpha = 120, 3.0
    theta = SF.detector_axes(n_b, d_alpha)
    dp = np.exp(-(theta / 60.0) ** 2)[None]              # a plausible forward-peaked pattern
    b = SF.dose_budget(dp, theta, [0.0, 50.0, 75.0, 95.0], 100.0, 180.0)
    fr = [b["hsa"][f"{h:.1f}"]["eels_fraction"] for h in (0.0, 50.0, 75.0, 95.0)]
    assert fr == sorted(fr), f"EELS fraction not monotonic in the hole size: {fr}"
    for h in (0.0, 50.0, 75.0, 95.0):
        e = b["hsa"][f"{h:.1f}"]
        assert abs(e["eels_fraction"] + e["ptycho_fraction"] - 1.0) < 1e-12


# ---------------------------------------------------------------- the ptychography channel
@test
def test_hypotheses_are_distinguishable():
    """The up and down structures must differ by more than a rounding error in profile space."""
    _, truth = T.build(8, 4)
    z = np.arange(0.0, 6 * float(truth["c"]), 0.1)
    for k in (0, 2):
        d = truth["deltas"][k]
        up, dn = AF.hypothesis_pair(d[:2], abs(d[2]), truth, z, 0.0)
        assert np.linalg.norm(up - dn) > 0.3, "the two hypotheses are nearly identical"


@test
def test_sign_readout_recovers_every_domain():
    """End to end on a synthetic reconstruction: all four signs, from depth sectioning alone."""
    assert AF.selftest(verbose=False), "the sign readout failed on the synthetic reconstruction"


@test
def test_sign_readout_survives_noise():
    """The per-domain decision must hold at 100% profile noise, which is far past a real recon."""
    _, truth = T.build(12, 5)
    budget = {"scan_window_A": 24.0, "scan_center_A": [float(truth["box_A"]) / 2] * 2,
              "beam_thickness_A": float(truth["box_z_A"]), "dx_object_A": 0.1}
    cell, dom = AF._score(truth, budget, 2.0, 1.0, 0)
    assert dom == 1.0, f"only {dom*100:.0f}% of domains correct at 100% noise"


@test
def test_fusion_is_exact_on_ground_truth():
    """Given the true sign and the true theta, the fused vector must be the true vector."""
    _, truth = T.build(8, 4)
    th = AF.eels_theta_from_contrast(truth, 100.0, 75.0)
    obs = np.array([np.sign(truth["deltas"][k][2]) for k in range(4)])
    dom = np.array([str(n) for n in truth["names"]])
    res = AF.fuse(obs, dom, truth, th,
                  {str(n): truth["deltas"][k][:2] for k, n in enumerate(truth["names"])})
    for n, r in res.items():
        assert np.allclose(r["delta"], r["truth"], atol=2e-3), \
            f"domain {n}: fused {np.round(r['delta'],4)} vs truth {np.round(r['truth'],4)}"


@test
def test_known_object_frame_round_trip():
    """known_object.register must recover the exact frame cut() produced, in both orientations.

    If the two disagree on an offset, the true object is written misaligned even though the
    registration reports success -- and the known-object diagnosis would be void.
    """
    import known_object as K
    dx = 0.049211
    Tr, box, _ = K.true_transmission(os.path.join(HERE, "sample", "toy_membrane.vasp"), 4, dx)
    n_full = 864; roi = int(round(0.55 * n_full)); o = (n_full - roi) // 2
    anchor = int(round((box / 2) / dx)) - n_full // 2 + o
    for orient in ("identity", "transpose"):
        s = anchor + 7
        fake = K.cut(Tr, orient, s, s + 3, n_full, roi)
        got = K.register(Tr, fake, roi, dx, anchor, T.REF.a_tet)
        assert got[:3] == (orient, s % Tr.shape[1], (s + 3) % Tr.shape[1]), f"{orient}: {got}"


# ---------------------------------------------------------------- the working constraint
@test
def test_labyrinth_code_untouched():
    """sim/simulate_4dstem.py and the PTO/STO MATLAB driver carry the validated labyrinth geometry.

    This work imports and parameterises them; it must never edit them. Skipped outside a git repo.
    """
    repo = os.path.abspath(os.path.join(HERE, ".."))
    guarded = ["sim/simulate_4dstem.py", "ptycho/run_synthetic_recon_ML.m"]
    try:
        r = subprocess.run(["git", "-C", repo, "diff", "--name-only", "HEAD", "--", *guarded],
                           capture_output=True, text=True, timeout=30)
    except Exception:
        print("    (skipped: no git)")
        return
    if r.returncode != 0:
        print("    (skipped: not a git repo)")
        return
    changed = [x for x in r.stdout.split() if x]
    assert not changed, f"MODIFIED protected files: {changed} -- the labyrinth pipeline must stay intact"


@test
def test_hollow_driver_exists_and_masks():
    """The MHP driver is present and really wires the hole into the engine's mask1."""
    p = os.path.join(HERE, "..", "ptycho", "run_fusion_hollow.m")
    assert os.path.exists(p), "ptycho/run_fusion_hollow.m is missing"
    src = open(p).read()
    for token in ("global mask1", "mask_hsa%.2f.mat", "getenv('HSA')", "meta.convergence_mrad"):
        assert token in src, f"run_fusion_hollow.m does not contain {token!r}"


def main() -> int:
    print(f"running {len(RESULTS)} fusion tests\n")
    failed = 0
    for fn in RESULTS:
        try:
            fn()
            print(f"  [ok  ] {fn.__name__}")
        except Exception as exc:
            failed += 1
            print(f"  [FAIL] {fn.__name__}: {exc}")
            if os.environ.get("VERBOSE"):
                traceback.print_exc()
    print(f"\n{len(RESULTS)-failed}/{len(RESULTS)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
