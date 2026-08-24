#!/usr/bin/env python
"""@file test_atomfind.py
@brief Fast, data-free test suite: validates an install without the 91 MB volume.

Every test here runs on SYNTHETIC data in a few seconds, so a peer reproducing the published
result can confirm their environment is sane before committing to the full pipeline, and can
tell an environment problem apart from a genuine disagreement with our numbers.

The end-to-end check on the real volume is a separate, slower thing: run
`python atomfind/run_atomfind.py --preset NL70_coherent` and compare against PEER.md.

    python atomfind/test_atomfind.py            # all tests
    python atomfind/test_atomfind.py -v         # with tracebacks
"""
from __future__ import annotations
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from atomfind import config, align, psf as psfmod, deconv, find, validate, uncertainty

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- helpers
def _cfg(**kw):
    """A small synthetic-volume config: 1 A/layer, 0.05 A/px, tiny field."""
    c = config.preset("NL70_coherent")
    c.dz = 1.0
    c.trim_z_A = (0.0, 40.0)
    c.zmax_show_A = 40.0
    c.bulk_z_A = (2.0, 38.0)
    c.exit_band_z_A = 38.0
    for k, v in kw.items():
        setattr(c, k, v)
    return c


def _gauss_kernel(nz=9, nxy=9, sz=1.6, sxy=1.0):
    """A separable Gaussian standing in for the measured single-atom response."""
    z = np.arange(nz) - nz // 2
    y = np.arange(nxy) - nxy // 2
    k = (np.exp(-0.5 * (z / sz) ** 2)[:, None, None]
         * np.exp(-0.5 * (y / sxy) ** 2)[None, :, None]
         * np.exp(-0.5 * (y / sxy) ** 2)[None, None, :])
    return k / k.max()


def _planted_volume(sites, shape=(40, 60, 60), amp=1.0, noise=0.0, seed=0):
    """Volume built by stamping a known kernel at known (layer,row,col) sites."""
    V = np.zeros(shape)
    K = _gauss_kernel()
    hz, hxy = K.shape[0] // 2, K.shape[1] // 2
    for (l, r, c) in sites:
        zs, ze = l - hz, l + hz + 1
        ys, ye = r - hxy, r + hxy + 1
        xs, xe = c - hxy, c + hxy + 1
        if zs < 0 or ys < 0 or xs < 0 or ze > shape[0] or ye > shape[1] or xe > shape[2]:
            continue
        V[zs:ze, ys:ye, xs:xe] += amp * K
    if noise:
        V = V + np.random.default_rng(seed).normal(0, noise, shape)
    return V


# ---------------------------------------------------------------- configuration / portability
def test_data_path_search_order():
    """data_path resolves by BASENAME through $ATOMFIND_DATA -> package data/ -> repo -> Desktop,
    and an absolute path passes through untouched. This is what makes the package portable."""
    with tempfile.TemporaryDirectory() as d:
        probe = os.path.join(d, "zz_probe_file.npy")
        np.save(probe, np.zeros(3))
        old = os.environ.get("ATOMFIND_DATA")
        try:
            os.environ["ATOMFIND_DATA"] = d
            assert os.path.samefile(config.data_path("zz_probe_file.npy"), probe)
            assert config.data_path(probe) == probe, "absolute paths must pass through"
        finally:
            os.environ.pop("ATOMFIND_DATA", None)
            if old:
                os.environ["ATOMFIND_DATA"] = old


def test_missing_data_raises_by_name():
    """A missing file must name itself and the search path, not raise a bare ENOENT."""
    try:
        config.data_path("zz_definitely_absent.npy", required=True)
    except FileNotFoundError as e:
        assert "zz_definitely_absent.npy" in str(e)
        assert "ATOMFIND_DATA" in str(e), "the error should say how to fix it"
    else:
        raise AssertionError("required=True must raise for a missing file")


def test_gt_cache_present_and_sane():
    """The shipped prepared-GT cache is what lets a peer run without abtem."""
    p = config.data_path(align.GT_CACHE)
    assert os.path.exists(p), f"{align.GT_CACHE} missing -- python -m atomfind.make_gt_cache"
    d = np.load(p)
    pos, Z = d["pos"], d["Z"]
    assert pos.shape[0] == Z.shape[0] and pos.shape[1] == 3
    assert len(Z) > 10_000, f"only {len(Z)} atoms in the cache"
    assert set(np.unique(Z)) <= {8, 22, 38, 82}, f"unexpected species {set(np.unique(Z))}"
    assert pos.min() >= -1e-6, "prepared cell should sit in the positive octant"


def test_load_object_accepts_npy_unchanged():
    """The .npy path must be exactly np.load -- the dispatch wrapper must not perturb the
    published result, which was produced through it."""
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "v.npy")
        a = (np.random.default_rng(0).normal(size=(4, 5, 5))
             + 1j * np.random.default_rng(1).normal(size=(4, 5, 5))).astype(np.complex64)
        np.save(p, a)
        assert np.array_equal(align.load_object(p), np.load(p))


def test_load_object_rejects_unknown_format_clearly():
    """A peer pointing this at the wrong file should be told what it wants, by name."""
    try:
        align.load_object("something.h5")
    except ValueError as e:
        assert ".mat" in str(e) and ".npy" in str(e), f"unhelpful message: {e}"
    else:
        raise AssertionError("an unknown extension must raise")


def test_phase_is_per_layer_median_subtracted():
    """Each layer carries its own arbitrary phase offset from the reconstruction; leaving it
    in would put a layer-dependent pedestal under every amplitude."""
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "v.npy")
        nz = 6
        phase = np.zeros((nz, 8, 8))
        for l in range(nz):                       # a different offset per layer
            phase[l] = 0.3 * l
        phase[:, 4, 4] += 0.5                     # one bright voxel per layer
        np.save(p, np.exp(1j * phase).astype(np.complex64))
        cfg = _cfg(recon_vol=p, dx=0.05)
        V, _ = align.load_phase(cfg)
        med = np.median(V, axis=(1, 2))
        assert np.abs(med).max() < 1e-9, f"per-layer offset survived (median {med})"


# ---------------------------------------------------------------- geometry
def test_alignment_affine_roundtrip():
    """site_to_index is the ONLY recon<->model map; a scale/offset error here silently
    corrupts every recall number, so pin it against an explicit inverse."""
    al = align.Alignment(dx=0.05, dz=1.0, X0=30.0, Y0=10.0, SGN=1, OFF=0.4,
                         CAL_X=0.0, CAL_Y=0.0, corr_depth=0.9,
                         mX=0.005, bX=0.47, mY=0.006, bY=0.47)
    X, Y, Zd = 33.0, 12.5, 20.0
    row, col, layer = al.site_to_index(X, Y, Zd)
    col_back = ((col - al.bX) / (1 + al.mX)) * al.dx + al.X0
    row_back = ((row - al.bY) / (1 + al.mY)) * al.dx + al.Y0
    z_back = ((layer / (1 + al.mZ) + 0.5) * al.dz - al.OFF) * al.SGN
    assert abs(col_back - X) < 1e-9 and abs(row_back - Y) < 1e-9, "in-plane map not invertible"
    assert abs(z_back - Zd) < 1e-9, "depth map not invertible"


def test_alignment_affine_is_not_constant():
    """The in-plane map must be AFFINE, not a constant offset: dx=window/N differs from the
    physical pixel by ~0.6%, which accumulates to ~2 px across the field (METHODS 11.4)."""
    al = align.Alignment(dx=0.05, dz=1.0, X0=0.0, Y0=0.0, SGN=1, OFF=0.0, CAL_X=0, CAL_Y=0,
                         corr_depth=0.9, mX=0.006, bX=0.0, mY=0.006, bY=0.0)
    _, c0, _ = al.site_to_index(0.0, 0.0, 0.0)
    _, c1, _ = al.site_to_index(20.0, 0.0, 0.0)
    drift = (c1 - c0) - (20.0 / al.dx)
    assert drift > 1.0, f"a 0.6% scale over a 400 px field must drift >1 px, got {drift:.2f}"


# ---------------------------------------------------------------- detection
def test_peaks3d_finds_planted_atoms():
    """The peak-picking baseline must recover well-separated planted atoms."""
    sites = [(10, 20, 20), (10, 20, 40), (20, 40, 20), (30, 30, 30)]
    V = _planted_volume(sites)
    got = find.peaks3d(V, _cfg(), 0.05, rel_floor=0.05)
    assert len(got) >= len(sites), f"found {len(got)} of {len(sites)} planted atoms"
    for (l, r, c) in sites:
        d = np.hypot(got["row"] - r, got["col"] - c) + np.abs(got["layer"] - l)
        assert d.min() < 2.0, f"planted atom at {(l,r,c)} not recovered (nearest {d.min():.1f})"


def test_peaks3d_respects_amplitude_floor():
    """A relative floor must actually suppress weak detections, or 'best-effort threshold'
    comparisons against the baselines are meaningless."""
    V = _planted_volume([(10, 20, 20)], amp=1.0)
    V += _planted_volume([(20, 40, 40)], amp=0.01)
    many = len(find.peaks3d(V, _cfg(), 0.05, rel_floor=0.001))
    few = len(find.peaks3d(V, _cfg(), 0.05, rel_floor=0.5))
    assert few < many, f"floor had no effect ({few} vs {many})"
    assert few >= 1, "the bright atom must survive any sane floor"


def test_deconv_is_non_negative_and_guarded():
    """Richardson-Lucy must stay non-negative and must not run away on noise."""
    V = _planted_volume([(10, 20, 20), (20, 30, 30)], noise=0.02, seed=1)
    K = _gauss_kernel()
    dec, info = deconv.richardson_lucy_3d(V, K, _cfg())
    assert dec.min() >= -1e-9, "RL output must be non-negative"
    assert np.isfinite(dec).all(), "RL produced non-finite values"
    assert info["iters_done"] >= 1
    assert dec.shape == V.shape


# ---------------------------------------------------------------- uncertainty (the deliverable)
def _fake_found_and_match(n=600, seed=0):
    """Synthetic found-atom table with KNOWN Gaussian errors, for calibration tests."""
    rng = np.random.default_rng(seed)
    found = np.zeros(n, dtype=[("row", "f8"), ("col", "f8"), ("layer", "f8"), ("z_A", "f8"),
                               ("amp", "f8"), ("sx_A", "f8"), ("sy_A", "f8"), ("sz_A", "f8"),
                               ("samp", "f8"), ("quality", "f8"), ("species", "i4"),
                               ("col_id", "i4"), ("guided", "i4")])
    found["species"] = rng.choice([82, 22, 8], n)
    found["z_A"] = rng.uniform(5, 35, n)
    found["amp"] = rng.uniform(0.5, 1.5, n)
    found["quality"] = 0.9
    found["col_id"] = rng.integers(0, 50, n)
    for a in "xyz":
        found[f"s{a}_A"] = 0.05
    true_sigma = 0.10                      # model sigma is deliberately WRONG by 2x
    match = {"match_gi": np.arange(n),
             "match_dx": rng.normal(0, true_sigma, n),
             "match_dy": rng.normal(0, true_sigma, n),
             "match_dz": rng.normal(0, true_sigma, n)}
    return found, match, true_sigma


def test_conformal_coverage_holds_at_nominal():
    """THE central claim of the uncertainty machinery: conformal intervals hit their nominal
    coverage even when the model sigma is badly wrong. Here sigma is understated 2x."""
    found, match, _ = _fake_found_and_match()
    cfg = _cfg()
    qtab = uncertainty.calibrate(found, match, cfg, alphas=(0.32, 0.05), min_n=20)
    for alpha, target in ((0.32, 0.68), (0.05, 0.95)):
        _rows, cov = uncertainty.coverage_table(found, match, cfg, qtab, alpha)
        for ax in "xyz":
            assert abs(cov[ax] - target) < 0.06, \
                f"alpha={alpha}: {ax} coverage {cov[ax]:.2f}, target {target}"


def test_conformal_scales_with_wrong_model_sigma():
    """The conformal quantile must ABSORB a mis-scaled model sigma: double sigma, and q must
    roughly halve, leaving the interval width unchanged. That is what makes it transfer."""
    cfg = _cfg()
    found, match, _ = _fake_found_and_match()
    q1 = uncertainty.calibrate(found, match, cfg, alphas=(0.05,), min_n=20)
    w1 = uncertainty.apply(found, q1, cfg, 0.05)["z"]
    found2 = found.copy()
    for a in "xyz":
        found2[f"s{a}_A"] *= 2.0
    q2 = uncertainty.calibrate(found2, match, cfg, alphas=(0.05,), min_n=20)
    w2 = uncertainty.apply(found2, q2, cfg, 0.05)["z"]
    rel = abs(np.median(w2) - np.median(w1)) / np.median(w1)
    assert rel < 0.05, f"interval width moved {rel:.1%} when only the model sigma was rescaled"


def test_intervals_widen_with_confidence():
    """A 95% interval must be wider than a 68% one, everywhere."""
    cfg = _cfg()
    found, match, _ = _fake_found_and_match()
    qtab = uncertainty.calibrate(found, match, cfg, alphas=(0.32, 0.05), min_n=20)
    w68 = uncertainty.apply(found, qtab, cfg, 0.32)["z"]
    w95 = uncertainty.apply(found, qtab, cfg, 0.05)["z"]
    assert (w95 >= w68 - 1e-12).all(), "95% interval narrower than 68% somewhere"


# ---------------------------------------------------------------- scoring
def test_finder_report_recall_precision_arithmetic():
    """Recall and precision must be computed against a known matching, not merely be plausible."""
    cfg = _cfg()
    al = align.Alignment(dx=1.0, dz=1.0, X0=0.0, Y0=0.0, SGN=1, OFF=0.0,
                         CAL_X=0, CAL_Y=0, corr_depth=0.9)
    pos = np.array([[38.0, 18.0, 10.0], [40.0, 18.0, 10.0], [42.0, 18.0, 10.0]])
    Z = np.array([82, 82, 82])
    n = len(pos)
    found = np.zeros(n, dtype=[("row", "f8"), ("col", "f8"), ("layer", "f8"), ("z_A", "f8"),
                               ("amp", "f8"), ("sx_A", "f8"), ("sy_A", "f8"), ("sz_A", "f8"),
                               ("samp", "f8"), ("quality", "f8"), ("species", "i4"),
                               ("col_id", "i4"), ("guided", "i4")])
    r, c, l = al.site_to_index(pos[:, 0], pos[:, 1], pos[:, 2])
    found["row"], found["col"], found["layer"] = r, c, l
    found["species"] = 82
    found["amp"] = 1.0
    rep, _m = validate.finder_report(found, pos, Z, al, cfg)
    assert rep["n_found"] == n
    assert rep["precision"] > 0.99, f"perfect detections scored precision {rep['precision']}"
    assert rep["Pb"]["recall"] > 0.99, f"perfect detections scored recall {rep['Pb']['recall']}"


def test_confusion_rate_counts_species_errors():
    """Confusion is the diagnostic the results section leans on; it must respond to a real
    mislabelling and not to a correct one."""
    cfg = _cfg()
    al = align.Alignment(dx=1.0, dz=1.0, X0=0.0, Y0=0.0, SGN=1, OFF=0.0,
                         CAL_X=0, CAL_Y=0, corr_depth=0.9)
    pos = np.array([[38.0, 18.0, 10.0], [42.0, 18.0, 10.0]])
    Z = np.array([22, 8])
    dt = [("row", "f8"), ("col", "f8"), ("layer", "f8"), ("z_A", "f8"), ("amp", "f8"),
          ("sx_A", "f8"), ("sy_A", "f8"), ("sz_A", "f8"), ("samp", "f8"), ("quality", "f8"),
          ("species", "i4"), ("col_id", "i4"), ("guided", "i4")]
    r, c, l = al.site_to_index(pos[:, 0], pos[:, 1], pos[:, 2])

    right = np.zeros(2, dtype=dt); right["row"], right["col"], right["layer"] = r, c, l
    right["species"] = [22, 8]; right["amp"] = 1.0
    swapped = right.copy(); swapped["species"] = [8, 22]

    rep_ok, _ = validate.finder_report(right, pos, Z, al, cfg)
    rep_bad, _ = validate.finder_report(swapped, pos, Z, al, cfg)
    assert validate.confusion_rate(rep_ok) < 1e-9, "correct labels scored non-zero confusion"
    assert validate.confusion_rate(rep_bad) > 0.9, "a full Ti/O swap must read as confusion"


def test_health_warnings_fire_on_bad_input():
    """The health warnings are the safety net on a new dataset; silence on a broken run is
    the failure mode that matters."""
    cf = {f"{a}->{b}": 0 for a in (82, 22, 8) for b in (82, 22, 8)}
    cf["82->82"] = 20            # only the bright Pb are placed correctly
    cf["22->8"], cf["8->22"] = 40, 40
    bad = {"precision": 0.95, "n_found": 100,
           "Pb": {"recall": 0.9, "recall_bulk": 0.9}, "Ti": {"recall": 0.0, "recall_bulk": 0.0},
           "O": {"recall": 0.0, "recall_bulk": 0.0}, "confusion": cf}
    w = validate.health_warnings(bad)
    assert w, "a run with zero Ti/O recall and 80% confusion produced no warning"
    joined = " ".join(w).lower()
    assert "confusion" in joined, "the confusion canary did not fire"
    assert any(s in joined for s in ("ti recall", "o recall")), "the recall canary did not fire"

    good = {"precision": 0.97, "n_found": 1834,
            "Pb": {"recall": 0.88}, "Ti": {"recall": 0.89}, "O": {"recall": 0.89},
            "confusion": {**{f"{a}->{b}": 0 for a in (82, 22, 8) for b in (82, 22, 8)},
                          "82->82": 600, "22->22": 600, "8->8": 600, "22->8": 7}}
    assert not validate.health_warnings(good), \
        f"a healthy run produced spurious warnings: {validate.health_warnings(good)}"


# ---------------------------------------------------------------- runner
def _run(verbose=False):
    tests = sorted(k for k, v in globals().items() if k.startswith("test_") and callable(v))
    npass = 0
    for name in tests:
        try:
            globals()[name]()
            print(f"  PASS  {name}")
            npass += 1
        except Exception as e:
            print(f"  FAIL  {name}: {type(e).__name__}: {e}")
            if verbose:
                import traceback
                traceback.print_exc()
    print(f"== {npass}/{len(tests)} passed ==")
    return npass == len(tests)


if __name__ == "__main__":
    raise SystemExit(0 if _run("-v" in sys.argv) else 1)
