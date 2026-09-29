#!/usr/bin/env python
"""@file test_haider.py
@brief The CEOS (Haider) -> abTEM conversion of the ARM200F tableau, checked three ways.

aberration_waves.haider_to_abtem() turns a CEOS tableau entry (magnitude, angle of the complex coefficient) into abTEM's
(C_nm, phi_nm). It is checked against the Haider formula itself (aberration_waves.haider_chi, evaluated straight from
the complex notation), then abTEM's own aberration function and the probe abTEM builds are checked against the same
formula, so the tableau the simulation receives is the tableau the corrector measured. A wrong factor or angle sign
shows up as hundreds of waves at 80 mrad (the naive reading, C_nm = CEOS magnitude and phi = +angle/m, is off by ~400).

    ~/hyperspy-bundle/bin/python campaign/test_haider.py        # or: pytest campaign/test_haider.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aberration_waves as aw          # noqa: E402

LAM = aw.LAMBDA_A
C3, C5 = -2.0e5, 4.0e7                 # a C3/C5 pair of the ARM plans' size


def _tableau():
    tab = {k: (v[0], v[1]) for k, v in aw.arm_measured().items()}
    tab["C3"], tab["C5"] = (C3, None), (C5, None)
    return tab


def _krivanek(d, th, ph):
    w = np.zeros_like(th)
    for k, v in d.items():
        if k.startswith("C"):
            n, m = int(k[1]), int(k[2])
            w += v * th ** (n + 1) / (n + 1) * np.cos(m * (ph - d.get("phi" + k[1:], 0.0)))
    return w


def test_each_term_matches_haider():
    t = np.linspace(-0.08, 0.08, 161)
    tx, ty = np.meshgrid(t, t, indexing="xy")
    th, ph = np.hypot(tx, ty), np.arctan2(ty, tx)
    for name, val in _tableau().items():
        d = aw.haider_to_abtem(name, *val)
        err = np.abs(_krivanek(d, th, ph) - aw.haider_chi({name: val}, tx, ty)).max() / LAM
        assert err < 1e-9, (name, err)


def test_mixed_terms_carry_their_factor():
    assert aw.haider_to_abtem("B2", 1.0, 0.0)["C21"] == 3.0
    assert aw.haider_to_abtem("S3", 1.0, 0.0)["C32"] == 4.0
    assert aw.haider_to_abtem("B4", 1.0, 0.0)["C41"] == 5.0
    assert aw.haider_to_abtem("D4", 1.0, 0.0)["C43"] == 5.0
    assert aw.haider_to_abtem("A5", 1.0, 0.0)["C56"] == 1.0


def test_abtem_aberration_function_matches_haider():
    import abtem
    t = np.linspace(-0.08, 0.08, 321)
    tx, ty = np.meshgrid(t, t, indexing="xy")
    alpha, phi = np.hypot(tx, ty), np.arctan2(ty, tx)
    A = abtem.transfer.Aberrations(aberration_coefficients=aw.arm_tableau(C3, C5), energy=300e3)
    E = np.asarray(A._evaluate_from_angular_grid(alpha, phi))            # exp(-i chi)
    ref = -2 * np.pi / LAM * aw.haider_chi(_tableau(), tx, ty)
    d = np.angle(E * np.exp(-1j * ref))[alpha <= 0.08]
    assert np.abs(d).max() < 0.02, np.abs(d).max()                       # float32 + rounding, of ~200 waves


def test_abtem_probe_carries_the_tableau():
    import abtem
    try:
        abtem.config.set({"local_diagnostics.progress_bar": False})
    except Exception:
        pass
    al, c1, L, N = 60, 50.0, 200.0, 2048
    P = np.asarray(abtem.Probe(energy=300e3, semiangle_cutoff=al, extent=L, gpts=N, defocus=c1,
                               aberrations=aw.arm_tableau(C3, C5)).build(lazy=False).array)
    F = np.fft.fft2(np.fft.ifftshift(P))            # the probe sits at the array centre: undo that before transforming
    k = np.fft.fftfreq(N, L / N)
    KX, KY = np.meshgrid(k, k, indexing="ij")       # abTEM arrays are [x, y]
    TX, TY = KX * LAM, KY * LAM
    th = np.hypot(TX, TY)
    m = th < 0.9 * al * 1e-3
    ref = -2 * np.pi / LAM * (aw.haider_chi(_tableau(), TX, TY) - c1 * th ** 2 / 2)
    d = np.angle(F * np.exp(-1j * ref))[m]
    d = np.angle(np.exp(1j * (d - np.angle(np.mean(np.exp(1j * d))))))
    assert np.abs(d).max() < 0.05, np.abs(d).max()


if __name__ == "__main__":
    for f in (test_each_term_matches_haider, test_mixed_terms_carry_their_factor,
              test_abtem_aberration_function_matches_haider, test_abtem_probe_carries_the_tableau):
        f()
        print("PASS", f.__name__)
