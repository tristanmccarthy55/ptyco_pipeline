#!/usr/bin/env python
"""
@file test_phonon_sigma.py
@brief Check that frozen phonons move each atom by the per-axis rms the simulator means, on the installed abTEM.

abTEM's FrozenPhonons `sigmas` is the 3-D rms on 1.0.5 and the per-axis std on 1.0.9. simulate_4dstem.abtem_sigmas()
measures which and converts; this draws 2000 configurations of one atom of each species through it, at room temperature
and at the cooled temperatures the campaign uses, and compares the per-axis std with u = sqrt(B / 8 pi^2) from the
simulator's own B table and Debye model. Exits non-zero if any species is off by more than 4 %.

    "$SHARE/phucrh/envs/abtem/bin/python" sim/test_phonon_sigma.py     # on Blythe, before any phonon sim
"""
from __future__ import annotations

import os
import sys

import numpy as np
from ase import Atoms

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import simulate_4dstem as sim  # noqa: E402

N, TOL = 2000, 0.04


def per_axis(symbol, u):
    a = Atoms(symbol, positions=[[5.0, 5.0, 5.0]], cell=[10.0, 10.0, 10.0])
    d = np.array([c.positions[0] - 5.0 for c in sim.abtem.FrozenPhonons(a, num_configs=N,
                                                                        sigmas=sim.abtem_sigmas({symbol: u}), seed=3)])
    return d.std(0)


def main():
    print(f"abTEM {sim.abtem.__version__}: sigmas = 1 -> per-axis {sim.abtem_per_axis():.4f}")
    bad = 0
    for T in (None, 100.0, 20.0):
        B = dict(sim._B_BY_SPECIES) if T is None else sim.B_at_temperature(T)
        for k, b in B.items():
            u = float(np.sqrt(b / (8 * np.pi ** 2)))
            s = per_axis(k, u)
            ok = np.all(np.abs(s / u - 1) < TOL)
            bad += not ok
            print(f"  {'RT' if T is None else f'{T:g} K':>6s} {k:2s}  u {u:.4f} A  got {np.round(s, 4)}  {'ok' if ok else 'OFF'}")
    print("PASS" if not bad else f"FAIL: {bad} species off by more than {TOL:.0%}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
