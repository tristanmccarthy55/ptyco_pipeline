#!/usr/bin/env python
"""@file make_gt_cache.py
@brief Precompute the prepared ground-truth frame so a reproduction run does not need abtem.

`align.load_gt` prepares the reference structure exactly as the simulation did (rotate to
beam=+z, orthogonalize, pad square, centre, vacuum). Only the orthogonalize step needs abtem,
which is a heavy dependency to install for one call, so the result is cached as a small .npz
that ships in the peer-reproduction tarball.

    python -m atomfind.make_gt_cache            # writes <package>/data/gt_prepared.npz
    python -m atomfind.make_gt_cache --check    # rebuild and compare against the cache
"""
from __future__ import annotations
import argparse
import os

import numpy as np

from . import align, config


def build():
    """(pos, Z) from the .vasp via the abtem path -- the definition the cache must match."""
    return align._prepare_gt(config.data_path(config.VASP_NAME, required=True))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                  "data", align.GT_CACHE))
    ap.add_argument("--check", action="store_true",
                    help="rebuild from the .vasp and verify the existing cache matches")
    a = ap.parse_args()

    pos, Z = build()
    if a.check:
        d = np.load(config.data_path(align.GT_CACHE, required=True))
        assert d["pos"].shape == pos.shape, f"shape {d['pos'].shape} != {pos.shape}"
        worst = float(np.abs(d["pos"] - pos).max())
        assert (d["Z"] == Z).all(), "atomic numbers differ"
        assert worst == 0.0, f"positions differ by up to {worst:.3e} A"
        print(f"cache matches the abtem path exactly ({len(Z)} atoms, max |dpos| = {worst:g})")
        return
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    np.savez_compressed(a.out, pos=pos, Z=Z)
    print(f"wrote {a.out}  ({len(Z)} atoms, {os.path.getsize(a.out)/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
