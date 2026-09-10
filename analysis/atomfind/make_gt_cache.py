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


def build(thin_cells: int = 0, z_vacuum: float = 4.0):
    """(pos, Z) from the .vasp via the abtem path -- the definition the cache must match.

    thin_cells > 0 builds the THIN aberration-campaign slab instead of the full 18-cell box
    (see align._prepare_gt_thin); the two are different structures, so a thin cache must be
    kept out of the package data dir and pointed at with --data-dir.
    """
    vasp = config.data_path(config.VASP_NAME, required=True)
    if thin_cells > 0:
        return align._prepare_gt_thin(vasp, thin_cells, z_vacuum)
    return align._prepare_gt(vasp)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                  "data", align.GT_CACHE))
    ap.add_argument("--check", action="store_true",
                    help="rebuild from the .vasp and verify the existing cache matches")
    ap.add_argument("--thin-cells", type=int, default=0,
                    help="[thin-ab] build the THIN slab GT of this many unit cells instead of "
                         "the full box (the aberration campaign uses 5); write it OUTSIDE the "
                         "package data dir and pass that dir to atomfind via --data-dir")
    ap.add_argument("--z-vacuum", type=float, default=4.0,
                    help="[thin-ab] vacuum padding each side along the beam (sim Z_VACUUM)")
    a = ap.parse_args()

    pos, Z = build(a.thin_cells, a.z_vacuum)
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
