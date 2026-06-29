#!/usr/bin/env python
"""Merge scan-tile sims (simulate_4dstem.py --scan-tile I/N) into one full dataset.

Each tile dir holds data_dp.hdf5 + data_position.hdf5 + tile_range.npy = [g0, g1, total].
Tiles are contiguous y-fastest blocks of the global scan, so we order by g0, verify
full coverage with no gaps/overlaps, then STREAM each tile into its slice of the output
(never materialising the whole stack — peak memory is one tile). probe_initial.mat and
sim_meta.mat are scan-independent, copied from the first tile.

  python merge_tiles.py --tiles-dir sim_out_<tag>/tiles --out-dir sim_out_<tag>
"""
import argparse, glob, shutil
from pathlib import Path
import numpy as np
import h5py


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tiles-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--scan", default="01")
    args = ap.parse_args(argv)

    tiles = []
    for dp in glob.glob(str(args.tiles_dir / "*" / args.scan / "data_dp.hdf5")):
        d = Path(dp).parent
        g0, g1, total = (int(v) for v in np.load(d / "tile_range.npy"))
        tiles.append((g0, g1, total, d))
    if not tiles:
        raise SystemExit(f"no tiles found under {args.tiles_dir}/*/{args.scan}/")
    tiles.sort(key=lambda t: t[0])
    total = tiles[0][2]

    # verify contiguous, complete, non-overlapping coverage
    cur = 0
    for g0, g1, tot, d in tiles:
        if tot != total or g0 != cur:
            raise SystemExit(f"coverage broken at {d}: expected start {cur}, got {g0} "
                             f"(total {tot} vs {total})")
        cur = g1
    if cur != total:
        raise SystemExit(f"tiles cover only {cur}/{total} positions — a tile is missing")
    print(f"[merge] {len(tiles)} tiles cover 0:{total} contiguously")

    # detector size from the first tile
    with h5py.File(tiles[0][3] / "data_dp.hdf5", "r") as f:
        _, nx_b, ny_b = f["dp"].shape
        dt = f["dp"].dtype

    outd = args.out_dir / args.scan
    outd.mkdir(parents=True, exist_ok=True)
    pos_out = np.zeros((2, total), np.float32)
    with h5py.File(outd / "data_dp.hdf5", "w") as fout:
        dset = fout.create_dataset("dp", shape=(total, nx_b, ny_b), dtype=dt)
        for g0, g1, tot, d in tiles:                       # stream tile -> its slice
            with h5py.File(d / "data_dp.hdf5", "r") as f:
                blk = f["dp"][...]
                assert blk.shape[0] == g1 - g0, f"{d}: {blk.shape[0]} frames != {g1-g0}"
                dset[g0:g1] = blk
            with h5py.File(d / "data_position.hdf5", "r") as f:
                pos_out[:, g0:g1] = f["probe_positions_0"][...]
            print(f"[merge]   tile {g0:>8d}:{g1:<8d}  ok")
    with h5py.File(outd / "data_position.hdf5", "w") as f:
        f.create_dataset("probe_positions_0", data=pos_out)
    print(f"[merge] wrote data_dp ({total},{nx_b},{ny_b}) + positions (2,{total}) -> {outd}")

    for name in ("probe_initial.mat", "sim_meta.mat"):
        src = tiles[0][3] / name
        if src.exists():
            shutil.copy2(src, outd / name)
    print(f"[merge] copied probe_initial.mat + sim_meta.mat from {tiles[0][3].parent.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
