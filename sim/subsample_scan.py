#!/usr/bin/env python
"""Downsample a 4D-STEM dataset's SCAN to a coarser regular grid (keep every f-th probe
position in each direction) -> a smaller dataset the recon can hold in RAM.

The full 0.05 A / 160k-position set is ~76 GB, and the recon needs ~2x that (patterns +
amplitudes) -> OOM on a 128 GB node. Subsample:
  f=2 -> 0.10 A / 200x200 = 40k  (~19 GB; more overlap than the 0.15 A baseline)
  f=3 -> 0.15 A / 133x133 = ~18k (~9 GB; the validated baseline size, safest)

Streams the big data_dp sequentially (low memory). Positions stay y-fastest, so the
result is a normal coarser dataset the recon launcher reconstructs unchanged. probe +
sim_meta are geometry-only (per-pattern) and copied as-is; the dose step derives the
new step from the positions, so nothing else needs editing.

  python subsample_scan.py --in-dir sim_out_step0.05_slice0.5_ph16s0.08 --factor 2
  -> sim_out_step0.05_slice0.5_ph16s0.08_sub2/01/
"""
import argparse, shutil
from pathlib import Path
import numpy as np
import h5py


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", required=True, type=Path, help="dataset dir (contains 01/)")
    ap.add_argument("--factor", required=True, type=int, help="keep every f-th position each way")
    ap.add_argument("--out-dir", type=Path, default=None, help="default: <in-dir>_sub<f>")
    ap.add_argument("--scan", default="01")
    args = ap.parse_args(argv)
    f = args.factor
    ind = args.in_dir / args.scan
    for name in ("data_dp.hdf5", "data_position.hdf5"):
        if not (ind / name).exists():
            raise SystemExit(f"ERROR: missing {ind/name}")

    # infer the scan grid (positions are y-fastest: global k = ix*ny + iy)
    with h5py.File(ind / "data_position.hdf5", "r") as fp:
        pos = fp["probe_positions_0"][...]                 # (2, Npos)
    npos = pos.shape[1]
    ny = int(len(np.unique(np.round(pos[1], 4))))
    nx = npos // ny
    if nx * ny != npos:
        raise SystemExit(f"grid inference failed: nx*ny={nx*ny} != Npos={npos}")
    k = np.arange(npos)
    mask = ((k // ny) % f == 0) & ((k % ny) % f == 0)      # keep every f-th row & col
    keep = np.where(mask)[0]
    nxk, nyk = len(np.unique(k[keep] // ny)), len(np.unique(k[keep] % ny))
    print(f"[sub] grid {nx}x{ny}={npos} -> every {f} -> {nxk}x{nyk}={len(keep)} positions")

    out = args.out_dir or args.in_dir.parent / f"{args.in_dir.name}_sub{f}"
    outd = out / args.scan; outd.mkdir(parents=True, exist_ok=True)

    # stream data_dp sequentially, filtering to kept rows (peak memory = one chunk)
    with h5py.File(ind / "data_dp.hdf5", "r") as fi, h5py.File(outd / "data_dp.hdf5", "w") as fo:
        dpi = fi["dp"]; N, nxb, nyb = dpi.shape
        dpo = fo.create_dataset("dp", shape=(len(keep), nxb, nyb), dtype=dpi.dtype)
        wp = 0; CH = 2000
        for i in range(0, N, CH):
            blk = dpi[i:i+CH]
            sub = blk[mask[i:i+CH]]
            dpo[wp:wp+sub.shape[0]] = sub
            wp += sub.shape[0]
        assert wp == len(keep), f"wrote {wp} != {len(keep)}"
    with h5py.File(outd / "data_position.hdf5", "w") as fo:
        fo.create_dataset("probe_positions_0", data=pos[:, keep])
    for name in ("probe_initial.mat", "sim_meta.mat"):
        if (ind / name).exists():
            shutil.copy2(ind / name, outd / name)
    print(f"[sub] wrote {outd}  ({len(keep)} positions, ~{len(keep)*nxb*nyb*4/1e9:.0f} GB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
