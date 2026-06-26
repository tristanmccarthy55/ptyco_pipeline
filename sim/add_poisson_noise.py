#!/usr/bin/env python
"""Apply Poisson (shot) noise to a noiseless sim at a chosen electron dose.

POST-processing step, on purpose: noise depends only on dose, so we can sweep dose
WITHOUT re-running the (expensive) multislice sim. Reads a sim output directory,
scales the noiseless patterns to the incident electron count (dose * scan_step^2,
preserving the relative per-pattern totals so TDS/ADF contrast survives), then
Poisson-samples. Writes a NEW tagged dataset; the original sim is untouched.

  python add_poisson_noise.py --in-dir ../sim_out_step0.1_ph8 --dose 1e5
    -> ../sim_out_step0.1_ph8_dose1e5/01/  (data_dp noisy; positions/probe/meta linked)

Notes
-----
* dose is the incident electron fluence in e/A^2. Per pattern the incident count is
  dose * step^2 (step from the saved positions). Typical EM ptychography: 1e4-1e6.
* the recon's load_from_p needs avg count >~ 1e-4 e/pixel; at step 0.1 A that means
  dose >~ a few x 1e3 e/A^2 (it prints the achieved mean count so you can check).
"""
import argparse, os
from pathlib import Path
import numpy as np
import h5py


def derive_step(pos_path: Path) -> float:
    with h5py.File(pos_path, "r") as f:
        pos = f["probe_positions_0"][...]          # (2, Npos), [x; y], y-fastest
    dy = np.diff(pos[1].astype(np.float64)); dy = np.abs(dy[np.abs(dy) > 1e-6])
    return float(np.median(dy))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", required=True, type=Path,
                    help="sim output directory (the one CONTAINING the scan folder, e.g. sim_out_*/)")
    ap.add_argument("--dose", required=True, type=float, help="incident electron dose [e/A^2]")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="default: <in-dir>_dose<D> alongside the input")
    ap.add_argument("--scan", default="01", help="scan subfolder name (default 01)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    ind = args.in_dir / args.scan
    dp_in, pos_in = ind / "data_dp.hdf5", ind / "data_position.hdf5"
    for p in (dp_in, pos_in):
        if not p.exists():
            raise SystemExit(f"ERROR: missing {p} — point --in-dir at a finished sim dir.")

    step = derive_step(pos_in)
    epp = args.dose * step**2                       # incident electrons per pattern
    dose_tag = f"{args.dose:.0e}".replace("e+0", "e").replace("e+", "e").replace("e-0", "e-")
    out = args.out_dir or args.in_dir.parent / f"{args.in_dir.name}_dose{dose_tag}"
    outd = out / args.scan; outd.mkdir(parents=True, exist_ok=True)

    with h5py.File(dp_in, "r") as f:
        H = f["dp"][...].astype(np.float64)         # (Npos, Nx, Ny), noiseless (rel. intensities)
    npos = H.shape[0]
    # global scale so the MEAN pattern carries `epp` electrons -> preserves the
    # relative per-pattern totals (TDS/ADF contrast), unlike per-pattern renorm.
    mean_total = float(H.reshape(npos, -1).sum(1).mean())
    counts = H * (epp / max(mean_total, 1e-30))
    rng = np.random.default_rng(args.seed)
    noisy = rng.poisson(counts).astype(np.float32)

    with h5py.File(outd / "data_dp.hdf5", "w") as f:
        f.create_dataset("dp", data=noisy)
    print(f"[poisson] dose {args.dose:.0e} e/A^2, step {step:.3f} A -> {epp:.1f} e/pattern")
    print(f"[poisson] mean count: {noisy.mean():.3g} e/pixel  (peak {noisy.max():.0f}); "
          f"avg/pixel {'OK' if noisy.mean() > 1e-4 else 'LOW (<1e-4 — raise dose)'}")

    # link the unchanged inputs so the new dir is a complete, recon-ready dataset
    for name in ("data_position.hdf5", "probe_initial.mat", "sim_meta.mat"):
        src = ind / name
        if src.exists():
            dst = outd / name
            if dst.exists() or dst.is_symlink(): dst.unlink()
            os.symlink(os.path.relpath(src, outd), dst)
    print(f"[poisson] wrote {outd}  (positions/probe/meta symlinked from the noiseless sim)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
