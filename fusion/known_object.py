#!/usr/bin/env python
"""@file known_object.py
@brief Hand PtychoShelves the TRUE object: is the depth failure the optimiser, or the forward model?

Every reconstruction so far resolves the lattice in projection and returns no depth structure, at
24, 20 and 16 layers and at 200 and 600 iterations, on noiseless data with the exact probe and the
exact positions. Two explanations remain and they call for opposite fixes:

  basin     the engine's forward model agrees with the one that made the data, the true object is
            a far better fit, and a random start simply never finds it  -> fix the initialisation;
  mismatch  the engine CANNOT reproduce these data even given the truth (slicing, sampling,
            propagator or orientation convention differs from abTEM's), so the depth structure sits
            below its own model error and no setting will ever recover it -> fix the model.

This writes the true multislice object in exactly the frame the engine uses, so a reconstruction can
start from it. Run once with the object FROZEN (OBJECT_START=inf) to read the residual of the truth
under PtychoShelves' own model, and once with it free to see whether the engine keeps or destroys it.

The frame is not assumed. The true projected phase is cross-correlated against a finished
reconstruction of the same data over every periodic shift and both orientations (the detector goes
in transposed, custom_data_flip = [0 0 1]); the object is then cut from the periodic true potential
at that registration. A clear winner between the two orientations is required.

    ~/hyperspy-bundle/bin/python known_object.py --ref <a finished Niter*.mat of the same data> \\
        --nlayers 16 --out known_object_NL16.mat
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def read_full_object(path: str):
    """@brief The reconstruction's FULL object (all layers), in MATLAB index order [row, col].

    h5py hands back a MATLAB v7.3 array transposed, so each layer is transposed once on the way in:
    element [i, j] of the result is MATLAB object(i+1, j+1). The same convention is used when the
    true object is written back, which is what makes the round trip exact.
    """
    import h5py
    with h5py.File(path, "r") as f:
        refs = f["outputs"]["object"][:, 0]
        layers = []
        for r in refs:
            a = f[r][:]
            if a.dtype.names:
                a = a["real"] + 1j * a["imag"]
            layers.append(np.asarray(a).T)
        dx = float(np.asarray(f["outputs"]["pixel_size"]).ravel()[0]) * 1e10      # m -> A
        zd = np.asarray(f["outputs"]["z_distance"]).ravel()
    zd = zd[np.isfinite(zd)] * 1e10
    return np.asarray(layers), dx, zd


def true_transmission(vasp: str, n_layers: int, dx: float, energy_eV: float = 300e3):
    """@brief exp(i sigma V_j) for each of the engine's layers, on the periodic specimen box.

    Same parametrisation and projection as the simulation that made the data (Lobato, infinite);
    what differs is only the slicing and sampling, which are the engine's -- that difference IS part
    of what the frozen-object run measures.
    @return (T[layer, x, y], box_xy_A, box_z_A); abTEM indexes arrays [x, y].
    """
    import abtem
    import ase.io
    atoms = ase.io.read(vasp)
    Lx, Ly, Lz = atoms.cell.lengths()
    gx = int(round(Lx / dx))
    if abs(gx * dx - Lx) > 1e-3:
        raise SystemExit(f"specimen box {Lx:.4f} A is not a whole number of {dx:.5f} A pixels")
    pot = abtem.Potential(atoms, gpts=(gx, gx), slice_thickness=[Lz / n_layers] * n_layers,
                          parametrization="lobato", projection="infinite", device="cpu")
    V = np.asarray(pot.build().compute().array)                   # (n_layers, gx, gx)
    if V.shape[0] != n_layers:
        raise SystemExit(f"abTEM made {V.shape[0]} slices, the engine uses {n_layers}")
    sigma = float(abtem.core.energy.energy2sigma(energy_eV))
    return np.exp(1j * sigma * V), float(Lx), float(Lz)


def _highpass(P, dx: float, scale_A: float = 1.2):
    """@brief Keep the atomic-column structure, drop the slow background the recon carries."""
    from scipy.ndimage import gaussian_filter
    return P - gaussian_filter(P, scale_A / dx, mode="wrap")


def register(T, R, roi: int, dx: float, anchor: int, a: float):
    """@brief Find the orientation and sub-cell offset of the true object in the engine's frame.

    Translation is NOT searched globally: a periodic lattice matches almost equally at every cell
    translation, so a global peak is close to arbitrary. The absolute position is fixed instead by
    geometry -- PtychoShelves' positions are object-centred, so the object's centre pixel is the scan
    centre -- and only a residual within half a cell of that anchor is searched. Both projected
    phases are high-passed first so the atomic columns, not the recon's slow background, drive it.
    What is left to decide is the orientation, which the in-plane polarisation pattern settles: a
    transpose swaps an x-displaced domain for a y-displaced one.
    @return (orientation, s0, s1, corr, corr_other)
    """
    Pt = _highpass(np.angle(T).sum(0), dx)                        # true projected phase [x, y]
    Pr = _highpass(np.angle(R).sum(0), dx)                        # recon projected phase [row, col]
    n = R.shape[1]
    o = (n - roi) // 2
    W = Pr[o:o + roi, o:o + roi]
    W = W - W.mean()
    half = int(round(0.5 * a / dx))
    best = {}
    for name, Q in (("identity", Pt), ("transpose", Pt.T)):
        Qz = Q - Q.mean()
        top = (-9.0, 0, 0)
        for s0 in range(anchor - half, anchor + half + 1):
            rows = np.take(Qz, np.arange(s0, s0 + roi), axis=0, mode="wrap")
            for s1 in range(anchor - half, anchor + half + 1, 2):         # coarse, refined below
                patch = np.take(rows, np.arange(s1, s1 + roi), axis=1, mode="wrap")
                c = float((patch * W).sum() / np.sqrt((patch ** 2).sum() * (W ** 2).sum()))
                if c > top[0]:
                    top = (c, s0, s1)
        c0, b0, b1 = top
        for s0 in range(b0 - 1, b0 + 2):                               # refine to 1 px
            rows = np.take(Qz, np.arange(s0, s0 + roi), axis=0, mode="wrap")
            for s1 in range(b1 - 2, b1 + 3):
                patch = np.take(rows, np.arange(s1, s1 + roi), axis=1, mode="wrap")
                c = float((patch * W).sum() / np.sqrt((patch ** 2).sum() * (W ** 2).sum()))
                if c > c0:
                    c0, b0, b1 = c, s0, s1
        best[name] = (b0 % Qz.shape[0], b1 % Qz.shape[1], c0)
    win = max(best, key=lambda k: best[k][2])
    other = [k for k in best if k != win][0]
    s0, s1, ncc = best[win]
    return win, s0, s1, ncc, best[other][2]


def cut(T, orient: str, s0: int, s1: int, n_full: int, roi: int):
    """@brief The true object on the engine's full n_full x n_full grid, at the registration found."""
    Q = T if orient == "identity" else np.transpose(T, (0, 2, 1))
    o = (n_full - roi) // 2                       # the correlation window started o pixels in
    N = Q.shape[1]
    ii = (s0 - o + np.arange(n_full)) % N
    jj = (s1 - o + np.arange(n_full)) % N
    return Q[:, ii][:, :, jj]                     # (n_layers, n_full, n_full), [layer, row, col]


def depth_ratio(O, a: float, c: float, dz: float, dx: float):
    """@brief Sanity check on the object itself: the true one MUST show the lattice comb along z."""
    P = np.angle(O)
    N = P.shape[1]
    step = int(round(a / dx))
    prof = P[:, N // 2 - 2 * step:N // 2 + 2 * step, N // 2 - 2 * step:N // 2 + 2 * step]
    prof = prof.reshape(P.shape[0], -1)
    k = np.argsort(prof.std(0))[-50:]             # the strongest (column) pixels
    z = (np.arange(P.shape[0]) + 0.5) * dz
    Rz = prof[:, k] - prof[:, k].mean(0)
    f = np.fft.rfftfreq(len(z), d=dz)
    S = (np.abs(np.fft.rfft(Rz, axis=0)) ** 2).mean(1)
    kk = int(np.argmin(np.abs(f - 1.0 / c)))
    band = (f > 0.05) & (f < 0.5 / dz)
    return float(S[kk] / S[band].mean())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref", required=True, help="a finished Niter*.mat of the SAME data (for the frame)")
    ap.add_argument("--vasp", default=os.path.join(HERE, "sample", "toy_membrane.vasp"))
    ap.add_argument("--nlayers", type=int, required=True, help="layers of the run that will use it")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    R, dx, zd = read_full_object(args.ref)
    n_full = R.shape[1]
    roi = int(round(0.55 * n_full))                # central, well-constrained part of the recon
    print(f"reference recon : {R.shape[0]} layers, {n_full}^2 px at {dx:.5f} A")
    T, box, box_z = true_transmission(args.vasp, args.nlayers, dx)
    dz = box_z / args.nlayers
    print(f"true object     : {T.shape[0]} layers x {dz:.4f} A on the {T.shape[1]}^2 periodic box"
          + (f"  (reference used {len(zd)} x {zd[0]:.4f} A)" if len(zd) else ""))

    # geometric anchor: object centre pixel <-> scan centre (the box centre, the four-domain corner)
    o = (n_full - roi) // 2
    anchor = int(round((box / 2.0) / dx)) - n_full // 2 + o
    import toy_sample as TS0
    orient, s0, s1, ncc, ncc_other = register(T, R, roi, dx, anchor, TS0.REF.a_tet)
    print(f"geometric anchor: shift {anchor} px (object centre = scan centre)")
    print(f"registration    : {orient}, shift ({s0}, {s1}) px, correlation {ncc:.3f} "
          f"(other orientation {ncc_other:.3f})")
    if ncc < 0.5 or ncc - ncc_other < 0.1:
        raise SystemExit("registration is not unambiguous -- refusing to write an object in an "
                         "unknown frame")

    O = cut(T, orient, s0, s1, n_full, roi)
    import toy_sample as TS
    r_true = depth_ratio(O, TS.REF.a_tet, TS.REF.c_tet, dz, dx)
    print(f"sanity          : the true object's depth-comb ratio is {r_true:.1f} "
          f"(the reconstructions scored 0.4-0.9; this must be >> 1)")
    if r_true < 3:
        raise SystemExit("the injected object shows no depth comb itself -- the test would be void")

    from scipy.io import savemat
    obj = np.transpose(O, (1, 2, 0))[:, :, None, :].astype(np.complex64)   # (row, col, 1, layer)
    savemat(args.out, {"object": obj}, do_compression=False)
    print(f"wrote {args.out}  object {obj.shape} complex64 "
          f"({os.path.getsize(args.out)/1e6:.0f} MB) -- layers in dim 4, as load_from_p expects")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
