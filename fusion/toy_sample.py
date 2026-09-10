#!/usr/bin/env python
"""@file toy_sample.py
@brief Thin PbTiO3 polar membrane with four engineered domains -- the data-fusion testbed.

The proof-of-concept specimen: a free-standing c-axis-clamped PbTiO3 membrane, `n_z` unit cells
thick, viewed down the polar axis (beam = z), divided into four quadrant domains whose Ti
off-centring delta = r_Ti - centroid(O6) has the SAME magnitude everywhere but different
directions. Every column is homogeneous along the beam (all n_z cells identical), so the EELS
signal needs no depth deconvolution -- it reads the polarisation where the probe lands.

The domains are chosen so that neither probe alone is sufficient (see README.md):
  A vs B  identical |delta_z| and identical in-plane delta -> EELS-degenerate AND
          projection-degenerate; ONLY depth sectioning separates them (the headline pair).
  C vs D  identical |delta_z|, different in-plane azimuth, opposite sign.
  AB vs CD  different |delta_z| -> the EELS magnitude channel.

Lattice: tetragonal a x a x c with c || z in EVERY domain (an in-plane-clamped film, which is
what a PTO/STO superlattice membrane is). Only the polar DISPLACEMENT pattern rotates. That
keeps the box commensurate and removes every confound (thickness, strain, composition) so any
measured contrast is polarisation and nothing else.

    ~/hyperspy-bundle/bin/python toy_sample.py            # build + validate + write
    ~/hyperspy-bundle/bin/python toy_sample.py --preview  # + a ground-truth figure
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "eels"))
import config as C            # noqa: E402  (RefCrystal: a_tet, c_tet, z_Ti, z_O_ap, z_O_eq)

REF = C.REF
Z_OF = {"Pb": 82, "Ti": 22, "O": 8}


# ---------------------------------------------------------------- the polar displacement mode
def polar_mode() -> tuple[np.ndarray, float]:
    """@brief The P4mm polar distortion as a rigid two-sublattice shift, net-translation removed.

    In the standard PbTiO3 setting the whole O sublattice moves by +0.4874 A and Ti by +0.1565 A
    along c relative to Pb. Subtracting the mean over the five atoms leaves a pattern with ZERO
    net cell translation, so it can be rotated to any direction and tiled next to a differently
    oriented domain without a rigid step at the wall.

    @return (w, d_ti) -- w = per-sublattice weights [Pb, Ti, O] in A for a unit direction,
            d_ti = the resulting |delta_Ti| = |r_Ti - centroid(O6)| in A.
    """
    c = REF.c_tet
    w = np.array([0.0,                                    # Pb (origin of the standard setting)
                  (REF.z_Ti - 0.5) * c,                   # Ti
                  (REF.z_O_ap - 0.0) * c])                # O  (equatorial moves by the same)
    assert abs((REF.z_O_eq - 0.5) * c - w[2]) < 1e-9, "apical and equatorial O must shift alike"
    w = w - w.repeat([1, 1, 3]).mean()                    # 1 Pb + 1 Ti + 3 O -> zero net shift
    return w, float(w[2] - w[1])                          # delta_Ti = r_Ti - r_O  (negative)


# ---------------------------------------------------------------- domain table
def default_domains(d_ti: float) -> list[dict]:
    """@brief The four quadrant domains, |delta| fixed, direction engineered (see module doc).

    `theta` is the angle of delta from the beam; `azimuth` its in-plane direction; `sign` the
    sign of delta_z. Quadrant keys are (ix >= n/2, iy >= n/2).
    """
    def d(theta_deg, azimuth_deg, sign):
        t, p = np.radians(theta_deg), np.radians(azimuth_deg)
        v = np.array([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), sign * np.cos(t)])
        return v * abs(d_ti)

    return [
        dict(name="A", quad=(0, 0), theta_deg=20.0, azimuth_deg=0.0,  sign=+1, delta=d(20, 0,  +1)),
        dict(name="B", quad=(1, 0), theta_deg=20.0, azimuth_deg=0.0,  sign=-1, delta=d(20, 0,  -1)),
        dict(name="C", quad=(0, 1), theta_deg=60.0, azimuth_deg=0.0,  sign=+1, delta=d(60, 0,  +1)),
        dict(name="D", quad=(1, 1), theta_deg=60.0, azimuth_deg=90.0, sign=-1, delta=d(60, 90, -1)),
    ]


def uniform_domains(delta) -> list[dict]:
    """@brief One polarisation everywhere -- a controlled specimen for the sign-encoding test.

    The four-domain membrane is the demonstration; this is the measurement. Two uniform membranes
    that differ ONLY in sign(delta_z) can be scanned at IDENTICAL probe positions, which removes
    the intra-cell sampling variance that makes a domain-to-domain comparison meaningless.
    """
    delta = np.asarray(delta, float)
    th = np.degrees(np.arccos(abs(delta[2]) / np.linalg.norm(delta)))
    az = np.degrees(np.arctan2(delta[1], delta[0]))
    return [dict(name="U", quad=q, theta_deg=th, azimuth_deg=az,
                 sign=int(np.sign(delta[2])), delta=delta)
            for q in ((0, 0), (1, 0), (0, 1), (1, 1))]


# ---------------------------------------------------------------- builder
def build(n_lat: int = 12, n_z: int = 5, vacuum_A: float = 2.0, domains=None, cap: bool = True):
    """@brief Tile the domains into one periodic membrane.

    @param n_lat  in-plane cells per side (must be even: the quadrants split it in half)
    @param n_z    cells along the beam; every cell in a column is IDENTICAL (homogeneous depth)
    @param vacuum_A  vacuum padding each side along the beam only (x,y stay periodic)
    @param cap    close the slab with a second PbO plane so BOTH surfaces are PbO-terminated.
                  Without it the stack runs PbO...TiO2, an up/down asymmetry in composition that
                  would sit on top of the polar one and confound the sign-of-delta_z test.
    @return (atoms, truth) -- ASE Atoms and a dict of per-cell ground truth.
    """
    from ase import Atoms
    assert n_lat % 2 == 0, "n_lat must be even so the quadrant split lands on a cell boundary"
    a, c = REF.a_tet, REF.c_tet
    w, d_ti = polar_mode()
    domains = domains or default_domains(d_ti)

    # centrosymmetric basis (build_cells.build_perovskite at s=0), fractional
    base_frac = np.array([[0.0, 0.0, 0.0],      # Pb
                          [0.5, 0.5, 0.5],      # Ti
                          [0.5, 0.5, 0.0],      # O apical   (Ti-O-Ti chain || z)
                          [0.5, 0.0, 0.5],      # O equatorial (chain || x)
                          [0.0, 0.5, 0.5]])     # O equatorial (chain || y)
    species = ["Pb", "Ti", "O", "O", "O"]
    chain = np.array([-1, -1, 2, 0, 1])         # Ti-O-Ti chain axis per site (-1 = cation)
    base_cart = base_frac * np.array([a, a, c])

    half = n_lat // 2
    dom_of = {tuple(d["quad"]): d for d in domains}
    pos, sym, cell_id, site_id, chain_id, dom_id = [], [], [], [], [], []
    grid_dom = np.empty((n_lat, n_lat), dtype="<U1")
    grid_delta = np.zeros((n_lat, n_lat, 3))

    for ix in range(n_lat):
        for iy in range(n_lat):
            dom = dom_of[(int(ix >= half), int(iy >= half))]
            grid_dom[ix, iy] = dom["name"]
            grid_delta[ix, iy] = dom["delta"]
            u = -dom["delta"] / abs(d_ti)                       # unit dir of the CATION shift
            shift = np.outer(w.repeat([1, 1, 3]), u)            # (5,3) per-site displacement
            for iz in range(n_z):
                origin = np.array([ix * a, iy * a, iz * c])
                pos.append(base_cart + shift + origin)
                sym.extend(species)
                cell_id.extend([(ix, iy, iz)] * 5)
                site_id.extend(range(5))
                chain_id.extend(chain)
                dom_id.extend([dom["name"]] * 5)
            if cap:                                             # closing PbO plane at z = n_z*c
                origin = np.array([ix * a, iy * a, n_z * c])
                pos.append(base_cart[[0, 2]] + shift[[0, 2]] + origin)
                sym.extend(["Pb", "O"])
                cell_id.extend([(ix, iy, n_z)] * 2)
                site_id.extend([0, 2])
                chain_id.extend([-1, 2])
                dom_id.extend([dom["name"]] * 2)

    pos = np.concatenate(pos, axis=0)
    atoms = Atoms(symbols=sym, positions=pos,
                  cell=[n_lat * a, n_lat * a, n_z * c + (c if cap else 0)], pbc=True)
    atoms.center(axis=2, vacuum=vacuum_A)                       # free-standing along the beam

    zp = atoms.get_positions()[:, 2]
    truth = dict(
        n_lat=n_lat, n_z=n_z, a=a, c=c, vacuum_A=vacuum_A, cap=int(cap),
        box_A=float(atoms.cell.lengths()[0]),
        box_z_A=float(atoms.cell.lengths()[2]),
        slab_thickness_A=float(zp.max() - zp.min()),
        delta_Ti_A=abs(d_ti),
        domain_grid=grid_dom, delta_grid=grid_delta,
        cell_id=np.array(cell_id), site_id=np.array(site_id),
        chain_axis=np.array(chain_id), atom_domain=np.array(dom_id),
        names=np.array([d["name"] for d in domains]),
        theta_deg=np.array([d["theta_deg"] for d in domains]),
        azimuth_deg=np.array([d["azimuth_deg"] for d in domains]),
        sign=np.array([d["sign"] for d in domains]),
        deltas=np.array([d["delta"] for d in domains]),
        quads=np.array([d["quad"] for d in domains]),
    )
    return atoms, truth


# ---------------------------------------------------------------- validation (the build gate)
def domain_index(atoms, truth):
    """@brief Per-Ti domain label, cell index and a BULK mask (cage entirely inside one domain).

    A Ti's O6 cage draws two of its equatorial O from the +x and +y neighbouring cells, so a Ti in
    the last cell before a domain wall genuinely has a mixed cage. Those columns are real physics
    (the wall), not a build error, but they must be excluded from a gate that asks whether the
    build is exact. `bulk` keeps the Ti whose whole cage sits in one domain.
    """
    n_lat, a = truth["n_lat"], truth["a"]
    ti = atoms.get_positions()[atoms.get_atomic_numbers() == 22]
    ix = np.floor(ti[:, 0] / a).astype(int) % n_lat
    iy = np.floor(ti[:, 1] / a).astype(int) % n_lat
    dom = truth["domain_grid"][ix, iy]
    same = (truth["domain_grid"][(ix + 1) % n_lat, iy] == dom) & \
           (truth["domain_grid"][ix, (iy + 1) % n_lat] == dom)
    return dom, ix, iy, same


def validate(atoms, truth, verbose: bool = True) -> bool:
    """@brief Assert the membrane really carries the designed polarisation field.

    Gates: atom count; both surfaces PbO-terminated (no composition asymmetry to confound the
    sign test); delta_Ti in each domain's BULK matches the design vector exactly; theta from the
    beam and sign(delta_z) as designed; every column homogeneous along the beam; and one
    undistorted lattice throughout, so any measured contrast is polarisation and nothing else.
    Wall columns (mixed cage) are reported, not gated -- they are the physical domain wall.
    """
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "eels"))
    import build_cells as B
    ok = True
    n_lat, n_z, a, c = truth["n_lat"], truth["n_z"], truth["a"], truth["c"]

    n_expect = n_lat * n_lat * (n_z * 5 + 2 * int(truth["cap"]))
    ok &= _gate(len(atoms) == n_expect, f"atom count {len(atoms)} == {n_expect}", verbose)

    Z, zp = atoms.get_atomic_numbers(), atoms.get_positions()[:, 2]
    bot = set(np.unique(Z[zp < zp.min() + 1.0]))
    top = set(np.unique(Z[zp > zp.max() - 1.0]))
    ok &= _gate(bot == top == {8, 82},
                f"both surfaces PbO-terminated (bottom {sorted(bot)}, top {sorted(top)}) "
                f"-> no composition up/down asymmetry", verbose)

    d_all = B.offcentering(atoms, 22, 6)                        # (N_Ti, 3), full PBC
    dom, ix, iy, bulk = domain_index(atoms, truth)
    for k, name in enumerate(truth["names"]):
        want = truth["deltas"][k]
        m = (dom == name) & bulk
        got = d_all[m]
        err = np.abs(got - want).max()
        th = np.degrees(np.arccos(np.clip(abs(got[:, 2]) / np.linalg.norm(got, axis=1), -1, 1)))
        ok &= _gate(err < 1e-6, f"domain {name} bulk (n={m.sum()}): delta = "
                                f"{np.round(got.mean(0), 4)} (target {np.round(want, 4)}, "
                                f"max err {err:.1e} A)", verbose)
        ok &= _gate(abs(th.mean() - truth["theta_deg"][k]) < 1e-6 and
                    np.all(np.sign(got[:, 2]) == truth["sign"][k]),
                    f"domain {name} bulk: theta = {th.mean():.2f} deg "
                    f"(target {truth['theta_deg'][k]:.1f}), sign(dz) = "
                    f"{truth['sign'][k]:+d} throughout", verbose)

    key = ix * n_lat + iy                                       # column homogeneity along beam
    spread = max(np.abs(d_all[(key == k) & bulk] - d_all[(key == k) & bulk].mean(0)).max()
                 for k in np.unique(key[bulk]))
    ok &= _gate(spread < 1e-9,
                f"bulk columns homogeneous along beam over all {n_z} cells "
                f"(max spread {spread:.1e} A) -> EELS needs no depth deconvolution", verbose)

    # one undistorted lattice: every Pb must sit on its ideal node PLUS exactly its own domain's
    # polar shift -- so a wall carries no rigid translation, only the change of polar direction.
    w, d_ti = polar_mode()
    pb = atoms.get_positions()[Z == 82]
    px = np.floor(pb[:, 0] / a).astype(int) % n_lat
    py = np.floor(pb[:, 1] / a).astype(int) % n_lat
    u = -truth["delta_grid"][px, py] / abs(d_ti)                 # per-Pb cation-shift direction
    ideal = np.stack([px * a, py * a], 1) + w[0] * u[:, :2]
    res = np.abs(((pb[:, :2] - ideal + a / 2) % a) - a / 2).max()
    ok &= _gate(res < 1e-9, f"single lattice: every Pb on its node + its own polar shift "
                            f"(max residual {res:.1e} A) -> no rigid step at a wall", verbose)

    nwall = int((~bulk).sum() / max(n_z, 1))
    if verbose:
        print(f"  [info] {nwall} of {n_lat*n_lat} columns are wall columns (mixed O6 cage); "
              f"the analysis scores domain interiors")
        print(f"[toy] {'ALL GATES PASSED' if ok else 'GATES FAILED'}")
    return bool(ok)


def _gate(cond, msg, verbose):
    if verbose:
        print(f"  [{'ok ' if cond else 'FAIL'}] {msg}")
    return bool(cond)


# ---------------------------------------------------------------- I/O
def write(atoms, truth, out_dir: str):
    """@brief Write the POSCAR the simulator reads plus the ground truth the analysis scores against."""
    import ase.io
    os.makedirs(out_dir, exist_ok=True)
    vasp = os.path.join(out_dir, "toy_membrane.vasp")
    ase.io.write(vasp, atoms, format="vasp", direct=True, sort=True)
    npz = os.path.join(out_dir, "toy_truth.npz")
    np.savez(npz, **truth)
    print(f"[toy] wrote {vasp}  ({len(atoms)} atoms, box "
          f"{np.round(atoms.cell.lengths(), 3)} A)")
    print(f"[toy] wrote {npz}")
    return vasp, npz


def summary(truth):
    """@brief One table: what each domain is, and which probe can see it."""
    d = truth["deltas"]
    print(f"\n{'dom':>4} {'delta (A)':>26} {'theta':>7} {'|dz|/dmax':>10} {'in-plane (A)':>14}")
    for k, n in enumerate(truth["names"]):
        sz = abs(d[k, 2]) / truth["delta_Ti_A"]
        print(f"{n:>4} {str(np.round(d[k], 3)):>26} {truth['theta_deg'][k]:>6.1f} "
              f"{sz:>10.3f} {str(np.round(d[k, :2], 3)):>14}")
    print(f"\n  A vs B : same |dz|, same in-plane  -> EELS-blind AND projection-blind; "
          f"depth sectioning only")
    print(f"  C vs D : same |dz|, azimuth x vs y  -> in-plane ptychography; opposite sign")
    print(f"  AB/CD  : |dz|/dmax {abs(d[0,2])/truth['delta_Ti_A']:.2f} vs "
          f"{abs(d[2,2])/truth['delta_Ti_A']:.2f} -> the EELS magnitude channel")


def preview(atoms, truth, path: str):
    """@brief Ground-truth figure: the designed delta field and its two projections."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    g, dg, a = truth["domain_grid"], truth["delta_grid"], truth["a"]
    n = truth["n_lat"]
    x = (np.arange(n) + 0.5) * a
    X, Y = np.meshgrid(x, x, indexing="ij")

    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    im = ax[0].pcolormesh(X, Y, dg[:, :, 2], cmap="RdBu_r",
                          vmin=-truth["delta_Ti_A"], vmax=truth["delta_Ti_A"], shading="nearest")
    ax[0].set_title(r"ground truth  $\delta_z$  (sign = ptychography's job)")
    fig.colorbar(im, ax=ax[0], label=r"$\delta_z$ (Å)")

    im = ax[1].pcolormesh(X, Y, np.abs(dg[:, :, 2]), cmap="viridis",
                          vmin=0, vmax=truth["delta_Ti_A"], shading="nearest")
    ax[1].set_title(r"$|\delta_z|$  (magnitude = EELS's job)")
    fig.colorbar(im, ax=ax[1], label=r"$|\delta_z|$ (Å)")

    s = max(1, n // 12)
    ax[2].quiver(X[::s, ::s], Y[::s, ::s], dg[::s, ::s, 0], dg[::s, ::s, 1],
                 np.abs(dg[::s, ::s, 2]), cmap="viridis", scale=3.0, width=6e-3)
    ax[2].set_title(r"in-plane $\delta_{xy}$ (arrows), $|\delta_z|$ (colour)")
    for k, nm in enumerate(truth["names"]):
        qx, qy = truth["quads"][k]
        ax[2].text((qx + 0.5) * n * a / 2, (qy + 0.5) * n * a / 2, nm,
                   ha="center", va="center", fontsize=22, color="w",
                   path_effects=None, alpha=0.75)
    for A_ in ax:
        A_.set_aspect("equal"); A_.set_xlabel("x (Å)")
    ax[0].set_ylabel("y (Å)")
    fig.suptitle(f"toy polar membrane: {n}x{n}x{truth['n_z']} cells, "
                 f"{truth['slab_thickness_A']:.1f} Å thick, |δ| = {truth['delta_Ti_A']:.3f} Å")
    fig.savefig(path, dpi=150)
    print(f"[toy] wrote {path}")


def main(argv=None) -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-lat", type=int, default=12, help="in-plane cells per side (even)")
    ap.add_argument("--n-z", type=int, default=5, help="cells along the beam (column-homogeneous)")
    ap.add_argument("--vacuum", type=float, default=2.0, help="vacuum each side along the beam [A]")
    ap.add_argument("--out-dir", default=os.path.join(here, "sample"))
    ap.add_argument("--preview", action="store_true", help="also write the ground-truth figure")
    args = ap.parse_args(argv)

    atoms, truth = build(args.n_lat, args.n_z, args.vacuum)
    print(f"[toy] {len(atoms)} atoms | {args.n_lat}x{args.n_lat}x{args.n_z} cells | "
          f"box {np.round(atoms.cell.lengths(), 3)} A | slab {truth['slab_thickness_A']:.2f} A")
    summary(truth)
    print()
    if not validate(atoms, truth):
        return 1
    write(atoms, truth, args.out_dir)
    if args.preview:
        preview(atoms, truth, os.path.join(args.out_dir, "toy_truth.png"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
