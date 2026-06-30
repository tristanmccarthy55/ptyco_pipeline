"""
simulate_4dstem.py — minimal abTEM 4D-STEM generator for PtychoShelves.

Simulates a 4D-STEM dataset with abTEM (multislice) and writes it straight into
the PtychoShelves worked-example format, ready for the MultiHollowPtycho pipeline
(hollow angle = 0 baseline) with no conversion step.

Strategy (full box, no cropping):
  - keep the FULL real-space cell (cropping would alias the broadened exit wave),
  - collect the full detector (~200 mrad),
  - return the measurement as a LAZY Dask array,
  - bin the detector 4x4 with dask.array.coarsen, then .compute() (~20 GB, fits RAM),
  - save a single data_dp.hdf5 via h5py (no Zarr).

Outputs (into OUT_DIR, default ./sim_out/01/):
  - data_dp.hdf5         /dp                  binned diffraction intensities
  - data_position.hdf5   /probe_positions_0   scan positions (Å)
  - probe_initial.mat    probe, p             initial probe (at the BINNED size)

Library: abTEM 1.0.5 (new API).  See DATA_FORMAT.md (ptyco repo) for the contract.

    python simulate_4dstem.py --test        # 3x3 local sanity check
    python simulate_4dstem.py               # full production scan (HPC)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import h5py
import dask.array as da
import ase.io
from scipy.io import savemat

import abtem


# ======================================================================
# CONFIG  (edit these for your data)
# ======================================================================
PROJECT_ROOT = Path(__file__).resolve().parent
POSCAR_PATH  = PROJECT_ROOT / "PTO6_STO6_18_18_labyrinthPoscar.vasp"
OUT_DIR      = PROJECT_ROOT / "sim_out" / "01"

# --- physics / optics ---
ENERGY_EV          = 300e3     # beam energy [eV]
CONVERGENCE_MRAD   = 100.0     # probe convergence semi-angle [mrad]
OVERFOCUS_A        = 20.0      # overfocus MAGNITUDE [Å] (2 nm). Sign handled below.

# --- detector / sampling / binning ---
DETECTOR_MAX_ANGLE_MRAD = 200.0   # full detector outer angle [mrad]
SLICE_THICKNESS_A       = 2.0     # multislice slice thickness [Å]
BIN_FACTOR              = 4       # NxN detector binning applied to the lazy array

# --- dose ---
# abTEM flux-normalises each pattern to total ≈ 1, so the per-pixel "photon count"
# is ~1e-5 and PtychoShelves rejects it (load_from_p: avg count must be >= 1e-4).
# Scale every pattern to DOSE_E electrons (integrated counts). No Poisson noise is
# added, so the data is noiseless; a high dose just keeps counts well above the
# threshold (avg ~8e4 e/pixel at 1e10, comfortable in float32).
DOSE_E = 1e10                     # electrons per diffraction pattern (noiseless)

# --- scan (production) ---
# Coordinates are in the PREPARED cell frame (after rotate+orthogonalize+pad+centre),
# NOT the raw POSCAR frame. The run prints the material bbox + exit-wave margins so
# you can place the scan on a feature with a safe buffer from the vacuum edges.
# (40, 20) targets the frustrated-vortex region with PTO/surface on the right.
SCAN_CENTER_X_A = 40.0
SCAN_CENTER_Y_A = 20.0
SCAN_WINDOW_A   = 20.0
SCAN_STEP_A     = 0.1

# --- thermal diffuse scattering (frozen phonons) ---
# OFF by default (coherent) so sim_out/ stays the validated coherent baseline. When
# N_PHONONS > 0 the multislice is run over N_PHONONS randomly displaced atomic
# configurations and the DIFFRACTION INTENSITIES are averaged (incoherent) -> TDS.
# Cost scales ~linearly with N_PHONONS. sigma is the rms 1-D displacement (~0.08 Å
# is a sensible room-temperature value for these elements).
N_PHONONS      = 0        # 0 = coherent; 8-16 for a realistic TDS sim
PHONON_SIGMA_A = 0.08     # rms thermal displacement [Å]
PHONON_SEED    = 1

# --- device ---
DEVICE = "gpu"   # "gpu" on the HPC L40; "cpu" for a laptop test

# --- structure prep ---
ROTATE_DEG_Y = -90.0     # rotation about y to set the beam (z) axis
Z_VACUUM_A   = 2.0       # center(axis=2, vacuum=...) padding each side along beam


def wavelength_a() -> float:
    """Electron wavelength [Å] at ENERGY_EV."""
    return float(abtem.core.energy.energy2wavelength(ENERGY_EV))


def potential_sampling_a() -> float:
    """Real-space sampling [Å] to band-limit beyond DETECTOR_MAX_ANGLE_MRAD.

    sampling <= lambda / (2 sin(theta_max)); abTEM antialiases to ~2/3 of the
    array Nyquist, so we sample a little finer to keep theta_max 'valid'.
    """
    theta = DETECTOR_MAX_ANGLE_MRAD * 1e-3
    return wavelength_a() / (2.0 * np.sin(theta)) / 1.5


# ======================================================================
# 1. STRUCTURE
# ======================================================================
def build_phantom_atoms():
    """Thin, ASYMMETRIC, CHIRAL phantom for the orientation (8-DOF) test.

    A handful of Pb atoms arranged as the letter 'F' in a single z-plane. 'F' has
    NO symmetry (C1), so all 8 dihedral-group orientations of the detector produce
    a visibly different reconstruction -> the correct custom_data_flip is the one
    that renders the F upright (vertical bar on the left, arms pointing +x, top arm
    higher than the middle arm). Single-slice is exact for this thin object, so the
    orientation result is NOT confounded by sample thickness.

    Box x,y are kept at the same 70.008 Å as the real sim so Ndpx/d_alpha (and thus
    the winning custom_data_flip) transfer directly to the real multislice recon.
    """
    from ase import Atoms
    side = 70.008                      # match the real sim's in-plane box -> same Ndpx
    zc = 4.0                           # atoms in one plane; box z below
    # 'F' around the scan centre (40, 20). x→ right, y→ up.
    coords = [(38, 14), (38, 17), (38, 20), (38, 23), (38, 26),   # vertical bar (left)
              (41, 26), (44, 26),                                  # top arm (long, +x)
              (41, 20)]                                            # middle arm (short)
    pos = [(x, y, zc) for (x, y) in coords]
    atoms = Atoms("Pb" * len(coords), positions=pos,
                  cell=[side, side, 8.0], pbc=True)
    print(f"[atoms] PHANTOM: {len(atoms)} Pb atoms as 'F' around (40,20); "
          f"box {side:.1f}×{side:.1f}×8.0 Å (thin)")
    return atoms, float(side)


def load_and_prepare_atoms():
    """Load POSCAR, orient for the beam, make the in-plane box square, add vacuum.

    1. read POSCAR
    2. rotate -90 deg about y   (sets which axis the beam looks down)
    3. orthogonalize_cell()     (required: multislice needs an orthogonal box)
    4. pad X and Y so the in-plane bounding box is square (NO cropping)
    5. center(axis=2, vacuum=2.0)
    """
    atoms = ase.io.read(str(POSCAR_PATH))
    print(f"[atoms] loaded {len(atoms)} atoms; cell "
          f"{np.round(atoms.cell.lengths(), 2)} Å")

    atoms.rotate(ROTATE_DEG_Y, "y", rotate_cell=True)
    atoms = abtem.orthogonalize_cell(atoms)
    Lx, Ly, Lz = atoms.cell.lengths()
    print(f"[atoms] after rotate+orthogonalize: {Lx:.2f} × {Ly:.2f} × {Lz:.2f} Å "
          f"(beam = z)")

    # pad X and Y to a square in-plane box (square box -> square DP sampling)
    side = max(Lx, Ly)
    atoms.cell[0, 0] = side
    atoms.cell[1, 1] = side
    atoms.center(axis=0)
    atoms.center(axis=1)
    atoms.center(axis=2, vacuum=Z_VACUUM_A)

    bx, by, bz = atoms.cell.lengths()
    print(f"[atoms] final box: {bx:.2f} × {by:.2f} × {bz:.2f} Å  "
          f"(beam path ≈ {bz:.1f} Å)")
    assert abs(bx - by) < 1e-6, f"in-plane box not square: {bx} vs {by}"
    return atoms, float(bx)


# ======================================================================
# 2. POTENTIAL & PROBES
# ======================================================================
def build_potential(atoms, announce=False):
    """Potential for ONE set of atom positions (a single frozen-phonon config, or the
    coherent structure). Frozen phonons are handled by LOOPING this in run_scan_binned
    (one config in memory at a time) rather than wrapping FrozenPhonons here, which
    would force abTEM to hold all configs at once -> OOM at 16 configs."""
    sampling = potential_sampling_a()
    pot = abtem.Potential(
        atoms,
        sampling=sampling,
        slice_thickness=SLICE_THICKNESS_A,
        parametrization="lobato",
        projection="infinite",
        device=DEVICE,
    )
    if announce:
        print(f"[potential] sampling = {sampling:.4f} Å  ->  gpts = {pot.gpts}  "
              f"({pot.num_slices} slices)")
    return pot


def _defocus():
    # Defocus sign derived from abTEM 1.0.5 source (not assumed):
    #   chi = (2*pi/lambda)*(1/2 * alpha^2 * C10 + ...), T(k)=exp(-i*chi),
    #   and abTEM defines defocus = -C10  =>  T(k) = exp(+i*pi*lambda*defocus*k^2).
    #   abTEM's Fresnel propagator (forward, dz>0) is exp(-i*pi*lambda*dz*k^2),
    #   so T = P(dz = -defocus). A NEGATIVE defocus therefore builds the entrance
    #   wave as the in-focus probe propagated FORWARD by |defocus| -> the crossover
    #   is |defocus| ABOVE the entrance surface (in vacuum) = OVERFOCUS.
    # => overfocus of magnitude OVERFOCUS_A is defocus = -OVERFOCUS_A.
    return -OVERFOCUS_A


def build_probe(potential):
    """Probe on the full simulation grid, used for the scan (accurate multislice)."""
    probe = abtem.Probe(
        energy=ENERGY_EV,
        semiangle_cutoff=CONVERGENCE_MRAD,
        defocus=_defocus(),
        device=DEVICE,
    )
    probe.grid.match(potential)
    print(f"[probe] semiangle = {CONVERGENCE_MRAD:.0f} mrad, "
          f"abTEM defocus = {_defocus():+.1f} Å (overfocus {OVERFOCUS_A:.1f} Å: "
          f"crossover before entrance surface)")
    return probe


def build_initial_probe(n_b: int, box_a: float):
    """Probe on the BINNED grid for probe_initial.mat.

    The binned detector has reciprocal pixel dk_b = BIN_FACTOR / box and N_b pixels.
    A probe with gpts=N_b and real-space extent = box / BIN_FACTOR has exactly that
    dk_b, the same object pixel dx = box/(N_b*BIN_FACTOR), and the same outer angle
    as the binned DP — so it is geometrically consistent with the saved data.
    """
    extent = box_a / BIN_FACTOR
    probe = abtem.Probe(
        energy=ENERGY_EV,
        semiangle_cutoff=CONVERGENCE_MRAD,
        defocus=_defocus(),
        gpts=(n_b, n_b),
        extent=(extent, extent),
        device="cpu",   # tiny; keep off the GPU
    )
    return np.asarray(probe.build().compute().array).astype(np.complex64)


# ======================================================================
# 3. SCAN + LAZY 4x4 BINNING
# ======================================================================
def make_scan(test_mode: bool, tile=None):
    """Return (scan, positions_xy (M,2) y-fastest, ny, (g0, g1, total)).

    tile=(I, N): take x-band I of N from the FULL grid's positions and scan them as a
    CustomScan -> bit-exact tiling (no grid re-derivation). The band is a contiguous
    block [g0:g1] of the global y-fastest order (x is the slow axis), so concatenating
    tiles 0..N-1 in order reconstructs the full dataset exactly. merge_tiles.py does that.
    """
    half = SCAN_WINDOW_A / 2.0
    start = (SCAN_CENTER_X_A - half, SCAN_CENTER_Y_A - half)
    end   = (SCAN_CENTER_X_A + half, SCAN_CENTER_Y_A + half)
    if test_mode:
        start = (SCAN_CENTER_X_A - 0.2, SCAN_CENTER_Y_A - 0.2)
        end   = (SCAN_CENTER_X_A + 0.2, SCAN_CENTER_Y_A + 0.2)
        full = abtem.GridScan(start=start, end=end, gpts=(3, 3))
    else:
        full = abtem.GridScan(start=start, end=end, sampling=SCAN_STEP_A)
    nx, ny = int(full.gpts[0]), int(full.gpts[1])
    full_pos = np.asarray(full.get_positions())          # (nx, ny, 2), y-fastest on C-flatten
    total = nx * ny
    if tile is None:
        print(f"[scan] {nx} × {ny} = {total} positions, "
              f"{tuple(np.round(start,2))} -> {tuple(np.round(end,2))} Å")
        return full, full_pos.reshape(-1, 2), ny, (0, total, total)
    I, N = tile
    assert 0 <= I < N <= nx, f"tile {I}/{N} invalid for nx={nx}"
    x0, x1 = I * nx // N, (I + 1) * nx // N
    band = full_pos[x0:x1].reshape(-1, 2)
    g0, g1 = x0 * ny, x1 * ny
    print(f"[scan] TILE {I+1}/{N}: x-rows [{x0}:{x1}] of {nx} -> {band.shape[0]} positions "
          f"= global [{g0}:{g1}] of {total}")
    return abtem.CustomScan(band), band, ny, (g0, g1, total)


def report_scan_geometry(atoms, beam_thickness_a: float):
    """Print where the scan sits relative to the material + the exit-wave halo,
    so a vacuum-padding shift can't silently mis-place the scan again.

    x is vacuum-padded (non-periodic) -> the halo must stay clear of the box edge
    there or multislice wraps/aliases. y is periodic (wrap is real material)."""
    pos = atoms.get_positions()
    bx, by, _ = atoms.cell.lengths()
    mx0, mx1 = pos[:, 0].min(), pos[:, 0].max()
    my0, my1 = pos[:, 1].min(), pos[:, 1].max()
    half = SCAN_WINDOW_A / 2.0
    # exit-wave radius ≈ geometric beam spread through the sample + probe radius
    ew_r = beam_thickness_a * np.tan(CONVERGENCE_MRAD * 1e-3) + OVERFOCUS_A * np.tan(CONVERGENCE_MRAD * 1e-3)
    hx0, hx1 = SCAN_CENTER_X_A - half - ew_r, SCAN_CENTER_X_A + half + ew_r
    print(f"[geom] material x[{mx0:.1f},{mx1:.1f}] y[{my0:.1f},{my1:.1f}]  box {bx:.1f}×{by:.1f} Å")
    print(f"[geom] scan centre ({SCAN_CENTER_X_A},{SCAN_CENTER_Y_A}) window {SCAN_WINDOW_A} Å; "
          f"exit-wave radius ≈ {ew_r:.1f} Å")
    print(f"[geom] x halo [{hx0:.1f},{hx1:.1f}] -> margin to box edges "
          f"{hx0:.1f} / {bx-hx1:.1f} Å  (x is VACUUM-padded; keep >~5 Å)")
    if hx0 < 5 or (bx - hx1) < 5:
        print("[geom] WARNING: exit-wave halo is close to the box edge in x -> "
              "risk of multislice wraparound; move the scan centre inward.")


def _crop_center_to_multiple(n_u: int, factor: int) -> tuple[int, int]:
    """Largest N_c <= N_u that is a multiple of 2*factor (keeps DC centred after
    binning), and the symmetric start index. Returns (start, N_c)."""
    n_c = (n_u // (2 * factor)) * (2 * factor)
    start = (n_u - n_c) // 2
    return start, n_c


def _scan_one_config(probe, potential, scan, detector):
    """Scan one (single-config) potential -> binned (M, N_b, N_b) NumPy, y-fastest."""
    meas = probe.scan(potential, scan=scan, detectors=detector, lazy=True)
    lazy = da.asarray(meas.array)                    # (..., N_u, N_u); no phonon axis here
    n_u = int(lazy.shape[-1])
    s, n_c = _crop_center_to_multiple(n_u, BIN_FACTOR)
    binned = da.coarsen(np.sum, lazy[..., s:s + n_c, s:s + n_c],
                        {lazy.ndim - 2: BIN_FACTOR, lazy.ndim - 1: BIN_FACTOR})
    n_b = n_c // BIN_FACTOR
    arr = np.asarray(binned.compute()).astype(np.float64).reshape(-1, n_b, n_b)
    return arr, n_b, n_u, n_c


def run_scan_binned(probe, atoms, scan):
    """Scan -> binned (M, N_b, N_b). Frozen phonons are processed ONE CONFIG AT A TIME
    (build that config's potential, scan, bin, accumulate the incoherent average), so
    peak memory is that of a single coherent sim regardless of phonon count — the fix
    for the 16-config OOM. Same total work as the ensemble path; just memory-bounded."""
    detector = abtem.PixelatedDetector(max_angle=DETECTOR_MAX_ANGLE_MRAD)
    if N_PHONONS and N_PHONONS > 0:
        print(f"[phonons] {N_PHONONS} configs, sigma={PHONON_SIGMA_A} Å — one at a time "
              f"(bounded memory; incoherent TDS average)")
        configs = list(abtem.FrozenPhonons(atoms, num_configs=N_PHONONS,
                                           sigmas=PHONON_SIGMA_A, seed=PHONON_SEED))
        acc = None
        for i, cfg in enumerate(configs):
            arr_i, n_b, n_u, n_c = _scan_one_config(
                probe, build_potential(cfg, announce=(i == 0)), scan, detector)
            acc = arr_i if acc is None else acc + arr_i
            print(f"[phonons] config {i+1}/{N_PHONONS} done")
        arr = (acc / N_PHONONS).astype(np.float32)
    else:
        print("[phonons] OFF (coherent — no TDS)")
        arr_i, n_b, n_u, n_c = _scan_one_config(
            probe, build_potential(atoms, announce=True), scan, detector)
        arr = arr_i.astype(np.float32)
    print(f"[bin] detector {n_u} -> crop {n_c} -> {BIN_FACTOR}×{BIN_FACTOR} bin -> N_b = {n_b}")
    print(f"[bin] binned measurement: {arr.shape[0]} positions × {n_b}×{n_b}")
    return arr


# ======================================================================
# 4. SAVE  (PtychoShelves contract + C-vs-Fortran handling)
# ======================================================================
def save_outputs(arr, pos_xy, out_dir: Path, box_a: float):
    """arr: (M, N_b, N_b) binned intensities (y-fastest).  pos_xy: (M, 2) [x, y] Å."""
    out_dir.mkdir(parents=True, exist_ok=True)
    npos, n_b, n_bx = arr.shape
    assert n_b == n_bx, f"binned DP not square: {n_b}x{n_bx}"
    A = arr.astype(np.float64)

    # FIXED dose scale (NOT normalised by this run's own mean) so independently
    # simulated scan tiles share one consistent scale and merge seamlessly. abTEM
    # flux-normalises each pattern (~1), so xDOSE_E gives ~DOSE_E e/pattern. The Poisson
    # step renormalises this away anyway; it only has to be consistent across tiles and
    # above the recon's count floor.
    A *= DOSE_E
    A = A.astype(np.float32)
    print(f"[dose] ×{DOSE_E:.0e} (fixed)  "
          f"(avg {A.reshape(npos,-1).sum(1).mean()/(n_b*n_b):.3g} e/pixel)")

    # --- data_dp.hdf5 -------------------------------------------------
    # MATLAB h5read reverses axes: HDF5 (s0,s1,s2) -> MATLAB [s2,s1,s0] with
    # M(p,q,r)=H[r-1,q-1,p-1]. For MATLAB dp[N_b,N_b,Npos] with dp(dy,dx,k)=A[k,dy,dx]
    # the dataset must be (Npos, N_b_x, N_b_y) = A.transpose(0,2,1).
    H_dp = np.ascontiguousarray(A.transpose(0, 2, 1))
    dp_path = out_dir / "data_dp.hdf5"        # NOTE: .hdf5 (PtychoShelves loader name)
    with h5py.File(dp_path, "w") as f:
        f.create_dataset("dp", data=H_dp)
    print(f"[save] {dp_path}  (HDF5 /dp {H_dp.shape} -> MATLAB [{n_b},{n_b},{npos}])")

    # --- data_position.hdf5 ------------------------------------------
    pos_xy = np.asarray(pos_xy).astype(np.float32)   # (Npos, 2) [x, y]
    assert pos_xy.shape == (npos, 2), f"positions {pos_xy.shape} != ({npos}, 2)"
    pos_path = out_dir / "data_position.hdf5"
    with h5py.File(pos_path, "w") as f:
        f.create_dataset("probe_positions_0", data=pos_xy.T)   # (2, Npos)
    print(f"[save] {pos_path}  (HDF5 /probe_positions_0 {pos_xy.T.shape} "
          f"-> MATLAB [{npos},2])")

    # --- probe_initial.mat (binned grid; scipy for MATLAB load()+complex) ----
    probe_wave = build_initial_probe(n_b, box_a)
    itot = float(A.reshape(npos, -1).sum(axis=1).mean())
    probe_wave *= np.sqrt(itot) / (np.linalg.norm(probe_wave) + 1e-30)
    probe_wave = probe_wave.astype(np.complex64)
    assert probe_wave.shape == (n_b, n_b), \
        f"probe {probe_wave.shape} != DP size {(n_b, n_b)}"
    mat_path = out_dir / "probe_initial.mat"
    savemat(str(mat_path), {
        "probe": probe_wave,
        "p": {"binning": False, "detector": {"binning": False}},
    })
    print(f"[save] {mat_path}  (probe {probe_wave.shape} complex)")

    return dp_path, pos_path, A


# ======================================================================
# 5. SELF-TESTS + DRIVER GEOMETRY
# ======================================================================
def selftest_ordering(dp_path, pos_path, A, pos_xy, ny):
    """Round-trip the written files under MATLAB's axis-reversal convention to
    catch the C-vs-Fortran transpose landmine without needing MATLAB. Works for the
    full grid AND a single tile band (which is also a y-fastest contiguous block)."""
    npos, n_b, n_bx = A.shape

    with h5py.File(dp_path, "r") as f:
        H = f["dp"][...]                       # (Npos, N_b_x, N_b_y)
    M_dp = np.transpose(H, (2, 1, 0))          # MATLAB [N_b, N_b, Npos]
    assert M_dp.shape == (n_b, n_bx, npos), M_dp.shape
    assert np.array_equal(M_dp, np.transpose(A, (1, 2, 0))), \
        "dp ordering mismatch (detector/scan transpose)"

    with h5py.File(pos_path, "r") as f:
        Hp = f["probe_positions_0"][...]       # (2, Npos)
    assert np.array_equal(np.transpose(Hp, (1, 0)), pos_xy), "position ordering mismatch"

    assert n_b == n_bx, "DP not square"
    assert npos == pos_xy.shape[0], f"{npos} frames != {pos_xy.shape[0]} positions"

    # y-fastest: first ny positions share one x, y increases uniformly (float32 tol)
    xs = pos_xy[:ny, 0].astype(np.float64)
    ys = pos_xy[:ny, 1].astype(np.float64)
    assert np.allclose(xs, xs[0]), "x should be constant over the first ny points"
    if ny > 1:
        dystep = np.diff(ys)
        assert np.all(dystep > 0), "y should increase monotonically (y-fastest)"
        assert np.allclose(dystep, dystep[0], rtol=1e-3), "y step should be uniform"

    # DC at the array centre (centre-of-MASS, robust to a strong Bragg disk)
    mean_dp = A.mean(axis=0).astype(np.float64)
    tot = mean_dp.sum()
    yy, xx = np.indices(mean_dp.shape)
    cy = float((yy * mean_dp).sum() / tot)
    cx = float((xx * mean_dp).sum() / tot)
    exp = n_b // 2
    tol = max(3.0, 0.02 * n_b)
    assert abs(cy - exp) <= tol and abs(cx - exp) <= tol, \
        f"diffraction COM at ({cy:.1f},{cx:.1f}), expected near ({exp},{exp})"

    print(f"[selftest] OK: dp {M_dp.shape}, pos {pos_xy.shape}, "
          f"y-fastest verified, DC-centred (COM ({cy:.1f},{cx:.1f}) of {n_b})")


def write_driver_geometry(n_b: int, box_a: float, beam_thickness_a: float,
                          out_dir: Path):
    """Print AND save (sim_meta.mat) the binned geometry the MATLAB driver needs,
    so the driver reads it instead of hardcoding (and drifting from) the sim."""
    lam = wavelength_a()
    dk_b = BIN_FACTOR / box_a                       # binned recip pixel [1/Å]
    d_alpha_rad = dk_b * lam                         # angle per binned pixel [rad]
    d_alpha_mrad = d_alpha_rad * 1e3
    rbf = CONVERGENCE_MRAD / d_alpha_mrad           # BF disk radius [binned px]
    dx = 1.0 / (n_b * dk_b)                         # object pixel [Å]

    meta = {
        "Ndpx": int(n_b),
        "d_alpha_rad": float(d_alpha_rad),
        "d_alpha_mrad": float(d_alpha_mrad),
        "rbf": float(rbf),
        "dx_object_A": float(dx),
        "energy_kev": float(ENERGY_EV / 1e3),
        "box_A": float(box_a),
        "bin_factor": int(BIN_FACTOR),
        "beam_thickness_A": float(beam_thickness_a),
        "convergence_mrad": float(CONVERGENCE_MRAD),
        "overfocus_A": float(OVERFOCUS_A),
        "scan_step_A": float(SCAN_STEP_A),
        "n_phonons": int(N_PHONONS),
        "phonon_sigma_A": float(PHONON_SIGMA_A),
        "ADU": 1.0,
    }
    savemat(str(out_dir / "sim_meta.mat"), {"meta": meta})

    print("\n" + "=" * 64)
    print("DRIVER GEOMETRY (saved to sim_meta.mat; the MATLAB driver reads it):")
    print(f"  Ndpx (p.asize)      = {n_b}")
    print(f"  d_alpha             = {d_alpha_mrad:.4f} mrad/pixel")
    print(f"  rbf (BF disk radius)= {rbf:.1f} binned pixels")
    print(f"  object pixel dx     = {dx:.4f} Å")
    print(f"  HT                  = {ENERGY_EV / 1e3:.0f} keV")
    print(f"  beam thickness      = {beam_thickness_a:.1f} Å (sample along beam)")
    print(f"  ADU                 = 1   (synthetic intensities; engine rescales probe)")
    print("=" * 64)


# ======================================================================
# MAIN
# ======================================================================
def main(argv=None) -> int:
    global DEVICE, SLICE_THICKNESS_A, SCAN_STEP_A, DOSE_E, N_PHONONS, PHONON_SIGMA_A, PHONON_SEED
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--test", action="store_true",
                    help="Tiny 3x3 scan for fast local shape/geometry validation.")
    ap.add_argument("--device", default=DEVICE, choices=["gpu", "cpu"])
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--slice-thickness", type=float, default=SLICE_THICKNESS_A,
                    help="Multislice slice thickness [Å]. Larger = fewer slices "
                         "= faster (use for the geometry test campaign).")
    ap.add_argument("--scan-step", type=float, default=SCAN_STEP_A,
                    help="Scan step [Å]. Larger = fewer positions = faster "
                         "(~0.4 for the test campaign).")
    ap.add_argument("--dose", type=float, default=DOSE_E,
                    help="Electrons per diffraction pattern (noiseless scaling).")
    ap.add_argument("--phantom", action="store_true",
                    help="Thin asymmetric 'F' Pb phantom for the orientation (8-DOF) "
                         "test (single-slice exact; same Ndpx as the real sim).")
    ap.add_argument("--phonons", type=int, default=N_PHONONS,
                    help="Frozen-phonon configs (thermal diffuse scattering). "
                         "0 = coherent (default); 8-16 for a realistic sim.")
    ap.add_argument("--phonon-sigma", type=float, default=PHONON_SIGMA_A,
                    help="RMS thermal displacement [Å] for frozen phonons (~0.08 ≈ RT).")
    ap.add_argument("--phonon-seed", type=int, default=PHONON_SEED)
    ap.add_argument("--scan-tile", default=None,
                    help='HPC tiling: run only scan x-band I of N, as "I/N" (0-indexed). '
                         'Each tile is an independent job; sim/merge_tiles.py reassembles '
                         'them into the full dataset (bit-exact, consistent dose scale).')
    args = ap.parse_args(argv)
    DEVICE = args.device
    SLICE_THICKNESS_A = args.slice_thickness
    SCAN_STEP_A = args.scan_step
    DOSE_E = args.dose
    N_PHONONS = args.phonons
    PHONON_SIGMA_A = args.phonon_sigma
    PHONON_SEED = args.phonon_seed

    tile = None
    if args.scan_tile:
        I, N = (int(v) for v in args.scan_tile.split("/"))
        tile = (I, N)

    print("=" * 64)
    print(f"4D-STEM simulation  |  device={DEVICE}  |  test={args.test}  |  "
          f"phantom={args.phantom}  |  tile={args.scan_tile or 'none'}")
    print(f"lambda = {wavelength_a():.5f} Å  |  bin = {BIN_FACTOR}×{BIN_FACTOR}  |  "
          f"slice = {SLICE_THICKNESS_A} Å  |  step = {SCAN_STEP_A} Å  |  "
          f"phonons = {N_PHONONS or 'off (coherent)'}")
    print("=" * 64)

    atoms, box_a = build_phantom_atoms() if args.phantom else load_and_prepare_atoms()
    beam_thickness_a = float(atoms.cell.lengths()[2] - 2 * Z_VACUUM_A)
    report_scan_geometry(atoms, beam_thickness_a)
    probe = build_probe(build_potential(atoms, announce=True))   # ref potential -> probe grid
    scan, pos_xy, ny, grange = make_scan(args.test, tile)

    arr = run_scan_binned(probe, atoms, scan)      # (M, N_b, N_b); builds per-config potentials
    n_b = arr.shape[-1]

    dp_path, pos_path, A = save_outputs(arr, pos_xy, args.out_dir, box_a)
    selftest_ordering(dp_path, pos_path, A, pos_xy, ny)
    write_driver_geometry(n_b, box_a, beam_thickness_a, args.out_dir)
    if tile is not None:
        # record this band's global position range so merge can order + verify
        np.save(args.out_dir / "tile_range.npy", np.array(grange, dtype=np.int64))
        print(f"[tile] global position range {grange[0]}:{grange[1]} of {grange[2]}")

    print("\nDone. Wrote files to", args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
