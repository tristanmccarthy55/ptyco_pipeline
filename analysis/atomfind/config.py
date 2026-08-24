#!/usr/bin/env python
"""@file config.py
@brief Central configuration -- every per-volume knob for the atom-finding pipeline.

Everything that changes between reconstructions lives here, so the same code runs unchanged
across volumes: edit one Config or pass --preset (no code edits). Calibration constants are
kept identical to the validated overlay (analysis/column_cross_section_overlay.py), so the
alignment reproduces the tested overlay exactly.

Data files (the reconstructed volume, the measured kernels, the reference structure) are
resolved by NAME through data_path(), never by an absolute path, so the pipeline runs on a
machine that has never seen this author's Desktop. Search order:

  1. $ATOMFIND_DATA            -- an explicit directory (also settable with --data-dir)
  2. <package>/data/           -- what the peer-reproduction tarball unpacks into
  3. <repo>/sim/               -- the reference structure as committed in this repo
  4. ~/Desktop/                -- the original development location, kept last so existing
                                 local runs keep working unchanged

Outputs go to $ATOMFIND_OUT, else ./atomfind_out.

Requires numpy, scipy, h5py, matplotlib and ase. abtem is optional: it is used only to
prepare the ground-truth frame, and only when the precomputed data/gt_prepared.npz cache is
absent (see align.load_gt).
"""
from __future__ import annotations
from dataclasses import dataclass, field, replace
import os

_PKG = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_PKG, "..", ".."))

#: the ground-truth structure the simulation consumed
VASP_NAME = "PTO6_STO6_18_18_labyrinthPoscar.vasp"


def data_dirs() -> list[str]:
    """Directories searched for data files, most specific first (see the module docstring)."""
    dirs = []
    env = os.environ.get("ATOMFIND_DATA")
    if env:
        dirs.extend(os.path.abspath(os.path.expanduser(d)) for d in env.split(os.pathsep) if d)
    dirs += [os.path.join(_PKG, "data"),
             os.path.join(_REPO, "sim"),
             os.path.expanduser("~/Desktop")]
    return dirs


def data_path(name: str, required: bool = False) -> str:
    """@brief Resolve a data file by BASENAME across data_dirs().
    @param name basename, e.g. "NL70_new_vol.npy"; an absolute path is returned unchanged.
    @param required raise FileNotFoundError instead of returning the first candidate.
    @return the first existing match, else the preferred location (so error messages name it).
    """
    if os.path.isabs(name):
        return name
    for d in data_dirs():
        p = os.path.join(d, name)
        if os.path.exists(p):
            return os.path.abspath(p)
    if required:
        raise FileNotFoundError(
            f"{name!r} not found in any of:\n  " + "\n  ".join(data_dirs()) +
            "\nSet $ATOMFIND_DATA (or pass --data-dir) to the directory holding the data files.")
    return os.path.join(data_dirs()[0] if os.environ.get("ATOMFIND_DATA")
                        else os.path.join(_PKG, "data"), name)


def set_data_dir(path: str) -> None:
    """Prepend a directory to the search path for the rest of the process (used by --data-dir)."""
    path = os.path.abspath(os.path.expanduser(path))
    cur = os.environ.get("ATOMFIND_DATA")
    os.environ["ATOMFIND_DATA"] = path + (os.pathsep + cur if cur else "")


def default_out_dir() -> str:
    """@brief Where results are written: $ATOMFIND_OUT, else ./atomfind_out.

    Exception, mirroring the data search path: if ./atomfind_out does not exist but the
    original development location ~/Desktop/atomfind_out does, that one is kept, so the
    figure scripts that read cfg.out_dir (fig_check, paper/make_*_fig) are not silently
    repointed at an empty directory on this machine. A fresh checkout has neither and gets
    ./atomfind_out. --out overrides in every entry point.
    """
    env = os.environ.get("ATOMFIND_OUT")
    if env:
        return os.path.abspath(os.path.expanduser(env))
    local = os.path.abspath(os.path.join(os.getcwd(), "atomfind_out"))
    legacy = os.path.expanduser("~/Desktop/atomfind_out")
    if not os.path.isdir(local) and os.path.isdir(legacy):
        return legacy
    return local


@dataclass
class Config:
    # ---- which reconstructed volume -------------------------------------
    name: str = "NL70_coherent"
    recon_vol: str = "NL70_new_vol.npy"       # complex64 (nL, Ny, Nx); resolved by data_path()

    # ---- calibration ----------------------------------------------------
    # In-plane object pixel. The VALIDATED overlay uses scan_window / Nx
    # (20 A / 404 px); the physics object pixel is 0.0492 A. They agree to <1%
    # and the sub-pixel CAL below absorbs the residual. dx=None => derive as
    # scan_window_A / Nx at load time (generalises to any ROI size).
    dx: float | None = None
    dz: float = 0.999                 # depth spacing (A/layer). NL70=0.999, NL42=1.665
    scan_center_xy: tuple = (40.0, 20.0)   # sim scan centre (A), prepared-cell frame
    scan_window_A: float = 20.0            # sim scan window (A) -> field of view

    # ---- recon <-> ground-truth map (VALIDATED, section 3.5 / 11.4) -----
    # recon (row r, col c) -> GT physical  X = X0 + c*dx,  Y = Y0 + r*dx
    X0: float = 30.0
    Y0: float = 10.0
    # depth registration search: (sign, off_lo, off_hi) branches, fitted data-drivenly
    depth_branches: tuple = ((+1, -8.0, 4.0), (-1, 66.0, 78.0))

    # ---- physics (for the synthetic PSF + sanity), section 3.1 / 11.6 ---
    energy_keV: float = 300.0
    wavelength_A: float = 0.01969
    convergence_mrad: float = 100.0
    defocus_A: float = -20.0          # overfocus 20 A
    # Optical axial resolution ~ lambda / NA^2 ; NA = sin(conv). Used only as a
    # physically-motivated default for the synthetic PSF z-width.
    # (For 100 mrad, lambda/NA^2 ~ 2 A -> matches the handover's ~2 A axial limit.)

    # ---- dose (drives RL iteration cap / regularisation) ----------------
    # None = noiseless/coherent (current NL70). Set e/A^2 for the dosed data.
    dose_e_per_A2: float | None = None

    # ---- synthetic PSF parameters (physics stand-in) --------------------
    # In-plane: ~ lambda/(2*NA) = 0.0197/(2*sin 100mrad) ~ 0.1 A (huge aperture). We
    # set 0.4 A (a touch broader than the diffraction limit) so the synthetic and the
    # measured data PSF BRACKET a plausible range for the sensitivity test.
    psf_inplane_fwhm_A: float = 0.4   # tight in-plane blob (probe/aperture limited)
    psf_axial_fwhm_A: float = 2.2     # stretched along beam (missing cone / axial NA); data measures ~3 A
    psf_half_xy_A: float = 1.4        # kernel half-extent in-plane
    psf_half_z_A: float = 3.5         # kernel half-extent along beam
    psf_missingcone_tail: float = 0.35  # extra exponential z-tail weight (asymmetry)

    # ---- data-derived PSF (average isolated Pb blobs) -------------------
    psf_data_species: int = 82        # Pb: heaviest/cleanest reference atom
    psf_data_zhalf_A: float = 1.9     # z half-window < Pb-Pb period/2 (3.9/2) so
                                      # neighbouring in-column Pb do NOT leak in
    psf_data_xyhalf_A: float = 1.2
    psf_data_min_sep_A: float = 2.2   # in-plane isolation from any OTHER heavy column

    # ---- deconvolution --------------------------------------------------
    rl_iters: int = 12                # dose-scaled cap (see effective_rl_iters)
    rl_iters_noiseless: int = 15      # cap for the noiseless/coherent NL70 (dose=None)
    rl_filter_epsilon: float = 1e-3   # skimage RL regulariser (avoids /0 noise blow-up)
    rl_blowup_ratio: float = 12.0     # divergence guard: stop if est peak/median exceeds
    #                                   this x the input's (noise-amplification runaway)

    # ---- model-based fit (the O detector) -------------------------------
    fit_tile_A: float = 3.0           # in-plane tile edge (A) for the windowed NNLS
    fit_ridge: float = 0.0            # optional ridge on NNLS (0 = pure NNLS)

    # ---- blind atom finder (find.py) ------------------------------------
    # in-plane column detection
    find_min_sep_A: float = 1.4       # min in-plane separation of column peaks (< a/2)
    find_col_pct: float = 90.0        # depth-mean percentile threshold for a column peak
    find_track_win_px: int = 3        # per-layer centroid tracking half-window (px)
    # v1 axial spike deconvolution (kept as comparison baseline)
    find_profile_navg: int = 3        # in-plane half-box (px) for column z-profile extraction;
    #                                   MUST match psf.axial_kernel navg (validated: gold-kernel
    #                                   comb @ GT z reproduces real column profile, corr 0.94)
    spike_grid_A: float = 0.5         # z candidate spacing (2x oversampled vs dz)
    spike_merge_A: float = 1.0        # merge nonzero spikes closer than this -> one atom
    spike_min_frac: float = 0.05      # drop spikes below this fraction of the column max
    # v2: 3-D tube CLEAN (matching pursuit with the measured 3-D kernel) -- the default.
    # Rationale (measured): 1-D profile reduction fails on dim B-O columns (lean-tracking
    # wanders 0.5-0.65 A); the 3-D tube fit gives each atom its own xy -> Ti recall 65%->~86%+.
    tube_halfwidth_px: int = 10       # tube half-width (~0.5 A): contains the column lean
    # CLEAN stop floor. PORTABILITY: an absolute floor does NOT transfer -- the matched-
    # filter peak scales with the volume's phase amplitude (NL70 p99.9 0.333 vs dose-series
    # 0.188 at comparable noise), so a fixed floor cuts ~2x deeper there and decapitates
    # the weak species. Default is noise-relative: floor = clean_floor_k * sigma_MF,
    # sigma_MF estimated per tube by MAD over quiet voxels (see find.clean_floor_for).
    clean_floor_relative: bool = True
    # k swept on BOTH volumes (RESULTS sec.7): recall is monotone in k on each. k=2.0 is the
    # only setting meeting the dose1e10 targets (Ti 74%, confusion 13.7%) with ONE config;
    # it costs NL70 precision 0.97->0.95 while bulk O recall rises 92->96% -- a
    # precision/recall trade, not a degradation (the extra NL70 atoms are filterable via the
    # exported `quality` column). k=3.5 restores NL70 precision 0.97 but drops dose Ti to 63%.
    clean_floor_k: float = 2.0        # k in floor = k * sigma_matched-filter
    clean_floor: float = 0.4          # absolute fallback (used if relative=False or the
    #                                   sigma estimate degenerates)
    clean_max_atoms: int = 45         # per-tube atom cap (~70 A / 1.95 A spacing + margin)
    clean_nms_z_A: float = 1.5        # non-max suppression half-extent along z, in ANGSTROM
    #                                   (was 'layers' -- silently changed physical meaning
    #                                   when dz went 0.999 -> 0.666)
    clean_nms_xy_px: int = 4
    refine_sweeps: int = 2            # joint Gauss-Newton refinement passes per tube
    quality_min_corr: float = 0.5     # junk cut: min normalised patch-vs-kernel correlation
    # error-bar resolution floors (added in quadrature with the formal CRB): on noiseless
    # data the CRB is systematic-limited, not counting-limited. Floors are CALIBRATED to
    # 68% COVERAGE (fraction of atoms with |true error| <= 1 sigma) on the NL70 validation
    # volume -- the honest metric. (Median-ratio calibration was tried and rejected: it
    # gave ratio=1.0 but only ~50% coverage because the error tails are heavy.)
    # Per-species z floors = measured blind p68(|err_z|): Pb 0.24 / Ti 0.27 / O 0.31 A.
    # xy floor: after the fiducial affine refinement of the recon<->GT map the measured
    # p68(|err_xy|) is ~0.01 A; 0.015 keeps a margin for map transfer to new volumes.
    # PORTABILITY: these are NL70-calibrated. They do NOT transfer as constants (dose1e10
    # z-coverage ~20-28% with the raw NL70 numbers). The systematic localisation floor
    # scales with the volume's actual BLUR, so by default we rescale them by the measured
    # kernel FWHM relative to the reference kernel they were calibrated on (blind: uses the
    # PSF, not GT). Residual shortfall on poorly-registered volumes is REAL and is surfaced
    # by the sigma-coverage health warning rather than hidden.
    sigma_floor_xy_A: float = 0.015
    sigma_floor_z_A: float = 0.30     # fallback for unknown species
    sigma_floor_z_species: dict = field(default_factory=lambda: {82: 0.24, 22: 0.27, 8: 0.31})
    sigma_floor_scale_with_psf: bool = True
    sigma_floor_ref_fwhm_z_A: float = 1.0    # NL70 Pb kernel axial FWHM (calibration ref)
    sigma_floor_ref_fwhm_xy_A: float = 0.1   # NL70 Pb kernel in-plane FWHM
    sigma_floor_scale_cap: float = 4.0       # sanity clamp on the rescale factor
    # ---- uncertainty (see uncertainty.py) --------------------------------
    # Model sigma = joint-CRLB (+) kernel mismatch, then split-conformal Mondrian
    # calibration. The tuned per-species floors above are RETIRED (kept only as the
    # legacy path / documentation of what was replaced).
    kernel_mismatch_on: bool = True
    uq_alphas: tuple = (0.32, 0.05)   # 68% and 95% prediction intervals
    uq_min_stratum: int = 20          # below this, fall back to the pooled quantile
    uq_default_alpha: float = 0.05    # the level exported as the DEFAULT interval
    guided_sigma_scale: float = 1.4   # guided atoms: measured p68 ratio guided/blind ~1.35
    # exit-band inflation: the last ~10 A before the exit surface carries reconstruction
    # artifacts; measured coverage there under-ran (z 63%, x 2-sigma 85%) -> widen bars.
    exit_band_z_A: float = 56.0
    exit_sigma_scale: float = 1.4

    # ---- v3: preprocessing + lattice-aware species + guided re-detection ----
    preprocess_bg: bool = True        # subtract a smooth per-layer background (depth haze)
    bg_smooth_A: float = 3.0          # background blur scale (>> atom width, ~ column pitch)
    # column typing: k-means (k=3) on per-column amp_p75 -> A / B-O / pure-O bands
    # (measured zero-overlap: A 3.8-4.5, B-O 1.4-1.8, pure-O 0.56-0.87)
    comb_period_A: float = 3.9        # nominal plane period (fallback; fitted per column)
    slot_empty_A: float = 0.8         # comb slot is "empty" if no same-comb atom within this
    guided_species: tuple = (8,)      # guided re-detection for O ONLY: heavy atoms are found
    #                                   blind at ~97% bulk, and the few guided Pb/Ti were
    #                                   MISLOCATED with overconfident bars (guided-Ti z
    #                                   coverage measured at 9% -- confidently wrong)
    guided_min_corr: float = 0.35     # guided fits: lower quality bar (position prior pays)
    guided_gate_z_A: float = 0.7      # guided fit must stay within this of the predicted slot
    guided_gate_xy_A: float = 0.35    #   ... and this in-plane (from the column lean)
    guided_dedup_A: float = 1.2       # post-fit occupancy guard: reject a guided atom landing
    #                                   within this 3-D distance of one already accepted in the
    #                                   tube (blind or guided). GT-set: NO two real atoms are
    #                                   closer than 1.70 A (1st pctile 1.72), so any closer pair
    #                                   is non-physical; 1.2 sits above the duplicate scale
    #                                   (<=0.8) yet below where two DISTINCT atoms 1.7 A apart can
    #                                   land once each is localised with ~0.37 A z-scatter. (The
    #                                   pre-fit slot_empty check + the 0.7 A fit gate otherwise
    #                                   let an atom land ~0.1 A from an existing one -> O doubles.)
    #                                   Transfer-checked on the held-out NL105 dose1e10 volume:
    #                                   confusion 13.7->12.4%, Ti recall held (74% bulk) -- it
    #                                   helps there too, so it is not an NL70-specific tweak.
    # GT matching tolerances (validation only)
    match_tol_xy_A: float = 0.6
    match_tol_z_A: float = 2.0
    bulk_z_A: tuple = (10.0, 56.0)    # "bulk" depth band for the headline recall metric

    # ---- per-species empirical kernels ----------------------------------
    # Ti single-atom kernel (clean per docs/history/PSF_SIM_RESPONSE.md): broader axially than Pb
    # (weaker channeling). Used for species classification + light-atom refinement.
    # O's kernel is noise (do NOT use) -> O gets the Ti shape.
    ti_kernel_vol: str | None = None

    # ---- analysis window / trimming -------------------------------------
    zmax_show_A: float = 66.0         # trim exit-surface dumping-ground artifact
    trim_z_A: tuple = (2.0, 66.0)     # interior kept for deconv/detection (entrance/exit
    #                                   dumping-ground layers behave like noise under RL)
    reference_columns: tuple = ()     # (row,col) seeds; empty => auto from GT Pb columns

    # ---- output ---------------------------------------------------------
    out_dir: str = field(default_factory=default_out_dir)

    # ---- external empirical PSF (drops in when the sim thread delivers) --
    # Path to a reconstructed single-isolated-atom volume (.npy, same format).
    # When set, psf.empirical_psf() uses THIS instead of the data/synth stand-in.
    single_atom_vol: str | None = None
    single_atom_species: int = 82

    def vasp(self) -> str:
        """Absolute path to the ground-truth structure."""
        return data_path(VASP_NAME)

    def resolve(self) -> "Config":
        """Return a copy with every data field turned into an absolute path. Called once at
        load time so a missing file is reported by name, with the full search path."""
        out = replace(self)
        out.recon_vol = data_path(self.recon_vol, required=True)
        for f in ("single_atom_vol", "ti_kernel_vol"):
            v = getattr(self, f)
            if v:
                setattr(out, f, data_path(v, required=True))
        out.out_dir = os.path.abspath(os.path.expanduser(self.out_dir))
        return out

    def effective_rl_iters(self) -> int:
        """Cap RL iterations more aggressively at low dose (deconv amplifies noise)."""
        if self.dose_e_per_A2 is None:
            return self.rl_iters_noiseless        # noiseless -> can afford more iterations
        if self.dose_e_per_A2 >= 1e8:
            return self.rl_iters
        if self.dose_e_per_A2 >= 1e6:
            return max(4, self.rl_iters // 2)
        return max(3, self.rl_iters // 4)


# ---------------------------------------------------------------- presets
def preset(name: str) -> Config:
    """Named volumes. Add the better data here when it lands; run_atomfind picks by --preset."""
    presets = {
        # current data (this is what the prototype runs on). Gold empirical PSF =
        # the measured single-Pb system kernel (docs/history/PSF_SIM_RESPONSE.md sec.1). Per the sims
        # thread we use the Pb SHAPE for every species (it's the imaging-system response;
        # the element only sets amplitude), so this is the default kernel for all atoms.
        "NL70_coherent": Config(name="NL70_coherent",
                                recon_vol="NL70_new_vol.npy",
                                dz=0.999, dose_e_per_A2=None,
                                single_atom_vol="psf_Pb_NL70_vol.npy",
                                single_atom_species=82,
                                ti_kernel_vol="psf_Ti_NL70_vol.npy"),
        "NL42_coherent": Config(name="NL42_coherent",
                                recon_vol="NL42_new_vol.npy",
                                dz=1.665, dose_e_per_A2=None,
                                single_atom_vol="psf_Pb_NL70_vol.npy",
                                single_atom_species=82,
                                ti_kernel_vol="psf_Ti_NL70_vol.npy"),
        # TEMPLATE for the reviewer-2 data (0.1 A-binned, 16-phonon, dosed). Set recon_vol
        # + dose when it lands; dose-MATCH the kernel (d1e8 kernel <-> 1e8 recon). dz=0.666
        # (NL105), per docs/history/PSF_SIM_RESPONSE.md sec.4.
        "reviewer2": Config(name="reviewer2",
                            recon_vol="REVIEWER2_vol.npy",
                            dz=0.666, dose_e_per_A2=1e8,
                            single_atom_vol="psf_Pb_rev2_d1e8_vol.npy",
                            single_atom_species=82),
    }
    if name not in presets:
        raise KeyError(f"unknown preset {name!r}; have {list(presets)}")
    return presets[name]
