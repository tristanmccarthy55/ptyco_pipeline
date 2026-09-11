% MULTISLICE HOLLOW ptychography (MHP) of the fusion toy membrane.
%
% Identical to run_synthetic_recon_ML.m (which is left untouched, along with everything
% specific to the PTO/STO labyrinth) except for ONE thing: the detector has a HOLE in it.
% The global `mask1` marks the pixels the engine may use (1 = use, 0 = ignore); the hole is
% the low-angle disc the EELS spectrometer takes, so those pixels are excluded from the
% phase retrieval. See +engines/+GPU/private/modulus_constraint.m, which forms the
% engine's excluded set as mask = -(mask1-1).
%
% Reference: Yu Lei and Peng Wang, "Hollow multi-slice electron ptychography for
% simultaneous 3D structural imaging and EELS in 4D-STEM", arXiv:2506.22352.
%
% Env:
%   SIM_BASE   dataset root (contains 01/)                 [../fusion/runs/fusion/]
%   HSA        hollow semi-angle as a fraction of alpha    [0.75]  0 = full detector
%   NLAYERS    depth slices                                [24]
%   plus every env the parent driver reads (NITER, GROUPING, REGLAYER, PROBE_MODES, ...)
%
% Run from:  ptychoshelves-clean/ptycho/

%%
clear variables
addpath(pwd)
addpath(fullfile(pwd,'utils'))
addpath(fullfile(pwd,'utils_EM'))
addpath(core.find_base_package)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scan_string_format = '%02d';
idx_scan  = 1;
base_path = getenv('SIM_BASE');
if isempty(base_path)
    base_path = '../fusion/runs/fusion/';
end
scan_dir = fullfile(base_path, sprintf(scan_string_format, idx_scan));

%%%%%%%%%%%%%%%%%%%%%%%%%%%% geometry from sim_meta.mat %%%%%%%%%%%%%%%%%%%%%%%%%
meta_file = fullfile(scan_dir, 'sim_meta.mat');
if ~exist(meta_file, 'file')
    error('sim_meta.mat not found at %s — run simulate_4dstem.py first.', meta_file);
end
S = load(meta_file); meta = S.meta;
Ndpx    = double(meta.Ndpx);
d_alpha = double(meta.d_alpha_rad);
HT      = double(meta.energy_kev);
thick   = double(meta.beam_thickness_A);    % sample thickness along beam [Å]
fprintf('sim_meta: Ndpx=%d, d_alpha=%.4f mrad, HT=%.0f keV, thick=%.1f A\n', ...
        Ndpx, d_alpha*1e3, HT, thick);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% HOLLOW DETECTOR MASK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mask1: 1 = pixel used for ptychography, 0 = pixel inside the hole (sent to EELS).
% simulate_fusion.py writes one mask per hollow semi-angle next to the data, so the SAME
% simulated dataset serves the whole HSA sweep -- only the mask changes.
hsa_env = getenv('HSA');
if ~isempty(hsa_env); hsa = str2double(hsa_env); else; hsa = 0.75; end
global mask1
if hsa <= 0
    mask1 = ones(Ndpx, Ndpx, 'single');
    fprintf('hollow semi-angle = 0 (full detector, MHP baseline)\n');
else
    mask_file = fullfile(scan_dir, sprintf('mask_hsa%.2f.mat', hsa));
    if ~exist(mask_file, 'file')
        error(['mask %s not found -- simulate_fusion.py must be run with this HSA ' ...
               'in --hsa (it writes one mask per value).'], mask_file);
    end
    mask1 = single(load(mask_file).mask);
    if size(mask1,1) ~= Ndpx
        error('mask is %dx%d but Ndpx = %d', size(mask1,1), size(mask1,2), Ndpx);
    end
    alpha_mrad = double(meta.convergence_mrad);
    fprintf(['hollow semi-angle = %.2f alpha = %.1f mrad : %d of %d detector pixels ' ...
             '(%.1f%%) sent to EELS\n'], hsa, hsa*alpha_mrad, nnz(mask1==0), numel(mask1), ...
            100*nnz(mask1==0)/numel(mask1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% multislice reconstruction parameters %%%%%%%%%%%%%%%
% Probe modes model partial coherence (real data). Our sim is a single fully-
% coherent probe, so 1 mode is correct; >1 just adds spurious incoherent background.
% Default 1; override with PROBE_MODES (e.g. for real, partially-coherent data).
pm_env = getenv('PROBE_MODES');
if ~isempty(pm_env); Nprobe = round(str2double(pm_env)); else; Nprobe = 1; end
fprintf('probe_modes = %d\n', Nprobe);
% Nlayers: set via NLAYERS env (e.g. 41 for the 4 A sublattice, 74/82 for the 2 A
% oxygen sublattice — chosen INCOMMENSURATE with the ~1.95 A atomic planes to avoid
% depth plane-locking). Falls back to ~10 A slices if NLAYERS is unset.
nl_env = getenv('NLAYERS');
if ~isempty(nl_env)
    Nlayers = max(1, round(str2double(nl_env)));
else
    % thin membrane: depth sectioning is the whole point, so default to ~1 A slices
    % (the PbO and TiO2 planes sit 1.8-2.3 A apart and their ORDER is the sign of P_z).
    Nlayers = max(1, round(thick / 1.0));
end
delta_z = thick / Nlayers;                      % actual layer spacing [Å]
fprintf('multislice: Nlayers=%d, delta_z=%.3f A (planes ~1.95 A; ratio %.2f)\n', ...
        Nlayers, delta_z, 1.95/delta_z);

% two-engine schedule (coarse presolve -> full), modelled on the proven baseline.
% grouping = diffraction patterns per GPU batch; GPU memory ~ grouping * Ndp^2, so large Ndp
% (BIN=1 -> 1426) OOMs at 32 (needs 52 GB > 47 GB L40). GROUPING env overrides the full-engine
% value ("g_pre g_full" or a single g_full); default 32.
grouping                  = [64,  32];
gr_env = getenv('GROUPING');
if ~isempty(gr_env)
    gv = str2num(gr_env); %#ok<ST2NM>
    if numel(gv) >= 2; grouping = [round(gv(1)), round(gv(2))];
    else;              grouping = [grouping(1), round(gv(1))]; end
end
fprintf('grouping (per engine) = [%d %d]\n', grouping(1), grouping(2));
% iterations per engine; the lattice converges well before 200, so NITER lets the
% heavy deep/fine runs fit walltime (e.g. NITER=120 for 70-layer ~1 A slices).
ni_env = getenv('NITER');
if ~isempty(ni_env); Niter = [round(str2double(ni_env)), round(str2double(ni_env))]; else; Niter = [200, 200]; end
% Probe update start. CONFIRMED: the fixed-probe 7-layer run reproduced the refined
% run's lattice to-a-tee, so the simulated probe we hand in is CORRECT. For synthetic
% data we therefore DON'T refine it — fixing the (true) probe is both accurate and
% removes the probe-update instability that NaN'd the deep refined runs at iter 20.
% (The earlier 41-layer noise was under-constraint, not the probe.) Default: fixed.
% For REAL data (unknown probe) re-enable refinement via PROBE_START (e.g. 20, or 60
% to delay so the object settles first on deep runs).
% BLIND-FIT STABILITY: the full-resolution (2nd) engine's probe update diverges to NaNs on the
% thin weak-phase slab, while the coarse PRESOLVE (1st) engine updates the probe stably. So we
% release the probe in the presolve only and FIX it (inf) for the full engine — the recovered
% (presolve) probe is carried into the full-res object refinement. PROBE_START2 overrides the
% 2nd-engine start if you ever want to refine there too.
ps_env = getenv('PROBE_START'); ps2_env = getenv('PROBE_START2');
if ~isempty(ps_env)
    ps2 = inf; if ~isempty(ps2_env); ps2 = str2double(ps2_env); end
    Nst_probe = [str2double(ps_env), ps2];      % release in presolve; FIX in full engine (NaN-safe)
else
    Nst_probe = [inf, inf];                     % synthetic: known/true probe -> fixed
end
fprintf('probe_change_start (per engine) = [%g %g]\n', Nst_probe(1), Nst_probe(2));

% TEM aperture constraint (Zhen Chen): constrain the recovered probe to the aperture in Fourier
% space, where the mask = abs(fft2(probe_initial)) thresholded (load_from_p). OFF by default (a
% FIXED true probe needs no constraint). ESSENTIAL for BLIND retrieval (PROBE_START set) — without
% it the probe update absorbs high-frequency aliasing into a junk probe. PROBE_SUPPORT_FFT=1 turns
% it on. NOTE: this drives p.probe_support_tem — do NOT set eng.probe_support_fft to a scalar (it is
% a MASK array; a scalar makes the engine crop_pad it to a single centre pixel and zero the probe).
psf_env = getenv('PROBE_SUPPORT_FFT'); probe_support_tem = ~isempty(psf_env) && str2double(psf_env)==1;
psr_env = getenv('PROBE_SUPPORT_RADIUS');
if ~isempty(psr_env) && str2double(psr_env)>0; probe_support_radius = str2double(psr_env); else; probe_support_radius = []; end
fprintf('probe_support_tem = %d ; probe_support_radius = %s\n', probe_support_tem, mat2str(probe_support_radius));
Npos_st                   = [inf, inf];     % positions are EXACT (from sim) -> fixed

% --- released / variable probe for noisy, partially-coherent data (TDS + dose) ---
% With >1 probe mode the mutually-incoherent modes must actually be updated together
% (apply_multimodal_update) or the extra modes just sit there. Auto-on when Nprobe>1;
% MULTIMODAL=0/1 overrides. VARIABLE_PROBE=<#modes> turns on orthogonal probe
% relaxation (the probe varies across the scan) — 0/unset = off. These let the recon
% absorb the incoherent TDS background instead of corrupting the object with it.
mm_env = getenv('MULTIMODAL');
if ~isempty(mm_env); multimodal = logical(str2double(mm_env)); else; multimodal = (Nprobe > 1); end
vp_env = getenv('VARIABLE_PROBE');
if ~isempty(vp_env); n_varprobe = max(0, round(str2double(vp_env))); else; n_varprobe = 0; end
fprintf('apply_multimodal_update = %d ; variable_probe_modes = %d\n', multimodal, n_varprobe);
% Depth (multilayer) regularizer: regulation_multilayers.m is a missing-cone low-pass in
% kz (W = 1-atan((R*|kz|/k_xy)^2)/(pi/2)) -> it BLURS DEPTH. Here the depth profile IS the
% measurement -- the sign of P_z is read from the stacking order of the Pb / Ti-O / O columns
% along the beam -- so unlike the parent driver this one defaults it OFF. The stabilising job
% it used to do is covered by the fixed (true) probe and the very high scan overlap.
% REGLAYER=0.05 is the fallback if a deep run diverges; anything larger destroys the signal.
rl_env = getenv('REGLAYER');
if ~isempty(rl_env); reglayer = [str2double(rl_env), str2double(rl_env)]; else; reglayer = [0, 0]; end
fprintf('regularize_layers (per engine) = [%g %g]\n', reglayer(1), reglayer(2));
if any(reglayer > 0.1)
    warning(['REGLAYER = %g will low-pass the depth axis, which is the quantity this ' ...
             'experiment measures. Use 0 (or <=0.05 only to stabilise a divergent run).'], reglayer(1));
end
Np_presolve               = [2*floor(Ndpx/4), Ndpx]; % half-Ndp, forced EVEN (the GPU engine
%   uses even FFT sizes; Ndpx=1426 -> round(/2)=713 is ODD -> 713/712 size clash). 356->178.
Niter_save_results        = [50,  50];
Niter_save_exit_wave      = [200, 200];
strcustom0                = sprintf('hollow_hsa%.2f', hsa);

% LSQ step size: deep multilayer recons are ill-conditioned and diverge to NaN at
% the proven-baseline 0.5 (esp. when the probe update starts). Yu's 30-layer hollow
% script used 0.1. Default 0.1; override via BETA_LSQ env (drop to 0.05 if 82 layers
% still NaNs).
blsq_env = getenv('BETA_LSQ');
if ~isempty(blsq_env); beta_LSQ_val = str2double(blsq_env); else; beta_LSQ_val = 0.1; end
fprintf('beta_LSQ = %.3g\n', beta_LSQ_val);

% -------- chained-restart (continue a long recon across the 2-day wall) ----------
% RESTART_DIR set  -> init the object from the previous chain segment's native
%   *_recons.h5 (io.load_ptycho_recons -> S.object) and run ONLY the full-resolution
%   engine for NITER more iterations (presolve is skipped — the object is full-res).
%   NB: PtychoShelves writes _recons.h5 only when a segment COMPLETES; its intermediate
%   Niter*.mat keep the object under outputs.object, which prepare_initial_object does
%   NOT read (that was the original chain bug). So segments must FINISH (budget
%   SEG_ITERS < wall) and the chain launcher chains with afterok.
% RESTART_DIR unset -> fresh run (coarse presolve + full).
% SAVE_EVERY (default 50): how often to write intermediate Niter*.mat (for monitoring /
%   a partial object); NOT used for resume — that's the _recons.h5.
se_env = getenv('SAVE_EVERY');
if ~isempty(se_env); save_every = round(str2double(se_env)); else; save_every = 50; end
Niter_save_results = [save_every, save_every];
restart_dir = getenv('RESTART_DIR');
% KNOWN_OBJECT: start from the TRUE object (fusion/known_object.py writes it). Diagnoses whether
% a failure to recover depth is the optimiser (the truth is a far better fit that a random start
% never finds) or the forward model (the engine cannot reproduce the data even from the truth).
% It reuses the restart path: collapse to the full-resolution engine, load the object from file.
known_obj   = getenv('KNOWN_OBJECT');
if ~isempty(restart_dir) && ~isempty(known_obj)
    error('set RESTART_DIR or KNOWN_OBJECT, not both');
end
do_restart  = ~isempty(restart_dir) || ~isempty(known_obj);
restart_obj = '';
if ~isempty(known_obj)
    if ~exist(known_obj, 'file')
        error('KNOWN_OBJECT not found: %s', known_obj);
    end
    restart_obj = known_obj;
    fprintf('KNOWN OBJECT: starting from the TRUE object %s (full-res engine only)\n', known_obj);
elseif do_restart
    dd = dir(fullfile(restart_dir, '**', '*_recons.h5'));   % PtychoShelves native recon file
    if isempty(dd)
        error('RESTART_DIR has no *_recons.h5 (did the previous segment complete?): %s', restart_dir);
    end
    [~, imax] = max([dd.datenum]);                          % most recently written
    restart_obj = fullfile(dd(imax).folder, dd(imax).name);
    fprintf('RESTART: continuing from %s (full-res engine only)\n', restart_obj);
end
if do_restart
    % collapse the two-engine schedule to the full-resolution (2nd) engine
    Niter = Niter(end); grouping = grouping(end); Nst_probe = Nst_probe(end);
    Npos_st = Npos_st(end); reglayer = reglayer(end); Np_presolve = Np_presolve(end);
    Niter_save_results = Niter_save_results(end); Niter_save_exit_wave = Niter_save_exit_wave(end);
else
    fprintf('FRESH run (coarse presolve + full)\n');
end

os_env = getenv('OBJECT_START');
if ~isempty(os_env); obj_start = str2double(os_env); else; obj_start = 1; end
if isempty(known_obj); obj_src = 'random'; else; obj_src = 'TRUE (known object)'; end
fprintf('initial object = %s ; object_change_start = %g%s\n', obj_src, obj_start, ...
        repmat(' (FROZEN: evaluating the object, not updating it)', 1, double(isinf(obj_start))));

% ---- the settings this experiment lives or dies on, echoed together ------------
fprintf(['\n== fusion recon preflight ==\n' ...
         '   hollow semi-angle    = %.2f alpha\n' ...
         '   Nlayers x delta_z    = %d x %.3f A  (over %.2f A)\n' ...
         '   regularize_layers    = [%g %g]   <- MUST be 0: depth is the measurement\n' ...
         '   probe_modes          = %d        <- 1 for a coherent synthetic probe\n' ...
         '   beta_LSQ             = %.3g\n' ...
         '   custom_data_flip     = [0 0 1]   (transpose; set at the engine below)\n' ...
         '   probe / positions    = FIXED (known from the simulation)\n\n'], ...
        hsa, Nlayers, delta_z, thick, reglayer(1), reglayer(end), Nprobe, beta_LSQ_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% p struct %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear p
p = struct();
p.   verbose_level = 2;
p.   use_display   = false;
p.   scan_number   = idx_scan;

p.   z            = 1 / d_alpha;
p.   asize        = [Ndpx, Ndpx];
p.   ctr          = [fix(Ndpx/2)+1, fix(Ndpx/2)+1];
p.   prop_regime  = 'farfield';
p.   energy       = HT;
p.   electron     = true;
p.   affine_matrix = [1, 0; 0, 1];

p.   src_metadata   = 'none';
p.   queue.lockfile = false;

p.   detector.name           = 'empad';
p.   detector.check_2_detpos = [];
p.   detector.data_prefix    = '';
p.   detector.binning        = false;
p.   detector.upsampling     = false;
p.   detector.burst_frames   = 1;
p.   prepare.data_preparator        = 'matlab_aps';
p.   prepare.auto_prepare_data      = true;
p.   prepare.force_preparation_data = true;
p.   prepare.store_prepared_data    = false;
p.   prepare.prepare_data_function  = '';
p.   prepare.auto_center_data       = false;

p.   src_positions   = 'hdf5_pos';
p.   positions_file  = '';
p.   scan.type       = 'custom';
p.   scan.roi_label  = [];
p.   scan.format     = scan_string_format;
p.   scan.custom_positions_source = '';
p.   scan.custom_params           = [];

p.   prefix              = '';
p.   suffix              = 'ML_recon';
p.   scan_string_format  = scan_string_format;
p.   base_path           = base_path;
p.   specfile            = '';
p.   ptycho_matlab_path  = '';
p.   cSAXS_matlab_path   = '';
p.   raw_data_path{1}    = '';
p.   prepare_data_path   = '';
p.   prepare_data_filename = [];
p.   save_path{1}        = '';
p.   io.default_mask_file = '';
p.   io.default_mask_type = 'binary';
p.   io.file_compression  = 0;
p.   io.data_compression  = 3;
p.   io.load_prep_pos     = false;

if do_restart
    p.   model_object  = false;                 % continue a previous chain segment
    p.   initial_iterate_object_file{1} = restart_obj;
else
    p.   model_object  = true;
    p.   model.object_type = 'rand';
    p.   initial_iterate_object_file{1} = '';
end
p.   model_probe   = false;
% TEM aperture constraint is a TOP-LEVEL p option: the solver checks par.p.probe_support_tem and
% load_from_p builds the mask from abs(fft2(probe_initial)). eng.* is nested in p.engines{} and
% would NOT be seen here, so set it on p directly (see the PROBE_SUPPORT_FFT env block above).
p.   probe_support_tem = probe_support_tem;
p.   model.probe_is_focused            = true;
p.   model.probe_central_stop          = true;
p.   model.probe_diameter              = 170e-6;
p.   model.probe_central_stop_diameter = 50e-6;
p.   model.probe_zone_plate_diameter   = 170e-6;
p.   model.probe_outer_zone_width      = [];
p.   model.probe_propagation_dist      = 3e-3;
p.   model.probe_focal_length          = 51e-3;
p.   model.probe_upsample              = 10;
p.   initial_probe_file     = fullfile(p.base_path, sprintf(p.scan.format, p.scan_number), 'probe_initial.mat');
p.   probe_file_propagation = 0.0e-3;
p.   share_probe   = 0;
p.   share_object  = 0;

p.   probe_modes    = Nprobe;
p.   object_modes   = 1;
p.   mode_start_pow = 0.02;
p.   mode_start     = 'herm';
p.   ortho_probes   = true;
p.   object_regular = 0;

p.   plot.prepared_data             = false;
p.   save.external                  = false;
p.   save.store_images              = false;
p.   save.store_images_intermediate = false;
p.   save.store_images_ids          = 1:4;
p.   save.store_images_format       = 'png';
p.   save.store_images_dpi          = 150;
p.   save.exclude                   = {'fmag','fmask','illum_sum'};
p.   save.save_reconstructions_intermediate = true;
p.   save.save_reconstructions      = true;
p.   save.output_file               = 'h5';

%%%%%%%%%%%%%%%%%%%%%%%%%%%% engines (multislice) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ieng = 1:length(Niter)
    eng = struct();
    eng. name           = 'GPU';
    eng. use_gpu        = true;
    eng. keep_on_gpu    = false;
    eng. compress_data  = false;
    eng. gpu_id         = [];
    eng. check_gpu_load = true;

    eng. number_iterations   = Niter(ieng);
    eng. asize_presolve      = [Np_presolve(ieng), Np_presolve(ieng)];
    eng. method              = 'MLs';
    eng. opt_errmetric       = 'L1';
    eng. grouping            = grouping(ieng);
    eng. probe_modes         = p.probe_modes;
    eng. object_change_start = obj_start;   % OBJECT_START env; inf = frozen (evaluate only)
    eng. probe_change_start  = Nst_probe(ieng);

    eng. reg_mu                       = 0;
    eng. delta                        = 0;
    eng. positivity_constraint_object = 0;
    eng. apply_multimodal_update      = multimodal;
    eng. probe_backpropagate          = 0;
    eng. probe_support_radius         = probe_support_radius;
    eng. probe_support_fft            = [];                   % MASK (built by the engine), never a scalar
    eng. probe_support_tem            = probe_support_tem;    % TEM aperture constraint (mask from probe_initial)

    eng. beta_object = 1;
    eng. beta_probe  = 1;
    eng. delta_p     = 0.1;
    eng. momentum    = 0;
    eng. beta_LSQ    = beta_LSQ_val;
    eng. accelerated_gradients_start = inf;

    eng. apply_subpix_shift             = true;
    eng. probe_position_search          = Npos_st(ieng);   % inf => fixed positions
    eng. probe_geometry_model           = {};
    eng. probe_position_error_max       = inf;
    eng. apply_relaxed_position_constraint = false;

    % multilayer extension: Nlayers slices of delta_z (Å -> m)
    eng. delta_z           = delta_z * ones(Nlayers,1) * 1e-10;
    eng. regularize_layers = reglayer(ieng);
    eng. preshift_ML_probe = false;

    eng. background       = 0;
    eng. background_width = inf;
    eng. clean_residua    = false;

    eng. probe_fourier_shift_search = inf;
    eng. estimate_NF_distance       = inf;
    eng. detector_rotation_search   = inf;
    eng. detector_scale_search      = inf;
    eng. variable_probe             = (n_varprobe > 0);
    eng. variable_probe_modes       = max(1, n_varprobe);
    eng. variable_probe_smooth      = 0;
    eng. variable_intensity         = false;

    eng. get_fsc_score     = false;
    eng. mirror_objects    = false;
    eng. auto_center_data  = false;
    eng. auto_center_probe = false;
    eng. custom_data_flip  = [0,0,1];   % [fliplr flipud TRANSPOSE]: the orientation
                                        % sweep (run_orientation_sweep) proved the
                                        % synthetic detector needs a transpose to bind
                                        % correctly to the scan. Only this resolves atoms.
    eng. apply_tilted_plane_correction = '';

    eng. plot_results_every           = inf;
    eng. save_results_every           = Niter_save_results(ieng);
    eng. save_results_every_exit_wave = Niter_save_exit_wave(ieng);
    eng. save_phase_image             = true;
    eng. save_probe_mag               = true;

    resultDir = strcat(p.base_path, sprintf(p.scan.format, p.scan_number));
    strcustom = strcat('_Npbst', num2str(Nst_probe(ieng)), '_', strcustom0);
    eng.fout  = generateResultDir(eng, resultDir, ...
        strcat(strcustom, '_Ndp', num2str(Np_presolve(ieng)), '_step', num2str(ieng,'%02d')));
    disp(eng.fout);
    mkdir(eng.fout);
    copyfile(strcat(mfilename('fullpath'), '.m'), eng.fout);

    [p, ~] = core.append_engine(p, eng);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out = core.ptycho_recons(p);
fprintf('Multislice HOLLOW reconstruction finished (HSA = %.2f alpha).\n', hsa);
