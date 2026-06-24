% MULTISLICE (depth-resolved) reconstruction of abTEM-synthetic 4D-STEM data
% through Yu's MultiHollowPtycho engine, hollow-angle = 0. The sample is ~70 Å of
% PbTiO3 at 100 mrad -> strongly dynamical, so single-slice ptychography fails;
% this reconstructs Nlayers depth slices along the beam (and can section a vortex).
%
% Geometry/calibration is read from sim_meta.mat so it never drifts from the sim.
% Uses the SAME sim_out/01/ data as the single-slice check (no re-simulation).
%
% Run from:  ptychoshelves-clean/ptycho/   (MATLAB pwd must be this directory)

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
    base_path = '../sim_out/';
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%% hollow inner mask (= ZERO) %%%%%%%%%%%%%%%%%%%%%%%%%
global mask1
mask1 = ones(Ndpx, Ndpx, 'single');

%%%%%%%%%%%%%%%%%%%%%%%%%%%% multislice reconstruction parameters %%%%%%%%%%%%%%%
Nprobe    = 8;
DELTA_Z_A = 10;                         % target recon slice thickness [Å]
Nlayers   = max(1, round(thick / DELTA_Z_A));   % ~7 layers through ~70 Å
delta_z   = thick / Nlayers;            % actual layer spacing [Å]
fprintf('multislice: Nlayers=%d, delta_z=%.1f A\n', Nlayers, delta_z);

% two-engine schedule (coarse presolve -> full), modelled on the proven baseline
grouping                  = [64,  32];
Niter                     = [200, 200];
Nst_probe                 = [20,  10];
Npos_st                   = [inf, inf];     % positions are EXACT (from sim) -> fixed
reglayer                  = [1,   0.5];     % symmetrise info between layers
Np_presolve               = [round(Ndpx/2), Ndpx];   % 178 -> 356
Niter_save_results        = [50,  50];
Niter_save_exit_wave      = [200, 200];
strcustom0                = 'synthetic_ML';

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

p.   model_object  = true;
p.   model.object_type = 'rand';
p.   initial_iterate_object_file{1} = '';
p.   model_probe   = false;
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
    eng. object_change_start = 1;
    eng. probe_change_start  = Nst_probe(ieng);

    eng. reg_mu                       = 0;
    eng. delta                        = 0;
    eng. positivity_constraint_object = 0;
    eng. apply_multimodal_update      = false;
    eng. probe_backpropagate          = 0;
    eng. probe_support_radius         = [];
    eng. probe_support_fft            = false;
    eng. probe_support_tem            = false;

    eng. beta_object = 1;
    eng. beta_probe  = 1;
    eng. delta_p     = 0.1;
    eng. momentum    = 0;
    eng. beta_LSQ    = 0.5;
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
    eng. variable_probe             = false;
    eng. variable_probe_modes       = 1;
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
fprintf('Multislice synthetic reconstruction finished.\n');
