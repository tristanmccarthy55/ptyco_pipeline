% ORIENTATION (8-DOF) TEST: reconstruct the thin 'F' phantom 8 times, once per
% detector orientation (the dihedral group D4 = all [fliplr, flipud, transpose]
% combinations, set via eng.custom_data_flip). The phantom is thin, so single-slice
% is exact and the result is NOT confounded by thickness. Exactly one orientation
% renders the F upright -> that is the correct custom_data_flip to lock into the
% real multislice recon.
%
% Run from:  ptychoshelves-clean/ptycho/
% Data:      <SIM_BASE>/01/  (point SIM_BASE at the PHANTOM sim output)

%%
clear variables
addpath(pwd)
addpath(fullfile(pwd,'utils'))
addpath(fullfile(pwd,'utils_EM'))
addpath(core.find_base_package)

scan_string_format = '%02d';
idx_scan  = 1;
base_path = getenv('SIM_BASE');
if isempty(base_path); base_path = '../sim_phantom/'; end
scan_dir  = fullfile(base_path, sprintf(scan_string_format, idx_scan));

meta_file = fullfile(scan_dir, 'sim_meta.mat');
if ~exist(meta_file,'file')
    error('sim_meta.mat not found at %s — run the phantom sim first.', meta_file);
end
S = load(meta_file); meta = S.meta;
Ndpx    = double(meta.Ndpx);
d_alpha = double(meta.d_alpha_rad);
HT      = double(meta.energy_kev);
fprintf('sim_meta: Ndpx=%d, d_alpha=%.4f mrad, HT=%.0f keV\n', Ndpx, d_alpha*1e3, HT);

global mask1
mask1 = ones(Ndpx, Ndpx, 'single');

% the 8 detector orientations: [fliplr, flipud, transpose]
flips = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1];
Niter = 30;     % short: a correct binding shows the F lattice within ~20 iters

for ii = 1:size(flips,1)
    flip = flips(ii,:);
    tag  = sprintf('orient_%d%d%d', flip(1), flip(2), flip(3));
    fprintf('\n================ %s  [fliplr flipud transpose] = [%d %d %d] ================\n', ...
            tag, flip(1), flip(2), flip(3));

    clear p
    p = struct();
    p.   verbose_level = 1;
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
    p.   probe_modes    = 8;
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

    eng = struct();
    eng. name           = 'GPU';
    eng. use_gpu        = true;
    eng. keep_on_gpu    = false;
    eng. compress_data  = false;
    eng. gpu_id         = [];
    eng. check_gpu_load = true;
    eng. number_iterations   = Niter;
    eng. asize_presolve      = [Ndpx, Ndpx];
    eng. method              = 'MLs';
    eng. opt_errmetric       = 'L1';
    eng. grouping            = 64;
    eng. probe_modes         = p.probe_modes;
    eng. object_change_start = 1;
    eng. probe_change_start  = 10;
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
    eng. probe_position_search          = inf;     % fixed exact positions
    eng. probe_geometry_model           = {};
    eng. probe_position_error_max       = inf;
    eng. apply_relaxed_position_constraint = false;
    eng. delta_z           = [];                    % single slice (thin phantom)
    eng. regularize_layers = 0;
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
    eng. custom_data_flip  = flip;                  % <<< the orientation under test
    eng. apply_tilted_plane_correction = '';
    eng. plot_results_every           = inf;
    eng. save_results_every           = Niter;
    eng. save_results_every_exit_wave = Niter;
    eng. save_phase_image             = true;
    eng. save_probe_mag               = true;

    resultDir = strcat(p.base_path, sprintf(p.scan.format, p.scan_number));
    eng.fout  = generateResultDir(eng, resultDir, strcat('_', tag));
    disp(eng.fout); mkdir(eng.fout);

    [p, ~] = core.append_engine(p, eng);
    try
        out = core.ptycho_recons(p);
    catch ME
        fprintf('  %s FAILED: %s\n', tag, ME.message);
    end
end

fprintf('\nOrientation sweep finished. Compare the 8 O_phase_roi images;\n');
fprintf('the correct custom_data_flip renders the F upright.\n');
