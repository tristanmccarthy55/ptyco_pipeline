% main demo script for multislice hollow ptychograpy
% By Yu Lei @ University of Warwick, 7/18/2025
% Reference:
%    Yu Lei and Peng Wang, Hollow multi-slice electron ptychography for simultaneous 3D structural imaging and EELS in 4D-STEM, arXiv:2506.22352 (https://arxiv.org/abs/2506.22352)
%%

% ### prepare your data first by using ./ptycho/util_EM/prepare_data_electron_exp.m

%%
clear variables
addpath(pwd)
addpath(fullfile(pwd,'utils'))
addpath(fullfile(pwd,'utils_EM'))
addpath(core.find_base_package)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I/O parameters
scan_string_format = '%02d';

% data param
%base_path= fullfile(pwd,'exampleData');
base_path= ['./Hollow_Data_exp_30nm/'] ; % change to your data path
if isempty(base_path)
    error('change the path of the data directory and run again');
end
range=[2 12 26 40 50];
for i=range

    idx_scan=i; % scan area index
    d_alpha = 21.4/26; % mrad, angle of each pixel in diffraction pattern.
    d_alpha=d_alpha*1e-3; %rad,
    HT = 300; % energy in keV
    Ndpx=256;  % size of diffraction

    %recon param
    Nprobe = 8; % # of probe modes
    thick=300; % sample thickness in angstrom
    Nlayers = 30; % # of slices for multi-slice, 1 for single-slice
    % delta_z = [10] ; % thickness of each slice in Angstrom
    % strcustom0 = 'z1nm_reg1'; % custom defined for output dir
    %params for different engines
    grouping = [80,20];
    Niter = [200, 200];
    Nst_probe= [20,10]; % start probe update
    Npos_st=[50,50];  % Iteration starting position search
    reglayer = [1,0.5]; %regularize_layers
    Np_presolve=[128,256];
    Niter_save_results=[20,20];
    Niter_save_results_exit_wave=[200,200];

    %load mask
    if i>1
    global mask1;
    mask = load([base_path,num2str(idx_scan,'%02d'),'/mask.mat']).mask;
    mask1 = single(mask);
    end
    %
    delta_z=thick/Nlayers;
    strcustom0='test';
    %% General
    clear p
    p = struct();
    p.   verbose_level = 2;                            % verbosity for standard output (0-1 for loops, 2-3 for testing and adjustments, >= 4 for debugging)
    p.   use_display = false; %true;                                      % global switch for display, if [] then true for verbose > 1
    p.   scan_number = [idx_scan];                                    % Multiple scan numbers for shared scans

    % Geometry
    p.   z = 1/d_alpha ;  % 1/delta_angle in rad^-1, for electron ptychography
    p.   asize = [Ndpx,Ndpx];                                     % Diffr. patt. array size
    p.   ctr = [fix(Ndpx/2)+1, fix(Ndpx/2)+1];                                       % Diffr. patt. center coordinates (y,x) (empty means middle of the array); e.g. [100 207;100+20 207+10];
    p.   prop_regime = 'farfield';                              % propagation regime: nearfield, farfield (default), !! nearfield is supported only by GPU engines
    p.   energy = HT ;                                      % Energy (in keV), leave empty to use spec entry mokev
    p.   electron = true;  % for electron ptychography, added by ZC

    p.   affine_matrix = [1 , 0; 0, 1] ; % Applies affine transformation (e.g. rotation, stretching) to the positions (ignore by = []). Convention [yn;xn] = M*[y;x]. For flOMNI we found in June 2019: = [1 , 0.0003583 ; 5.811e-05 , 1 ]; for OMNY we found in October 2018: = [1 0;tan(0.4*pi/180) 1]; laMNI in June 2018  [1,0.0154;-0.0017,1.01]; laMNI in August [1.01 0.0031; -0.0018 1.00]
    % Scan meta data
    p.   src_metadata = 'none';                                 % source of the meta data, following options are supported: 'spec', 'none' , 'artificial' - or add new to +scan/+meta/
    p.   queue.lockfile = false;                                % If true writes a lock file, if lock file exists skips recontruction

    % Data preparation
    p.   detector.name = 'empad';                           % see +detectors/ folder
    p.   detector.check_2_detpos = [];                          % = []; (ignores)   = 270; compares to dettrx to see if p.ctr should be reversed (for OMNY shared scans 1221122), make equal to the middle point of dettrx between the 2 detector positions
    p.   detector.data_prefix = '';                             % Default using current eaccount e.g. e14169_1_
    p.   detector.binning = false;                              % = true to perform 2x2 binning of detector pixels, for binning = N do 2^Nx2^N binning
    p.   detector.upsampling = false;                           % upsample the measured data by 2^data_upsampling, (transposed operator to the binning), it can be used for superresolution in nearfield ptychography or to account for undersampling in a far-field dataset
    p.   detector.burst_frames = 1;                             % number of frames collected per scan position

    p.   prepare.data_preparator = 'matlab_aps';                    % data preparator; 'python' or 'matlab' or 'matlab_aps'
    p.   prepare.auto_prepare_data = true;                      % if true: prepare dataset from raw measurements if the prepared data does not exist
    p.   prepare.force_preparation_data = true;                 % Prepare dataset even if it exists, it will overwrite the file % Default: @prepare_data_2d
    p.   prepare.store_prepared_data = false;                    % store the loaded data to h5 even for non-external engines (i.e. other than c_solver)
    p.   prepare.prepare_data_function = '';                    % (used only if data should be prepared) custom data preparation function handle;
    p.   prepare.auto_center_data = false;                      % if matlab data preparator is used, try to automatically center the diffraction pattern to keep center of mass in center of diffraction

    % Scan positions
    p.   src_positions = 'hdf5_pos';                           % 'spec', 'orchestra', 'load_from_file', 'matlab_pos' (scan params are defined below) or add new position loaders to +scan/+positions/
    p.   positions_file = [''];    %Filename pattern for position files, Example: ['../../specES1/scan_positions/scan_%05d.dat']; (the scan number will be automatically filled in)

    % scan parameters for option src_positions = 'matlab_pos';
    p.   scan.type = 'custom';                                  % {'round', 'raster', 'round_roi', 'custom'}
    p.   scan.roi_label = [];                            % For APS data
    p.   scan.format = scan_string_format;                      % For APS data format for scan directory generation
    p.   scan.custom_positions_source = '';                     % custom: a string name of a function that defines the positions; also accepts mat file with entry 'pos', see +scans/+positions/+mat_pos.m
    p.   scan.custom_params = [];                               % custom: the parameters to feed to the custom position function.

    % I/O
    p.   prefix = '';                                              % For automatic output filenames. If empty: scan number
    p.   suffix = strcat('ML_recon');              % Optional suffix for reconstruction
    p.   scan_string_format = scan_string_format;                  % format for scan string generation, it is used e.g for plotting and data saving
    p.   base_path = base_path;     % base path : used for automatic generation of other paths

    p.   specfile = '';                                         % Name of spec file to get motor positions and check end of scan, defaut is p.spec_file == p.base_path;
    p.   ptycho_matlab_path = '';                               % cSAXS ptycho package path
    p.   cSAXS_matlab_path = '';                                % cSAXS base package path
    p.   raw_data_path{1} = '';                                 % Default using compile_x12sa_filename, used only if data should be prepared automatically
    p.   prepare_data_path = '';                                % Default: base_path + 'analysis'. Other example: '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/'; also supports %u to insert the scan number at a later point (e.g. '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/S%.5u')
    p.   prepare_data_filename = [];                            % Leave empty for default file name generation, otherwise use [sprintf('S%05d_data_%03dx%03d',p.scan_number(1), p.asize(1), p.asize(2)) p.prep_data_suffix '.h5'] as default
    p.   save_path{1} = '';                                     % Default: base_path + 'analysis'. Other example: '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/'; also supports %u to insert the scan number at a later point (e.g. '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/S%.5u')
    p.   io.default_mask_file = '';                             % load detector mask defined in this file instead of the mask in the detector packages, (used only if data should be prepared)
    p.   io.default_mask_type = 'binary';                       % (used only if data should be prepared) ['binary', 'indices']. Default: 'binary'
    p.   io.file_compression = 0;                               % reconstruction file compression for HDF5 files; 0 for no compression
    p.   io.data_compression = 3;                               % prepared data file compression for HDF5 files; 0 for no compression
    p.   io.load_prep_pos = false;                              % load positions from prepared data file and ignore positions provided by metadata
    %% Reconstruction
    % Initial iterate object
    p.   model_object = true;                                   % Use model object, if false load it from file
    p.   model.object_type = 'rand';                            % specify how the object shall be created; use 'rand' for a random initial guess; use 'amplitude' for an initial guess based on the prepared data
    p.   initial_iterate_object_file{1} = '';                   %  use this mat-file as initial guess of object, it is possible to use wild characters and pattern filling, example: '../analysis/S%05i/wrap_*_1024x1024_1_recons*'

    % Initial iterate probe
    p.   model_probe = false;                                    % Use model probe, if false load it from file
    p.   model.probe_is_focused = true;                         % Model probe is focused (false: just a pinhole)
    p.   model.probe_central_stop = true;                       % Model central stop
    p.   model.probe_diameter = 170e-6;                         % Model probe pupil diameter
    p.   model.probe_central_stop_diameter = 50e-6;             % Model central stop diameter
    p.   model.probe_zone_plate_diameter = 170e-6;              % Model probe zone plate diameter
    p.   model.probe_outer_zone_width = [];                     % Model probe zone plate outermost zone width (not used if not a focused probe)
    p.   model.probe_propagation_dist = 3e-3;                 % Model probe propagation distance (pinhole <-> sample for unfocused, focal-plane <-> sample for focused)
    p.   model.probe_focal_length = 51e-3;                      % Model probe focal length (used only if model_is_focused is true
    %   AND model_outer_zone_width is empty)
    p.   model.probe_upsample = 10;                             % Model probe upsample factor (for focused probes)

    %Use probe from this mat-file (not used if model_probe is true)
    p.   initial_probe_file = fullfile(p.base_path,sprintf(p.scan.format, p.scan_number),'probe_initial.mat');
    p.   probe_file_propagation = 0.0e-3;                            % Distance for propagating the probe from file in meters, = 0 to ignore

    % Shared scans - Currently working only for sharing probe and object
    p.   share_probe  = 0;                                      % Share probe between scans. Can be either a number/boolean or a list of numbers, specifying the probe index; e.g. [1 2 2] to share the probes between the second and third scan.
    p.   share_object = 0;                                      % Share object between scans. Can be either a number/boolean or a list of numbers, specifying the object index; e.g. [1 2 2] to share the objects between the second and third scan.

    % Modes
    p.   probe_modes  = Nprobe;                                 % Number of coherent modes for probe
    p.   object_modes = 1;                                      % Number of coherent modes for object
    % Mode starting guess
    p.   mode_start_pow = [0.02];                               % Normalized intensity on probe modes > 1. Can be a number (all higher modes equal) or a vector
    p.   mode_start = 'herm';                                   % (for probe) = 'rand', = 'herm' (Hermitian-like base), = 'hermver' (vertical modes only), = 'hermhor' (horizontal modes only)
    p.   ortho_probes = true;                                   % orthogonalize probes after each engine
    p.   object_regular = 0;                                  % should be smaller than 1/8, smooth object amplitude, usefull for many layers

    %% Plot, save and analyze
    p.   plot.prepared_data = false;                         % plot prepared data
    p.   save.external = true;                             % Use a new Matlab session to run save final figures (saves ~6s per reconstruction). Please be aware that this might lead to an accumulation of Matlab sessions if your single reconstruction is very fast.
    p.   save.store_images = false;                              % Write preview images containing the final reconstructions in [p.base_path,'analysis/online/ptycho/'] if p.use_display = 0 then the figures are opened invisible in order to create the nice layout. It writes images in analysis/online/ptycho
    p.   save.store_images_intermediate = false;                % save images to disk after each engine
    p.   save.store_images_ids = 1:4;                           % identifiers  of the figure to be stored, 1=obj. amplitude, 2=obj. phase, 3=probes, 4=errors, 5=probes spectrum, 6=object spectrum
    p.   save.store_images_format = 'png';                      % data type of the stored images jpg or png
    p.   save.store_images_dpi = 150;                           % DPI of the stored bitmap images
    p.   save.exclude = {'fmag', 'fmask', 'illum_sum'};         % exclude variables to reduce the file size on disk
    p.   save.save_reconstructions_intermediate = true;        % save final object and probes after each engine
    p.   save.save_reconstructions = true;                      % save reconstructed object and probe when full reconstruction is finished
    p.   save.output_file = 'h5';                               % data type of reconstruction file; 'h5' or 'mat'

    %% ENGINES
    for ieng=1:length(Niter)
        % --------- GPU engines  -------------   See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
        eng = struct();                        % reset settings for this engine
        eng. name = 'GPU';
        eng. use_gpu = true;                   % if false, run CPU code, but it will get very slow
        eng. keep_on_gpu = false;               % keep data + projections on GPU, false is useful for large data if DM is used
        eng. compress_data = false;             % use automatic online memory compression to limit need of GPU memory
        eng. gpu_id = 1;                      % default GPU id, [] means choosen by matlab
        eng. check_gpu_load = true;            % check available GPU memory before starting GPU engines

        % general
        eng. number_iterations = Niter(ieng);          % number of iterations for selected method
        eng. asize_presolve = [Np_presolve(ieng),Np_presolve(ieng)];

        eng. method = 'MLs';                   % choose GPU solver: DM, ePIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
        eng. opt_errmetric = 'L1';            % optimization likelihood - poisson, L1
        eng. grouping = grouping(ieng);                    % size of processed blocks, larger blocks need more memory but they use GPU more effeciently, !!! grouping == inf means use as large as possible to fit into memory
        % * for hPIE, ePIE, MLs methods smaller blocks lead to faster convergence,
        % * for MLc the convergence is similar
        % * for DM is has no effect on convergence
        eng. probe_modes  = p.probe_modes;                % Number of coherent modes for probe
        eng. object_change_start = 1;          % Start updating object at this iteration number
        eng. probe_change_start = Nst_probe(ieng);           % Start updating probe at this iteration number

        % regularizations
        eng. reg_mu = 0;                       % Regularization (smooting) constant ( reg_mu = 0 for no regularization)
        eng. delta = 0;                        % press values to zero out of the illumination area in th object, usually 1e-2 is enough
        eng. positivity_constraint_object = 0; % enforce weak (relaxed) positivity in object, ie O = O*(1-a)+a*|O|, usually a=1e-2 is already enough. Useful in conbination with OPRP or probe_fourier_shift_search

        eng. apply_multimodal_update = false; % apply all incoherent modes to object, it can cause isses if the modes collect some crap
        eng. probe_backpropagate = 0;         % backpropagation distance the probe mask, 0 == apply in the object plane. Useful for pinhole imaging where the support can be applied  at the pinhole plane
        eng. probe_support_radius = [];       % Normalized radius of circular support, = 1 for radius touching the window
        eng. probe_support_fft = false;       % assume that there is not illumination intensity out of the central FZP cone and enforce this contraint. Useful for imaging with focusing optics. Helps to remove issues from the gaps between detector modules.
        eng. probe_support_tem = false;        % assume a binary mask aperture for TEM, generated from initial probe, by Zhen Chen

        % basic recontruction parameters
        % PIE / ML methods                    % See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
        eng. beta_object = 1;                 % object step size, larger == faster convergence, smaller == more robust, should not exceed 1
        eng. beta_probe = 1;                  % probe step size, larger == faster convergence, smaller == more robust, should not exceed 1
        eng. delta_p = 0.1;                   % LSQ dumping constant, 0 == no preconditioner, 0.1 is usually safe, Preconditioner accelerates convergence and ML methods become approximations of the second order solvers
        eng. momentum = 0;                    % add momentum acceleration term to the MLc method, useful if the probe guess is very poor or for acceleration of multilayer solver, but it is quite computationally expensive to be used in conventional ptycho without any refinement.
        eng. beta_LSQ = 0.1;                                      % The momentum method works usually well even with the accelerated_gradients option.  eng.momentum = multiplication gain for velocity, eng.momentum == 0 -> no acceleration, eng.momentum == 0.5 is a good value
        % momentum is enabled only when par.Niter < par.accelerated_gradients_start;
        eng. accelerated_gradients_start = inf; % iteration number from which the Nesterov gradient acceleration should be applied, this option is supported only for MLc method. It is very computationally cheap way of convergence acceleration.

        % ADVANCED OPTIONS                     See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
        % position refinement
        eng. apply_subpix_shift = true;       % apply FFT-based subpixel shift, it is automatically allowed for position refinement
        eng. probe_position_search = Npos_st(ieng);      % iteration number from which the engine will reconstruct probe positions, from iteration == probe_position_search, assume they have to match geometry model with error less than probe_position_error_max

        eng. probe_geometry_model = {};  % list of free parameters in the geometry model, choose from: {'scale', 'asymmetry', 'rotation', 'shear'}
        eng. probe_position_error_max = inf; % in meters, maximal expected random position errors, probe prositions are confined in a circle with radius defined by probe_position_error_max and with center defined by original positions scaled by probe_geometry_model
        eng. apply_relaxed_position_constraint = false;  % added by YJ
        % multilayer extension
        eng. delta_z = delta_z*ones(Nlayers,1)*1e-10;                     % if not empty, use multilayer ptycho extension , see ML_MS code for example of use, [] == common single layer ptychography , note that delta_z provides only relative propagation distance from the previous layer, ie delta_z can be either positive or negative. If preshift_ML_probe == false, the first layer is defined by position of initial probe plane. It is useful to use eng.momentum for convergence acceleration
        eng. regularize_layers = reglayer(ieng);            % multilayer extension: 0<R<<1 -> apply regularization on the reconstructed object layers, 0 == no regularization, 0.01 == weak regularization that will slowly symmetrize information content between layers
        eng. preshift_ML_probe = false;         % multilayer extension: if true, assume that the provided probe is reconstructed in center of the sample and the layers are centered around this position
        % other extensions
        eng. background = 0;               % average background scattering level, for OMNI values around 0.3 for 100ms, for flOMNI <0.1 per 100ms exposure, see for more details: Odstrcil, M., et al., Optics letters 40.23 (2015): 5574-5577.
        eng. background_width = inf;           % width of the background function in pixels,  inf == flat background, background function is then convolved with the average diffraction pattern in order to account for beam diversion
        eng. clean_residua = false;            % remove phase residua from reconstruction by iterative unwrapping, it will result in low spatial freq. artefacts -> object can be used as an residua-free initial guess for netx engine

        % wavefront & camera geometry refinement     See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
        eng. probe_fourier_shift_search = inf; % iteration number from which the engine will: refine farfield position of the beam (ie angle) from iteration == probe_fourier_shift_search
        eng. estimate_NF_distance = inf;       % iteration number from which the engine will: try to estimate the nearfield propagation distance using gradient descent optimization
        eng. detector_rotation_search = inf;   % iteration number from which the engine will: search for optimal detector rotation, preferably use with option mirror_scan = true , rotation of the detector axis with respect to the sample axis, similar as rotation option in the position refinement geometry model but works also for 0/180deg rotation shared scans
        eng. detector_scale_search = inf;      % iteration number from which the engine will: refine pixel scale of the detector, can be used to refine propagation distance in ptycho
        eng. variable_probe = true;           % Use SVD to account for variable illumination during a single (coupled) scan, see for more details:  Odstrcil, M. et al. Optics express 24.8 (2016): 8360-8369.
        eng. variable_probe_modes = 1;         % OPRP settings , number of SVD modes using to describe the probe evolution.
        eng. variable_probe_smooth = 0;        % OPRP settings , enforce of smooth evolution of the OPRP modes -> N is order of polynomial fit used for smoothing, 0 == do not apply any smoothing. Smoothing is useful if only a smooth drift is assumed during the ptycho acquisition
        eng. variable_intensity = false;       % account to changes in probe intensity

        % extra analysis
        eng. get_fsc_score = false;            % measure evolution of the Fourier ring correlation during convergence
        eng. mirror_objects = false;           % mirror objects, useful for 0/180deg scan sharing -> geometry refinement for tomography, works only if 2 scans are provided

        % custom data adjustments, useful for offaxis ptychography
        eng.auto_center_data = false;           % autoestimate the center of mass from data and shift the diffraction patterns so that the average center of mass corresponds to center of mass of the provided probe
        eng.auto_center_probe = false;          % center the probe position in real space before reconstruction is started
        eng.custom_data_flip = [0,0,0];         % apply custom flip of the data [fliplr, flipud, transpose]  - can be used for quick testing of reconstruction with various flips or for reflection ptychography
        eng.apply_tilted_plane_correction = ''; % if any(p.sample_rotation_angles([1,2]) ~= 0),  this option will apply tilted plane correction. (a) 'diffraction' apply correction into the data, note that it is valid only for "low NA" illumination  Gardner, D. et al., Optics express 20.17 (2012): 19050-19059. (b) 'propagation' - use tilted plane propagation, (c) '' - will not apply any correction

        eng.plot_results_every = inf;
        eng.save_results_every = Niter_save_results(ieng);
        eng.save_results_every_exit_wave = Niter_save_results_exit_wave(ieng);
        eng.save_phase_image = true;
        eng.save_probe_mag = true;

        resultDir = strcat(p.base_path,sprintf(p.scan.format, p.scan_number));
        strcustom=strcat('_Npbst',num2str(Nst_probe(ieng)),'_',strcustom0);
        eng.fout =  generateResultDir(eng, resultDir,strcat(strcustom,'_Ndp',num2str(Np_presolve(ieng)),'_step',num2str(ieng,'%02d')));
        disp(eng.fout)

        mkdir(eng.fout);

        copyfile(strcat(mfilename('fullpath'),'.m'),eng.fout);

        %% add engine
        [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
    end

    %% Run the reconstruction
    out = core.ptycho_recons(p);
end

%%
%{
%%
  Academic License Agreement
================================

 Introduction 
 •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").

 Terms and Conditions of the LICENSE
 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
       hereinafter set out and until termination of this license as set forth below.
 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
       in the commercial use, application or exploitation of works similar to the PROGRAM.
 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
       another computing language:
       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
       Scherrer Institut, Switzerland."

    Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
 P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
   (doi: 10.1126/science.1158573),
 for mixed coherent modes:
 P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
 for LSQ-ML method 
 M. Odstrcil, A. Menzel, M. Guizar-Sicairos, Iterative least-squares solver for generalized maximum-likelihood ptychography, Opt. Express, in press (2018). (doi: ).
 for OPRP method 
 M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  Ptychographic coherent diffractive imaging with orthogonal probe relaxation, 
   Opt. Express 24, 8360 (2016). (doi: 10.1364/OE.24.008360).
 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
       in particular of patent with international application number PCT/GB2005/001464 and US9401042B2. The LICENSOR explicitly declares not to indemnify the users of the software 
       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
 10. This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
       the courts of Zürich, Switzerland. 

%}
