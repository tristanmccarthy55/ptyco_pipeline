% INITIALIZE generate list of default parameters 
% [param] = initialize 
% 
%
% returns: 
% ++ param       structure containing parameters for the engines 

% Academic License Agreement
% 
% Source Code
% 
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
%   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
% 
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
% 
% Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% for LSQ-ML method 
% M. Odstrcil, A. Menzel, M.G. Sicairos,  Iterative least-squares solver for generalized maximum-likelihood ptychography, Optics Express, 2018
% for OPRP method 
%  M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  "Ptychographic coherent diffractive imaging with orthogonal probe relaxation." Optics express 24.8 (2016): 8360-8369
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Zürich, Switzerland. 
% 
%   


function [param] = get_defaults

    %%%%%%%%%%%%%% GPU SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%
    param.use_gpu = true;        % use GPU if possible 
    param.keep_on_gpu = true;    % keep the data all the time on GPU
    param.compress_data = true;  % apply online compress on the GPU data 
    param.gpu_id = []; % default GPU id, [] means choosen by matlab
    param.check_gpu_load = true;
    
    %% basic recontruction parameters 
    %% PIE 
    param.beta_object = 1;
    param.beta_probe = 1;  % step size, faster convergence , more instable ?? 
    %% DM
    param.pfft_relaxation = 0.1; 
    param.probe_inertia = 0.3; % add inertia to the probe reconstruction to avoid oscilations 
    %% general 
    param.share_probe = true;
    param.share_object = false;
    param.delta = 0;  % press values to zero out of the probe area !!  illim < max*delta is removed 
    param.relax_noise = 0.0;  % relaxation for noise, lower => slower convergence, more robust 
    param.positivity_constraint_object = 0; % enforce weak positivity in object 
    param.Nmodes = 1;  %  number of multi apertures , always better to start wih one !! 
    param.probe_modes = 1; % number of probes 
    param.object_modes = 1;  %  number of multi apertures , always better to start wih one !! 
    param.probe_change_start = 1;  % iteration when the probe reconstruction is started
    param.object_change_start = 1;% iteration when the object reconstruction is started
    param.apply_relaxed_position_constraint = false; % by ZC for geometry constrain
       
    param.number_iterations = 300 ; 
    param.grouping = inf;
    param.method = 'MLs';
    param.likelihood = 'L1' ; % l1 or poisson,   - choose which likelihood should be used for solver, poisson is suported only for PIE 
    param.verbose_level = 1;
    param.plot_results_every = 50;

    param.remove_residues = false; % autodetect and remove phase residua 
    param.extension = ''; 

    %% multislice by ZC
     param.rmvac = true; % by ZC, remove last vacuum layer for multilayer
     param.Nlayers=1; % default single-slice
      
     param.layer4pos=[];
     param.background = 0;
    %% data handling 
    param.upsampling_data_factor = 0;           % assume that the data were created by upsampling using function utils.unbinning 

    param.damped_mask = 5e-3;  % if damped_mask = 0 -> do nothing, if 1>x>0  ->  push masked regions weakly towards measured magnitude value in each iteration
    
    param.background_detection = false; 
    param.background_width = inf;

    
    %% ADVANCED OPTIONS   
    
    param.object_regular =  [0, 0]; %  enforce smoothness !!!, use between [0-0.1 ]
    param.remove_object_ambiguity = true;    % remove intensity ambiguity between the object and the probes 
    param. variable_probe = false;           % Use SVD to account for variable illumination during a single (coupled) scan
    param. apply_subpix_shift = false;       % apply FFT-based subpixel shift, important for good position refinement but it is slow

    param.probe_geometry_model = {'scale', 'asymmetry', 'rotation', 'shear'};  % list of free parameters in the geometry model
    param.probe_position_search = inf;
    param.probe_fourier_shift_search = inf; 
    param.estimate_NF_distance = inf;
    param.detector_rotation_search = inf;   % rotation of the detector axis with respect to the sample axis, similar as rotation option in the position refinement geometry model but works also for 0/180deg rotation shared scans 
    param.detector_scale_search = inf;      % pixel scale of the detector, can be used to refine propagation distance in ptycho 

    param.apply_multimodal_update = false; % use thibault modes to get higher signal, it can cause isses, not real gain  if blur method is used 
    param.probe_backpropagate = 0; 
    param.beta_LSQ = 0.9;       % use predictive step length
    param.delta_p = 0.1;     % LSQ damping constant 
    param.variable_probe_modes = 1; % OPRP settings 
    param.variable_probe_smooth = 0;% OPRP settings 
    param. variable_intensity = false; % account fort variable intensity
    param.relaxed_object_constrain = 0; % enforce known object (inputs.object_orig)
    param.probe_position_error_max = 10e-9; % max expected error of the stages 
    param.probe_fourier_shift_search = inf; 
    param.momentum = 0;             % use mementume accelerated gradient decsent method 
    
    param.regularize_layers = 0;    % 0<R<1 -> apply regularization on the reconstructed layers 
    param.preshift_ML_probe = true; % multilayer ptycho extension: if true, assume that the provided probe is reconstructed in center of the sample. 
    
    param. initial_probe_rescaling = true;  % find the optimal scaling correction for the provided probe guess in the initial iteration 
    param. accelerated_gradients_start = inf;  % use accelerated gradients to speed up the convergence
    param. align_shared_objects = false;      % align multiple objects from various scans 

    
    % extra analysis
    param. get_fsc_score = false;         % measure evolution of the Fourier ring correlation during convergence 
    param. mirror_objects = false;        % mirror objects, useful for 0/180deg scan sharing 
    param.align_shared_objects = false;   % align the objects before sharing them onto single one 

    % fly scans 
    param.flyscan_offset = 0; 
    param.flyscan_dutycycle = 1;
    rng('default');
    rng('shuffle');


end
