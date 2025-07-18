% PTYCHO_SOLVER  the main loop of ptychography. Calls the selected engine 
% apply additional constraints, and tries to remove ambiguities
% 
% [outputs, fourier_error, fsc_score] = ptycho_solver(self, par, cache)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
% ** cache     structure with precalculated values to avoid unnecessary overhead
%
% returns:
% ++ outputs        self-like structure with final reconstruction
% ++ fourier_error  array [Npos,1] containing evolution of reconstruction error 
% ++ fsc_score      [] or a structure with outputs from online estimation of FSC curve
%
%


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
%       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.�?
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

function [outputs, fourier_error, fsc_score] = ptycho_solver(self, par, cache)

import engines.GPU.analysis.*
import engines.GPU.shared.*
import engines.GPU.initialize.*
import engines.GPU.GPU_wrapper.*
import math.*
import utils.*

% precalculate the parallel block sizes and sets 
[cache, par] = get_parallel_blocks(self, par, cache); 
 
global gpu use_gpu
verbose( par.verbose_level )    

if  (nargout == 1 && verbose()  == 0 && isinf(par.plot_results_every))
    par.get_error = false;
elseif ~isfield(par, 'get_error')
    par.get_error = true;
end

if use_gpu 
    verbose(struct('prefix',['GPU-',num2str(gpu.Index),'_', par.method]))%modified by YJ to print out more info
    %verbose(struct('prefix',['GPU-', par.method]))
    verbose(0,'Started solver using %s method on GPU %i', par.method, gpu.Index )
    
else
    verbose(struct('prefix',['CPU-', par.method]))
    verbose(0,'Started solver using %s method on CPU', par.method)
end


lastwarn('') 
fsc_score = cell(1,0);
par.Nscans = length(self.reconstruct_ind);

%% move everything on GPU if needed; 
if par.use_gpu
    % if not sparse solvers as ePIE, hPIE, MLs are used presplit data into
    % bunches (allow larger data to be processed )
    %%%%split_data = is_method(par, {'MLc', 'DM'});
    split_data = false; % modified by YJ, seems to avoid some errors from GPU
    [self, cache] =  move_to_gpu(self,cache, par.keep_on_gpu, split_data);
end 
if par.share_object && par.object_modes == 1
    % enforce only a single object 
    self.object = self.object(1,:);
    cache.illum_sum_0 = cache.illum_sum_0(1);
end

%% allocate memory 
fourier_error = Garray( nan(par.number_iterations, self.Npos));
if is_method(par, {'DM'})
    psi_dash = cell(max(par.probe_modes,par.object_modes), length(cache.preloaded_indices_simple{1}.indices));
end
if is_method(par, {'PIE', 'ML'})
    if par.beta_LSQ
        cache.beta_object = ones(self.Npos,par.Nlayers,'single')*par.beta_object;
        cache.beta_probe  = ones(self.Npos,par.Nlayers,'single')*par.beta_probe;
    else
        cache.beta_object = single(par.beta_object);
        cache.beta_probe  = single(par.beta_probe);
    end
    switch lower(par.likelihood)
        case 'l1', cache.beta_xi = 1; % optimal step for gauss
        case 'poisson', cache.beta_xi = 0.5*ones(1,1,self.Npos,'single'); % for poisson it will be further refined 
    end
end

%% in case of multilayer extension assume that the provided probe is positioned in middle of the sample -> shift it at the beginning
if par.Nlayers > 1 && par.preshift_ML_probe
   probe_offset = -sum(self.z_distance(1:end-1))/2; 
   for ii = 1:par.probe_modes
        self.probe{ii} = utils.prop_free_nf(self.probe{ii}, self.lambda , probe_offset ,self.pixel_size) ;
   end
end

%% in case of tilted plane ptychography, tilt the provided probe 
if any(par.p.sample_rotation_angles(1:2)) && check_option(par.p, 'apply_tilted_plane_correction', 'propagation') 
    % apply propagators to the tilted plane 
    for ii = 1:par.probe_modes % bugs! because modes are layers
        self.probe{ii} = self.modes{ii}.tilted_plane_propagate_fwd(self.probe{ii});
    end
end
    

global pprev;
pprev = -1;

mode_id = 1;  % main mode (assume single most important mode for approchimations) 

t0 = tic;
t_start = tic;

% object averaging for DM code 
for ll = 1:length(self.object)
    object_avg{ll} = 0; 
end
N_object_avg = 0;
%par.initial_probe_rescaling = true or false
for iter =  (1-par.initial_probe_rescaling):par.number_iterations
    
    if iter > 1
        %{
        if verbose() == 0
            progressbar(iter, par.number_iterations, max(20,round(sqrt(par.number_iterations))))
        else
            verbose(1,'Iteration %s: %i / %i  (time %3.3g  avg:%3.3g)', par.method, iter, par.number_iterations, toc(t_start), toc(t0)/(iter-1))
        end
        %}
        %modified by YJ to print out more details

        %verbose(0,'Iteration %s: %i / %i  (time %3.3g  avg:%3.3g)', par.method, iter, par.number_iterations, toc(t_start), toc(t0)/(iter-1))
        %avgTimePerIter = toc(t0)/(iter-1);
        %timeLeft = (par.number_iterations-iter+1)*avgTimePerIter;
        %verbose(0, 'Method: %s, GPU id: %i',par.method, gpu.Index)
        %if timeLeft>3600
         %   verbose(0, 'Iteration: %i / %i  (Time left:%3.3g hour. avg:%3.3g sec)', iter, par.number_iterations, timeLeft/3600, avgTimePerIter)
        %elseif timeLeft>60
         %   verbose(0,'Iteration: %i / %i  (Time left:%3.3g min. avg:%3.3g sec)', iter, par.number_iterations, timeLeft/60, avgTimePerIter)
        %else
        %    verbose(0,'Iteration: %i / %i  (Time left:%3.3g sec. avg:%3.3g sec)', iter, par.number_iterations, timeLeft, avgTimePerIter)
        %end
        avgTimePerIter = toc(t0)/(iter-1);
        timeLeft = (par.number_iterations - iter + 1) * avgTimePerIter;
        iterTime = toc(iter_t_start);  % 当前 iteration 的耗时

        if timeLeft > 3600
            verbose(0, 'Iteration: %i / %i  (This iter: %.3g sec, Time left: %.3g hour, avg: %.3g sec)', ...
                iter, par.number_iterations, iterTime, timeLeft/3600, avgTimePerIter);
        elseif timeLeft > 60
            verbose(0, 'Iteration: %i / %i  (This iter: %.3g sec, Time left: %.3g min, avg: %.3g sec)', ...
                iter, par.number_iterations, iterTime, timeLeft/60, avgTimePerIter);
        else
            verbose(0, 'Iteration: %i / %i  (This iter: %.3g sec, Time left: %.3g sec, avg: %.3g sec)', ...
                iter, par.number_iterations, iterTime, timeLeft, avgTimePerIter);
        end

    end
    iter_t_start = tic;  % 开始计时当前 iteration
    t_start = tic;
    if  iter > 0.9*par.number_iterations && is_method(par, 'DM')
        for ll = 1:length(self.object)
            object_avg{ll} = object_avg{ll}  + self.object{ll}; 
        end
        N_object_avg = N_object_avg +1 ; 
        verbose(1,'==== Averaging DM result ======')
    end
     
    %% GEOMETRICAL CORRECTIONS        
    if (iter > par.probe_position_search || iter > par.detector_rotation_search) && is_method(par, {'PIE', 'ML'})
        self = find_geom_correction(self,cache,par,iter,mode_id);
    end  
    %% remove extra degree of freedom for OPRP and other optimizations
    if iter > par.probe_fourier_shift_search
        for kk = 1:par.Nscans
            ind = self.reconstruct_ind{kk};
            self.modes{1}.probe_fourier_shift(ind,:) = self.modes{1}.probe_fourier_shift(ind,:) - mean(self.modes{1}.probe_fourier_shift(ind,:));
        end    
    end
        
    %% remove ambiguity related to the variable probe 
    if par.variable_probe && iter > par.probe_change_start && is_method(par, 'ML')
         self = remove_variable_probe_ambiguities(self,par); 
    end
         
   %% remove the ambiguity in the probe / object reconstruction => keep average object transmission around 1
   if  mod(iter,10)==1 &&  par.remove_object_ambiguity  && ~is_used(par, {'fly_scan'}) &&  ~is_method(par, {'DM', 'PIE'})  % too slow for variable probe 
        self = remove_object_ambiguity(self, cache, par) ; 
   end
   
    if (mod(iter, 10) == 1 || iter  < 5) && check_option(par, 'get_fsc_score')   && ...
       (((par.Nscans > 1 ) && size(self.object,1) == par.Nscans) || ... 
       ( check_option(self, 'object_orig') ))
            
        %% Fourier ring correlation between two scans with independend objects 
        aux = online_FSC_estimate(self, par, cache, fsc_score(end,:), iter); 
        fsc_score(end+1,1:size(aux,1), 1:size(aux,2)) = aux; 
    end   
    
    %% ADVANCED FLY SCAN 
    if  is_used(par, 'fly_scan')
        if iter == 1
           disp(['== AVG fly scan step ', num2str( median(sqrt(sum(diff(self.probe_positions_0).^2,1)))  )]) 
        end
        self = prepare_flyscan_positions(self, par); 
    end
    
    
    %% update current probe positions (views)
    if iter <= 1 || iter >= par.probe_position_search
        %%%%%%%%  crop only ROI of the full image for subsequent calculations  %%%%%%%%%% 
        for ll = 1:par.Nmodes
            if ~is_used(par, 'fly_scan') %  && isempty(self.modes{1}.ASM_factor)
                %% conventional farfield/nearfield ptycho -> update view coordinates and keep subpixels shift < 1 px 
                [cache.oROI_s{ll},cache.oROI{ll},sub_px_shift] = find_reconstruction_ROI( self.modes{1}.probe_positions,self.Np_o, self.Np_p); 
                self.modes{1}.sub_px_shift = sub_px_shift; 
            else % if is_used(par, 'fly_scan')   
                %% flyscan farfield ptycho -> update view coordinates and keep subpixels shift < 1 px 
                [cache.oROI_s{ll},cache.oROI{ll},sub_px_shift] = find_reconstruction_ROI( self.modes{1}.probe_positions,self.Np_o, self.Np_p); 
                %% use the fftshift  for much larger corrections in the case of the fly scan 
                self.modes{ll}.sub_px_shift = self.modes{ll}.probe_positions -  self.modes{1}.probe_positions_0 + sub_px_shift;
%             else
%                 %% nearfield ptycho -> keep view coordinates and update subpixels shift only -> assume that positon correction was only minor
%                 if iter <= 1
%                     [cache.oROI_s{ll},cache.oROI{ll},sub_px_shift] = find_reconstruction_ROI( self.modes{1}.probe_positions_0,self.Np_o, self.Np_p); 
%                 end
%                 self.modes{ll}.sub_px_shift = self.modes{ll}.probe_positions -  self.modes{1}.probe_positions_0;               
            end
        end
        
    end
    
    %% update probe fft support window  if mode.probe_scale_upd(end) ~= 0, important to avoid issues during subpixel probe rescaling (is pixel scale search)
    if self.modes{1}.probe_scale_upd(end) > 0
        self.modes{1}.probe_scale_window = get_window(self.Np_p, 1+self.modes{1}.probe_scale_upd(end), 1) .* get_window(self.Np_p, 1+self.modes{1}.probe_scale_upd(end), 2); 
    elseif self.modes{1}.probe_scale_upd(end) < 0
        self.modes{1}.probe_scale_window = fftshift(get_window(self.Np_p, 1-self.modes{1}.probe_scale_upd(end), 1) .* get_window(self.Np_p, 1-self.modes{1}.probe_scale_upd(end), 2)) ; 
    else
        self.modes{1}.probe_scale_window = [];
    end    
   
    %% updated illumination
    if iter <= 1 || ( iter > par.probe_change_start && (mod(iter, 10) == 1 || iter < par.probe_change_start+10 ))
        aprobe2 = abs(self.probe{1}(:,:,1)).^2; 
        for ll = 1:size(self.object,1)
            if par.share_object
                ind = [self.reconstruct_ind{:}];
            else
                ind = self.reconstruct_ind{ll};
            end
            % avoid oscilations by adding momentum term 
            cache.illum_sum_0{ll} = set_views(cache.illum_sum_0{ll}, Garray(aprobe2), 1,1,ind, cache)/2;
            cache.illum_norm(ll) = norm2(cache.illum_sum_0{ll});
            cache.MAX_ILLUM(ll) = max2(cache.illum_sum_0{ll});
        end
    end
    
    %% improve convergence speed by gradient acceleration 
    if is_method(par, 'MLc') && iter >= par.accelerated_gradients_start
        [self, cache] = accelerate_gradients(self, par, cache, iter); 
    end
    
        %% suppress large amplitude of object, quick try by Zhen Chen
%     if iter >= par.object_change_start && par.Nlayers > 1
%         for ll=1:par.Nlayers
%             temp = abs(self.object{ll});
%             if any(gather(temp(:))>10)
%                 keyboard;
%             end
%             temp (temp> 10) =1;
%             self.object{ll} = temp.* exp(1i* angle(self.object{ll}));
%         end
%     end
    %% %%%%%%%%%%%%%%%%%%%%%%%%  PERFORM ONE ITERATION OF THE SELECTED METHOD %%%%%%%%%%%%%%%%%%%%%

    switch  lower(par.method)
        case {'epie', 'hpie'}
            [self, cache, fourier_error] = engines.GPU.PIE(self,par,cache,fourier_error,iter);
        case { 'mls','mlc'}
            [self, cache, fourier_error] = engines.GPU.LSQML(self,par,cache,fourier_error,iter);
        case 'dm'
            [self, cache,psi_dash,fourier_error] =  engines.GPU.DM(self,par,cache,psi_dash,fourier_error,iter);
        otherwise
            error('Not implemented method')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if iter == 0; continue; end  % interation 0 is used only to calibrate iinitial probe intensity
    
    if verbose() > 0  && any(~isnan(fourier_error(iter,:)))       
        switch lower(par.likelihood)
            case  'l1', verbose(1,'=====  Fourier error = %3.4g ', nanmedian(fourier_error(iter,:)) ); 
            case 'poisson'
                err = fourier_error(iter,:) - fourier_error(1,:);
                verbose(1,'=====  Log likelihood = %3.5g ', nanmedian(err)); 
        end
    end
 
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%% CORRECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % apply probe constraints 
    if iter >= par.probe_change_start
        if ~ (check_option(par.p,'probe_support_tem') && check_option(par.p,'probe_support_tem_Nend') && iter > par.p.probe_support_tem_Nend)  % quick implement for TEM aperture constrain by Zhen Chen, bugs for probe rescale.
                    
            % propagate probe on the detector 
            probe_temp = fft2_safe(self.probe{1});  

            if ~isempty(self.modes{1}.probe_support_fft)
                probe_temp = probe_temp .* self.modes{1}.probe_support_fft;
            end

            % propagate probe back to the sample plane 
            self.probe{1} = ifft2_safe(probe_temp);  
        elseif ~ check_option(par.p,'probe_support_tem') % force not use near field propagation if probe_support_tem exists, by Zhen Chen
            self.probe{mode_id} = apply_probe_contraints(self.probe{mode_id}, self.modes{mode_id});
        end
    end
  
    % push low illum regions of object to zero  
    if par.delta > 0 
        if  iter > par.object_change_start
            for ll = 1:par.Nlayers % max(par.object_modes, par.Nscans) % For multilayer by Zhen Chen
                % push everywhere to zero, even out of the object region
                self.object{ll} = Gfun(@regular_object_out_regions, self.object{ll}, cache.illum_sum_0{1}, cache.MAX_ILLUM(1),par.delta); 
            end
        end
        if  iter > par.probe_change_start && ~par.variable_probe
            for ll = 1:par.probe_modes
                % push everywhere to zero, even out of the object region
                self.probe{ll} = self.probe{ll}  .* (1-par.delta);
            end
        end
    end
             
    
    %% suppress uncontrained values !! 
    if iter > par.object_change_start && par.object_regular(1) > 0
        for ll = 1:par.Nlayers % max(par.object_modes, par.Nscans)     % quick fix for multilayer, but bugs for multiple scans or modes, by Zhen Chen      
            self.object{ll} = apply_smoothness_constraint(self.object{ll},par.object_regular(1)); % blur only intensity, not phase 
        end
    end
    
%     %% suppress large amplitude of object, quick try by Zhen Chen
    if iter >= par.object_change_start && par.Nlayers > 1   && iter < 100
        for ll=1:par.Nlayers
            temp = abs(self.object{ll});
            temp (temp> 1.5) =1;
            self.object{ll} = temp.* exp(1i* angle(self.object{ll}));
        end
    end

%% reverse amplitude larger than 1 , try by Zhen Chen
%     if iter >= par.object_change_start && par.Nlayers > 1 && par.p.positivity_absorption_constraint
%         for ll=1:par.Nlayers
%             temp = abs(self.object{ll});
%             temp_ang=angle(self.object{ll});
%             temp_ang(temp>1)=-temp_ang(temp>1);
%             temp (temp> 1) = 2-temp(temp>1);            
%             self.object{ll} = temp.* exp(1i* temp_ang);
%         end
%     end
    %% weak positivity object 
    if iter > par.object_change_start  && any(par.positivity_constraint_object)
         for ll = 1:par.Nlayers % par.object_modes %  For multilayer by Zhen Chen 
             self.object{ll}(cache.object_ROI{:}) = ...
                 Gfun(@positivity_constraint_object,self.object{ll}(cache.object_ROI{:}), par.positivity_constraint_object);
         end
    end
    
    %% probe orthogonalization 
    if par.probe_modes > par.Nrec && (~is_method(par, 'DM') || iter == par.number_iterations)
        %  orthogonalization of incoherent probe modes 
        if is_used(par, 'fly_scan')
            probes = self.probe;
            % orthogonalize the modes with all the other shifted modes 
            for i = 1:par.Nrec
                dx = mean(self.modes{i}.sub_px_shift - self.modes{1}.sub_px_shift); 
                % apply average shift 
                probes{i} = imshift_fft(probes{i}, dx); 
            end
            probes = ortho_modes(probes);  % perform othogonalization 
            % update only the incoherent orthogonal modes 
            self.probe(1+par.Nrec:par.probe_modes) = probes(1+par.Nrec:par.probe_modes);
        else
            %% orthogonalize the incoherent probe modes 
            for ii = 1:par.probe_modes
                P(:,:,ii) = self.probe{ii}(:,:,1); 
            end
            P = core.probe_modes_ortho(P);
            for ii = 1:par.probe_modes
                self.probe{ii}(:,:,1) = P(:,:,ii); 
            end
        end
    end

  %% regularize multilayer reconstruction 
  if par.regularize_layers > 0 && par.Nlayers > 1 %  && mod(iter, 2) == 1, ## change from Nlayers > 1 to Nlayers > 2 by Zhen Chen
        self = regulation_multilayers(self, par, cache);
  end
  
  % constraint periodic along propagation, by Zhen Chen.
        % regularize_layers works better (set regularize_layers > 0), then this is unnecessary.
    if par.Nlayers > 1 && isfield(par.p,'forcelayer') && par.p.forcelayer && (iter >= par.p.Nst_forcelayer && iter <= par.p.Nend_forcelayer)
        object_avg=0;
        for layer=1:par.Nlayers 
            object_avg = object_avg + self.object{layer};
        end
        object_avg = object_avg ./(par.Nlayers);
%             object_avg = (object_avg .* exp(1i*angle(self.object{par.Nlayers})))./(par.Nlayers-1);
        for layer=1:par.Nlayers
            self.object{layer} = object_avg ;
        end  
    end
    
    
  %% PLOTTING 
    %%%% plot  results %%%%%%%%%%%%
    if mod(iter, par.plot_results_every ) == 0 &&  par.plot_results_every ~=0   
          try
              if verbose()  <= 0
                  % use cSAXS plorring rutines 
                  ptycho_plot_wrapper(self, par, fourier_error)
              else
                %   use more detailed plotting rutines 
                if (par.probe_modes > 1 )
                    %% probe incoherent modes 
                    plot_probe_modes(self,par);                        
                end
                if ( par.Nlayers > 1) || (par.Nscans > 1 && ~par.share_object)
                    %% object incoherent modes 
                    plot_object_modes(self, cache)
                end

                plot_results(self,cache, par, Ggather(fourier_error), ...
                self.modes{mode_id}.probe_positions)
              
              end

                % show variable modes 
                if (par.variable_probe  && par.variable_probe_modes > 0)
                    plot_variable_probe(self, par)
                end
                
              
                % show position correction in the fourier plane 
                if iter > par.probe_fourier_shift_search
                    plotting.smart_figure(24654)
                    clf
                    hold all
                    for ll = 1:par.Nmodes    
                        plot(self.modes{ll}.probe_fourier_shift)
                    end    
                    hold off
                    grid on 
                    axis tight
                    xlabel('Position #')
                    ylabel('Corrected probe shift in Fourier plane [px]')
                    title('Fourier space probe shift')
                end                
                
   
                % show position correction 
                if iter > min([par.probe_position_search, par.estimate_NF_distance, par.detector_rotation_search, par.detector_scale_search]) ...
                        && is_method(par, {'PIE', 'ML'}) 
                         plot_geom_corrections(self, self.modes{1}, Ggather(self.object{1}),iter, par, cache)
                    if iter > min(par.detector_rotation_search, par.detector_scale_search)
                    for i = 1:max(1,par.Nrec)
                        verbose(1,sprintf(' Reconstruction id: %i =============  Detector pixel scale: %0.5g  Detector rotation: %0.5g deg', ...
                                            i,  1-self.modes{i}.probe_scale_upd(end) , self.modes{i}.probe_rotation(end,1)))
                    end
                    end
                end
                    
            if par.get_fsc_score && ~isempty(fsc_score)
                % plot score estimated by the fourier ring correlation 
                plot_frc_analysis(fsc_score, par)
            end
        
            drawnow
   
        catch err
            warning(err.message)
            if verbose()  > 1
                keyboard
            end
        end
   
    end
    %% save intermediate images, added by YJ
    if isfield(par,'save_results_every') && (mod(iter, par.save_results_every ) == 0 &&  par.save_results_every ~=0) || iter == par.number_iterations
         if ~exist(par.fout, 'dir')
            mkdir(par.fout);
         end
        % save temporary outputs
        if iter == par.number_iterations
            outputs = self;
            outputs.diffraction = [];
        else
            outputs = struct();
        end
        for ll = 1:par.probe_modes
            outputs.probe{ll} = Ggather(self.probe{ll});
        end
        for ll = 1:max(par.object_modes,par.Nlayers)
            outputs.object{ll} = Ggather(self.object{ll});
            object_roi = Ggather(self.object{ll});
            outputs.object_roi{ll} = object_roi(cache.object_ROI{:});
        end
        if strcmp(par.likelihood, 'poisson')
            fourier_error_out = fourier_error- fourier_error(1,:);
            fourier_error_out = Ggather(fourier_error_out);
        else
            fourier_error_out = Ggather(fourier_error);
        end
        fourier_error_out = mean(fourier_error_out,2,'omitnan'); % fix nan bug by ZC
        outputs.probe_positions = self.modes{1}.probe_positions;

        %if par.variable_probe
        %    outputs.probe_evolution = self.probe_evolution;
        %end
        outputs.pixel_size = self.pixel_size;
        %outputs.probe_positions_0 = self.probe_positions_0;
        if iter == par.number_iterations
            par.p.fmag = []; % delete large but not usefull params
            par.p.fmask = [];
            par.p.scanidxs=[];
            par.p.probes=[];
%             for ieng=1:length(par.p.engines)
%                 if isfield(par.p.engines{ieng},'probes_final')
%                     par.p.engines{ieng}.probes_final=[];
%                 end
%                 if isfield(par.p.engines{ieng},'object_final')
%                     par.p.engines{ieng}.object_final=[];
%                 end
%             end
            save(strcat(par.fout,'Niter',num2str(iter),'.mat'),'outputs','fourier_error_out','par','-v7.3');
        else
            save(strcat(par.fout,'Niter',num2str(iter),'.mat'),'outputs','fourier_error_out','-v7.3');
        end

        %save(strcat(par.fout,'Niter',num2str(iter),'.mat'),'outputs','fourier_error_out','avgTimePerIter');
        %% save object phase
        if isfield(par,'save_phase_image') && par.save_phase_image
            for ll=1:par.Nlayers
                O_phase_roi = phase_unwrap(angle(outputs.object_roi{ll}));
                %O_phase_roi = angle(outputs.object_roi);
                saveName = strcat('O_phase_roi_','Niter',num2str(iter),'_Layer',num2str(ll),'.tiff');
                saveDir = strcat(par.fout,'/O_phase_roi/');
                if ~exist(saveDir, 'dir')
                    mkdir(saveDir)
                end
                %imwrite(mat2gray(O_phase_roi),strcat(saveDir,fileName),'tiff');

                % save as single tiff
                tagstruct.ImageLength     = size(O_phase_roi,1);
                tagstruct.ImageWidth      = size(O_phase_roi,2);
                tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
                tagstruct.BitsPerSample   = 16;
                tagstruct.SamplesPerPixel = 1;
                tagstruct.RowsPerStrip    = 16;
                tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                tagstruct.Software        = 'MATLAB';
                t = Tiff(strcat(saveDir,saveName),'w');
                t.setTag(tagstruct)
                t.write(uint16(mat2gray(O_phase_roi)*2^16));
                t.close();
            end
        end
        %% save probe mage
        if isfield(par,'save_probe_mag') && par.save_probe_mag
            for ii=1:size(outputs.probe{1},3)
                probe_mag = zeros(size(outputs.probe{1}(:,:,ii,1),1),size(outputs.probe{1}(:,:,ii,1),2)*length(outputs.probe));
                for jj=1:length(outputs.probe)
                    x_lb = (jj-1)*size(outputs.probe{1}(:,:,ii,1),2)+1;
                    x_ub = jj* size(outputs.probe{1}(:,:,ii,1),2);
                    probe_mag(:,x_lb:x_ub) = abs(outputs.probe{jj}(:,:,ii,1));
                end

                saveName = strcat('probe_mag_Niter',num2str(iter),'_',num2str(ii),'th_probe','.tiff');
                saveDir = strcat(par.fout,'/probe_mag/');
                if ~exist(saveDir, 'dir')
                    mkdir(saveDir)
                end
                imwrite(mat2gray(probe_mag)*64,parula,strcat(saveDir,saveName),'tiff')
            end
        end
    end
    
end

    %% return results 
    if is_method(par, 'DM')
        % return average from last 10% iterations 
        for ll = 1:length(self.object)
             self.object{ll} = object_avg{ll} /N_object_avg; 
        end
    end


    iter_time = toc(t0) / par.number_iterations; 
    verbose(1,' ====  Time per one iteration %3.3fs', iter_time)    
    verbose(1,' ====  Total time %3.2fs', toc(t0))    

    
    % clip outliers from the low illum regions 
    MAX_OBJ = 0;
    for ii = 1:length(self.object)
        MAX_OBJ = max(MAX_OBJ, max2(abs(self.object{ii}(cache.object_ROI{:}))));
    end
    for ii = 1:length(self.object)
        aobj = abs(self.object{ii});
        self.object{ii} = min(MAX_OBJ, aobj) .* self.object{ii} ./ max(aobj, 1e-3);
    end
    

    %% in case of tilted plane ptychography, back-tilt to the detector probe 
    if any(par.p.sample_rotation_angles(1:2)) && check_option(par.p, 'apply_tilted_plane_correction', 'propagation') 
        % apply propagators to the tilted plane 
        for ii = 1:par.probe_modes
            self.probe{ii} = self.modes{ii}.tilted_plane_propagate_back(self.probe{ii});
        end
    end

    %% in case of multilayer extension assume return the reconstructed probe in middle of the sample !! 
    if par.Nlayers > 1 && par.preshift_ML_probe
       probe_offset = sum(self.z_distance(1:end-1))/2; 
       for ii = 1:par.probe_modes
            self.probe{ii} = utils.prop_free_nf(self.probe{ii}, self.lambda , probe_offset ,self.pixel_size) ;
       end
    end
    
    outputs = self;
    
    % avoid duplication in memory 
    self.noise = [];
    self.diffraction= [];
    self.mask = [];
    
    % store useful parameters back to the main structure only for the 1st mode 
    outputs.relative_pixel_scale = self.modes{1}.scales(end,:);
    outputs.rotation =  self.modes{1}.rotation(end,:);
    outputs.shear =   self.modes{1}.shear(end,:);
    outputs.z_distance = self.modes{1}.distances(end,:);
    outputs.shift_scans = self.modes{1}.shift_scans;
    outputs.probe_fourier_shift = self.modes{1}.probe_fourier_shift; 
    outputs.probe_positions = self.modes{1}.probe_positions;
    outputs.detector_rotation = self.modes{1}.probe_rotation(end,:);
    outputs.detector_scale = 1+self.modes{1}.probe_scale_upd(end);

    if strcmp(par.likelihood, 'poisson')
        fourier_error = fourier_error- fourier_error(1,:);
    end
    fourier_error = Ggather(fourier_error);


    outputs.illum_sum = cache.illum_sum_0; 

    outputs.diffraction = [];
    outputs.noise = [];
    outputs.mask = [];

    %% move everything back to RAM from GPU 
    if par.use_gpu
        outputs =  move_from_gpu(outputs);
    end
    
    %% report results
    try
        verbose(0,'==== REPORT ==== \n SNR %4.3g  RES %3.2g (%3.3gnm) AuC: %3.3g\n\n', fsc_score{end,1}.SNR_avg,...
                    fsc_score{end,1}.resolution,mean(self.pixel_size)*1e9/fsc_score{end,1}.resolution, fsc_score{end,1}.AUC)
    end
    iter_t_start = tic;  % 开始计时当前 iteration
end

function x = positivity_constraint_object(x, relax)
   x = relax.*abs(x) + (1-relax).*x;
end
function object = regular_object_out_regions(object, illum, max_illum, delta)
    % push everywhere to zero, even out of the object region
    W = illum / max_illum;
    W = W ./ (0.1+W);
    object = object .* (W + (1-W).*(1-delta));
end

function win = get_window(Np_p, scale, ax)
    % aux function 
    % apodize window for img to prevent periodic boundary errors 
    import engines.GPU.GPU_wrapper.Garray
    win = ones(floor(Np_p(ax)/scale/2-2)*2); 
    win = utils.crop_pad(win, [Np_p(ax),1]);
    win = shiftdim(win, 1-ax);
    win = Garray(win);
end


 
