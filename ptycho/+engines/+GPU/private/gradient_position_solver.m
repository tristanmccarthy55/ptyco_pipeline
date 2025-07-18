% GRADIENT_POSITION_SOLVER solve position errors in the real space 
%
% [pos_update, cache] = gradient_position_solver(self,xi,O,P,ind, iter, cache)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** xi        exit wave update vector 
% ** O         object views 
% ** P         probe or probes 
% ** ind       indices of of processed position 
% ++ iter      current iteration 
% ** cache     structure with precalculated values to avoid unnecessary overhead
%
% returns:
% ++ pos_update     position updates for each of the indices   
% ++ cache          updated  structure with precalculated values 
%
% see also: engines.GPU.LSQML, engines.GPU.PIE 


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


function [pos_update, probe_rotation,probe_scale,cache] = gradient_position_solver(self,xi,O,P,ind, iter, cache, par)
    import engines.GPU.GPU_wrapper.*
    import math.*
    import utils.*
    % use gradinent solver for position correction  

    low_mem_errs = {'id:parallel:gpu:array:OOMForOperation',...
     'id:MATLAB:LowGPUMem','MATLAB:LowGPUMem',...
     'parallel:gpu:array:OOM',...
     'parallel:gpu:device:UnknownCUDAError', ...
     'parallel:gpu:array:OOMForOperation', ...
     'parallel:gpu:array:FFTInternalError'};
    
    % wrapper around get_img_grad, in case of low memory it will try to repeat
    % Ntimes before giving up 
    pos_update = 0; probe_rotation = 0; probe_scale = 0; 
    N = 5; 
    for ii = 1:N
        try
            % reuse dx_O, dy_O to save memory !! 
            [dx_O,dy_O]=get_img_grad(O);
            
            if iter >= par.detector_rotation_search
                %% estimate detector rotation 
                xgrid  = Garray(linspace(-1,1,self.Np_p(1))'); 
                ygrid  = Garray(-linspace(-1,1,self.Np_p(2))); 
                [nom, denom] = Gfun(@get_coefs_mixed,xi, P, dx_O, dy_O, xgrid, ygrid);
                probe_rotation =  gather(sum2(nom)./ sum2(denom));
            end
            
            if iter >= par.detector_scale_search
                %% estimate detector scale (ie pixel scale error in farfield mode) 
                xgrid  = Garray(-linspace(-1,1,self.Np_p(2)) .* tukeywin(self.Np_p(2), 0.1)'); 
                ygrid  = Garray(-linspace(-1,1,self.Np_p(1))'.* tukeywin(self.Np_p(1), 0.1)); 
                [nom, denom] = Gfun(@get_coefs_mixed,xi, P, dx_O, dy_O, xgrid, ygrid);
                probe_scale =  gather(sum2(nom)./ sum2(denom));
                probe_scale = 0.5*mean(probe_scale) / mean(self.Np_p);
            end

            if iter >= par.probe_position_search
                %% estimate sample shift 
                [dx_O, denom_dx, dy_O, denom_dy] = Gfun(@get_coefs_shift,xi,P,dx_O, dy_O);
                dx =  sum2(dx_O)./ sum2(denom_dx);
                dy =  sum2(dy_O)./ sum2(denom_dy);
            end
            break
        catch ME
            warning('Low memory')
            if ~any(strcmpi(ME.identifier, low_mem_errs))
                rethrow(ME)
            end
            pause(1)
        end
    end
    if ii == 5
        rethrow(ME) 
    end

    if iter < par.probe_position_search
        return
    end
    
    shift = squeeze(Ggather(cat(4,dx, dy)));
    
    %disable limit by YJ. maybe better for small probes
    % prevent outliers and too rapid shifts 
    max_shift = min(0.1, 10*mad(shift)); 
    shift = min(abs(shift), max_shift) .* sign(shift); % avoid too fast jumps, <0.5px/iter is enough 
    %old code
    %shift = min(abs(shift), 0.2) .* sign(shift); % avoid too fast jumps, <0.5px/iter is enough 

    pos_update = reshape(shift,[],2); 
    
    if ~isfield(cache, 'velocity_map_positions')
        cache.velocity_map_positions = zeros(self.Npos,2,'single');
    end
    if ~isfield(cache, 'position_update_memory')
        cache.position_update_memory = {};
    end
    
    cache.position_update_memory{iter}(ind,:) = pos_update;
   
    %% USE MOMENTUM ACCELERATION TO MAKE THE CONVERGENCE FASTER
%     try
    ACC = 0; 
%     momentum_memory = 5; % remember 5 iterations 
%     
%     % only in case far field ptychography 
%     if isinf(self.z_distance) && sum(cellfun(@length, cache.position_update_memory) > 0) > momentum_memory 
%         for ii = 1:momentum_memory
%             corr_level(ii) = mean(diag(corr(cache.position_update_memory{end}(ind,:), cache.position_update_memory{end-ii}(ind,:))));
%         end
%         if all(corr_level > 0 )
%             %estimate optimal friction from previous steps 
%             poly_fit = polyfit(0:momentum_memory,log([1,corr_level]),1); 
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             gain = 0.5;                           % smaller -> lower relative speed (less momentum)
%             friction =  0.1*max(-poly_fit(1),0);   % smaller -> longer memory, more momentum 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         else
%            gain = 0; friction = 0.5; 
%         end
%         
%         cache.velocity_map_positions(ind,:) = cache.velocity_map_positions(ind,:)*(1-friction) + pos_update;
%         % apply the velocity to the refined positions , if the postition updated are sufficiently small 
% 
%         if max(abs(pos_update)) < 0.1
%             ACC = norm2(pos_update + gain*cache.velocity_map_positions(ind,:)) / norm2(pos_update); 
%             pos_update = pos_update + gain*cache.velocity_map_positions(ind,:); 
%         end
% 
%     end
%     catch
%         keyboard
%     end
    if any(ind==1)
        verbose(1,'Grad pos corr -- AVG step  %3.3g px , acceleration = %4.1f', max(abs(pos_update(:))), ACC)
    end
    
end
    
function [nom1, denom1, nom2, denom2] = get_coefs_shift(xi, P, dx_O, dy_O)

    dx_OP = dx_O.*P;
    nom1 = real(conj(dx_OP) .* xi);
    denom1 = abs(dx_OP).^2; 

    dy_OP = dy_O.*P;
    nom2 = real(conj(dy_OP) .* xi);
    denom2 = abs(dy_OP).^2; 

end

    
function [nom, denom] = get_coefs_mixed(xi, P, dx_O, dy_O, xgrid, ygrid)

    dm_O = dx_O .* xgrid + dy_O .* ygrid; 

    dm_OP = dm_O.*P;
    nom = real(conj(dm_OP) .* xi);
    denom = abs(dm_OP).^2; 

end
