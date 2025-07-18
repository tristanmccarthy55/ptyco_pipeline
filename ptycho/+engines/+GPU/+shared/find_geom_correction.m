% FIND_GEOM_CORRECTION use current probe positions estimates to update geometry model and
% improve the  new probe positions 
%
% [self] = find_geom_correction(self,cache, par, iter,best_mode_id)
% 
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** cache     structure with precalculated values to avoid unnecessary overhead
% ** par       structure containing parameters for the engines 
% ** iter      current iteration 
% ** best_mode_id   strongest mode id
%
% returns:
% ++ self        self-like structure with final reconstruction
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


function [self] = find_geom_correction(self,cache, par, iter,best_mode_id)

    import engines.GPU.GPU_wrapper.*
    import engines.GPU.shared.*
    import utils.*
    import math.*
    
    mode = self.modes{best_mode_id};
   
    %% constrain the detector rotation 
     
    % store only the single update per scan 
    if iter > par.detector_rotation_search
        for ii = 1:length(self.reconstruct_ind)
            ind = self.reconstruct_ind{ii};
            mrot(ii) = mean(mode.probe_rotation_all(ind)); 
        end
        if par.mirror_objects
            % enforce zero average rotation if two mirror scans are provided
            mrot = mrot - mean(mrot);
        end
        for ii = 1:length(self.reconstruct_ind)
            mode.probe_rotation_all(ind) = mrot(ii) ;
        end
        mode.probe_rotation(end+1,:) = mrot;
    end
    
    
    if iter <= par.probe_position_search
        self.modes{best_mode_id} = mode;
        return
    end

    
    pos = mode.probe_positions;
    pos_0 = mode.probe_positions_0;
       
    
    if all(isnan(mode.probe_positions_weight(:))) || all(mode.probe_positions_weight(:)==0)
        %% EMPIRICAL ESTIMATION OF POSITION RELIABILITY
        verbose(1,'EMPIRICAL ESTIMATION OF POSITION RELIABILITY\n')
        illum = utils.crop_pad(abs(self.probe{1}(:,:,1)).^2, self.Np_p/2); 
        total_variation = zeros(self.Npos,2, 'single'); 
        
        for ii = 1:par.Nscans
%             best_layer = ceil(par.Nlayers/2); % change to middle layer by ZC
            best_layer = par.Nlayers;
            o_tmp =  self.object{min(end,ii), best_layer}; 
            o_tmp = o_tmp ./ max2(abs(o_tmp(cache.object_ROI{:})));
            % keep it more memory effecient (important for GPU !! )
            Npos = length(self.reconstruct_ind{ii});
            for jj = 1:ceil(Npos/par.grouping)
                ind = 1+(jj-1)*par.grouping:min(Npos, jj*par.grouping);
                obj_proj = get_views(o_tmp,[],1,1,self.reconstruct_ind{ii}(ind),cache);
                obj_proj = utils.crop_pad(obj_proj, self.Np_p/2); 

                [nx, ny,~]  = size(obj_proj); 
                [X,Y] = meshgrid(-ny/2:ny/2-1, -nx/2:nx/2-1);
                % suppress edge effects of the FFT derivatives 
                spatial_filter = exp(-(X.^16+Y.^16)/(min(nx,ny)/2.2)^16);
                obj_proj = obj_proj.* spatial_filter;
                [dX, dY] = get_img_grad(obj_proj);
                clear obj_proj 
                illum_proj = get_views(utils.imgaussfilt2_fft(cache.illum_sum_0{min(ii,end)},self.Np_p/10),[],1,1,self.reconstruct_ind{ii}(ind),cache);
                illum_proj = utils.crop_pad(illum_proj, self.Np_p/2); 

                dX = abs(dX) .* illum_proj.* illum; 
                dY = abs(dY) .* illum_proj.* illum;
                clear illum_proj
                total_variation(self.reconstruct_ind{ii}(ind),:) = Ggather(sqrt(squeeze([mean2(dX),mean2(dY)]))');
                clear dX dY
            end
        end
        mode.probe_positions_weight = total_variation.^4./mean(total_variation.^4); 
    end

    probe_positions_weight = double(mode.probe_positions_weight);

    jj = size(mode.scales,1)+1; 
        
    % find geometry for each scan separatelly 
    for ii = 1:par.Nscans
        ind = self.reconstruct_ind{ii};
        C0 = mode.affine_matrix(:,:,ii) - eye(2);
        C0 = C0(:);

        if par.Nscans > 1 && par.share_object 
            % it the case of multiple scans allow also freedom of coordinates shifts
            pos_fun = @(C)(( [1+C(1), C(2); C(3), 1+C(4)]*pos_0(ind,:)')' + C([5,6])' );
            if isfield(mode, 'shift_scans' ) && size(mode.shift_scans,2)>=ii
                C0(5:6) = mode.shift_scans(:,ii);
            else
                C0(5:6) = 0;
            end
        else
            pos_fun = @(C)(( [1+C(1), C(2); C(3), 1+C(4)]*pos_0(ind,:)')' );
        end
            
        err_fun = @(C)( probe_positions_weight(ind,:) .* (pos(ind,:) - pos_fun(C))); 

        options = optimoptions('lsqnonlin','Display','off');
        C(:,ii) = lsqnonlin( err_fun, C0,[],[],options) ; 
        

        %% restrict the geometry model only to the allowed degreed of freedom 
        % ===================================================================
        M{ii} = reshape(C(1:4,ii),2,2)+eye(2); 

        [scale, asymmetry, rotation, shear] = decompose_affine_matrix(M{ii}); 
        if ~ismember('scale', par.probe_geometry_model)
            scale = 1;
        end
        if ~ismember('asymmetry', par.probe_geometry_model)
            asymmetry = 0;
        end
        if ~ismember('rotation', par.probe_geometry_model)
            rotation = 0;
        end
        if ~ismember('shear', par.probe_geometry_model)
            shear = 0;
        end
        M{ii} = compose_affine_matrix(scale, asymmetry, rotation, shear);
        % ===================================================================

        
        mode.scales(jj,ii) =  scale;
        mode.asymmetry(jj,ii) = asymmetry;
        mode.rotation(jj,ii) =  rotation;
        mode.shear(jj,ii)    =  shear;
        if par.Nscans > 1 && par.share_object 
            mode.shift_scans(:,ii) =   C(5:6,ii);
        else
            mode.shift_scans(:,ii) =   [0,0];
        end
        
        % store initial guess 
        mode.affine_matrix(:,:,ii) = M{ii}; 

        % calculate ideal model positions 
        pos_model(ind,:) = pos_fun([reshape(M{ii} - eye(2), [],1); mode.shift_scans(:,ii)]); 

    end


    self.affine_matrix = M;    
    verbose(2,['-----  Geom. correction  ', repmat('%3.3g ', 1,length(C))], C)


    % use average 
    resid_pos= pos - pos_model; 
    
    % ignore errors in the global shift of the positions 
    for ii = 1:par.Nscans
        ind = self.reconstruct_ind{ii};
        resid_pos(ind,:)  = resid_pos(ind,:) - mean(resid_pos(ind,:)); 
    end
    
    err = abs(resid_pos); 

    max_err =  par.probe_position_error_max ./ self.pixel_size .* self.relative_pixel_scale; 

    verbose(1, '==== AVG position error %3.2g px MAX error %3.2g LIMIT %3.2g px ', mean(err(:)), max(err(:)),  max(max_err))
    
     
    %% apply only relaxed constrain on the probe positions  !!! 
    if par.apply_relaxed_position_constraint  % added by ZC
        relax = 0.1;
        % constrain more the probes in flat regions 
        W = relax*(1-  (probe_positions_weight./ (1+probe_positions_weight)));  
        % penalize positions that are further than max_err from origin 
        W = min(10*relax, W+max(0,err - max_err).^2 ./ max_err.^2  );  % avoid travel larger than max error
    else
        W=0;
    end
    
   % allow free movement in depenence on realibility and max allowed error 
    pos_new =  pos .*(1-W)+W.*pos_model; %disabled by YJ
%     pos_new = pos;
    mode.probe_positions = pos_new;
    mode.probe_positions_model = pos_model;

    
    if any(isnan(mode.probe_positions(:)))
        keyboard
    end


    
    
    
    self.modes{best_mode_id} = mode;

end

