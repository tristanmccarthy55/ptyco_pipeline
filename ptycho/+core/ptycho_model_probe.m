% pout = ptycho_model_probe(p)

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
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089).
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

function pout = ptycho_model_probe(p)
import utils.*
import io.*

% Define often-used variables
lambda = p.lambda;
asize = p.asize; % Diffr. patt. array size   
dx_spec = p.dx_spec;
a2 = prod(asize);

if check_option(p, 'prop_regime', 'nearfield')
    % for this task use original values of the pixel sizes
    dx_spec =  p.lambda*p.z*p.nearfield_magnification ./ (p.asize*p.ds);
end

% Prepare probe
if p.model_probe
    if p.model.probe_is_focused
        verbose(2, 'Using focused probe as initial model.');
        if asize(1) ~= asize(2)
            error('Focused probe modeling is only implemented for square arrays (please feel free to change that).');
        end

        if isempty(p.model.probe_zone_plate_diameter) || isempty(p.model.probe_outer_zone_width)
            zp_f = p.model.probe_focal_length;
            verbose(3, 'Using model.probe_focal_length for modeled probe.');
        else
            zp_f = p.model.probe_zone_plate_diameter * p.model.probe_outer_zone_width / lambda;
        end
        
        % The probe is generated in a larger array to avoid aliasing
        upsample = p.model.probe_upsample;

        defocus = p.model.probe_propagation_dist;
        Nprobe = upsample*asize(1);                       % Array dimension for the simulation
        dx = (zp_f+defocus)*lambda/(Nprobe*dx_spec(1));   % pixel size in the pupil plane
        r1_pix = p.model.probe_diameter / dx;              % size in pixels of first pinhole
        r2_pix = p.model.probe_central_stop_diameter / dx;       % size in pixels of central stop
        
        
        % Pupil
        [x,y] = meshgrid(-Nprobe/2:floor((Nprobe-1)/2),-Nprobe/2:floor((Nprobe-1)/2));
        r2 = x.^2 + y.^2;
%         w = (r2 < (r1_pix)^2);
        if upsample*asize(1) < round(r1_pix)-5
            error(sprintf('For this experimental parameters asize must be at least %d in order for the lens to fit in the window.',ceil((round(r1_pix)-5)./upsample+1)))
        end
        w = fftshift(filt2d_pad(upsample*asize(1), round(r1_pix)+5, round(r1_pix)-5, 'circ'));
        if p.model.probe_central_stop
            w = w .*(1-fftshift(filt2d_pad(upsample*asize(1), round(r2_pix)+2, round(r2_pix-2), 'circ')));
        end
        if isfield(p.model,'probe_structured_illum_power') && p.model.probe_structured_illum_power
            rng default
            r = utils.imgaussfilt2_fft(randn(upsample*p.asize),upsample*2); 
            r = r / math.norm2(r); 
            r = exp(1i*r*p.model.probe_structured_illum_power);
            w = imgaussfilt(w,upsample/2).*r; 
        end

        % Propagation
        probe_hr = prop_free_ff(w .* exp(-1i * pi * r2 * dx^2 / (lambda * zp_f)), lambda, zp_f + defocus, dx);

        % Cropping back to field of view
        probe = crop_pad(probe_hr, asize); 
        
        % prevent unreal sharp edges from the cropped tails in probe        
        [probe] = utils.apply_3D_apodization(probe, 0); 
                
        probe = probe .* sqrt(1e5/sum(sum(abs(probe).^2)));
        clear x y r2 w probe_hr  
        
    else
        verbose(2, 'Using circular pinhole as initial model.');
        [x1,x2] = ndgrid(-asize(1)/2:floor((asize(1)-1)/2),-asize(2)/2:floor((asize(2)-1)/2));
        probe = ( (x1 * dx_spec(1)).^2 + (x2 * dx_spec(2)).^2 < (p.model.probe_diameter/2)^2);
        probe = prop_free_nf(double(probe), lambda, p.model.probe_propagation_dist, dx_spec);
        clear x1 x2
    end
    verbose(3, 'Successfully generated model probe.');
else
    if ~isfield(p,'probe_file_propagation')
        p.probe_file_propagation = [];
    end
    verbose(2, 'Using previous run as initial probe.');
    
    % if string allows it, fill in the scan numbers     
     p.initial_probe_file = sprintf(replace(p.initial_probe_file,'\','\\'), p.scan_number(1)); 

    for searchpath = {'', p.ptycho_matlab_path}             
        fpath = dir(fullfile(searchpath{1},p.initial_probe_file)); 
        % check if only one unique file is found
        if length(fpath) > 1
            error('Too many paths corresponding to patterns %s were found', p.initial_probe_file)
        elseif length(fpath) == 1
            p.initial_probe_file = fullfile(fpath.folder, fpath.name);
            break
        end
    end
    if isempty(fpath)
        error(['Did not find initial probe file: ' p.initial_probe_file])
    end

    fileokflag = 0;
    while ~fileokflag
        try
            S = load_ptycho_recons(p.initial_probe_file, 'probe'); % avoid object loading when it is not needed 
            probe = S.probe;
            S = load_ptycho_recons(p.initial_probe_file, 'p');
            fileokflag = 1; 
            verbose(2, 'Loaded probe from: %s',p.initial_probe_file );

            %% check if the loaded probe was binned or no
            if isfield(S, 'p') && isfield(S.p, 'binning')
                binning = S.p.binning;
            elseif isfield(S, 'p') && isfield(S.p, 'detector') && isfield(S.p.detector, 'binning')
                binning = S.p.detector.binning;
            else
                if verbose() > 0
                    binning = []; 
                    while isempty(binning)
                        binning = str2num(input('Define binning factor 2^x for loaded initial probe (i.e. 0 for no binning):','s'));
                    end
                    % save provided binning option to the loaded probe file
                    S.p.detector.binning = binning; 
                    save(p.initial_probe_file, '-append', '-struct', 'S')
                else
                    % prevent stopping code if automatic reconstructions are running
                    verbose(0, 'Initial probe binning could not be determined, assuming no binning')
                    binning = 0; 
                end
            end
            % modify the loaded probe into a nonbinned version
            probe = crop_pad(probe, [size(probe,1),size(probe,2)]*2^binning); 
            
            
        catch err
            disp(['File corrupt: ' p.initial_probe_file])
            disp(err.message)
            disp('Retrying')
            
            pause(1)
        end
    end
    
    % ignore this to enable multiple probes by Zhen Chen
%     if ndims(probe)==3
%         sz_pr = size(probe);
%         probe = reshape(probe, [sz_pr(1) sz_pr(2) 1 sz_pr(3)]);
%     end
    verbose(3, 'File %s loaded successfully.', p.initial_probe_file);

    if ~all([size(probe,1) size(probe,2)] == asize)
        verbose(2,'Loaded probe has the wrong size.');
        verbose(2,'Interpolating probe in file %s, from (%d,%d) to (%d,%d).', p.initial_probe_file,size(probe,1),size(probe,2),asize(1),asize(2));
        probe = interpolateFT(probe,asize);
    end
    if ~isempty(p.probe_file_propagation) && any(p.probe_file_propagation ~= 0)
        verbose(2,'Propagating probe from file by %f mm',p.probe_file_propagation*1e3);
        probe= prop_free_nf(double(probe), lambda, p.probe_file_propagation, dx_spec);
    end    
    
end

    probe = probe .* sqrt(a2 ./ sum(sum(abs(probe).^2)));

pout = p;
pout.probe_initial = probe;
