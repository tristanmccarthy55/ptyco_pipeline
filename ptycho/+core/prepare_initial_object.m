%PREPARE_INITIAL_OBJECT 
% prepare an initial guess for the object reconstruction
%
% ** p      p structure
%
% returns:
% ++ p      p structure
%
% see also: core.prepare_initial_guess
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

function [ p ] = prepare_initial_object( p )
import utils.verbose 
import utils.crop_pad 
import utils.interpolateFT 


if p.model_object
    switch p.model.object_type
        case 'rand'
            verbose(2,'Using random object as initial guess.')
        case 'amplitude'
            % create dummy object which will be overwritten once the
            % prepared data is available
            assert(p.fourier_ptycho, 'An initial guess based on the prepared data is only suited for Fourier ptychography.')
    end
    for obnum = 1:p.numobjs 
        p.object{obnum} = (1+1i*1e-6*rand([p.object_size(obnum,:) p.object_modes])).*ones([p.object_size(obnum,:) p.object_modes]);
    end
else
    if isfield(p, 'initial_iterate_object')
        warning('Loading initial object guess from file given by p.initial_iterate_object_file.')
    end
    verbose(2,'Using loaded object as initial guess.')

    if numel(p.initial_iterate_object_file) ~= p.numobjs
        verbose(2,'Number of initial iterate files and number of objects does not match')
        for ii=numel(p.initial_iterate_object_file):p.numobjs
            p.initial_iterate_object_file{ii} = p.initial_iterate_object_file{end};
        end
    end
    
    % make a bit smarter the use of initial_iterate_object_file and allow
    % some automatic patten filling + file search 
    for obnum = unique(p.share_object_ID)
        % if string allows it, fill in the scan numbers 
        p.initial_iterate_object_file{obnum} = sprintf(p.initial_iterate_object_file{obnum}, p.scan_number(obnum)); 
        if contains(p.initial_iterate_object_file{obnum}, '*') % if string contains wild character *, try to find the file 
           fpath = dir(p.initial_iterate_object_file{obnum}) ; 
           if isempty(fpath)
               warning('No file corresponding to pattern %s was found, using random initial guess', p.initial_iterate_object_file{obnum})
               p.object{obnum} = (1+1i*1e-6*rand([p.object_size(obnum,:) p.object_modes])).*ones([p.object_size(obnum,:) p.object_modes]);
               p.initial_iterate_object_file{obnum} = []; 
               continue
           elseif length(fpath) > 1
               warning('Too many files corresponding to pattern %s were found, using the last', p.initial_iterate_object_file{obnum})
               fpath = fpath(end); 
           end
           p.initial_iterate_object_file{obnum} = [fpath.folder,'/',fpath.name]; 
        end
    end
    
    % load data from disk
    for ii = unique(p.share_object_ID) % avoid loading datasets twice
        if isempty(p.initial_iterate_object_file{ii})
            continue
        end
        if ~exist(p.initial_iterate_object_file{ii}, 'file')
            error(['Did not find initial iterate: ' p.initial_iterate_object_file{ii}])
        end
       
        verbose(2,'Loading object %d from: %s',ii,p.initial_iterate_object_file{ii})
        S = io.load_ptycho_recons(p.initial_iterate_object_file{ii});
        object = double(S.object);
        
        % reinterpolate to the right pixel size 
        if isfield(S, 'p') && any(S.p.dx_spec ~= p.dx_spec)
            verbose(3, 'Warning: Reinterpolate loaded object to new pixels size')
            object = interpolateFT(object, ceil(size(object(:,:,1)).*S.p.dx_spec./p.dx_spec));
        end

        %%% check the object size
        if ~isequal(size(squeeze(object(:,:,1))), squeeze(p.object_size(ii,:)))
            % if the loaded dataset does not have the expected object size,
            % crop/pad it to p.object_size
            verbose(3, 'Warning: Object taken from file %s does not have the expected size of %d x %d.', ...
            p.initial_iterate_object_file{ii}, p.object_size(ii,1), ...
            p.object_size(ii,2))
            p.object{ii} = crop_pad(object, p.object_size(ii,:));
        else
            % the the object sizes are the same, just copy everything
            % to p.object
            p.object{ii} = object;
        end

        % now let's check the object modes
        mode_diff = p.object_modes-size(object,3);
        if mode_diff > 0
            % add (random) object modes
            p.object{ii}(:,:,size(object,3)+1:p.object_modes,:) = (1+1i*1e-6*rand([p.object_size(ii,:) mode_diff])).*ones([p.object_size(ii,:) mode_diff]);
        elseif mode_diff < 0
            % remove object modes
            p.object{ii}(:,:,p.object_modes+1:size(object,3),:) = [];
        end

    end
end

end
