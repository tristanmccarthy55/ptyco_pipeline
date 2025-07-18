%OMNY Load positions from Orchestra scan file


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


function [ p ] = orchestra( p )
import beamline.*
import utils.*

if isempty(p.positions_file)
    error('OMNY positions file is not specified. Please check p.positions_file in your template.')
end

if ~isfield(p,'angular_correction_setup') || isempty(p.angular_correction_setup)
    error('p.angular_correction_setup is not specified. Please check p.angular_correction_setup in your template.')
end

if isfield(p,'omny_interferometer')
    error(' p.omny_interferometer is not supported, use p.angular_correction_setup')
end

if ~isfield(p.detector,'burst_frames')||isempty(p.detector.burst_frames)
    p.detector.burst_frames = 1;
end

switch lower(p.angular_correction_setup)
    case 'omny'
        p.   orchestra.laser_height=-10.0e-3;                                 % Height of horizontal laser beam on the sphere compared to pin tip (only for p.fromspec='opos_angle', 13.5e-3 for OMNI (not fully tested, better with opos than opos_angle), -10.0e-3 for OMNY)
        p.   orchestra.mirrdis=-9.0e-3;                                       % Distance mirror-pin tip (only for p.fromspec='opos_angle', 22.0e-3 for OMNI (not fully tested, better with opos than opos_angle), -9.0e-3 for OMNY)
        p.   orchestra.beam_separation=7.5e-3;                                % Distance mirror-pin tip (only for p.fromspec='opos_angle', 13.0e-3 for OMNI (not fully tested, better with opos than opos_angle), 7.5e-3 for OMNY)
        apply_correction = true; 
    case 'flomni'
        p.   orchestra.laser_height=-13.5e-3;                                  % Height of horizontal laser beam on the sphere compared to pin tip (only for p.fromspec='opos_angle', 13.5e-3 for OMNI (not fully tested, better with opos than opos_angle), -10.0e-3 for OMNY)
        p.   orchestra.mirrdis=-17.4e-3;                                       % Distance mirror-pin tip (only for p.fromspec='opos_angle', 22.0e-3 for OMNI (not fully tested, better with opos than opos_angle), -9.0e-3 for OMNY)
        p.   orchestra.beam_separation=-16e-3;                                 % Distance mirror-pin tip (only for p.fromspec='opos_angle', 13.0e-3 for OMNI (not fully tested, better with opos than opos_angle), 7.5e-3 for OMNY)
        apply_correction = true; 
    case {'lamni', 'none'}
        apply_correction = false; 
    otherwise
        error('Wrong  p.angular_correction_setup, choose from ''omny'', ''flomni'',''lamni'',''none'' ')
end



for ii = 1:length(p.scan_number)
    p.scan.is_cont = true;  % So that burst data is prepared normally rather than integrated
    if ~exist(sprintf(p.positions_file,p.scan_number(ii)), 'file' )
        error('Missing OMNY specs file  %s', sprintf(p.positions_file,p.scan_number(ii)))
    end
    out_orch = read_omny_pos(sprintf(p.positions_file,p.scan_number(ii)));
    if ~isfield(out_orch,'Average_y_st_fzp') ||  ~isfield(out_orch,'Average_x_st_fzp')
        out_orch.Average_y_st_fzp = out_orch.Average_y;
        out_orch.Average_x_st_fzp = out_orch.Average_x;
    end
    if ~isfield(out_orch, 'Average_rotz_st')
        apply_correction =false;
    end
    if ~apply_correction
        if isfield(out_orch, 'Average_y_st_fzp')
            positions_real = [out_orch.Average_y_st_fzp*1e-6 out_orch.Average_x_st_fzp*1e-6];
        else  % outdated position format
            positions_real = [out_orch.Average_y*1e-6 out_orch.Average_x*1e-6];
        end
    else
        deltax = p.orchestra.laser_height*out_orch.Average_rotz_st*1e-6/p.orchestra.beam_separation; % p.orchestra.beam_separation: separation between two laser beams for angular measurement
        % p.orchestra.laser_height: height of horizontal laser beam on the sphere compared to pin tip
        deltay = p.orchestra.mirrdis*out_orch.Average_rotz_st*1e-6/p.orchestra.beam_separation; % p.orchestra.beam_separation: separation between two laser beams for angular measurement
        % p.orchestra.mirrdis dist mirror-pin tip
        posx = out_orch.Average_x_st_fzp*1e-6 - deltax;
        posy = out_orch.Average_y_st_fzp*1e-6 - deltay;
        positions_real = [posy posx];
    end
    
    p.numpts(ii) = size(positions_real,1)*p.detector.burst_frames;
    
    positions_tmp = zeros(p.numpts(ii), 2);
    positions_tmp(:,1) = reshape(repmat(positions_real(:,1)',[p.detector.burst_frames 1]),[],1); 
    positions_tmp(:,2) = reshape(repmat(positions_real(:,2)',[p.detector.burst_frames 1]),[],1); 
    
    p.positions_real = [p.positions_real ; positions_tmp];

    
end

end

