%PTYCHO_EXIT 
%onCleanup function for core.ptycho_recons
% 
% ** p          p structure
%
% see also: core.ptycho_recons

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

function ptycho_exit(p)
import utils.*

if isfield(p.getReport, 'crashed') && p.getReport.crashed
    fprintf('\n\n\n###############################################################\n')
    fprintf('########################## PONG! ##############################\n')
    fprintf('###############################################################\n')
    disp([p.getReport.ME.getReport '\n\n\n']);

    fprintf('Reconstruction stopped! Set verbose level to >3 for debugging. \n')
    if ~p.getReport.completed
        if ~isempty(p.io.phone_number) && p.io.send_crashed_recon_SMS
            io.sendSMS(p.io.phone_number, sprintf('Failed to reconstruct scan %s.', num2str(p.scan_number)));
        end
        if isfield(p.queue, 'file_this_recons')
            try
                fprintf('Moving queue file file back to %s.\n', fullfile(p.queue.path, p.queue.file_this_recons))
                
                % get log file name
                [~, ~, fext] = fileparts(p.queue.file_this_recons);
                log_dir = fullfile(p.queue.path, 'failed', 'log');
                if ~exist(log_dir, 'dir')
                    mkdir(log_dir)
                end
                log_file = fullfile(log_dir, strrep(p.queue.file_this_recons, fext, '.log'));
                
                % check if log file extists and update its content; move
                % queue file back to in_progess 
                if exist(log_file, 'file')
                    fid = fopen(log_file);
                    log_line = fgetl(fid);
                    fclose(fid);
                    log_int = strtrim(strsplit(log_line, ':'));
                    log_int = log_int{end};
                    log_int = str2double(log_int);
                    if ~isfield(p, 'queue_max_attempts')
                        p.queue.max_attempts = 5;
                        fprintf('Code crashed before parsing p.queue.max_attempts.\n')
                    end
                    if log_int >= p.queue.max_attempts
                        io.movefile_fast(fullfile(p.queue.path,'in_progress', p.queue.file_this_recons),fullfile(p.queue.path, 'failed', p.queue.file_this_recons))
                        fid = fopen(log_file, 'w');
                        fprintf(fid, [p.getReport.ME.getReport '\n\n\n']);
                        fprintf('Failed more than %u times. Moving file to ''failed''.\n', p.queue.max_attempts);
                        fclose(fid);
                        if ~isempty(p.io.phone_number) && p.io.send_failed_scans_SMS
                            io.sendSMS(p.io.phone_number, sprintf('Failed to reconstruct scan %s. I will move it to "failed".', num2str(p.scan_number)), 'sleep', p.SMS_sleep, 'logfile', fullfile(log_dir, 'sendSMS.log'));
                        end
                    else    
                        io.movefile_fast(fullfile(p.queue.path,'in_progress', p.queue.file_this_recons),fullfile(p.queue.path, p.queue.file_this_recons));
                        fid = fopen(log_file, 'w');
                        fprintf(fid, 'failed attempts: %u\n\n', log_int+1);
                        fprintf(fid, [p.getReport.ME.getReport '\n\n\n']);
                        fclose(fid);
                    end
                else
                    fid = fopen(log_file, 'w');
                    fprintf(fid, 'failed attempts: 1');
                    fclose(fid);
                    io.movefile_fast(fullfile(p.queue.path,'in_progress', p.queue.file_this_recons),fullfile(p.queue.path, p.queue.file_this_recons))
                end
                    
            catch
                fprintf('Failed to move file back to queue search path.\n')
            end
        end
        if isfield(p.queue, 'lockfile') 
            if isempty(p.queue.lockfile)
                if verbose > 2
                    p.queue.lockfile = false;
                else
                    p.queue.lockfile = true;
                end
            end
            if p.queue.lockfile
                if isempty(p.save_path{1})
                    try
                        for ii = 1:length(p.scan_number)
                            p.scan_str{ii} = sprintf(p.scan_string_format, p.scan_number(ii));        % Scan string
                        end
                        p = core.ptycho_prepare_paths(p);
                    catch
                        fprintf('Could not find lock file. \n')
                    end
                end
                for ii=1:length(p.save_path)
                    lock_filename = [p.save_path{ii} '/' p.run_name '_lock'];
                    if exist(lock_filename, 'file')
                        try
                            unix(['rm ' lock_filename]);
                            fprintf('Removing lock file %s\n',lock_filename)
                        catch
                            fprintf('Removing lock file %s failed\n',lock_filename)
                        end
                    end
                end
            end
        end
        if isfield(p.queue, 'remote_recons') && p.queue.remote_recons
            keyboard
        end
        
    end
    
    fprintf('Pausing for 5 seconds.\n')
    fprintf('###############################################################\n\n')
    pause(5);

elseif ~p.getReport.completed
    if isfield(p, 'remote_failed') && p.remote_failed
        try
            io.movefile_fast(fullfile(p.queue.path,'in_progress', p.queue.file_this_recons),fullfile(p.queue.path, 'failed', p.queue.file_this_recons));
        catch
            fprintf('Failed to move file to failed.\n')
        end
    elseif isfield(p.queue, 'file_this_recons')
        try
            fprintf('Reconstruction stopped, moving file to %s.\n', fullfile(p.queue.path, p.queue.file_this_recons))
            io.movefile_fast(fullfile(p.queue.path,'in_progress', p.queue.file_this_recons),fullfile(p.queue.path, p.queue.file_this_recons))
        catch
            fprintf('Failed to move file back to queue search path.\n')
        end
    end
    if isfield(p.queue, 'lockfile') 
        if isempty(p.queue.lockfile)
            if verbose > 2
                p.queue.lockfile = false;
            else
                p.queue.lockfile = true;
            end
        end
        if p.queue.lockfile
            if isempty(p.save_path{1})
                try
                    for ii = 1:length(p.scan_number)
                        p.scan_str{ii} = sprintf(p.scan_string_format, p.scan_number(ii));        % Scan string
                    end
                    p = core.ptycho_prepare_paths(p);
                catch
                    fprintf('Could not find lock file. \n')
                end
            end
            if isfield(p, 'run_name')  && ~isempty(p.run_name)
            for ii=1:length(p.save_path)
                lock_filename = [p.save_path{ii} '/' p.run_name '_lock'];
                if exist(lock_filename, 'file')
                    try
                        delete(lock_filename);
                        fprintf('Removing lock file %s\n',lock_filename)
                    catch
                        fprintf('Removing lock file %s failed\n',lock_filename)
                    end
                end
            end
            end
        end
    end
    if isfield(p, 'remote_file_this_recons')
        try
            fprintf('Removing remote file.\n')
            if exist(p.queue.remote_file_this_recons, 'file')
                delete(p.queue.remote_file_this_recons)
            end
            
            [~, this_file] = fileparts(p.queue.remote_file_this_recons);
            if p.queue.isreplica
                system(['touch ' fullfile(p.queue.remote_path, [this_file '.crash'])]);
            end
            
            this_file = [this_file '.mat'];
            
            if exist(fullfile(p.queue.remote_path, 'in_progress', this_file), 'file')
                delete(fullfile(p.queue.remote_path, 'in_progress', this_file));
            end
            if exist(fullfile(p.queue.remote_path, 'done', this_file), 'file')
                delete(fullfile(p.queue.remote_path, 'done', this_file));
            end
            if exist(fullfile(p.queue.remote_path, 'done', this_file), 'file')
                delete(fullfile(p.queue.remote_path, 'done', this_file));
            end

        catch
            fprintf('Failed to remove remote file.\n')
        end
    end
    
    if ~isempty(p.io.phone_number) && p.io.send_crashed_recon_SMS
        io.sendSMS(p.io.phone_number, sprintf('Failed to reconstruct scan %s.', num2str(p.scan_number)));
    end
end

pid = [p.ptycho_matlab_path './utils/.tmp_procID/proc_' num2str(feature('getpid')) '.dat'];
if exist(pid, 'file')
    delete(pid)
end

if ~isempty(p.io.phone_number) && p.io.send_finished_recon_SMS && p.getReport.completed
    io.sendSMS(p.io.phone_number, sprintf('Finished reconstructing scan %s.', num2str(p.scan_number)));
end

try
    verbose(struct('prefix', {[]}))
catch
end

end
