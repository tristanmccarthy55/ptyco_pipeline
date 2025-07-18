function [output] = generateResultDir(param, resultDir, extra)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
output = strcat(resultDir,'/',param.method,'_',param.opt_errmetric,'_p',num2str(param.probe_modes),'_g',num2str(param.grouping));

if strcmp(param.method, 'MLc') && param.accelerated_gradients_start < param.number_iterations    
    output = strcat(output,'_accGrad',num2str(param.accelerated_gradients_start));
end

if strcmp(param.method, 'MLc') && param.momentum ~=0  
    output = strcat(output,'_mom',num2str(param.momentum));
end

if param.beta_object < 1
    output = strcat(output,'_betaObj',num2str(param.beta_object));
end
if param.beta_probe < 1
    output = strcat(output,'_betaProb',num2str(param.beta_probe));
end

if isfield(param,'beta_LSQ')
    output = strcat(output,'_betaLSQ',num2str(param.beta_LSQ));
end
if param.delta_p ~= 0.1
    output = strcat(output,'_LSQdamping',num2str(param.delta_p));
end

if param.probe_position_search < param.number_iterations
    output = strcat(output,'_pc_Niter',num2str(param.probe_position_search));
    if ~isempty(param.probe_geometry_model)
        output = strcat(output,'_geo_model');
        for i=1:length(param.probe_geometry_model)
            output = strcat(output,'_',param.probe_geometry_model{i});
        end
    end
    if param.probe_position_error_max < inf
        output = strcat(output,'_maxPosError',num2str(param.probe_position_error_max/1e-9),'nm');
    end
end
if param.detector_rotation_search < param.number_iterations
    output = strcat(output,'_detRotSearchNiter',num2str(param.detector_rotation_search));
end
if param.detector_scale_search < param.number_iterations
    output = strcat(output,'_detScaleSearchNiter',num2str(param.detector_scale_search));
end
if param.probe_fourier_shift_search < param.number_iterations
    output = strcat(output,'_fpc_Niter',num2str(param.probe_fourier_shift_search));
end

if param.background>0
    output = strcat(output,'_bg',num2str(param.background));
end

if param.reg_mu>0
    output = strcat(output,'_regSmooth',num2str(param.reg_mu));
end

if param.positivity_constraint_object>0
    output = strcat(output,'_posConObj',num2str(param.positivity_constraint_object));
end


if param.variable_probe
    output = strcat(output,'_vp',num2str(param.variable_probe_modes));
    if param.variable_probe_smooth>0
        output = strcat(output,'_smooth',num2str(param.variable_probe_smooth));
    end
end
if param.variable_intensity
    output = strcat(output,'_vi');
end
if param.apply_multimodal_update
    output = strcat(output,'_mm');
end
%if isfield(param, 'extension') && param.extension{} = {'fly_scan'}
if is_used(param, 'fly_scan')
    output = strcat(output,'_afs');
end

%{
if any(param.det_bad_pixels(:))
    output = strcat(output,'_badPixels');
end
%}
if nargin==3
    output = strcat(output,extra,'/');
else
    output = strcat(output,'/');
end

end

