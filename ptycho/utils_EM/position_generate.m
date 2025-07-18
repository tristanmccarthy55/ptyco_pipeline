function probe_positions_0=position_generate(npx,npy,scanStepSize_x,scanStepSize_y, rot_ang)
% generate raster scan positions from scan step
% Zhen Chen @ Cornell University, 3/27/2021
% inputs: scan points (npx,npy)
%         scan step size in angstrom
%         relative orientation angle between scan direction and diffraction
%         direction, usually 0 degree 

[xx,yy]=scan_position_rot(npx,npy,scanStepSize_x,scanStepSize_y,rot_ang);% in Angstrom
probe_positions_0=[xx(:),yy(:)];
probe_positions_0=single(probe_positions_0);

%%
function [ppX_rot,ppY_rot]=scan_position_rot(N_scan_x,N_scan_y,scanStepSize_x,scanStepSize_y,rot_ang)
% scan positions after rotation
    ppx = linspace(-floor(N_scan_x/2),ceil(N_scan_x/2)-1,N_scan_x)*scanStepSize_x;
    ppy = linspace(-floor(N_scan_y/2),ceil(N_scan_y/2)-1,N_scan_y)*scanStepSize_y;
    [ppX,ppY] = meshgrid(ppx,ppy);

    ppY_rot = ppX*(-sind(rot_ang)) + ppY*cosd(rot_ang);
    ppX_rot = ppX*cosd(rot_ang) + ppY*sind(rot_ang);

end
end
