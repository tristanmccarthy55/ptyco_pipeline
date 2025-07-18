%hdf5_pos loads scan positions from hdf5 files 
%Written by Yi Jiang & Zhen Chen

function [ p ] = hdf5_pos( p )

for ii = 1:p.numscans
    positions_real = zeros(0,2); 

    switch p.scan.type
        case 'custom'
            if isempty(p.scan.custom_positions_source) %guess the position file name from base path
                pos_file = strcat(p.base_path,sprintf(p.scan.format, p.scan_number(ii)),'/data_position.hdf5');                
            end
            if exist(pos_file,'file')
                
            	ps = h5read(pos_file,'/probe_positions_0');
                ppX = ps(:,1)*1e-10; % Angstrom to meter
                ppY = ps(:,2)*1e-10;
                positions_real = zeros(length(ppX),2); 

                positions_real(:,1) = -ppY(:);
                positions_real(:,2) = -ppX(:);                
            else
            	error('Could not find function or data file %s', pos_file);
            end
        otherwise
            error('Unknown scan type %s.', p.scan.type);
    end
    
    p.numpts(ii) = size(positions_real,1);
    p.positions_real = [p.positions_real ; positions_real]; %append position
end
    
end

