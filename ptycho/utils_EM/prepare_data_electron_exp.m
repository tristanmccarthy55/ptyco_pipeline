% This script prepares data for PtychoShelves EM Hollow from clean EMPAD dataset,
% *.mat file provided in paper:
%  Yu Lei and Peng Wang, Hollow multi-slice electron ptychography for simultaneous 3D structural imaging and EELS in 4D-STEM, arXiv:2506.22352 (https://arxiv.org/abs/2506.22352)
% Yu Lei @ University of Warwick, 7/18/2025

%% Step 1: download sample data (rawdata.mat) from the link provided in
% PARADIM website to release, stay tuned: https://data.paradim.org/
% to start a data_dp.mat (3 dimension) is required.

clear;
%% Step 2: load data
data_dir = '../Hollow_Data_exp_30nm/'; %change this
fmat='sample_data_PrScO3.mat';  % *.mat file
for i=[2 12 26 40 50]
    scan_number = i; %Ptychoshelves needs
    % initial lize mask
    mask = ones(256, 256);

    if i==1
        r=-1;
        % radius
    else
        r = i / 2;
    end
    % center coordinate
    center_x = (256+1)/2;
    center_y = (256+1)/2;

    % mesh grid
    [x, y] = meshgrid(1:256, 1:256);

    % calculate distance to center
    dist_from_center = sqrt((x - center_x).^2 + (y - center_y).^2);

    % mask
    mask(dist_from_center <= r) = 0;
    %%
    Np_p = 256; % size of diffraction patterns used during reconstruction. can also pad to 256
    ADU=580; % counts per electron
    alpha0=21.4; % mrad
    rbf=26; % radius of center disk in pixels
    voltage=300; % kev
    df=-200; % estimated defocuse
    scanstep=0.41; % scan step size in angstrom
    rot_ang=0; % relative angle between scan and diffraction
    %% Step 3: go back to .../fold_slice/ptycho and pre-process data
    Np_p=[Np_p,Np_p];
    fdata=fullfile(data_dir,fmat);
    
    if i==2
        load(fdata);
        dp_o=dp;
        dp_o=dp_o/ADU;
    end
    dp=dp_o.*mask;
    npx=64;
    npy=64;

    crop_idx=[1,64,1,64]; % start from smaller data

    dp=reshape(dp,Np_p(1),Np_p(2),[]);
    Itot=mean(squeeze(sum(sum(dp_o,1),2))); %need this for normalizting initial probe

    % calculate pxiel size (1/A) in diffraction plane
    lambda = 12.3986./sqrt((2*511.0+voltage).*voltage);
    dk=alpha0/1e3/rbf/lambda; %%% PtychoShelves script needs this %%%

    %% Step 4: save CBED in a .hdf5 file (needed by Ptychoshelves)
    save_dir = fullfile(data_dir,num2str(scan_number,'%02d'));
    mkdir(save_dir)

    save(fullfile(save_dir,'data_dp.mat'),'dp','-v7.3');
    save(fullfile(save_dir,'mask.mat'),'mask','-v7.3');

    %% Step 5: prepare initial probe
    dx=1/Np_p(1)/dk; %% pixel size in real space (angstrom)
    cs = 0;
    probe=generateProbeFunction(dx,Np_p(1),0,0,df,cs,1,voltage,alpha0,0);
    probe=probe/sqrt(sum(sum(abs(probe.^2))))*sqrt(Itot)/sqrt(Np_p(1)*Np_p(2));
    probe=single(probe);
    % add parameters for PtychoShelves_electron
    p = {};
    p.binning = false;
    p.detector.binning = false;

    save(fullfile(save_dir,'probe_initial.mat'),'probe','p')
    copyfile(strcat(mfilename('fullpath'),'.m'),save_dir);

    %% prepare probe position
    probe_positions_0=position_generate(npx,npy,scanstep,scanstep, rot_ang);
    hdf5write(fullfile(save_dir,'data_position.hdf5'),'/probe_positions_0',probe_positions_0);
end
