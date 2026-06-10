clear;
clc;
close all;
pixels = 1024; %#ok<NASGU>



min_dist = 4;       % voxel distance for NMS

label_map = {'Pr', 'Sc', 'O'};


load('Niter60.mat');

obj = imresize3(obj, [size(obj,1)*2, size(obj,2)*2, size(obj,3)*2]);
x_step = 0.5*0.061524;      % Å
z_step = 0.5*1.0123;     % Å

obj_ref=max(obj,[],3);


x_len = size(obj,1);
y_len = size(obj,2);

xy_tol_ref = 0.3;   


BW = imregionalmax(obj_ref);

[row,col] = find(BW);

peak_val = obj_ref(BW);

label = kmeans(peak_val,4,'Replicates',20);

mean_val = zeros(4,1);
for i = 1:4
    mean_val(i) = mean(peak_val(label==i));
end

[~,order] = sort(mean_val,'descend');

Pr_cluster = order(1);
Sc_cluster = order(2);
O_cluster  = order(3);

Pr_pos = [row(label==Pr_cluster), col(label==Pr_cluster)];
Sc_pos = [row(label==Sc_cluster), col(label==Sc_cluster)];
O_pos  = [row(label==O_cluster),  col(label==O_cluster)];

Pr_pos_ref = Pr_pos;
Sc_pos_ref = Sc_pos;
O_pos_ref  = O_pos;

D = pdist2(O_pos,Pr_pos);

[minD,idxPr] = min(D,[],2);

O_val = obj_ref(sub2ind(size(obj_ref),...
    O_pos(:,1),O_pos(:,2)));

Pr_val_local = obj_ref(sub2ind(size(obj_ref),...
    Pr_pos(idxPr,1),Pr_pos(idxPr,2)));

keep = ~( minD<20 & O_val < 0.5*Pr_val_local );

O_pos = O_pos(keep,:);
D = pdist2(O_pos,Sc_pos);

[minD,idxSc] = min(D,[],2);
O_val = obj_ref(sub2ind(size(obj_ref),...
    O_pos(:,1),O_pos(:,2)));

Sc_val_local = obj_ref(sub2ind(size(obj_ref),...
    Sc_pos(idxSc,1),Sc_pos(idxSc,2)));

keep = ~( minD<20 & O_val < 0.6*Sc_val_local );

O_pos = O_pos(keep,:);
search_radius.Pr = 6;   % pixels
search_radius.Sc = 6;
search_radius.O  = 6;

Pr_pos = move_to_local_brightest(obj_ref, Pr_pos, search_radius.Pr);
Sc_pos = move_to_local_brightest(obj_ref, Sc_pos, search_radius.Sc);
O_pos  = move_to_local_brightest(obj_ref, O_pos,  search_radius.O);
Pr_val = sample_obj_ref_at_xy(obj_ref, Pr_pos);
Sc_val = sample_obj_ref_at_xy(obj_ref, Sc_pos);
O_val  = sample_obj_ref_at_xy(obj_ref, O_pos);
thresholds.Pr = 0.5*mean(Pr_val, 'omitnan');
thresholds.Sc =0.7*mean(Sc_val, 'omitnan');
thresholds.O  = 0.4*mean(O_val,  'omitnan');
fprintf('\nDynamic thresholds:\n');
fprintf('Pr = %.4f\n', thresholds.Pr);
fprintf('Sc = %.4f\n', thresholds.Sc);
fprintf('O  = %.4f\n\n', thresholds.O);

%% ---- Detect atoms ----
[detected_positions_centered, detected_symbols, detected_values, accepted_atoms] = ...
    detect_atoms_no_reference(obj, x_step, z_step, thresholds, min_dist, label_map, ...
    Pr_pos, Sc_pos, O_pos);


%% ===== Draw detected atoms =====

detected_positions = detected_positions_centered;

fig = figure('Color','w');
hold on;

target_symbols  = ["Pr","Sc","O"];
display_symbols = ["Pr","Sc","O"];

colors = lines(numel(target_symbols));
colors = [colors(2,:);colors(3,:);colors(1,:)];

for kk = 1:numel(target_symbols)

    idx = detected_symbols == target_symbols(kk);

    if ~any(idx)
        continue;
    end

    scatter3( ...
        detected_positions(idx,1), ...
        detected_positions(idx,2), ...
        detected_positions(idx,3), ...
        50, ...
        colors(kk,:), ...
        'filled', ...
        'DisplayName', char(display_symbols(kk)));

end
xlim([-4 4])
ylim([-4 4])
zlim([-40 40])
xlabel('X (Å)');
ylabel('Y (Å)');
zlabel('Z (Å)');

box on;
grid off;
legend('Location','best','FontSize',14);

view(3);
daspect([1 1 10]);

saveas(fig,'First_detected_slice_Pr_Sc_O.png');
%% ============================================================
%  Local functions
%% ============================================================

function [detected_positions_centered, detected_symbols, detected_values, accepted_atoms] = ...
    detect_atoms_no_reference(obj, x_step, z_step, thresholds, min_dist, label_map, ...
    Pr_pos, Sc_pos, O_pos)
[x_len, y_len, z_len] = size(obj);

% fprintf('Detecting 3D local maxima...\n');
h = 0.01;
maxima_mask = imregionalmax(imhmax(obj,h));
[x_idx, y_idx, z_idx] = ind2sub(size(obj), find(maxima_mask));
values = obj(maxima_mask);

[values_sorted, idx_sort] = sort(values, 'descend');
coords_sorted = [x_idx(idx_sort), y_idx(idx_sort), z_idx(idx_sort)];

% fprintf('Applying thresholding and 3D NMS...\n');
accepted_atoms = [];  % [x_idx, y_idx, z_idx, value, label]

for jj = 1:size(coords_sorted, 1)
    this_pos = coords_sorted(jj, :);
    v = values_sorted(jj);

    if v > thresholds.Pr
        label = 1;      % Pr
    elseif v > thresholds.Sc
        label = 2;      % Sc
    elseif v > thresholds.O
        label = 3;      % O
    else
        continue;
    end

    if isempty(accepted_atoms)
        accepted_atoms = [this_pos, v, label]; %#ok<AGROW>
    else
        dists = vecnorm(accepted_atoms(:, 1:3) - this_pos, 2, 2);
        if all(dists > min_dist)
            accepted_atoms(end+1, :) = [this_pos, v, label]; %#ok<AGROW>
        end
    end
end
%% ===== Remove atoms inconsistent with column type =====

column_tol = 3;   % pixel

col_xy = [Pr_pos; Sc_pos; O_pos];

col_label = [ ...
    ones(size(Pr_pos,1),1); ...
    2*ones(size(Sc_pos,1),1); ...
    3*ones(size(O_pos,1),1)];

if ~isempty(col_xy) && ~isempty(accepted_atoms)

    atom_xy = accepted_atoms(:,1:2);

    D = pdist2(atom_xy, col_xy);
    [minD, idxCol] = min(D, [], 2);

    nearest_col_label = col_label(idxCol);

    keep = (minD <= column_tol) & ...
        (accepted_atoms(:,5) == nearest_col_label);

    accepted_atoms = accepted_atoms(keep,:);

end
if isempty(accepted_atoms)
    detected_positions_centered = zeros(0, 3);
    detected_symbols = strings(0, 1);
    detected_values = zeros(0, 1);
    return;
end

%
x_phys = (accepted_atoms(:,1) - (x_len+1)/2) * x_step;
y_phys = (accepted_atoms(:,2) - (y_len+1)/2) * x_step;
z_phys = (accepted_atoms(:,3) - (z_len+1)/2) * z_step;

% detected_positions_centered = [x_phys, y_phys, z_phys];
detected_symbols = strings(size(accepted_atoms,1), 1);
for jj = 1:size(accepted_atoms,1)
    detected_symbols(jj) = string(label_map{accepted_atoms(jj,5)});
end
detected_values = accepted_atoms(:,4);
%% ===== Refine detected positions by local center of mass =====

refine_ranges.Pr = [0.2, 0.2, 2];   % [xy_x, xy_y, z] in Å
refine_ranges.Sc = [0.2, 0.2, 2];
refine_ranges.O  = [0.2, 0.2, 2];

n_iter = 3;

detected_positions_centered = refine_positions_center_of_mass( ...
    obj, accepted_atoms(:,1:3), detected_symbols, ...
    x_step, z_step, refine_ranges, n_iter);
end


function refined_positions_centered = ...
    refine_positions_center_of_mass( ...
    obj, peak_idx, symbols, ...
    x_step, z_step, refine_ranges, n_iter)

if nargin < 7
    n_iter = 1;
end

current_idx = peak_idx;

for iter = 1:n_iter

    refined_positions_centered = ...
        refine_positions_center_of_mass_once( ...
        obj, current_idx, symbols, ...
        x_step, z_step, refine_ranges);

    current_idx(:,1) = ...
        refined_positions_centered(:,1)/x_step + (size(obj,1)+1)/2;

    current_idx(:,2) = ...
        refined_positions_centered(:,2)/x_step + (size(obj,2)+1)/2;

    current_idx(:,3) = ...
        refined_positions_centered(:,3)/z_step + (size(obj,3)+1)/2;

end
end
function refined_positions_centered = refine_positions_center_of_mass_once( ...
    obj, peak_idx, symbols, x_step, z_step, refine_ranges)

[x_len, y_len, z_len] = size(obj);
n_atoms = size(peak_idx,1);

refined_positions_centered = zeros(n_atoms,3);

for ii = 1:n_atoms

    sym = symbols(ii);

    if isfield(refine_ranges, char(sym))
        range_A = refine_ranges.(char(sym));
    else
        range_A = [0.5, 0.5, 0.6];
    end

    rx = max(1, round(range_A(1) / x_step));
    ry = max(1, round(range_A(2) / x_step));
    rz = max(1, round(range_A(3) / z_step));

    cx = round(peak_idx(ii,1));
    cy = round(peak_idx(ii,2));
    cz = round(peak_idx(ii,3));

    x1 = max(1, cx-rx); x2 = min(x_len, cx+rx);
    y1 = max(1, cy-ry); y2 = min(y_len, cy+ry);
    z1 = max(1, cz-rz); z2 = min(z_len, cz+rz);

    local_obj = obj(x1:x2, y1:y2, z1:z2);

    
    local_obj = local_obj - min(local_obj(:));
    local_obj(local_obj < 0) = 0;

    if sum(local_obj(:)) <= 0
        refined_idx = [cx, cy, cz];
    else
        [X,Y,Z] = ndgrid(x1:x2, y1:y2, z1:z2);

        w = local_obj(:);

        x_com = sum(X(:).*w) / sum(w);
        y_com = sum(Y(:).*w) / sum(w);
        z_com = sum(Z(:).*w) / sum(w);

        refined_idx = [x_com, y_com, z_com];
    end

    x_phys = (refined_idx(1) - (x_len+1)/2) * x_step;
    y_phys = (refined_idx(2) - (y_len+1)/2) * x_step;
    z_phys = (refined_idx(3) - (z_len+1)/2) * z_step;

    refined_positions_centered(ii,:) = [x_phys, y_phys, z_phys];

end

end


function vals = sample_obj_ref_at_xy(obj_ref, pix_xy)

if isempty(pix_xy)
    vals = nan;
    return;
end

idx = sub2ind(size(obj_ref), pix_xy(:,1), pix_xy(:,2));
vals = obj_ref(idx);

end
function pos_new = move_to_local_brightest(img, pos, radius)

pos_new = pos;

[nx, ny] = size(img);

for ii = 1:size(pos,1)

    cx = pos(ii,1);
    cy = pos(ii,2);

    x1 = max(1, cx-radius);
    x2 = min(nx, cx+radius);
    y1 = max(1, cy-radius);
    y2 = min(ny, cy+radius);

    local_img = img(x1:x2, y1:y2);

    [~, ind] = max(local_img(:));
    [dx, dy] = ind2sub(size(local_img), ind);

    pos_new(ii,1) = x1 + dx - 1;
    pos_new(ii,2) = y1 + dy - 1;

end

end