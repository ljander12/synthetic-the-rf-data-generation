%% Direct SWS reconstruction from k-Wave ground-truth displacement

clearvars -except disp_gt shear_params prf array_length scene_depth ...
    prb_center prb_theta

% Load dataset if variables are not already in workspace
if ~exist('disp_gt','var')
    load('high_contrast_displacement_with_quality.mat', ...
        'disp_gt', ...
        'shear_params', ...
        'prf', ...
        'array_length', ...
        'scene_depth', ...
        'prb_center', ...
        'prb_theta');
end

%% Recreate phantom grid
[kgrid, medium] = create_phantom(shear_params);

%% Ground-truth image grid
x_gt = (-array_length/2 : kgrid.dx : array_length/2).';
z_gt = (5e-3 : kgrid.dz : scene_depth).';

fprintf('disp_gt size:       %s\n', mat2str(size(disp_gt)));
fprintf('Expected grid size: [%d %d]\n', numel(z_gt), numel(x_gt));

assert(size(disp_gt,1) == numel(z_gt));
assert(size(disp_gt,2) == numel(x_gt));

scan_gt = uff.linear_scan( ...
    'x_axis', x_gt, ...
    'z_axis', z_gt);

%% Build elastography object directly from k-Wave GT displacement
elasto_gt = uff.elastography();

elasto_gt.PRF = prf;
elasto_gt.target_frequency = shear_params.source_freq;
elasto_gt.scan = scan_gt;

elasto_gt.autocorrelation_field = uff.data_field();

elasto_gt.autocorrelation_field.set_component( ...
    'u', ...
    reshape( ...
        double(disp_gt), ...
        size(disp_gt,1), ...
        size(disp_gt,2), ...
        1, ...
        size(disp_gt,3)));

%% Directional filtering setup
num_angles = 12;

elasto_gt.directional_filtering_angles = ...
    2*pi*linspace(0, 1 - 1/num_angles, num_angles);

%% Extract 200 Hz component
extractor = postprocess.elastography.extract_freq_component(elasto_gt);

extractor.target_frequency = shear_params.source_freq;
extractor.component = 'u';

extractor.run();

%% Same radial filter as baseline analysis
%% Disable radial filtering
elasto_gt.radial_filter = uff.spatial_filter( ...
    'ones', ...
    'radial');

fprintf('Radial filter type: %s\n', ...
    elasto_gt.radial_filter.type);
elasto_gt.angular_filter = uff.spatial_filter( ...
    'gaussian_charite', ...
    'angular');

elasto_gt.angular_filter.deltaTheta = ...
    median(diff(elasto_gt.directional_filtering_angles));

%% Apply filter
filter_operator = postprocess.elastography.apply_filter(elasto_gt);
filter_operator.run();

%% Phase-gradient SWS
pgm = postprocess.elastography.phase_gradient_method(elasto_gt);
pgm.run();

sws_from_gt = squeeze(double(elasto_gt.sws_map.u));

%% Create exact material SWS map on same GT grid
[TF, ~] = coor_transformation(prb_center, prb_theta);

[X, Z] = meshgrid(x_gt, z_gt);

local_points = [ ...
    X(:).'; ...
    zeros(1,numel(X)); ...
    Z(:).'; ...
    ones(1,numel(X))];

global_points = TF * local_points;

F = griddedInterpolant( ...
    { ...
        double(kgrid.x_vec(:)), ...
        double(kgrid.y_vec(:)), ...
        double(kgrid.z_vec(:)) ...
    }, ...
    double(medium.sound_speed_shear), ...
    'nearest', ...
    'none');

sws_material_gt = reshape( ...
    F( ...
        global_points(1,:).', ...
        global_points(2,:).', ...
        global_points(3,:).'), ...
    size(X));

%% Simple regional statistics
incl_mask = ...
    abs(sws_material_gt - shear_params.c_shear_incl) < 1e-6;

bkg_mask = ...
    abs(sws_material_gt - shear_params.c_shear_bkg) < 1e-6;

valid = isfinite(sws_from_gt);

incl_values = sws_from_gt(incl_mask & valid);
bkg_values  = sws_from_gt(bkg_mask & valid);

fprintf('\n========================================\n');
fprintf('DIRECT K-WAVE GT -> PHASE GRADIENT\n');
fprintf('========================================\n');

fprintf('\nInclusion\n');
fprintf('GT:     %.4f m/s\n', shear_params.c_shear_incl);
fprintf('Mean:   %.4f m/s\n', mean(incl_values));
fprintf('Median: %.4f m/s\n', median(incl_values));
fprintf('SD:     %.4f m/s\n', std(incl_values));

fprintf('\nBackground\n');
fprintf('GT:     %.4f m/s\n', shear_params.c_shear_bkg);
fprintf('Mean:   %.4f m/s\n', mean(bkg_values));
fprintf('Median: %.4f m/s\n', median(bkg_values));
fprintf('SD:     %.4f m/s\n', std(bkg_values));

%% Visualize
figure;
imagesc(x_gt*1e3, z_gt*1e3, sws_material_gt);
axis image;
set(gca,'YDir','normal');
colorbar;
xlabel('Lateral position [mm]');
ylabel('Depth [mm]');
title('Material ground-truth SWS');

figure;
imagesc(x_gt*1e3, z_gt*1e3, sws_from_gt);
axis image;
set(gca,'YDir','normal');
colorbar;
xlabel('Lateral position [mm]');
ylabel('Depth [mm]');
title('SWS reconstructed directly from k-Wave displacement');

figure;
imagesc(x_gt*1e3, z_gt*1e3, sws_from_gt - sws_material_gt);
axis image;
set(gca,'YDir','normal');
colorbar;
xlabel('Lateral position [mm]');
ylabel('Depth [mm]');
title('Direct GT reconstruction error');

%% Raw 200 Hz phase-gradient diagnostic

U = squeeze(double(elasto_gt.unfiltered_phase_image_xz.u));

dx = mean(diff(x_gt));
dz = mean(diff(z_gt));

omega = 2*pi*shear_params.source_freq;

% Centered complex phase differences
kx = nan(size(U));
kz = nan(size(U));

kx(:,2:end-1) = ...
    angle(conj(U(:,1:end-2)) .* U(:,3:end)) / (2*dx);

kz(2:end-1,:) = ...
    angle(conj(U(1:end-2,:)) .* U(3:end,:)) / (2*dz);

k_raw = sqrt(kx.^2 + kz.^2);

sws_raw = omega ./ k_raw;

% Reject very low-amplitude pixels where phase is unstable
amp = abs(U);
amp_mask = amp > 0.1 * max(amp(:));

valid_raw = ...
    isfinite(sws_raw) & ...
    amp_mask & ...
    sws_raw > 0 & ...
    sws_raw < 10;

incl_raw = sws_raw(incl_mask & valid_raw);
bkg_raw  = sws_raw(bkg_mask & valid_raw);

fprintf('\n========================================\n');
fprintf('RAW 200 Hz PHASE GRADIENT\n');
fprintf('========================================\n');

fprintf('\nInclusion\n');
fprintf('GT:     %.4f m/s\n', shear_params.c_shear_incl);
fprintf('Mean:   %.4f m/s\n', mean(incl_raw));
fprintf('Median: %.4f m/s\n', median(incl_raw));
fprintf('SD:     %.4f m/s\n', std(incl_raw));

fprintf('\nBackground\n');
fprintf('GT:     %.4f m/s\n', shear_params.c_shear_bkg);
fprintf('Mean:   %.4f m/s\n', mean(bkg_raw));
fprintf('Median: %.4f m/s\n', median(bkg_raw));
fprintf('SD:     %.4f m/s\n', std(bkg_raw));

figure;
imagesc(x_gt*1e3, z_gt*1e3, sws_raw);
axis image;
set(gca,'YDir','normal');
colorbar;
clim([0 7]);
xlabel('Lateral position [mm]');
ylabel('Depth [mm]');
title('Raw 200 Hz phase-gradient apparent SWS');

%% Compare the two inclusions separately

CC = bwconncomp(incl_mask);

fprintf('\n========================================\n');
fprintf('RAW PHASE GRADIENT - EACH INCLUSION\n');
fprintf('========================================\n');

for n = 1:CC.NumObjects

    mask_n = false(size(incl_mask));
    mask_n(CC.PixelIdxList{n}) = true;

    vals = sws_raw(mask_n & valid_raw);

    fprintf('\nInclusion %d\n', n);
    fprintf('Pixels: %d\n', nnz(mask_n));
    fprintf('Mean:   %.4f m/s\n', mean(vals));
    fprintf('Median: %.4f m/s\n', median(vals));
    fprintf('SD:     %.4f m/s\n', std(vals));

end