clear; clc; close all;

%% Setup

repo_dir = fileparts(mfilename('fullpath'));
cd(repo_dir);

toolbox_dir = ...
    'C:\Users\louis\Documents\GitHub\MasterOppgave\elastography\ElastographyToolbox';

addpath(genpath(toolbox_dir));

%% Load homogeneous ground-truth simulation

load('homogeneous_ground_truth.mat');
%% Ground-truth displacement

gt_disp = cat(3, zeros(size(disp_gt(:,:,1))), disp_gt);

d_gt = 1 / shear_params.source_freq / 5;

x_gt = -array_length/2 + ...
    (0:size(gt_disp,2)-1) * d_gt;

z_gt = 5e-3 + ...
    (0:size(gt_disp,1)-1) * d_gt;

gt_scan = uff.linear_scan( ...
    'x_axis', x_gt, ...
    'z_axis', z_gt);

%% Create elastography object

elasto = uff.elastography();

elasto.PRF = prf;
elasto.target_frequency = shear_params.source_freq;
elasto.scan = gt_scan;

%% Add ground-truth displacement field

gt_disp = reshape(gt_disp, ...
    size(gt_disp,1), ...
    size(gt_disp,2), ...
    1, ...
    size(gt_disp,3));

elasto.autocorrelation_field = uff.data_field();
elasto.autocorrelation_field.set_component('u', gt_disp);

fprintf('PRF: %.1f Hz\n', elasto.PRF);
fprintf('Target frequency: %.1f Hz\n', elasto.target_frequency);
fprintf('GT displacement size: %s\n', mat2str(size(gt_disp)));
%% Check

fprintf('PRF: %.1f Hz\n', elasto.PRF);
fprintf('Target frequency: %.1f Hz\n', elasto.target_frequency);
fprintf('GT displacement size: %s\n', mat2str(size(gt_disp)));
%% Check

fprintf('PRF: %.1f Hz\n', elasto.PRF);
fprintf('Target frequency: %.1f Hz\n', elasto.target_frequency);
fprintf('GT displacement size: %s\n', mat2str(size(gt_disp)));

%% Extract target frequency component

extractor = postprocess.elastography.extract_freq_component(elasto);

extractor.target_frequency = elasto.target_frequency;
extractor.component = 'u';

extractor.run();

disp('Frequency extraction completed')

%% Spatial filtering

c0 = shear_params.c_shear_bkg;
wj = 2*pi*elasto.target_frequency;

elasto.radial_filter = uff.spatial_filter( ...
    'gaussian_charite','radial');

elasto.radial_filter.sigma = c0/wj;

numOfAngles = 12;

elasto.directional_filtering_angles = ...
    2*pi*linspace(0,1-1/numOfAngles,numOfAngles);

elasto.angular_filter = uff.spatial_filter( ...
    'gaussian_charite','angular');

elasto.angular_filter.deltaTheta = ...
    median(diff(elasto.directional_filtering_angles));

filterOp = postprocess.elastography.apply_filter(elasto);
filterOp.run();

disp('Spatial filtering completed')

%% Shear-wave speed estimation - phase gradient

elastoPG = uff.elastography();
elastoPG.copy(elasto);

pgm = postprocess.elastography.phase_gradient_method(elastoPG);
pgm.run();

disp('Phase-gradient SWS estimation completed')

%% Visualize SWS map

dispSWS = visualization.display();
dispSWS.source = elastoPG;
dispSWS.field_name = 'sws_map';
dispSWS.component = 'u';
dispSWS.colormap_name = 'parula';
dispSWS.clim = [2 5];



title('Phase gradient SWS');

dispSWS.clim = [2 5];
dispSWS.show();


%% Inspect SWS values

sws = elastoPG.sws_map.u;

p = prctile(sws(:), [1 5 25 50 75 95 99]);

fprintf('1%%:  %.3f m/s\n', p(1));
fprintf('5%%:  %.3f m/s\n', p(2));
fprintf('25%%: %.3f m/s\n', p(3));
fprintf('50%%: %.3f m/s\n', p(4));
fprintf('75%%: %.3f m/s\n', p(5));
fprintf('95%%: %.3f m/s\n', p(6));
fprintf('99%%: %.3f m/s\n', p(7));

%% Create ground-truth SWS map

sws = squeeze(elastoPG.sws_map.u);

[X, Z] = meshgrid( ...
    elastoPG.scan.x_axis, ...
    elastoPG.scan.z_axis);

% Only inclusion 1 intersects the current ultrasound imaging plane
mask_in1 = ...
    (X - shear_params.cx1).^2 + ...
    (Z - shear_params.cz1).^2 <= shear_params.cr1^2;

gt_sws = shear_params.c_shear_bkg * ones(size(X));

gt_sws(mask_in1) = ...
    shear_params.c_shear_incl;

%% Compare ground truth and estimated SWS

figure;

subplot(1,2,1)
imagesc( ...
    elastoPG.scan.x_axis*1e3, ...
    elastoPG.scan.z_axis*1e3, ...
    gt_sws);

axis image
set(gca,'YDir','normal')
clim([2 5])
colorbar
xlabel('x [mm]')
ylabel('z [mm]')
title('Ground truth SWS')

subplot(1,2,2)
imagesc( ...
    elastoPG.scan.x_axis*1e3, ...
    elastoPG.scan.z_axis*1e3, ...
    sws);

axis image
set(gca,'YDir','normal')
clim([2 5])
colorbar
xlabel('x [mm]')
ylabel('z [mm]')
title('Estimated SWS')

%% ROI quality metrics

%% Define valid evaluation region

edge_margin = 3e-3;   % 3 mm

valid_region = ...
    X > min(X(:)) + edge_margin & ...
    X < max(X(:)) - edge_margin & ...
    Z > min(Z(:)) + edge_margin & ...
    Z < max(Z(:)) - edge_margin;

mask_bkg = valid_region & ~mask_in1;

mask_in1_eval = mask_in1 & valid_region;


mean_in1 = mean(sws(mask_in1_eval));

mean_bkg = mean(sws(mask_bkg));

bias_in1 = mean_in1 - shear_params.c_shear_incl;

bias_bkg = mean_bkg - shear_params.c_shear_bkg;

std_in1 = std(sws(mask_in1_eval));

std_bkg = std(sws(mask_bkg));

fprintf('\nSWS ROI metrics:\n');
fprintf('Inclusion 1: mean = %.3f, bias = %.3f, std = %.3f m/s\n', ...
    mean_in1, bias_in1, std_in1);


fprintf('Background:  mean = %.3f, bias = %.3f, std = %.3f m/s\n', ...
    mean_bkg, bias_bkg, std_bkg);

%% Interior ROI masks

mask_in1_core = ...
    (X - shear_params.cx1).^2 + ...
    (Z - shear_params.cz1).^2 <= (0.7*shear_params.cr1)^2;


mask_in1_core = mask_in1_core & valid_region;


mean_in1_core = mean(sws(mask_in1_core));


std_in1_core = std(sws(mask_in1_core));


fprintf('\nInterior ROI metrics:\n');
fprintf('Inclusion 1 core: mean = %.3f, bias = %.3f, std = %.3f m/s\n', ...
    mean_in1_core, ...
    mean_in1_core - shear_params.c_shear_incl, ...
    std_in1_core);
%% Save SWS result for comparison

sws_from_gt_disp = sws;
x_gt_disp = elastoPG.scan.x_axis(:).';
z_gt_disp = elastoPG.scan.z_axis(:);

save(fullfile(tempdir, 'sws_from_gt_disp.mat'), ...
    'sws_from_gt_disp', ...
    'x_gt_disp', ...
    'z_gt_disp');

disp('Saved ground-truth displacement SWS result for comparison')

%% Homogeneous SWS metrics

sws = squeeze(elastoPG.sws_map.u);

edge_margin = 3e-3;

[X, Z] = meshgrid( ...
    elastoPG.scan.x_axis, ...
    elastoPG.scan.z_axis);

valid_region = ...
    X > min(X(:)) + edge_margin & ...
    X < max(X(:)) - edge_margin & ...
    Z > min(Z(:)) + edge_margin & ...
    Z < max(Z(:)) - edge_margin;

sws_valid = sws(valid_region);

fprintf('\nHomogeneous phantom SWS:\n');
fprintf('Ground truth: %.3f m/s\n', shear_params.c_shear_bkg);
fprintf('Mean:         %.3f m/s\n', mean(sws_valid,'omitnan'));
fprintf('Median:       %.3f m/s\n', median(sws_valid,'omitnan'));
fprintf('Bias:         %.3f m/s\n', ...
    mean(sws_valid,'omitnan') - shear_params.c_shear_bkg);
fprintf('Std:          %.3f m/s\n', std(sws_valid,'omitnan'));