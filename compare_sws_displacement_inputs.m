clear; clc; close all;

%% Load SWS results

auto_file = fullfile(tempdir, 'sws_from_autocorr.mat');
gt_file   = fullfile(tempdir, 'sws_from_gt_disp.mat');

load(auto_file);
load(gt_file);

fprintf('Autocorr SWS size: %s\n', mat2str(size(sws_from_autocorr)));
fprintf('GT-disp SWS size:  %s\n', mat2str(size(sws_from_gt_disp)));

%% Interpolate autocorrelation SWS onto GT-displacement grid

[X_gt, Z_gt] = meshgrid(x_gt_disp, z_gt_disp);

sws_autocorr_on_gt = interp2( ...
    x_autocorr, ...
    z_autocorr, ...
    double(sws_from_autocorr), ...
    X_gt, ...
    Z_gt, ...
    'linear', ...
    NaN);

%% Calculate difference

difference = sws_autocorr_on_gt - double(sws_from_gt_disp);

valid = ...
    isfinite(sws_autocorr_on_gt) & ...
    isfinite(sws_from_gt_disp);

err = difference(valid);

%% Metrics

RMSE = sqrt(mean(err.^2));
MAE  = mean(abs(err));
Bias = mean(err);

fprintf('\nSWS difference caused by displacement estimation:\n');
fprintf('RMSE: %.4f m/s\n', RMSE);
fprintf('MAE:  %.4f m/s\n', MAE);
fprintf('Bias: %.4f m/s\n', Bias);
fprintf('Valid pixels: %d\n', sum(valid(:)));
%% ROI comparison

repo_dir = fileparts(mfilename('fullpath'));

load(fullfile(repo_dir, ...
    'large_inclusions_displacement.mat'), ...
    'shear_params');

edge_margin = 3e-3;

valid_region = ...
    X_gt > min(X_gt(:)) + edge_margin & ...
    X_gt < max(X_gt(:)) - edge_margin & ...
    Z_gt > min(Z_gt(:)) + edge_margin & ...
    Z_gt < max(Z_gt(:)) - edge_margin;

mask_in1 = ...
    (X_gt - shear_params.cx1).^2 + ...
    (Z_gt - shear_params.cz1).^2 <= shear_params.cr1^2;

mask_in1 = mask_in1 & valid_region;
mask_bkg = valid_region & ~mask_in1;

%% Inclusion 1

err_in1 = difference(mask_in1);

RMSE_in1 = sqrt(mean(err_in1.^2));
MAE_in1  = mean(abs(err_in1));
Bias_in1 = mean(err_in1);

%% Background

err_bkg = difference(mask_bkg);

RMSE_bkg = sqrt(mean(err_bkg.^2));
MAE_bkg  = mean(abs(err_bkg));
Bias_bkg = mean(err_bkg);

fprintf('\nSWS difference inside valid ROIs:\n');

fprintf('Inclusion 1: RMSE = %.4f, MAE = %.4f, Bias = %.4f m/s\n', ...
    RMSE_in1, MAE_in1, Bias_in1);

fprintf('Background:  RMSE = %.4f, MAE = %.4f, Bias = %.4f m/s\n', ...
    RMSE_bkg, MAE_bkg, Bias_bkg);

%% Visual comparison

figure;

subplot(1,3,1)

imagesc( ...
    x_gt_disp*1e3, ...
    z_gt_disp*1e3, ...
    sws_from_gt_disp);

axis image
set(gca,'YDir','normal')
clim([2 5])
colorbar
xlabel('x [mm]')
ylabel('z [mm]')
title('SWS from GT displacement')

subplot(1,3,2)

imagesc( ...
    x_gt_disp*1e3, ...
    z_gt_disp*1e3, ...
    sws_autocorr_on_gt);

axis image
set(gca,'YDir','normal')
clim([2 5])
colorbar
xlabel('x [mm]')
ylabel('z [mm]')
title('SWS from autocorrelation')

subplot(1,3,3)

imagesc( ...
    x_gt_disp*1e3, ...
    z_gt_disp*1e3, ...
    difference);

axis image
set(gca,'YDir','normal')
colorbar
xlabel('x [mm]')
ylabel('z [mm]')
title('Autocorr - GT displacement')