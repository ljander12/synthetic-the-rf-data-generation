clear;
clc;
close all;

%% ========================================================================
% LARGE INCLUSIONS: DISPLACEMENT QUALITY METRICS
%
% Compares autocorrelation displacement estimates with k-Wave ground truth.
%
% Important:
% disp_gt already contains PRF-sampled displacement increments.
% Therefore, disp_gt must NOT be differentiated again.
% ========================================================================


%% -------------------------------------------------------------
% 1. Add ElastographyToolbox / USTB
% --------------------------------------------------------------

toolbox_dir = ...
    'C:\Users\louis\Documents\GitHub\MasterOppgave\elastography\ElastographyToolbox';

addpath(genpath(toolbox_dir));


%% -------------------------------------------------------------
% 2. Load data
% --------------------------------------------------------------

load('large_inclusions_displacement.mat');

fprintf('\n============================================================\n');
fprintf(' LARGE INCLUSIONS: DISPLACEMENT METRICS\n');
fprintf('============================================================\n');


%% -------------------------------------------------------------
% 3. Extract estimated displacement
% --------------------------------------------------------------

Nz_est = length(b_data.scan.z_axis);   % Expected: 1024
Nx_est = length(b_data.scan.x_axis);   % Expected: 512

Nt_est = size(displacement_estimation_autocorr.data, 4); % Expected: 59

disp_est = reshape( ...
    double(displacement_estimation_autocorr.data), ...
    Nz_est, ...
    Nx_est, ...
    Nt_est);

fprintf('\nEstimated displacement size: [%d %d %d]\n', ...
    size(disp_est));


%% -------------------------------------------------------------
% 4. Prepare ground-truth displacement
% -------------------------------------------------------------
%
% disp_gt already contains displacement increments sampled at the PRF.
%
% Each frame corresponds to displacement over one pulse interval.
%
% Autocorrelation uses packet_size = 4 RF frames, which corresponds
% to 3 adjacent frame-to-frame displacement intervals.
%
% Therefore, the correct GT comparison is the MEAN of three consecutive
% disp_gt frames.

packet_size = 4;

packet_intervals = packet_size - 1;

N_frames = size(disp_gt, 3) - packet_intervals + 1;

gt_packet = zeros( ...
    size(disp_gt,1), ...
    size(disp_gt,2), ...
    N_frames);

for k = 1:N_frames

    gt_packet(:,:,k) = mean( ...
        double(disp_gt(:,:,k:k+packet_intervals-1)), ...
        3);

end

fprintf('Ground-truth input size: [%d %d %d]\n', ...
    size(disp_gt));

fprintf('Aligned ground-truth size: [%d %d %d]\n', ...
    size(gt_packet));


%% -------------------------------------------------------------
% 5. Interpolate estimated displacement onto GT spatial grid
% -------------------------------------------------------------
%
% Estimated displacement:
%   1024 x 512 x 59
%
% Ground truth:
%   46 x 39 x 59
%
% The two grids cover approximately the same physical image region,
% so normalized spatial coordinates are used for interpolation.

Nz_gt = size(gt_packet,1);
Nx_gt = size(gt_packet,2);

x_est_norm = linspace(0, 1, Nx_est);
z_est_norm = linspace(0, 1, Nz_est);

x_gt_norm = linspace(0, 1, Nx_gt);
z_gt_norm = linspace(0, 1, Nz_gt);

[X_est, Z_est] = meshgrid( ...
    x_est_norm, ...
    z_est_norm);

[X_gt, Z_gt] = meshgrid( ...
    x_gt_norm, ...
    z_gt_norm);

disp_est_gt_grid = zeros( ...
    Nz_gt, ...
    Nx_gt, ...
    Nt_est);

for k = 1:Nt_est

    disp_est_gt_grid(:,:,k) = interp2( ...
        X_est, ...
        Z_est, ...
        disp_est(:,:,k), ...
        X_gt, ...
        Z_gt, ...
        'linear');

end

fprintf('Interpolated estimate size: [%d %d %d]\n', ...
    size(disp_est_gt_grid));


%% -------------------------------------------------------------
% 6. Check temporal dimensions
% --------------------------------------------------------------

if size(gt_packet,3) ~= size(disp_est_gt_grid,3)

    error( ...
        'Temporal mismatch: estimate has %d frames, GT has %d frames.', ...
        size(disp_est_gt_grid,3), ...
        size(gt_packet,3));

end


%% -------------------------------------------------------------
% 7. Calculate quality metrics
% --------------------------------------------------------------

metrics_disp = compute_quality_metrics( ...
    disp_est_gt_grid, ...
    gt_packet);


%% -------------------------------------------------------------
% 8. Display displacement metrics
% --------------------------------------------------------------

fprintf('\n------------------------------------------------------------\n');
fprintf(' DISPLACEMENT QUALITY METRICS\n');
fprintf('------------------------------------------------------------\n');

fprintf('N:           %d\n', ...
    metrics_disp.N);

fprintf('RMSE:        %.4f um\n', ...
    metrics_disp.RMSE * 1e6);

fprintf('NRMSE:       %.2f %%\n', ...
    metrics_disp.NRMSE_percent);

fprintf('MAE:         %.4f um\n', ...
    metrics_disp.MAE * 1e6);

fprintf('Bias:        %.4f um\n', ...
    metrics_disp.bias * 1e6);

fprintf('Error std:   %.4f um\n', ...
    metrics_disp.error_std * 1e6);

fprintf('Estimate std: %.4f um\n', ...
    metrics_disp.std * 1e6);

fprintf('Correlation: %.4f\n', ...
    metrics_disp.correlation);


%% -------------------------------------------------------------
% 9. Diagnostic amplitude statistics
% --------------------------------------------------------------

fprintf('\n------------------------------------------------------------\n');
fprintf(' AMPLITUDE DIAGNOSTICS\n');
fprintf('------------------------------------------------------------\n');

gt_rms = sqrt(mean(gt_packet(:).^2));

est_rms = sqrt(mean(disp_est_gt_grid(:).^2));

fprintf('GT RMS:        %.4f um\n', ...
    gt_rms * 1e6);

fprintf('Estimate RMS:  %.4f um\n', ...
    est_rms * 1e6);

fprintf('GT min/max:    %.4f / %.4f um\n', ...
    min(gt_packet(:))*1e6, ...
    max(gt_packet(:))*1e6);

fprintf('Est min/max:   %.4f / %.4f um\n', ...
    min(disp_est_gt_grid(:))*1e6, ...
    max(disp_est_gt_grid(:))*1e6);


%% -------------------------------------------------------------
% 10. Sign diagnostic
% -------------------------------------------------------------
%
% Check whether the estimator uses the opposite displacement sign
% convention relative to the ground truth.

metrics_flipped = compute_quality_metrics( ...
    -disp_est_gt_grid, ...
    gt_packet);

fprintf('\n------------------------------------------------------------\n');
fprintf(' SIGN-CONVENTION TEST\n');
fprintf('------------------------------------------------------------\n');

fprintf('Normal sign correlation:  %.4f\n', ...
    metrics_disp.correlation);

fprintf('Flipped sign correlation: %.4f\n', ...
    metrics_flipped.correlation);

fprintf('Normal sign RMSE:          %.4f um\n', ...
    metrics_disp.RMSE * 1e6);

fprintf('Flipped sign RMSE:         %.4f um\n', ...
    metrics_flipped.RMSE * 1e6);


%% -------------------------------------------------------------
% 11. Plot example frame
% -------------------------------------------------------------
%
% Visual comparison of one representative frame.

frame_idx = round(Nt_est / 2);

figure;

subplot(1,3,1);

imagesc( ...
    squeeze(gt_packet(:,:,frame_idx)) * 1e6);

axis image;
colorbar;

title(sprintf( ...
    'Ground truth, frame %d', ...
    frame_idx));

xlabel('Lateral sample');
ylabel('Depth sample');


subplot(1,3,2);

imagesc( ...
    squeeze(disp_est_gt_grid(:,:,frame_idx)) * 1e6);

axis image;
colorbar;

title(sprintf( ...
    'Autocorrelation estimate, frame %d', ...
    frame_idx));

xlabel('Lateral sample');
ylabel('Depth sample');


subplot(1,3,3);

imagesc( ...
    squeeze( ...
        disp_est_gt_grid(:,:,frame_idx) - ...
        gt_packet(:,:,frame_idx)) * 1e6);

axis image;
colorbar;

title('Error');

xlabel('Lateral sample');
ylabel('Depth sample');


%% -------------------------------------------------------------
% 12. Save result structure
% --------------------------------------------------------------

results = struct();

results.phantom = 'large_inclusions';

results.displacement.autocorrelation = metrics_disp;

results.displacement.gt_RMS = gt_rms;

results.displacement.estimate_RMS = est_rms;

results.displacement.packet_size = packet_size;

results.displacement.packet_intervals = packet_intervals;


save( ...
    'large_inclusions_quality_metrics.mat', ...
    'results');


fprintf('\n============================================================\n');
fprintf(' Metrics saved to:\n');
fprintf(' large_inclusions_quality_metrics.mat\n');
fprintf('============================================================\n');