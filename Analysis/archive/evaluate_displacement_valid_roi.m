%% =============================================================
% VALID ROI DISPLACEMENT EVALUATION
% =============================================================

% Correct sign convention
disp_est_corrected = -disp_est_gt_grid;

%% -------------------------------------------------------------
% 1. Define interior ROI
% -------------------------------------------------------------
%
% Exclude a small spatial margin to avoid estimator/window/boundary effects.

Nz = size(gt_packet,1);
Nx = size(gt_packet,2);

margin_z = 3;
margin_x = 3;

mask_2d = false(Nz, Nx);

mask_2d( ...
    1+margin_z:end-margin_z, ...
    1+margin_x:end-margin_x) = true;

% Extend ROI over all time frames
mask_3d = repmat( ...
    mask_2d, ...
    1, ...
    1, ...
    size(gt_packet,3));

%% -------------------------------------------------------------
% 2. Calculate metrics in valid ROI
% -------------------------------------------------------------

metrics_roi = compute_quality_metrics( ...
    disp_est_corrected, ...
    gt_packet, ...
    mask_3d);

fprintf('\n============================================================\n');
fprintf(' DISPLACEMENT METRICS - INTERIOR ROI\n');
fprintf('============================================================\n');

fprintf('N:           %d\n', ...
    metrics_roi.N);

fprintf('RMSE:        %.4f um\n', ...
    metrics_roi.RMSE * 1e6);

fprintf('NRMSE:       %.2f %%\n', ...
    metrics_roi.NRMSE_percent);

fprintf('MAE:         %.4f um\n', ...
    metrics_roi.MAE * 1e6);

fprintf('Bias:        %.4f um\n', ...
    metrics_roi.bias * 1e6);

fprintf('Error std:   %.4f um\n', ...
    metrics_roi.error_std * 1e6);

fprintf('Correlation: %.4f\n', ...
    metrics_roi.correlation);

%% -------------------------------------------------------------
% 3. Examine error distribution
% -------------------------------------------------------------

valid_est = disp_est_corrected(mask_3d);
valid_gt  = gt_packet(mask_3d);

err = valid_est - valid_gt;

fprintf('\n--- Error distribution ---\n');

fprintf('Estimate 1st/99th percentile: %.4f / %.4f um\n', ...
    prctile(valid_est,1)*1e6, ...
    prctile(valid_est,99)*1e6);

fprintf('Error 1st/99th percentile:    %.4f / %.4f um\n', ...
    prctile(err,1)*1e6, ...
    prctile(err,99)*1e6);

fprintf('Maximum absolute error:        %.4f um\n', ...
    max(abs(err))*1e6);

%% -------------------------------------------------------------
% 4. Error histogram
% -------------------------------------------------------------

figure;

histogram(err*1e6, 100);

xlabel('Displacement error [\mum]');
ylabel('Number of samples');
title('Autocorrelation displacement error');
grid on;

%% -------------------------------------------------------------
% 5. Robust error statistics
% -------------------------------------------------------------

abs_err_um = abs(err) * 1e6;

fprintf('\n--- Robust error statistics ---\n');

fprintf('Median absolute error: %.4f um\n', ...
    median(abs_err_um));

fprintf('95th percentile abs. error: %.4f um\n', ...
    prctile(abs_err_um,95));

fprintf('99th percentile abs. error: %.4f um\n', ...
    prctile(abs_err_um,99));

fprintf('Fraction |error| > 2 um: %.3f %%\n', ...
    100 * mean(abs_err_um > 2));

fprintf('Fraction |error| > 3 um: %.3f %%\n', ...
    100 * mean(abs_err_um > 3));

fprintf('Fraction |error| > 5 um: %.3f %%\n', ...
    100 * mean(abs_err_um > 5));


%% -------------------------------------------------------------
% 6. Maximum absolute error over time
% -------------------------------------------------------------

error_volume = ...
    disp_est_corrected - gt_packet;

max_abs_error_map = ...
    max(abs(error_volume), [], 3) * 1e6;

figure;

imagesc(max_abs_error_map);
axis image;
colorbar;

xlabel('Lateral sample');
ylabel('Depth sample');

title('Maximum absolute displacement error over time [\mum]');


%% -------------------------------------------------------------
% 7. RMSE per frame
% -------------------------------------------------------------

rmse_per_frame = zeros(size(error_volume,3),1);

for k = 1:size(error_volume,3)

    frame_error = error_volume(:,:,k);

    valid_frame_error = frame_error(mask_2d);

    rmse_per_frame(k) = ...
        sqrt(mean(valid_frame_error.^2)) * 1e6;

end

figure;

plot(1:length(rmse_per_frame), ...
     rmse_per_frame, ...
     '-o');

xlabel('Frame');
ylabel('RMSE [\mum]');
title('Displacement RMSE per frame');
grid on;


fprintf('\nWorst frame: %d\n', ...
    find(rmse_per_frame == max(rmse_per_frame),1));

fprintf('Worst-frame RMSE: %.4f um\n', ...
    max(rmse_per_frame));