clear;
clc;
close all;

%% ========================================================================
% AUTOCORRELATION QUALITY INDICATOR ANALYSIS
%
% Purpose:
%   Test whether the autocorrelation quality indicator C is associated
%   with displacement estimation error.
%
% Dataset:
%   large_inclusions_displacement_with_quality.mat
%
% The displacement estimate is compared with packet-averaged k-Wave
% ground-truth displacement.
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

load('large_inclusions_displacement_with_quality.mat');

fprintf('\n============================================================\n');
fprintf(' AUTOCORRELATION QUALITY INDICATOR ANALYSIS\n');
fprintf('============================================================\n');


%% -------------------------------------------------------------
% 3. Dimensions
% --------------------------------------------------------------

Nz_est = length(b_data.scan.z_axis);
Nx_est = length(b_data.scan.x_axis);

Nt_est = size(displacement_estimation_autocorr.data,4);

fprintf('\nEstimated grid: %d x %d x %d\n', ...
    Nz_est, Nx_est, Nt_est);


%% -------------------------------------------------------------
% 4. Estimated displacement
% --------------------------------------------------------------

disp_est = reshape( ...
    double(displacement_estimation_autocorr.data), ...
    Nz_est, ...
    Nx_est, ...
    Nt_est);

% Previous validation showed opposite sign convention relative to GT
disp_est = -disp_est;


%% -------------------------------------------------------------
% 5. Quality indicator
% --------------------------------------------------------------

C_est = reshape( ...
    double(quality_indicator_autocorr), ...
    Nz_est, ...
    Nx_est, ...
    Nt_est);

fprintf('Quality-indicator size: [%d %d %d]\n', ...
    size(C_est));


%% -------------------------------------------------------------
% 6. Construct packet-averaged ground truth
% -------------------------------------------------------------
%
% disp_gt already contains PRF-sampled displacement increments.
%
% packet_size = 4 RF frames corresponds to 3 pulse intervals.

packet_size = 4;
packet_intervals = packet_size - 1;

N_frames = size(disp_gt,3) - packet_intervals + 1;

gt_packet = zeros( ...
    size(disp_gt,1), ...
    size(disp_gt,2), ...
    N_frames);

for k = 1:N_frames

    gt_packet(:,:,k) = mean( ...
        double(disp_gt(:,:,k:k+packet_intervals-1)), ...
        3);

end

fprintf('Ground-truth size: [%d %d %d]\n', ...
    size(gt_packet));


%% -------------------------------------------------------------
% 7. Check temporal dimensions
% --------------------------------------------------------------

if size(gt_packet,3) ~= Nt_est

    error( ...
        'Temporal mismatch: estimate has %d frames and GT has %d frames.', ...
        Nt_est, ...
        size(gt_packet,3));

end


%% -------------------------------------------------------------
% 8. Interpolate displacement estimate and C onto GT grid
% --------------------------------------------------------------

Nz_gt = size(gt_packet,1);
Nx_gt = size(gt_packet,2);

x_est_norm = linspace(0,1,Nx_est);
z_est_norm = linspace(0,1,Nz_est);

x_gt_norm = linspace(0,1,Nx_gt);
z_gt_norm = linspace(0,1,Nz_gt);

[X_est,Z_est] = meshgrid( ...
    x_est_norm, ...
    z_est_norm);

[X_gt,Z_gt] = meshgrid( ...
    x_gt_norm, ...
    z_gt_norm);

disp_est_gt = zeros( ...
    Nz_gt, ...
    Nx_gt, ...
    Nt_est);

C_gt = nan( ...
    Nz_gt, ...
    Nx_gt, ...
    Nt_est);

for k = 1:Nt_est

    disp_est_gt(:,:,k) = interp2( ...
        X_est, ...
        Z_est, ...
        disp_est(:,:,k), ...
        X_gt, ...
        Z_gt, ...
        'linear');

    C_gt(:,:,k) = interp2( ...
        X_est, ...
        Z_est, ...
        C_est(:,:,k), ...
        X_gt, ...
        Z_gt, ...
        'linear');

end


%% -------------------------------------------------------------
% 9. Interior ROI
% -------------------------------------------------------------
%
% Exclude outer image margins.

margin_z = 3;
margin_x = 3;

mask_2d = false(Nz_gt,Nx_gt);

mask_2d( ...
    1+margin_z:end-margin_z, ...
    1+margin_x:end-margin_x) = true;

mask_3d = repmat( ...
    mask_2d, ...
    1, ...
    1, ...
    Nt_est);


%% -------------------------------------------------------------
% 10. Error volume
% --------------------------------------------------------------

error_volume = ...
    disp_est_gt - gt_packet;

abs_error_volume = ...
    abs(error_volume);


%% -------------------------------------------------------------
% 11. Select valid samples
% --------------------------------------------------------------

valid = ...
    mask_3d & ...
    isfinite(C_gt) & ...
    isfinite(abs_error_volume);

C_values = ...
    C_gt(valid);

error_values = ...
    abs_error_volume(valid);


%% -------------------------------------------------------------
% 12. Remove impossible / numerical C values
% -------------------------------------------------------------
%
% C should nominally represent normalized temporal coherence.
% Keep only finite non-negative values here.

valid_C = ...
    isfinite(C_values) & ...
    C_values >= 0;

C_values = ...
    C_values(valid_C);

error_values = ...
    error_values(valid_C);


%% -------------------------------------------------------------
% 13. Basic quality-indicator statistics
% --------------------------------------------------------------

fprintf('\n------------------------------------------------------------\n');
fprintf(' QUALITY INDICATOR STATISTICS\n');
fprintf('------------------------------------------------------------\n');

fprintf('Number of valid samples: %d\n', ...
    numel(C_values));

fprintf('C mean:              %.4f\n', ...
    mean(C_values));

fprintf('C median:            %.4f\n', ...
    median(C_values));

fprintf('C 1st percentile:    %.4f\n', ...
    prctile(C_values,1));

fprintf('C 5th percentile:    %.4f\n', ...
    prctile(C_values,5));

fprintf('C 95th percentile:   %.4f\n', ...
    prctile(C_values,95));

fprintf('C 99th percentile:   %.4f\n', ...
    prctile(C_values,99));


%% -------------------------------------------------------------
% 14. Error statistics
% --------------------------------------------------------------

fprintf('\n------------------------------------------------------------\n');
fprintf(' DISPLACEMENT ERROR STATISTICS\n');
fprintf('------------------------------------------------------------\n');

fprintf('Median abs. error:       %.4f um\n', ...
    median(error_values)*1e6);

fprintf('95th percentile abs.err: %.4f um\n', ...
    prctile(error_values,95)*1e6);

fprintf('99th percentile abs.err: %.4f um\n', ...
    prctile(error_values,99)*1e6);

fprintf('Maximum abs. error:      %.4f um\n', ...
    max(error_values)*1e6);


%% -------------------------------------------------------------
% 15. Correlation between C and absolute error
% --------------------------------------------------------------

R = corrcoef( ...
    C_values, ...
    error_values);

quality_error_corr = ...
    R(1,2);

fprintf('\n------------------------------------------------------------\n');
fprintf(' QUALITY VS ERROR\n');
fprintf('------------------------------------------------------------\n');

fprintf('Pearson correlation C vs |error|: %.4f\n', ...
    quality_error_corr);


%% -------------------------------------------------------------
% 16. Spearman correlation
% --------------------------------------------------------------

spearman_corr = corr( ...
    C_values, ...
    error_values, ...
    'Type', ...
    'Spearman');

fprintf('Spearman correlation C vs |error|: %.4f\n', ...
    spearman_corr);


%% -------------------------------------------------------------
% 17. Threshold analysis
% -------------------------------------------------------------

thresholds = ...
    0.50:0.05:0.95;

threshold_results = zeros( ...
    length(thresholds), ...
    5);

fprintf('\n------------------------------------------------------------\n');
fprintf(' QUALITY-THRESHOLD ANALYSIS\n');
fprintf('------------------------------------------------------------\n');

fprintf('\n');

for k = 1:length(thresholds)

    threshold = thresholds(k);

    keep = ...
        C_values >= threshold;

    retained_fraction = ...
        mean(keep);

    if any(keep)

        err_keep = ...
            error_values(keep);

        rmse_keep = ...
            sqrt(mean(err_keep.^2));

        mae_keep = ...
            mean(err_keep);

        median_keep = ...
            median(err_keep);

        p95_keep = ...
            prctile(err_keep,95);

    else

        rmse_keep = NaN;
        mae_keep = NaN;
        median_keep = NaN;
        p95_keep = NaN;

    end

    threshold_results(k,:) = [ ...
        retained_fraction*100, ...
        rmse_keep*1e6, ...
        mae_keep*1e6, ...
        median_keep*1e6, ...
        p95_keep*1e6];

    fprintf([ ...
        'C >= %.2f | retained = %6.2f %% | ' ...
        'RMSE = %.4f um | MAE = %.4f um | ' ...
        'median = %.4f um | P95 = %.4f um\n'], ...
        threshold, ...
        retained_fraction*100, ...
        rmse_keep*1e6, ...
        mae_keep*1e6, ...
        median_keep*1e6, ...
        p95_keep*1e6);

end


%% -------------------------------------------------------------
% 18. Error rates for large-error samples
% -------------------------------------------------------------

fprintf('\n------------------------------------------------------------\n');
fprintf(' LARGE-ERROR SAMPLE ANALYSIS\n');
fprintf('------------------------------------------------------------\n');

error_thresholds_um = [1 2 3 5];

for e_thr = error_thresholds_um

    bad = ...
        error_values*1e6 > e_thr;

    fprintf('|error| > %d um: %.3f %% of samples\n', ...
        e_thr, ...
        100*mean(bad));

    if any(bad)

        fprintf('  Median C for these samples: %.4f\n', ...
            median(C_values(bad)));

    end

end


%% -------------------------------------------------------------
% 19. Compare low and high C directly
% --------------------------------------------------------------

C_split = 0.8;

low_quality = ...
    C_values < C_split;

high_quality = ...
    C_values >= C_split;

fprintf('\n------------------------------------------------------------\n');
fprintf(' C < %.2f VERSUS C >= %.2f\n', ...
    C_split, C_split);
fprintf('------------------------------------------------------------\n');

if any(low_quality)

    fprintf('Low-C samples: %.2f %%\n', ...
        100*mean(low_quality));

    fprintf('Low-C median abs error: %.4f um\n', ...
        median(error_values(low_quality))*1e6);

    fprintf('Low-C RMSE:             %.4f um\n', ...
        sqrt(mean(error_values(low_quality).^2))*1e6);

end

if any(high_quality)

    fprintf('\nHigh-C samples: %.2f %%\n', ...
        100*mean(high_quality));

    fprintf('High-C median abs error: %.4f um\n', ...
        median(error_values(high_quality))*1e6);

    fprintf('High-C RMSE:             %.4f um\n', ...
        sqrt(mean(error_values(high_quality).^2))*1e6);

end


%% -------------------------------------------------------------
% 20. Scatter plot
% -------------------------------------------------------------

figure;

scatter( ...
    C_values, ...
    error_values*1e6, ...
    5, ...
    'filled');

xlabel('Autocorrelation quality indicator C');
ylabel('Absolute displacement error [\mum]');

title( ...
    sprintf( ...
        'Quality indicator vs displacement error, r = %.3f', ...
        quality_error_corr));

grid on;


%% -------------------------------------------------------------
% 21. Median error binned by C
% -------------------------------------------------------------

C_min = max(0, floor(min(C_values)*20)/20);
C_max = ceil(max(C_values)*20)/20;

edges = ...
    C_min:0.05:C_max;

if length(edges) < 2
    edges = [C_min C_min+0.05];
end

bin_centers = ...
    edges(1:end-1) + diff(edges)/2;

median_error_bin = ...
    nan(size(bin_centers));

p95_error_bin = ...
    nan(size(bin_centers));

N_bin = ...
    zeros(size(bin_centers));

for k = 1:length(bin_centers)

    idx = ...
        C_values >= edges(k) & ...
        C_values < edges(k+1);

    N_bin(k) = ...
        sum(idx);

    if any(idx)

        median_error_bin(k) = ...
            median(error_values(idx))*1e6;

        p95_error_bin(k) = ...
            prctile(error_values(idx),95)*1e6;

    end

end


figure;

plot( ...
    bin_centers, ...
    median_error_bin, ...
    '-o', ...
    'LineWidth', ...
    1.5);

xlabel('Autocorrelation quality indicator C');
ylabel('Median absolute displacement error [\mum]');

title('Median displacement error versus autocorrelation quality');

grid on;


figure;

plot( ...
    bin_centers, ...
    p95_error_bin, ...
    '-o', ...
    'LineWidth', ...
    1.5);

xlabel('Autocorrelation quality indicator C');
ylabel('95th percentile displacement error [\mum]');

title('95th percentile displacement error versus autocorrelation quality');

grid on;


%% -------------------------------------------------------------
% 22. Threshold-performance plots
% -------------------------------------------------------------

figure;

plot( ...
    thresholds, ...
    threshold_results(:,2), ...
    '-o', ...
    'LineWidth', ...
    1.5);

xlabel('C threshold');
ylabel('RMSE [\mum]');

title('RMSE after autocorrelation quality thresholding');

grid on;


figure;

plot( ...
    thresholds, ...
    threshold_results(:,1), ...
    '-o', ...
    'LineWidth', ...
    1.5);

xlabel('C threshold');
ylabel('Retained samples [%]');

title('Data retained after autocorrelation quality thresholding');

grid on;


%% -------------------------------------------------------------
% 23. Save results
% --------------------------------------------------------------

quality_results = struct();

quality_results.C_values = ...
    C_values;

quality_results.error_values = ...
    error_values;

quality_results.pearson_correlation = ...
    quality_error_corr;

quality_results.spearman_correlation = ...
    spearman_corr;

quality_results.thresholds = ...
    thresholds;

quality_results.threshold_results = ...
    threshold_results;

quality_results.C_split = ...
    C_split;

quality_results.bin_centers = ...
    bin_centers;

quality_results.median_error_bin = ...
    median_error_bin;

quality_results.p95_error_bin = ...
    p95_error_bin;

quality_results.N_bin = ...
    N_bin;


save( ...
    'large_inclusions_autocorr_quality_analysis.mat', ...
    'quality_results');


fprintf('\n============================================================\n');
fprintf(' QUALITY-INDICATOR ANALYSIS COMPLETE\n');
fprintf('============================================================\n');

fprintf('\nResults saved to:\n');
fprintf('large_inclusions_autocorr_quality_analysis.mat\n');