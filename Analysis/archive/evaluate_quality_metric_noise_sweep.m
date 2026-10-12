clear;
clc;
close all;

%% ========================================================================
% QUALITY-METRIC NOISE SWEEP
%
% Purpose:
%   Evaluate robustness of the autocorrelation displacement estimator and
%   investigate whether
%
%                       Q_C = |C - 1|
%
%   remains informative when additional noise is introduced.
%
% Noise is added directly to the existing beamformed data. k-Wave,
% Field II and beamforming are NOT rerun.
%
% Tested conditions:
%   Baseline (no added noise)
%   40 dB
%   30 dB
%   20 dB
%   10 dB
%
% IMPORTANT:
%   The stated SNR is the SNR of the ORIGINAL beamformed data relative to
%   the ADDED Gaussian noise. The original simulation may already contain
%   other estimation errors/noise.
% ========================================================================


%% ------------------------------------------------------------------------
% 1. Paths
% -------------------------------------------------------------------------

toolbox_dir = ...
    'C:\Users\louis\Documents\GitHub\MasterOppgave\elastography\ElastographyToolbox';

addpath(genpath(toolbox_dir));


%% ------------------------------------------------------------------------
% 2. Load existing data
% -------------------------------------------------------------------------

load( ...
    'large_inclusions_displacement_with_quality.mat', ...
    'b_data', ...
    'channel_data', ...
    'disp_gt');

fprintf('\n============================================================\n');
fprintf(' QUALITY-METRIC NOISE SWEEP\n');
fprintf('============================================================\n');


%% ------------------------------------------------------------------------
% 3. Reproducible random noise
% -------------------------------------------------------------------------

rng(1);


%% ------------------------------------------------------------------------
% 4. Noise levels
% -------------------------------------------------------------------------
%
% Inf = no additional noise.

snr_levels = [Inf 40 30 20 10];

condition_labels = { ...
    'Baseline', ...
    '40 dB', ...
    '30 dB', ...
    '20 dB', ...
    '10 dB'};

N_conditions = numel(snr_levels);


%% ------------------------------------------------------------------------
% 5. Autocorrelation parameters
% -------------------------------------------------------------------------

packet_size = 4;
z_gate = 4;
x_gate = 2;

packet_intervals = packet_size - 1;


%% ------------------------------------------------------------------------
% 6. Prepare correct packet-averaged ground truth
% -------------------------------------------------------------------------
%
% disp_gt already contains displacement increments over individual PRF
% intervals. Do NOT differentiate disp_gt again.
%
% Four RF frames contain three pulse intervals. The autocorrelation
% estimate is therefore compared with the mean displacement across
% three consecutive GT intervals.

N_frames_gt = ...
    size(disp_gt,3) - packet_intervals + 1;

gt_packet = zeros( ...
    size(disp_gt,1), ...
    size(disp_gt,2), ...
    N_frames_gt);

for k = 1:N_frames_gt

    gt_packet(:,:,k) = mean( ...
        double(disp_gt(:,:,k:k+packet_intervals-1)), ...
        3);

end

fprintf('\nGround-truth packet size: [%d %d %d]\n', ...
    size(gt_packet));


%% ------------------------------------------------------------------------
% 7. Spatial dimensions
% -------------------------------------------------------------------------

Nz_est = length(b_data.scan.z_axis);
Nx_est = length(b_data.scan.x_axis);

Nz_gt = size(gt_packet,1);
Nx_gt = size(gt_packet,2);


%% ------------------------------------------------------------------------
% 8. Interpolation grids
% -------------------------------------------------------------------------
%
% Previous geometry checks showed that b_data and disp_gt span the same
% physical imaging region.

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


%% ------------------------------------------------------------------------
% 9. Interior ROI
% -------------------------------------------------------------------------
%
% Remove a small spatial margin to avoid boundary/windowing effects.

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
    N_frames_gt);


%% ------------------------------------------------------------------------
% 10. Store original beamformed data
% -------------------------------------------------------------------------

original_data = b_data.data;

signal_power = ...
    mean(abs(double(original_data(:))).^2);

fprintf('\nOriginal beamformed-data RMS: %.6g\n', ...
    sqrt(signal_power));


%% ------------------------------------------------------------------------
% 11. Result buffers
% -------------------------------------------------------------------------

achieved_snr_db = nan(N_conditions,1);

RMSE_um = nan(N_conditions,1);
MAE_um = nan(N_conditions,1);
median_abs_error_um = nan(N_conditions,1);
correlation_GT = nan(N_conditions,1);

median_Q = nan(N_conditions,1);
mean_Q = nan(N_conditions,1);

pearson_Q_error = nan(N_conditions,1);
spearman_Q_error = nan(N_conditions,1);

threshold_Q = 0.05;

threshold_RMSE_um = nan(N_conditions,1);
threshold_MAE_um = nan(N_conditions,1);
threshold_retained_percent = nan(N_conditions,1);

error_gt_1um_percent = nan(N_conditions,1);
error_gt_2um_percent = nan(N_conditions,1);
error_gt_3um_percent = nan(N_conditions,1);


%% ========================================================================
% 12. NOISE SWEEP
% ========================================================================

for n = 1:N_conditions

    target_snr = snr_levels(n);

    fprintf('\n\n============================================================\n');
    fprintf(' CONDITION %d/%d: %s\n', ...
        n, ...
        N_conditions, ...
        condition_labels{n});
    fprintf('============================================================\n');


    %% --------------------------------------------------------------------
    % 12.1 Create independent copy of beamformed data
    % ---------------------------------------------------------------------

    b_test = uff.beamformed_data(b_data);


    %% --------------------------------------------------------------------
    % 12.2 Add Gaussian noise
    % ---------------------------------------------------------------------

    if isinf(target_snr)

        % Baseline: no added noise
        noisy_data = original_data;

        achieved_snr_db(n) = Inf;

        fprintf('No additional noise added.\n');

    else

        noise_power = ...
            signal_power / (10^(target_snr/10));

        % Generate noise with the same numeric type as the beamformed data
        if ~isreal(original_data)

            % Complex circular Gaussian noise
            noise_real = randn( ...
                size(original_data), ...
                'like', ...
                real(original_data));

            noise_imag = randn( ...
                size(original_data), ...
                'like', ...
                real(original_data));

            noise = ...
                sqrt(noise_power/2) .* ...
                complex(noise_real,noise_imag);

        else

            % Real Gaussian noise
            noise = ...
                sqrt(noise_power) .* ...
                randn( ...
                    size(original_data), ...
                    'like', ...
                    original_data);

        end


        noisy_data = ...
            original_data + noise;


        %% ----------------------------------------------------------------
        % Measure achieved added-noise SNR
        % -----------------------------------------------------------------

        actual_noise = ...
            double(noisy_data) - double(original_data);

        actual_noise_power = ...
            mean(abs(actual_noise(:)).^2);

        achieved_snr_db(n) = ...
            10*log10(signal_power/actual_noise_power);

        fprintf('Target added-noise SNR:   %.2f dB\n', ...
            target_snr);

        fprintf('Measured added-noise SNR: %.2f dB\n', ...
            achieved_snr_db(n));


        clear noise noise_real noise_imag actual_noise

    end


    %% --------------------------------------------------------------------
    % 12.3 Assign noisy beamformed data
    % ---------------------------------------------------------------------

    b_test.data = noisy_data;

    clear noisy_data


    %% --------------------------------------------------------------------
    % 12.4 Create autocorrelation estimator
    % ---------------------------------------------------------------------

    estimator = ...
        postprocess.autocorrelation_displacement_estimation();

    estimator.input = b_test;
    estimator.channel_data = channel_data;

    estimator.packet_size = packet_size;
    estimator.z_gate = z_gate;
    estimator.x_gate = x_gate;


    %% --------------------------------------------------------------------
    % 12.5 Run estimator
    % ---------------------------------------------------------------------

    fprintf('Running autocorrelation estimator...\n');

    displacement_test = ...
        estimator.go();

    C_test_flat = ...
        estimator.quality_indicator;


    %% --------------------------------------------------------------------
    % 12.6 Reshape displacement and quality indicator
    % ---------------------------------------------------------------------

    Nt_est = ...
        size(displacement_test.data,4);

    if Nt_est ~= N_frames_gt

        error( ...
            ['Temporal mismatch: displacement estimator returned %d ' ...
             'frames, while GT contains %d frames.'], ...
            Nt_est, ...
            N_frames_gt);

    end


    disp_est = reshape( ...
        double(displacement_test.data), ...
        Nz_est, ...
        Nx_est, ...
        Nt_est);

    C_est = reshape( ...
        double(C_test_flat), ...
        Nz_est, ...
        Nx_est, ...
        Nt_est);


    %% --------------------------------------------------------------------
    % 12.7 Correct displacement sign convention
    % ---------------------------------------------------------------------
    %
    % Previous GT validation showed that the autocorrelation estimator uses
    % the opposite sign convention relative to get_GT_kwave.

    disp_est = -disp_est;


    %% --------------------------------------------------------------------
    % 12.8 Interpolate displacement and C onto GT grid
    % ---------------------------------------------------------------------

    disp_gt_grid = zeros( ...
        Nz_gt, ...
        Nx_gt, ...
        Nt_est);

    C_gt_grid = nan( ...
        Nz_gt, ...
        Nx_gt, ...
        Nt_est);

    for k = 1:Nt_est

        disp_gt_grid(:,:,k) = interp2( ...
            X_est, ...
            Z_est, ...
            disp_est(:,:,k), ...
            X_gt, ...
            Z_gt, ...
            'linear');

        C_gt_grid(:,:,k) = interp2( ...
            X_est, ...
            Z_est, ...
            C_est(:,:,k), ...
            X_gt, ...
            Z_gt, ...
            'linear');

    end


    %% --------------------------------------------------------------------
    % 12.9 Valid samples
    % ---------------------------------------------------------------------

    valid = ...
        mask_3d & ...
        isfinite(disp_gt_grid) & ...
        isfinite(gt_packet) & ...
        isfinite(C_gt_grid);


    est_values = ...
        disp_gt_grid(valid);

    gt_values = ...
        gt_packet(valid);

    C_values = ...
        C_gt_grid(valid);


    %% --------------------------------------------------------------------
    % 12.10 Remove invalid negative quality values
    % ---------------------------------------------------------------------

    valid_C = ...
        isfinite(C_values) & ...
        C_values >= 0;

    est_values = est_values(valid_C);
    gt_values = gt_values(valid_C);
    C_values = C_values(valid_C);


    %% --------------------------------------------------------------------
    % 12.11 Displacement error
    % ---------------------------------------------------------------------

    err = ...
        est_values - gt_values;

    abs_err = ...
        abs(err);


    %% --------------------------------------------------------------------
    % 12.12 Standard displacement metrics
    % ---------------------------------------------------------------------

    RMSE_um(n) = ...
        sqrt(mean(err.^2))*1e6;

    MAE_um(n) = ...
        mean(abs_err)*1e6;

    median_abs_error_um(n) = ...
        median(abs_err)*1e6;


    if std(est_values) > eps && std(gt_values) > eps

        R_GT = corrcoef( ...
            est_values, ...
            gt_values);

        correlation_GT(n) = ...
            R_GT(1,2);

    end


    %% --------------------------------------------------------------------
    % 12.13 Candidate metric Q_C = |C - 1|
    % ---------------------------------------------------------------------

    Q = ...
        abs(C_values - 1);

    median_Q(n) = ...
        median(Q);

    mean_Q(n) = ...
        mean(Q);


    %% --------------------------------------------------------------------
    % 12.14 Correlation between Q and true displacement error
    % ---------------------------------------------------------------------

    if std(Q) > eps && std(abs_err) > eps

        R_Q = corrcoef( ...
            Q, ...
            abs_err);

        pearson_Q_error(n) = ...
            R_Q(1,2);

        spearman_Q_error(n) = corr( ...
            Q, ...
            abs_err, ...
            'Type', ...
            'Spearman');

    end


    %% --------------------------------------------------------------------
    % 12.15 Threshold Q <= 0.05
    % ---------------------------------------------------------------------

    keep = ...
        Q <= threshold_Q;

    threshold_retained_percent(n) = ...
        100*mean(keep);

    if any(keep)

        err_keep = ...
            err(keep);

        threshold_RMSE_um(n) = ...
            sqrt(mean(err_keep.^2))*1e6;

        threshold_MAE_um(n) = ...
            mean(abs(err_keep))*1e6;

    end


    %% --------------------------------------------------------------------
    % 12.16 Large-error fractions
    % ---------------------------------------------------------------------

    error_gt_1um_percent(n) = ...
        100*mean(abs_err*1e6 > 1);

    error_gt_2um_percent(n) = ...
        100*mean(abs_err*1e6 > 2);

    error_gt_3um_percent(n) = ...
        100*mean(abs_err*1e6 > 3);


    %% --------------------------------------------------------------------
    % 12.17 Print condition results
    % ---------------------------------------------------------------------

    fprintf('\n--- Displacement ---\n');

    fprintf('RMSE:                 %.4f um\n', ...
        RMSE_um(n));

    fprintf('MAE:                  %.4f um\n', ...
        MAE_um(n));

    fprintf('Median abs. error:    %.4f um\n', ...
        median_abs_error_um(n));

    fprintf('Correlation with GT:  %.4f\n', ...
        correlation_GT(n));


    fprintf('\n--- Candidate metric Q_C = |C-1| ---\n');

    fprintf('Median Q_C:           %.5f\n', ...
        median_Q(n));

    fprintf('Mean Q_C:             %.5f\n', ...
        mean_Q(n));

    fprintf('Pearson Q_C vs error: %.4f\n', ...
        pearson_Q_error(n));

    fprintf('Spearman Q_C vs error: %.4f\n', ...
        spearman_Q_error(n));


    fprintf('\n--- Q_C <= %.2f ---\n', ...
        threshold_Q);

    fprintf('Retained samples:     %.2f %%\n', ...
        threshold_retained_percent(n));

    fprintf('Thresholded RMSE:     %.4f um\n', ...
        threshold_RMSE_um(n));

    fprintf('Thresholded MAE:      %.4f um\n', ...
        threshold_MAE_um(n));


    fprintf('\n--- Large displacement errors ---\n');

    fprintf('|error| > 1 um:       %.3f %%\n', ...
        error_gt_1um_percent(n));

    fprintf('|error| > 2 um:       %.3f %%\n', ...
        error_gt_2um_percent(n));

    fprintf('|error| > 3 um:       %.3f %%\n', ...
        error_gt_3um_percent(n));


  %% --------------------------------------------------------------------
% 12.18 Free large temporary arrays before next condition
% ---------------------------------------------------------------------

clear b_test estimator displacement_test C_test_flat ...
      disp_est C_est disp_gt_grid C_gt_grid ...
      est_values gt_values C_values Q err abs_err err_keep
end


%% ========================================================================
% 13. Create results table
% ========================================================================

Condition = ...
    condition_labels(:);

TargetSNR_dB = ...
    snr_levels(:);

noise_sweep_results = table( ...
    Condition, ...
    TargetSNR_dB, ...
    achieved_snr_db, ...
    RMSE_um, ...
    MAE_um, ...
    median_abs_error_um, ...
    correlation_GT, ...
    median_Q, ...
    mean_Q, ...
    pearson_Q_error, ...
    spearman_Q_error, ...
    threshold_RMSE_um, ...
    threshold_MAE_um, ...
    threshold_retained_percent, ...
    error_gt_1um_percent, ...
    error_gt_2um_percent, ...
    error_gt_3um_percent);


%% ========================================================================
% 14. Display complete table
% ========================================================================

fprintf('\n\n============================================================\n');
fprintf(' NOISE-SWEEP SUMMARY\n');
fprintf('============================================================\n\n');

disp(noise_sweep_results);


%% ========================================================================
% 15. Plot displacement RMSE
% ========================================================================

x = 1:N_conditions;

figure;

plot( ...
    x, ...
    RMSE_um, ...
    '-o', ...
    'LineWidth', ...
    1.5);

hold on;

plot( ...
    x, ...
    threshold_RMSE_um, ...
    '-o', ...
    'LineWidth', ...
    1.5);

hold off;

xticks(x);
xticklabels(condition_labels);

xlabel('Noise condition');
ylabel('Displacement RMSE [\mum]');

title('Displacement accuracy versus added noise');

legend( ...
    'All valid samples', ...
    sprintf('|C-1| \\leq %.2f',threshold_Q), ...
    'Location', ...
    'best');

grid on;


%% ========================================================================
% 16. Plot candidate quality metric
% ========================================================================

figure;

plot( ...
    x, ...
    median_Q, ...
    '-o', ...
    'LineWidth', ...
    1.5);

xticks(x);
xticklabels(condition_labels);

xlabel('Noise condition');
ylabel('Median |C-1|');

title('Autocorrelation quality deviation versus added noise');

grid on;


%% ========================================================================
% 17. Plot correlation between Q and actual error
% ========================================================================

figure;

plot( ...
    x, ...
    pearson_Q_error, ...
    '-o', ...
    'LineWidth', ...
    1.5);

hold on;

plot( ...
    x, ...
    spearman_Q_error, ...
    '-o', ...
    'LineWidth', ...
    1.5);

hold off;

xticks(x);
xticklabels(condition_labels);

xlabel('Noise condition');
ylabel('Correlation coefficient');

title('|C-1| versus absolute displacement error');

legend( ...
    'Pearson', ...
    'Spearman', ...
    'Location', ...
    'best');

grid on;


%% ========================================================================
% 18. Plot retained samples
% ========================================================================

figure;

plot( ...
    x, ...
    threshold_retained_percent, ...
    '-o', ...
    'LineWidth', ...
    1.5);

xticks(x);
xticklabels(condition_labels);

xlabel('Noise condition');
ylabel('Retained samples [%]');

title(sprintf( ...
    'Samples retained with |C-1| \\leq %.2f', ...
    threshold_Q));

ylim([0 100]);

grid on;


%% ========================================================================
% 19. Save results
% ========================================================================

save( ...
    'large_inclusions_quality_metric_noise_sweep.mat', ...
    'noise_sweep_results', ...
    'snr_levels', ...
    'condition_labels', ...
    'threshold_Q');


writetable( ...
    noise_sweep_results, ...
    'large_inclusions_quality_metric_noise_sweep.csv');


fprintf('\n============================================================\n');
fprintf(' NOISE SWEEP COMPLETE\n');
fprintf('============================================================\n');

fprintf('\nResults saved to:\n');
fprintf('large_inclusions_quality_metric_noise_sweep.mat\n');
fprintf('large_inclusions_quality_metric_noise_sweep.csv\n');