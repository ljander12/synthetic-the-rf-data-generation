clear;
clc;
close all;

%% ========================================================================
% REPEATED QUALITY-METRIC NOISE SWEEP
%
% Purpose:
%   Repeat the noise-sweep experiment for multiple random noise
%   realizations to assess statistical robustness of:
%
%                       Q_C = |C - 1|
%
% Dataset:
%   large_inclusions_displacement_with_quality.mat
%
% Noise levels:
%   Baseline, 40, 30, 20, 10 dB
%
% Repetitions:
%   10 random seeds per non-baseline SNR condition
% ========================================================================


%% ------------------------------------------------------------------------
% 1. Paths
% -------------------------------------------------------------------------

toolbox_dir = ...
    'C:\Users\louis\Documents\GitHub\MasterOppgave\elastography\ElastographyToolbox';

addpath(genpath(toolbox_dir));


%% ------------------------------------------------------------------------
% 2. Load data
% -------------------------------------------------------------------------

load( ...
    'large_inclusions_displacement_with_quality.mat', ...
    'b_data', ...
    'channel_data', ...
    'disp_gt');

fprintf('\n============================================================\n');
fprintf(' REPEATED QUALITY-METRIC NOISE SWEEP\n');
fprintf('============================================================\n');


%% ------------------------------------------------------------------------
% 3. Experiment settings
% -------------------------------------------------------------------------

snr_levels = [Inf 40 30 20 10];

condition_labels = { ...
    'Baseline', ...
    '40 dB', ...
    '30 dB', ...
    '20 dB', ...
    '10 dB'};

N_conditions = numel(snr_levels);

N_repeats = 10;

packet_size = 4;
packet_intervals = packet_size - 1;

z_gate = 4;
x_gate = 2;

Q_threshold = 0.05;


%% ------------------------------------------------------------------------
% 4. Ground truth
% -------------------------------------------------------------------------

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


%% ------------------------------------------------------------------------
% 5. Dimensions
% -------------------------------------------------------------------------

Nz_est = length(b_data.scan.z_axis);
Nx_est = length(b_data.scan.x_axis);

Nz_gt = size(gt_packet,1);
Nx_gt = size(gt_packet,2);


%% ------------------------------------------------------------------------
% 6. Interpolation grids
% -------------------------------------------------------------------------

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
% 7. Interior ROI
% -------------------------------------------------------------------------

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
% 8. Original data and signal power
% -------------------------------------------------------------------------

original_data = b_data.data;

signal_power = ...
    mean(abs(double(original_data(:))).^2);


%% ------------------------------------------------------------------------
% 9. Allocate result matrices
%
% Rows:
%   SNR conditions
%
% Columns:
%   repeated noise realizations
% -------------------------------------------------------------------------

RMSE_um = nan(N_conditions,N_repeats);

MAE_um = nan(N_conditions,N_repeats);

median_abs_error_um = nan(N_conditions,N_repeats);

correlation_GT = nan(N_conditions,N_repeats);

median_Q = nan(N_conditions,N_repeats);

mean_Q = nan(N_conditions,N_repeats);

pearson_Q_error = nan(N_conditions,N_repeats);

spearman_Q_error = nan(N_conditions,N_repeats);

threshold_RMSE_um = nan(N_conditions,N_repeats);

threshold_MAE_um = nan(N_conditions,N_repeats);

retained_percent = nan(N_conditions,N_repeats);

error_gt_1um_percent = nan(N_conditions,N_repeats);

error_gt_2um_percent = nan(N_conditions,N_repeats);

error_gt_3um_percent = nan(N_conditions,N_repeats);

achieved_snr_db = nan(N_conditions,N_repeats);


%% ========================================================================
% 10. Main experiment
% ========================================================================

for n = 1:N_conditions

    target_snr = snr_levels(n);

    fprintf('\n\n============================================================\n');
    fprintf(' CONDITION: %s\n', condition_labels{n});
    fprintf('============================================================\n');


    %% Baseline does not need multiple identical runs
    if isinf(target_snr)

        repeat_list = 1;

    else

        repeat_list = 1:N_repeats;

    end


    for r = repeat_list

        fprintf('\nRepeat %d\n',r);


        %% ----------------------------------------------------------------
        % 10.1 Reproducible seed
        % -----------------------------------------------------------------

        rng(1000*n + r);


        %% ----------------------------------------------------------------
        % 10.2 Copy beamformed data
        % -----------------------------------------------------------------

        b_test = uff.beamformed_data(b_data);


        %% ----------------------------------------------------------------
        % 10.3 Add noise
        % -----------------------------------------------------------------

        if isinf(target_snr)

            noisy_data = original_data;

            achieved_snr_db(n,r) = Inf;

        else

            noise_power = ...
                signal_power / (10^(target_snr/10));


            if ~isreal(original_data)

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

                noise = ...
                    sqrt(noise_power) .* ...
                    randn( ...
                        size(original_data), ...
                        'like', ...
                        original_data);

            end


            noisy_data = ...
                original_data + noise;


            actual_noise = ...
                double(noisy_data) - ...
                double(original_data);

            actual_noise_power = ...
                mean(abs(actual_noise(:)).^2);

            achieved_snr_db(n,r) = ...
                10*log10(signal_power/actual_noise_power);

        end


        b_test.data = noisy_data;


        %% ----------------------------------------------------------------
        % 10.4 Run autocorrelation estimator
        % -----------------------------------------------------------------

        estimator = ...
            postprocess.autocorrelation_displacement_estimation();

        estimator.input = b_test;

        estimator.channel_data = channel_data;

        estimator.packet_size = packet_size;

        estimator.z_gate = z_gate;

        estimator.x_gate = x_gate;


        displacement_test = ...
            estimator.go();


        %% ----------------------------------------------------------------
        % 10.5 Reshape output
        % -----------------------------------------------------------------

        Nt_est = ...
            size(displacement_test.data,4);


        disp_est = reshape( ...
            double(displacement_test.data), ...
            Nz_est, ...
            Nx_est, ...
            Nt_est);


        C_est = reshape( ...
            double(estimator.quality_indicator), ...
            Nz_est, ...
            Nx_est, ...
            Nt_est);


        %% ----------------------------------------------------------------
        % 10.6 Correct sign convention
        % -----------------------------------------------------------------

        disp_est = -disp_est;


        %% ----------------------------------------------------------------
        % 10.7 Interpolate onto GT grid
        % -----------------------------------------------------------------

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


        %% ----------------------------------------------------------------
        % 10.8 Valid samples
        % -----------------------------------------------------------------

        valid = ...
            mask_3d & ...
            isfinite(disp_est_gt) & ...
            isfinite(gt_packet) & ...
            isfinite(C_gt);


        est_values = ...
            disp_est_gt(valid);

        gt_values = ...
            gt_packet(valid);

        C_values = ...
            C_gt(valid);


        valid_C = ...
            C_values >= 0 & ...
            isfinite(C_values);

        est_values = est_values(valid_C);
        gt_values = gt_values(valid_C);
        C_values = C_values(valid_C);


        %% ----------------------------------------------------------------
        % 10.9 Error
        % -----------------------------------------------------------------

        err = ...
            est_values - gt_values;

        abs_err = ...
            abs(err);


        %% ----------------------------------------------------------------
        % 10.10 Standard metrics
        % -----------------------------------------------------------------

        RMSE_um(n,r) = ...
            sqrt(mean(err.^2))*1e6;

        MAE_um(n,r) = ...
            mean(abs_err)*1e6;

        median_abs_error_um(n,r) = ...
            median(abs_err)*1e6;


        R_GT = corrcoef( ...
            est_values, ...
            gt_values);

        correlation_GT(n,r) = ...
            R_GT(1,2);


        %% ----------------------------------------------------------------
        % 10.11 Candidate metric
        % -----------------------------------------------------------------

        Q = ...
            abs(C_values - 1);

        median_Q(n,r) = ...
            median(Q);

        mean_Q(n,r) = ...
            mean(Q);


        %% ----------------------------------------------------------------
        % 10.12 Q versus true error
        % -----------------------------------------------------------------

        R_Q = corrcoef( ...
            Q, ...
            abs_err);

        pearson_Q_error(n,r) = ...
            R_Q(1,2);


        spearman_Q_error(n,r) = corr( ...
            Q, ...
            abs_err, ...
            'Type', ...
            'Spearman');


        %% ----------------------------------------------------------------
        % 10.13 Quality threshold
        % -----------------------------------------------------------------

        keep = ...
            Q <= Q_threshold;

        retained_percent(n,r) = ...
            100*mean(keep);


        err_keep = ...
            err(keep);

        threshold_RMSE_um(n,r) = ...
            sqrt(mean(err_keep.^2))*1e6;

        threshold_MAE_um(n,r) = ...
            mean(abs(err_keep))*1e6;


        %% ----------------------------------------------------------------
        % 10.14 Large-error fractions
        % -----------------------------------------------------------------

        error_gt_1um_percent(n,r) = ...
            100*mean(abs_err*1e6 > 1);

        error_gt_2um_percent(n,r) = ...
            100*mean(abs_err*1e6 > 2);

        error_gt_3um_percent(n,r) = ...
            100*mean(abs_err*1e6 > 3);


        %% ----------------------------------------------------------------
        % 10.15 Print repeat
        % -----------------------------------------------------------------

        fprintf([ ...
            'RMSE %.3f um | Qmed %.4f | ' ...
            'Spearman %.3f | retained %.1f %%\n'], ...
            RMSE_um(n,r), ...
            median_Q(n,r), ...
            spearman_Q_error(n,r), ...
            retained_percent(n,r));


        %% ----------------------------------------------------------------
        % 10.16 Clear large temporary arrays
        % -----------------------------------------------------------------

        clear b_test estimator displacement_test ...
              disp_est C_est disp_est_gt C_gt ...
              est_values gt_values C_values ...
              Q err abs_err err_keep ...
              noisy_data noise noise_real noise_imag actual_noise

    end


    %% --------------------------------------------------------------------
    % Copy baseline result across columns for easier statistics
    % ---------------------------------------------------------------------

    if isinf(target_snr)

        RMSE_um(n,:) = RMSE_um(n,1);

        MAE_um(n,:) = MAE_um(n,1);

        median_abs_error_um(n,:) = ...
            median_abs_error_um(n,1);

        correlation_GT(n,:) = ...
            correlation_GT(n,1);

        median_Q(n,:) = ...
            median_Q(n,1);

        mean_Q(n,:) = ...
            mean_Q(n,1);

        pearson_Q_error(n,:) = ...
            pearson_Q_error(n,1);

        spearman_Q_error(n,:) = ...
            spearman_Q_error(n,1);

        threshold_RMSE_um(n,:) = ...
            threshold_RMSE_um(n,1);

        threshold_MAE_um(n,:) = ...
            threshold_MAE_um(n,1);

        retained_percent(n,:) = ...
            retained_percent(n,1);

        error_gt_1um_percent(n,:) = ...
            error_gt_1um_percent(n,1);

        error_gt_2um_percent(n,:) = ...
            error_gt_2um_percent(n,1);

        error_gt_3um_percent(n,:) = ...
            error_gt_3um_percent(n,1);

        achieved_snr_db(n,:) = Inf;

    end

end


%% ========================================================================
% 11. Mean and standard deviation across repetitions
% ========================================================================

RMSE_mean = mean(RMSE_um,2);
RMSE_std  = std(RMSE_um,0,2);

median_Q_mean = mean(median_Q,2);
median_Q_std  = std(median_Q,0,2);

spearman_mean = mean(spearman_Q_error,2);
spearman_std  = std(spearman_Q_error,0,2);

pearson_mean = mean(pearson_Q_error,2);
pearson_std  = std(pearson_Q_error,0,2);

retained_mean = mean(retained_percent,2);
retained_std  = std(retained_percent,0,2);

threshold_RMSE_mean = mean(threshold_RMSE_um,2);
threshold_RMSE_std  = std(threshold_RMSE_um,0,2);

correlation_GT_mean = mean(correlation_GT,2);
correlation_GT_std  = std(correlation_GT,0,2);


%% ========================================================================
% 12. Summary table
% ========================================================================

Condition = condition_labels(:);

TargetSNR_dB = snr_levels(:);

summary_table = table( ...
    Condition, ...
    TargetSNR_dB, ...
    RMSE_mean, ...
    RMSE_std, ...
    correlation_GT_mean, ...
    correlation_GT_std, ...
    median_Q_mean, ...
    median_Q_std, ...
    pearson_mean, ...
    pearson_std, ...
    spearman_mean, ...
    spearman_std, ...
    threshold_RMSE_mean, ...
    threshold_RMSE_std, ...
    retained_mean, ...
    retained_std);


fprintf('\n\n============================================================\n');
fprintf(' REPEATED NOISE-SWEEP SUMMARY\n');
fprintf('============================================================\n\n');

disp(summary_table);


%% ========================================================================
% 13. RMSE plot with error bars
% ========================================================================

x = 1:N_conditions;

figure;

errorbar( ...
    x, ...
    RMSE_mean, ...
    RMSE_std, ...
    '-o', ...
    'LineWidth', ...
    1.5);

hold on;

errorbar( ...
    x, ...
    threshold_RMSE_mean, ...
    threshold_RMSE_std, ...
    '-o', ...
    'LineWidth', ...
    1.5);

hold off;

xticks(x);
xticklabels(condition_labels);

xlabel('Noise condition');
ylabel('Displacement RMSE [\mum]');

title('Displacement RMSE across repeated noise realizations');

legend( ...
    'All samples', ...
    sprintf('|C-1| \\leq %.2f',Q_threshold), ...
    'Location', ...
    'best');

grid on;


%% ========================================================================
% 14. Candidate metric plot
% ========================================================================

figure;

errorbar( ...
    x, ...
    median_Q_mean, ...
    median_Q_std, ...
    '-o', ...
    'LineWidth', ...
    1.5);

xticks(x);
xticklabels(condition_labels);

xlabel('Noise condition');
ylabel('Median |C-1|');

title('Quality metric versus added noise');

grid on;


%% ========================================================================
% 15. Spearman correlation plot
% ========================================================================

figure;

errorbar( ...
    x, ...
    spearman_mean, ...
    spearman_std, ...
    '-o', ...
    'LineWidth', ...
    1.5);

xticks(x);
xticklabels(condition_labels);

xlabel('Noise condition');
ylabel('Spearman correlation');

title('|C-1| versus absolute displacement error');

grid on;


%% ========================================================================
% 16. Retained fraction plot
% ========================================================================

figure;

errorbar( ...
    x, ...
    retained_mean, ...
    retained_std, ...
    '-o', ...
    'LineWidth', ...
    1.5);

xticks(x);
xticklabels(condition_labels);

xlabel('Noise condition');
ylabel('Retained samples [%]');

title(sprintf( ...
    'Data retained with |C-1| \\leq %.2f', ...
    Q_threshold));

ylim([0 100]);

grid on;


%% ========================================================================
% 17. Save
% ========================================================================

save( ...
    'large_inclusions_quality_metric_noise_sweep_repeated.mat', ...
    'summary_table', ...
    'RMSE_um', ...
    'MAE_um', ...
    'median_abs_error_um', ...
    'correlation_GT', ...
    'median_Q', ...
    'mean_Q', ...
    'pearson_Q_error', ...
    'spearman_Q_error', ...
    'threshold_RMSE_um', ...
    'threshold_MAE_um', ...
    'retained_percent', ...
    'error_gt_1um_percent', ...
    'error_gt_2um_percent', ...
    'error_gt_3um_percent', ...
    'achieved_snr_db', ...
    'snr_levels', ...
    'condition_labels', ...
    'N_repeats', ...
    'Q_threshold');


writetable( ...
    summary_table, ...
    'large_inclusions_quality_metric_noise_sweep_repeated.csv');


fprintf('\n============================================================\n');
fprintf(' REPEATED NOISE SWEEP COMPLETE\n');
fprintf('============================================================\n');