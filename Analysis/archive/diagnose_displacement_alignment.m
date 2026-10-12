clear;
clc;
close all;

%% ========================================================================
% DIAGNOSE DISPLACEMENT ALIGNMENT
%
% Purpose:
%   Investigate how the autocorrelation displacement estimate should be
%   compared with the k-Wave ground-truth displacement.
%
% Tests:
%   1. Mean displacement across autocorrelation packet
%   2. Total/summed displacement across autocorrelation packet
%   3. Positive and negative sign convention
%   4. Temporal shifts between estimated and ground-truth displacement
%
% Dataset:
%   large_inclusions_displacement.mat
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
fprintf(' LARGE INCLUSIONS: DISPLACEMENT ALIGNMENT DIAGNOSTICS\n');
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

fprintf('Estimated scan x-range: %.3f to %.3f mm\n', ...
    min(b_data.scan.x_axis)*1e3, ...
    max(b_data.scan.x_axis)*1e3);

fprintf('Estimated scan z-range: %.3f to %.3f mm\n', ...
    min(b_data.scan.z_axis)*1e3, ...
    max(b_data.scan.z_axis)*1e3);


%% -------------------------------------------------------------
% 4. Prepare ground-truth displacement
% -------------------------------------------------------------
%
% disp_gt contains 61 sampled absolute displacement frames.
%
% Add the initial zero-displacement frame so that frame-to-frame
% displacement can be calculated.

gt_abs = cat( ...
    3, ...
    zeros(size(disp_gt(:,:,1))), ...
    double(disp_gt));

% Frame-to-frame ground-truth displacement
gt_delta = diff(gt_abs, 1, 3);

fprintf('\nOriginal GT displacement size: [%d %d %d]\n', ...
    size(disp_gt));

fprintf('GT absolute displacement size: [%d %d %d]\n', ...
    size(gt_abs));

fprintf('GT frame-difference size:      [%d %d %d]\n', ...
    size(gt_delta));


%% -------------------------------------------------------------
% 5. Construct alternative packet interpretations
% -------------------------------------------------------------
%
% The autocorrelation displacement estimator used packet_size = 4.
%
% Four RF frames contain three frame-to-frame intervals:
%
%       frame 1 -> frame 2
%       frame 2 -> frame 3
%       frame 3 -> frame 4
%
% We test two possible interpretations:
%
%   MEAN:
%       Average displacement across these three intervals.
%
%   SUM:
%       Total displacement across these three intervals.
%
% The purpose is to determine which quantity corresponds to the output
% produced by the autocorrelation estimator.

packet_size = 4;

packet_intervals = packet_size - 1;

N_frames = size(gt_delta,3) - packet_size + 2;

gt_mean = zeros( ...
    size(gt_delta,1), ...
    size(gt_delta,2), ...
    N_frames);

gt_sum = zeros( ...
    size(gt_delta,1), ...
    size(gt_delta,2), ...
    N_frames);

for k = 1:N_frames

    packet = gt_delta(:,:,k:k+packet_intervals-1);

    % Mean displacement per frame interval
    gt_mean(:,:,k) = mean(packet, 3);

    % Total displacement across the packet
    gt_sum(:,:,k) = sum(packet, 3);

end

fprintf('\nPacket size:      %d frames\n', packet_size);
fprintf('Packet intervals: %d\n', packet_intervals);

fprintf('GT mean size: [%d %d %d]\n', size(gt_mean));
fprintf('GT sum size:  [%d %d %d]\n', size(gt_sum));


%% -------------------------------------------------------------
% 6. Interpolate estimated displacement onto GT spatial grid
% -------------------------------------------------------------
%
% Estimated displacement:
%       1024 x 512 x 59
%
% Ground truth:
%       46 x 39 x 59
%
% Previous checks showed that the estimated and GT grids cover the same
% imaging region:
%
%       x approximately -19 to +19 mm
%       z approximately   5 to 50 mm
%
% Therefore normalized coordinates can be used to map between the grids.

Nz_gt = size(gt_mean,1);
Nx_gt = size(gt_mean,2);

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

fprintf('\nInterpolated estimate size: [%d %d %d]\n', ...
    size(disp_est_gt_grid));


%% -------------------------------------------------------------
% 7. Basic RMS and amplitude diagnostics
% --------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' BASIC AMPLITUDE DIAGNOSTICS\n');
fprintf('============================================================\n');

est_rms = sqrt(mean(disp_est_gt_grid(:).^2));

gt_mean_rms = sqrt(mean(gt_mean(:).^2));

gt_sum_rms = sqrt(mean(gt_sum(:).^2));

fprintf('\nEstimate RMS: %.4f um\n', ...
    est_rms * 1e6);

fprintf('GT MEAN RMS: %.4f um\n', ...
    gt_mean_rms * 1e6);

fprintf('GT SUM RMS:  %.4f um\n', ...
    gt_sum_rms * 1e6);

fprintf('\nEstimate min/max: %.4f / %.4f um\n', ...
    min(disp_est_gt_grid(:))*1e6, ...
    max(disp_est_gt_grid(:))*1e6);

fprintf('GT MEAN min/max: %.4f / %.4f um\n', ...
    min(gt_mean(:))*1e6, ...
    max(gt_mean(:))*1e6);

fprintf('GT SUM min/max:  %.4f / %.4f um\n', ...
    min(gt_sum(:))*1e6, ...
    max(gt_sum(:))*1e6);


%% -------------------------------------------------------------
% 8. Temporal alignment, packet interpretation and sign test
% -------------------------------------------------------------
%
% Test:
%
%   MEAN or SUM
%   sign = +1 or -1
%   lag  = -5 ... +5 frames
%
% The preferred solution should ideally have:
%
%   - low RMSE
%   - low NRMSE
%   - high positive correlation

gt_candidates = { ...
    gt_mean, ...
    gt_sum};

candidate_names = { ...
    'MEAN', ...
    'SUM'};

sign_candidates = [1 -1];

lag_values = -5:5;

fprintf('\n============================================================\n');
fprintf(' TEMPORAL ALIGNMENT TEST\n');
fprintf('============================================================\n\n');


%% Storage for results

result_interpretation = {};
result_sign = [];
result_lag = [];

result_RMSE_um = [];
result_NRMSE_percent = [];
result_MAE_um = [];
result_bias_um = [];
result_correlation = [];
result_N = [];


result_index = 0;


%% -------------------------------------------------------------
% Test all combinations
% --------------------------------------------------------------

for c = 1:length(gt_candidates)

    gt_test = gt_candidates{c};

    for sign_factor = sign_candidates

        for lag = lag_values

            %% -------------------------------------------------
            % Temporal shift
            % --------------------------------------------------

            if lag > 0

                % Estimate shifted forward relative to GT

                est_test = ...
                    disp_est_gt_grid(:,:,1+lag:end);

                gt_lag = ...
                    gt_test(:,:,1:end-lag);


            elseif lag < 0

                shift = abs(lag);

                est_test = ...
                    disp_est_gt_grid(:,:,1:end-shift);

                gt_lag = ...
                    gt_test(:,:,1+shift:end);


            else

                % No temporal shift

                est_test = disp_est_gt_grid;

                gt_lag = gt_test;

            end


            %% -------------------------------------------------
            % Apply sign convention
            % --------------------------------------------------

            est_test = sign_factor * est_test;


            %% -------------------------------------------------
            % Calculate metrics
            % --------------------------------------------------

            m = compute_quality_metrics( ...
                est_test, ...
                gt_lag);


            %% -------------------------------------------------
            % Store results
            % --------------------------------------------------

            result_index = result_index + 1;

            result_interpretation{result_index,1} = ...
                candidate_names{c};

            result_sign(result_index,1) = ...
                sign_factor;

            result_lag(result_index,1) = ...
                lag;

            result_RMSE_um(result_index,1) = ...
                m.RMSE * 1e6;

            result_NRMSE_percent(result_index,1) = ...
                m.NRMSE_percent;

            result_MAE_um(result_index,1) = ...
                m.MAE * 1e6;

            result_bias_um(result_index,1) = ...
                m.bias * 1e6;

            result_correlation(result_index,1) = ...
                m.correlation;

            result_N(result_index,1) = ...
                m.N;


            %% -------------------------------------------------
            % Display
            % --------------------------------------------------

            fprintf([ ...
                '%-4s | sign %+d | lag %+d | ' ...
                'RMSE = %7.4f um | ' ...
                'NRMSE = %7.2f %% | ' ...
                'corr = %+7.4f\n'], ...
                candidate_names{c}, ...
                sign_factor, ...
                lag, ...
                m.RMSE*1e6, ...
                m.NRMSE_percent, ...
                m.correlation);

        end

    end

end


%% -------------------------------------------------------------
% 9. Create results table
% -------------------------------------------------------------

alignment_results = table( ...
    result_interpretation, ...
    result_sign, ...
    result_lag, ...
    result_RMSE_um, ...
    result_NRMSE_percent, ...
    result_MAE_um, ...
    result_bias_um, ...
    result_correlation, ...
    result_N, ...
    'VariableNames', { ...
        'Interpretation', ...
        'Sign', ...
        'Lag', ...
        'RMSE_um', ...
        'NRMSE_percent', ...
        'MAE_um', ...
        'Bias_um', ...
        'Correlation', ...
        'N'});


%% -------------------------------------------------------------
% 10. Find best result based on NRMSE
% -------------------------------------------------------------

alignment_sorted_NRMSE = sortrows( ...
    alignment_results, ...
    'NRMSE_percent', ...
    'ascend');


fprintf('\n============================================================\n');
fprintf(' BEST RESULTS BY NRMSE\n');
fprintf('============================================================\n\n');

disp(alignment_sorted_NRMSE( ...
    1:min(10,height(alignment_sorted_NRMSE)), :));


%% -------------------------------------------------------------
% 11. Find best result based on correlation
% -------------------------------------------------------------

alignment_sorted_corr = sortrows( ...
    alignment_results, ...
    'Correlation', ...
    'descend');


fprintf('\n============================================================\n');
fprintf(' BEST RESULTS BY CORRELATION\n');
fprintf('============================================================\n\n');

disp(alignment_sorted_corr( ...
    1:min(10,height(alignment_sorted_corr)), :));


%% -------------------------------------------------------------
% 12. Print best overall candidates
% -------------------------------------------------------------

best_NRMSE = alignment_sorted_NRMSE(1,:);

best_corr = alignment_sorted_corr(1,:);


fprintf('\n============================================================\n');
fprintf(' BEST NRMSE CANDIDATE\n');
fprintf('============================================================\n');

fprintf('Interpretation: %s\n', ...
    best_NRMSE.Interpretation{1});

fprintf('Sign:           %+d\n', ...
    best_NRMSE.Sign);

fprintf('Lag:            %+d frames\n', ...
    best_NRMSE.Lag);

fprintf('RMSE:           %.4f um\n', ...
    best_NRMSE.RMSE_um);

fprintf('NRMSE:          %.2f %%\n', ...
    best_NRMSE.NRMSE_percent);

fprintf('MAE:            %.4f um\n', ...
    best_NRMSE.MAE_um);

fprintf('Bias:           %.4f um\n', ...
    best_NRMSE.Bias_um);

fprintf('Correlation:    %.4f\n', ...
    best_NRMSE.Correlation);


fprintf('\n============================================================\n');
fprintf(' BEST CORRELATION CANDIDATE\n');
fprintf('============================================================\n');

fprintf('Interpretation: %s\n', ...
    best_corr.Interpretation{1});

fprintf('Sign:           %+d\n', ...
    best_corr.Sign);

fprintf('Lag:            %+d frames\n', ...
    best_corr.Lag);

fprintf('RMSE:           %.4f um\n', ...
    best_corr.RMSE_um);

fprintf('NRMSE:          %.2f %%\n', ...
    best_corr.NRMSE_percent);

fprintf('MAE:            %.4f um\n', ...
    best_corr.MAE_um);

fprintf('Bias:           %.4f um\n', ...
    best_corr.Bias_um);

fprintf('Correlation:    %.4f\n', ...
    best_corr.Correlation);


%% -------------------------------------------------------------
% 13. Plot correlation versus temporal lag
% -------------------------------------------------------------

figure;
hold on;

legend_entries = {};

for c = 1:length(candidate_names)

    for sign_factor = sign_candidates

        rows = ...
            strcmp(alignment_results.Interpretation, candidate_names{c}) & ...
            alignment_results.Sign == sign_factor;

        current = alignment_results(rows,:);

        current = sortrows(current, 'Lag');

        plot( ...
            current.Lag, ...
            current.Correlation, ...
            '-o', ...
            'LineWidth', 1.5);

        legend_entries{end+1} = sprintf( ...
            '%s, sign %+d', ...
            candidate_names{c}, ...
            sign_factor);

    end

end

xlabel('Temporal lag [frames]');
ylabel('Correlation coefficient');

title('Displacement temporal alignment');

legend( ...
    legend_entries, ...
    'Location', ...
    'best');

grid on;
hold off;


%% -------------------------------------------------------------
% 14. Plot NRMSE versus temporal lag
% -------------------------------------------------------------

figure;
hold on;

legend_entries = {};

for c = 1:length(candidate_names)

    for sign_factor = sign_candidates

        rows = ...
            strcmp(alignment_results.Interpretation, candidate_names{c}) & ...
            alignment_results.Sign == sign_factor;

        current = alignment_results(rows,:);

        current = sortrows(current, 'Lag');

        plot( ...
            current.Lag, ...
            current.NRMSE_percent, ...
            '-o', ...
            'LineWidth', 1.5);

        legend_entries{end+1} = sprintf( ...
            '%s, sign %+d', ...
            candidate_names{c}, ...
            sign_factor);

    end

end

xlabel('Temporal lag [frames]');
ylabel('NRMSE [%]');

title('Displacement error versus temporal alignment');

legend( ...
    legend_entries, ...
    'Location', ...
    'best');

grid on;
hold off;


%% -------------------------------------------------------------
% 15. Save diagnostic results
% -------------------------------------------------------------

save( ...
    'large_inclusions_alignment_diagnostics.mat', ...
    'alignment_results', ...
    'alignment_sorted_NRMSE', ...
    'alignment_sorted_corr', ...
    'best_NRMSE', ...
    'best_corr');

fprintf('\nDiagnostic results saved to:\n');
fprintf('large_inclusions_alignment_diagnostics.mat\n');

fprintf('\n============================================================\n');
fprintf(' ALIGNMENT DIAGNOSTICS COMPLETE\n');
fprintf('============================================================\n');