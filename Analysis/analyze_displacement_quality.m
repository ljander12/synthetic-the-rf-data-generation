function results = analyze_displacement_quality(config)
%ANALYZE_DISPLACEMENT_QUALITY
%
% Quantitative validation of autocorrelation displacement estimation against
% synthetic ground truth.
%
% Metrics:
%   - RMSE
%   - MAE
%   - median absolute error
%   - bias
%   - standard deviation of the error
%   - 95th / 99th percentile absolute error
%   - Pearson correlation with ground truth
%
% Optional:
%   - repeated AWGN noise sweep
%
% IMPORTANT:
%   disp_gt is already incremental displacement. Do NOT apply diff().
%

arguments
    config struct
end

fprintf('\n');
fprintf('========================================\n');
fprintf('DISPLACEMENT QUALITY ANALYSIS\n');
fprintf('========================================\n');
fprintf('Phantom: %s\n', config.phantom_case);
fprintf('Loading: %s\n\n', config.data_file);

%% 1. Load data
data = load(config.data_file);

required_fields = { ...
    'b_data', ...
    'disp_gt', ...
    'displacement_estimation_autocorr'};

for k = 1:numel(required_fields)
    if ~isfield(data, required_fields{k})
        error('Dataset is missing required variable: %s', required_fields{k});
    end
end

b_data = data.b_data;
disp_gt = double(data.disp_gt);

%% 2. Autocorrelation estimator parameters
if isfield(data, 'autocorr_params')
    autocorr_params = data.autocorr_params;
else
    warning(['autocorr_params not found in dataset. ' ...
        'Using legacy defaults z_gate=4, x_gate=2, packet_size=4.']);

    autocorr_params.z_gate = 4;
    autocorr_params.x_gate = 2;
    autocorr_params.packet_size = 4;
end

packet_size = autocorr_params.packet_size;
packet_intervals = packet_size - 1;

%% 3. Prepare ground truth
% disp_gt already contains incremental displacement.
% With packet_size = 4, one estimator output corresponds to the mean of
% three consecutive displacement intervals.

N_gt_frames = size(disp_gt, 3) - packet_intervals + 1;

if N_gt_frames < 1
    error('Not enough ground-truth frames for packet_size=%d.', packet_size);
end

gt_packet = zeros( ...
    size(disp_gt, 1), ...
    size(disp_gt, 2), ...
    N_gt_frames);

for k = 1:N_gt_frames
    gt_packet(:, :, k) = mean( ...
        disp_gt(:, :, k:(k + packet_intervals - 1)), ...
        3);
end

%% 4. Reshape estimated displacement
scan = b_data.scan;

x_est = double(scan.x_axis(:));
z_est = double(scan.z_axis(:));

Nx = numel(x_est);
Nz = numel(z_est);

est_raw = double(data.displacement_estimation_autocorr.data);
Nt_est = numel(est_raw) / (Nz * Nx);

if mod(Nt_est, 1) ~= 0
    error(['Estimated displacement cannot be reshaped to the beamformed ' ...
        'scan dimensions [%d x %d].'], Nz, Nx);
end

est_native = reshape(est_raw, Nz, Nx, Nt_est);

% Opposite sign convention relative to k-Wave ground truth.
est_native = -est_native;

N_frames = min(size(est_native, 3), size(gt_packet, 3));

est_native = est_native(:, :, 1:N_frames);
gt_packet = gt_packet(:, :, 1:N_frames);

%% 5. Interpolate estimate to ground-truth grid
x_gt = linspace(x_est(1), x_est(end), size(gt_packet, 2));
z_gt = linspace(z_est(1), z_est(end), size(gt_packet, 1));

[X_gt, Z_gt] = meshgrid(x_gt, z_gt);

est_gt_grid = nan(size(gt_packet));

for k = 1:N_frames
    est_gt_grid(:, :, k) = interp2( ...
        x_est, ...
        z_est, ...
        est_native(:, :, k), ...
        X_gt, ...
        Z_gt, ...
        'linear', ...
        NaN);
end

%% 6. Valid evaluation ROI
margin = config.edge_margin_samples;

spatial_mask = true(size(gt_packet, 1), size(gt_packet, 2));

if margin > 0
    if 2 * margin >= min(size(spatial_mask))
        error('edge_margin_samples is too large for the ground-truth grid.');
    end

    spatial_mask(1:margin, :) = false;
    spatial_mask(end-margin+1:end, :) = false;
    spatial_mask(:, 1:margin) = false;
    spatial_mask(:, end-margin+1:end) = false;
end

valid_mask = ...
    repmat(spatial_mask, 1, 1, N_frames) & ...
    isfinite(gt_packet) & ...
    isfinite(est_gt_grid);

%% 7. Baseline metrics
baseline = compute_displacement_metrics( ...
    est_gt_grid, ...
    gt_packet, ...
    valid_mask);

fprintf('BASELINE DISPLACEMENT METRICS\n');
fprintf('----------------------------------------\n');
fprintf('Valid samples:            %d\n', baseline.N);
fprintf('RMSE:                     %.5f um\n', baseline.RMSE_um);
fprintf('MAE:                      %.5f um\n', baseline.MAE_um);
fprintf('Median absolute error:    %.5f um\n', baseline.MedianAbsError_um);
fprintf('Bias:                     %.5f um\n', baseline.Bias_um);
fprintf('Error standard deviation: %.5f um\n', baseline.ErrorSD_um);
fprintf('P95 absolute error:       %.5f um\n', baseline.P95AbsError_um);
fprintf('P99 absolute error:       %.5f um\n', baseline.P99AbsError_um);
fprintf('Estimate-GT correlation:  %.5f\n', baseline.Correlation);
fprintf('----------------------------------------\n\n');

%% 8. Save baseline results
results = struct();
results.phantom_case = config.phantom_case;
results.baseline = baseline;
results.autocorr_params = autocorr_params;
results.n_frames = N_frames;

baseline_table = table( ...
    baseline.N, ...
    baseline.RMSE_um, ...
    baseline.MAE_um, ...
    baseline.MedianAbsError_um, ...
    baseline.Bias_um, ...
    baseline.ErrorSD_um, ...
    baseline.P95AbsError_um, ...
    baseline.P99AbsError_um, ...
    baseline.Correlation, ...
    'VariableNames', { ...
        'N', ...
        'RMSE_um', ...
        'MAE_um', ...
        'MedianAbsError_um', ...
        'Bias_um', ...
        'ErrorSD_um', ...
        'P95AbsError_um', ...
        'P99AbsError_um', ...
        'Correlation'});

baseline_csv = fullfile( ...
    config.results_dir, ...
    [config.phantom_case '_displacement_baseline.csv']);

writetable(baseline_table, baseline_csv);

%% 9. Optional repeated noise sweep
if config.run_noise_sweep

    if ~isfield(data, 'channel_data')
        error(['Noise sweep requires channel_data in the dataset because ' ...
            'the autocorrelation estimator uses it.']);
    end

    fprintf('\n');
    fprintf('========================================\n');
    fprintf('DISPLACEMENT NOISE SWEEP\n');
    fprintf('========================================\n');

    noise_levels_dB = config.noise_levels_dB(:);
    NR = config.n_noise_realizations;

    clean_b = b_data.data;
    signal_power = mean(abs(clean_b(:)).^2);

    n_levels = numel(noise_levels_dB);

    rmse_all = nan(n_levels, NR);
    mae_all = nan(n_levels, NR);
    bias_all = nan(n_levels, NR);
    corr_all = nan(n_levels, NR);

    for s = 1:n_levels

        snr_dB = noise_levels_dB(s);
        fprintf('\nSNR = %.1f dB\n', snr_dB);

        noise_power = signal_power / (10^(snr_dB / 10));

        for r = 1:NR

            rng(1000*s + r, 'twister');

            if isreal(clean_b)
                noise = sqrt(noise_power) .* randn(size(clean_b));
            else
                noise = sqrt(noise_power / 2) .* ...
                    (randn(size(clean_b)) + 1i * randn(size(clean_b)));
            end

            b_data.data = clean_b + noise;

            estimator = postprocess.autocorrelation_displacement_estimation();

            estimator.channel_data = data.channel_data;
            estimator.input = b_data;
            estimator.z_gate = autocorr_params.z_gate;
            estimator.x_gate = autocorr_params.x_gate;
            estimator.packet_size = autocorr_params.packet_size;

            displacement_noisy = estimator.go();

            noisy_raw = double(displacement_noisy.data);
            Nt_noisy = numel(noisy_raw) / (Nz * Nx);

            if mod(Nt_noisy, 1) ~= 0
                error('Unexpected noisy displacement dimensions.');
            end

            noisy_native = reshape(noisy_raw, Nz, Nx, Nt_noisy);
            noisy_native = -noisy_native;

            N_here = min(N_frames, size(noisy_native, 3));

            noisy_gt_grid = nan( ...
                size(gt_packet, 1), ...
                size(gt_packet, 2), ...
                N_here);

            for k = 1:N_here
                noisy_gt_grid(:, :, k) = interp2( ...
                    x_est, ...
                    z_est, ...
                    noisy_native(:, :, k), ...
                    X_gt, ...
                    Z_gt, ...
                    'linear', ...
                    NaN);
            end

            mask_here = ...
                repmat(spatial_mask, 1, 1, N_here) & ...
                isfinite(gt_packet(:, :, 1:N_here)) & ...
                isfinite(noisy_gt_grid);

            m = compute_displacement_metrics( ...
                noisy_gt_grid, ...
                gt_packet(:, :, 1:N_here), ...
                mask_here);

            rmse_all(s, r) = m.RMSE_um;
            mae_all(s, r) = m.MAE_um;
            bias_all(s, r) = m.Bias_um;
            corr_all(s, r) = m.Correlation;

            fprintf('  realization %2d/%2d: RMSE = %.4f um, corr = %.4f\n', ...
                r, NR, m.RMSE_um, m.Correlation);
        end
    end

    % Restore clean beamformed data.
    b_data.data = clean_b;

    noise_summary = table( ...
        noise_levels_dB, ...
        mean(rmse_all, 2, 'omitnan'), ...
        std(rmse_all, 0, 2, 'omitnan'), ...
        mean(mae_all, 2, 'omitnan'), ...
        std(mae_all, 0, 2, 'omitnan'), ...
        mean(bias_all, 2, 'omitnan'), ...
        std(bias_all, 0, 2, 'omitnan'), ...
        mean(corr_all, 2, 'omitnan'), ...
        std(corr_all, 0, 2, 'omitnan'), ...
        'VariableNames', { ...
            'SNR_dB', ...
            'RMSE_mean_um', ...
            'RMSE_std_um', ...
            'MAE_mean_um', ...
            'MAE_std_um', ...
            'Bias_mean_um', ...
            'Bias_std_um', ...
            'Correlation_mean', ...
            'Correlation_std'});

    fprintf('\n');
    disp(noise_summary);

    noise_csv = fullfile( ...
        config.results_dir, ...
        [config.phantom_case '_displacement_noise_summary.csv']);

    writetable(noise_summary, noise_csv);

    results.noise_levels_dB = noise_levels_dB;
    results.noise_summary = noise_summary;
    results.noise_raw.RMSE_um = rmse_all;
    results.noise_raw.MAE_um = mae_all;
    results.noise_raw.Bias_um = bias_all;
    results.noise_raw.Correlation = corr_all;
end

%% 10. Save MAT result
save_file = fullfile( ...
    config.results_dir, ...
    [config.phantom_case '_displacement_quality_results.mat']);

save(save_file, 'results');

fprintf('\nSaved displacement results:\n%s\n', save_file);
fprintf('Saved baseline table:\n%s\n', baseline_csv);

if config.run_noise_sweep
    fprintf('Saved noise summary:\n%s\n', noise_csv);
end

fprintf('\nDISPLACEMENT QUALITY ANALYSIS COMPLETE\n\n');

end

%% Local function: displacement metrics
function metrics = compute_displacement_metrics(est, gt, mask)

est_v = double(est(mask));
gt_v = double(gt(mask));

valid = isfinite(est_v) & isfinite(gt_v);

est_v = est_v(valid);
gt_v = gt_v(valid);

if isempty(est_v)
    error('No valid displacement samples remain for metric calculation.');
end

err = est_v - gt_v;
abs_err = abs(err);

metrics = struct();

metrics.N = numel(err);
metrics.RMSE_um = sqrt(mean(err.^2)) * 1e6;
metrics.MAE_um = mean(abs_err) * 1e6;
metrics.MedianAbsError_um = median(abs_err) * 1e6;
metrics.Bias_um = mean(err) * 1e6;
metrics.ErrorSD_um = std(err) * 1e6;
metrics.P95AbsError_um = percentile_linear(abs_err, 95) * 1e6;
metrics.P99AbsError_um = percentile_linear(abs_err, 99) * 1e6;

if numel(est_v) > 1
    R = corrcoef(est_v, gt_v);

    if all(size(R) == [2 2])
        metrics.Correlation = R(1, 2);
    else
        metrics.Correlation = NaN;
    end
else
    metrics.Correlation = NaN;
end

end

%% Local function: percentile without Statistics Toolbox
function q = percentile_linear(x, p)

x = sort(x(:));

if isempty(x)
    q = NaN;
    return;
end

if numel(x) == 1
    q = x;
    return;
end

position = 1 + (numel(x) - 1) * (p / 100);

lo = floor(position);
hi = ceil(position);

if lo == hi
    q = x(lo);
else
    alpha = position - lo;
    q = (1 - alpha) * x(lo) + alpha * x(hi);
end

end
