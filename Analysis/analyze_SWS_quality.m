function results = analyze_SWS_quality(config)

%ANALYZE_SWS_QUALITY

% Quantitative validation of phase-gradient SWS reconstruction against the

% known digital-phantom ground truth.

%

% Metrics:

%   - Inclusion/background mean, median, SD, bias, MAE, RMSE

%   - Global bias, MAE, RMSE

%   - CNR

%   - Edge accuracy (Dong et al., 2025)

%   - Edge sharpness (Dong et al., 2025)

%

% Dong et al. use MATLAB edge() with default settings. Their grid is

% isotropic, while the present beamformed scan is anisotropic. Therefore,

% the edge metrics are evaluated after resampling the SWS map and GT map to

% an isotropic grid. No Q_C metric or Q_C-based screening is used.

arguments

    config struct

end

fprintf('\n');

fprintf('========================================\n');

fprintf('SWS QUALITY ANALYSIS\n');

fprintf('========================================\n');

fprintf('Phantom: %s\n', config.phantom_case);

fprintf('Loading: %s\n\n', config.data_file);

%% 1. Dependencies

% k-Wave is not required when a precomputed ground-truth SWS map is
% already stored in the dataset. It is loaded only as a fallback in
% Section 6 when the ground-truth map must be created for the first time.

if isempty(which('edge')) || isempty(which('bwdist')) || isempty(which('bwperim'))

    error(['Edge accuracy and edge sharpness require Image Processing ' ...
        'Toolbox functions edge(), bwdist(), and bwperim().']);

end

%% 2. Load dataset

data = load(config.data_file);

required_fields = { ...
    'b_data', ...
    'displacement_estimation_autocorr', ...
    'shear_params'};

for k = 1:numel(required_fields)

    if ~isfield(data, required_fields{k})

        error('Dataset is missing required variable: %s', required_fields{k});

    end

end

b_data = data.b_data;

shear_params = data.shear_params;

%% ========================================================================

% Get pulse repetition frequency

% ========================================================================

if isfield(data, 'prf')

    PRF = ...
        double(data.prf);

elseif isfield(data, 'channel_data') && ...
        isprop(data.channel_data, 'PRF')

    PRF = ...
        double(data.channel_data.PRF);

else

    % Legacy synthetic datasets were generated with PRF = 3000 Hz.

    warning( ...
        ['PRF was not stored explicitly in this legacy dataset. ' ...
         'Using PRF = 3000 Hz.']);

    PRF = ...
        3000;

end

fprintf( ...
    'PRF: %.1f Hz\n', ...
    PRF);

target_frequency = double(shear_params.source_freq);

%% 3. Autocorrelation estimator parameters

if isfield(data, 'autocorr_params')

    autocorr_params = data.autocorr_params;

else

    warning(['autocorr_params not found in dataset. ' ...
        'Using legacy defaults z_gate=4, x_gate=2, packet_size=4.']);

    autocorr_params.z_gate = 4;

    autocorr_params.x_gate = 2;

    autocorr_params.packet_size = 4;

end

%% 4. Reshape displacement

scan = b_data.scan;

x_axis = double(scan.x_axis(:));

z_axis = double(scan.z_axis(:));

Nx = numel(x_axis);

Nz = numel(z_axis);

u_raw = double(data.displacement_estimation_autocorr.data);

Nt = numel(u_raw) / (Nz * Nx);

if mod(Nt, 1) ~= 0

    error(['Estimated displacement cannot be reshaped to scan dimensions ' ...
        '[%d x %d].'], Nz, Nx);

end

u_original = reshape(u_raw, Nz, Nx, Nt);

% Opposite sign convention relative to k-Wave GT.

u_original = -u_original;

%% 5. Reconstruct baseline SWS

fprintf('Running phase-gradient SWS reconstruction...\n');

sws_map = run_phase_gradient_pipeline( ...
    u_original, ...
    scan, ...
    PRF, ...
    target_frequency, ...
    config.sws_filter_speed);

%% 6. Load or create exact digital-phantom SWS map

% Prefer a ground-truth SWS map stored in the dataset. This makes the
% analysis portable: k-Wave is then needed only when generating a new
% phantom or when upgrading an older dataset that does not yet contain
% sws_gt.

if isfield(data, 'sws_gt')

    fprintf('Using stored ground-truth SWS map from dataset.\n');

    sws_gt = double(data.sws_gt);

    if ~isequal(size(sws_gt), [Nz Nx])
        error(['Stored sws_gt has size %s, but the beamformed scan is ' ...
            '[%d %d].'], mat2str(size(sws_gt)), Nz, Nx);
    end

    if isfield(data, 'sws_gt_info')
        gt_info = data.sws_gt_info;
    else
        gt_info = struct();
        gt_info.source = 'stored_sws_gt';
    end

else

    fprintf(['No stored ground-truth SWS map found.\n' ...
        'Creating it once using k-Wave...\n']);

    % k-Wave is needed only for this fallback path.
    if isempty(which('kWaveGrid'))

        candidate_kwave = { ...
            fullfile(getenv('HOME'), 'k-Wave'), ...
            fullfile(getenv('USERPROFILE'), 'Documents', 'MATLAB', ...
            'Add-Ons', 'Toolboxes', 'k-Wave')};

        for k = 1:numel(candidate_kwave)
            if isfolder(candidate_kwave{k})
                addpath(genpath(candidate_kwave{k}));
            end
        end

    end

    if isempty(which('kWaveGrid'))
        error(['k-Wave was not found and this dataset does not contain ' ...
            'a stored sws_gt map. Create sws_gt once on a machine with ' ...
            'k-Wave, then copy the updated MAT-file.']);
    end

    %% --------------------------------------------------------------------
    % Probe geometry
    % ---------------------------------------------------------------------

    if isfield(data, 'prb_center') && isfield(data, 'prb_theta')

        prb_center = data.prb_center;
        prb_theta = data.prb_theta;

    elseif strcmpi(config.phantom_case, 'large_inclusions')

        % Legacy large-inclusions dataset. This fallback is retained only
        % for exploratory legacy analysis. Regenerated datasets should save
        % prb_center and prb_theta explicitly.
        fprintf('\n');
        fprintf('Using documented legacy probe geometry.\n');
        fprintf('Source: Git commit 6c517d3\n');

        [kgrid_geometry, ~] = create_phantom(shear_params);

        prb_center = [ ...
            0, ...
            -20e-3, ...
            kgrid_geometry.z_vec(1) ...
        ];

        prb_theta = [ ...
            deg2rad(90), ...
            deg2rad(0) ...
        ];

    else

        error( ...
            ['Probe geometry is missing from dataset for %s. ' ...
             'No documented fallback geometry is available.'], ...
            config.phantom_case);

    end

    %% --------------------------------------------------------------------
    % Build and store exact digital-phantom SWS map
    % ---------------------------------------------------------------------

    [sws_gt, gt_info] = build_ground_truth_sws_map( ...
        shear_params, ...
        prb_center, ...
        prb_theta, ...
        x_axis, ...
        z_axis);

    sws_gt_info = gt_info;
    sws_gt_info.source = 'recreated_from_phantom';
    sws_gt_info.x_axis = x_axis;
    sws_gt_info.z_axis = z_axis;

    save( ...
        config.data_file, ...
        'sws_gt', ...
        'sws_gt_info', ...
        '-append');

    fprintf('Stored sws_gt in dataset for future analysis.\n');

end

%% 7. Build evaluation ROIs

[roi, roi_info] = build_sws_rois( ...
    sws_gt, ...
    x_axis, ...
    z_axis, ...
    shear_params.c_shear_bkg, ...
    shear_params.c_shear_incl, ...
    config.sws_roi_margin, ...
    config.sws_image_margin);

%% 8. Baseline SWS metrics

metrics = compute_sws_metrics( ...
    sws_map, ...
    sws_gt, ...
    roi);

%% 9. Dong et al. edge metrics

edge_metrics = compute_edge_metrics_dong( ...
    sws_map, ...
    sws_gt, ...
    x_axis, ...
    z_axis, ...
    shear_params.c_shear_bkg, ...
    shear_params.c_shear_incl);

metrics.edge_accuracy_mm = edge_metrics.edge_accuracy_mm;

metrics.edge_sharpness_mps_per_pixel = edge_metrics.edge_sharpness_mps_per_pixel;

metrics.edge_sharpness_mps_per_mm = edge_metrics.edge_sharpness_mps_per_mm;

metrics.detected_edge_pixels = edge_metrics.detected_edge_pixels;

%% 10. Print baseline metrics

fprintf('\n');

fprintf('BASELINE SWS METRICS\n');

fprintf('========================================\n');

fprintf('\nInclusion\n');

fprintf('  GT SWS:          %.4f m/s\n', ...
    shear_params.c_shear_incl);

fprintf('  Mean:            %.4f m/s\n', ...
    metrics.inclusion.mean);

fprintf('  Median:          %.4f m/s\n', ...
    metrics.inclusion.median);

fprintf('  SD:              %.4f m/s\n', ...
    metrics.inclusion.std);

fprintf('  Bias:            %.4f m/s\n', ...
    metrics.inclusion.bias);

fprintf('  Relative bias:   %.2f %%\n', ...
    metrics.inclusion.relative_bias_percent);

fprintf('  MAE:             %.4f m/s\n', ...
    metrics.inclusion.MAE);

fprintf('  RMSE:            %.4f m/s\n', ...
    metrics.inclusion.RMSE);
fprintf('\nBackground\n');

fprintf('  GT SWS:          %.4f m/s\n', ...
    shear_params.c_shear_bkg);

fprintf('  Mean:            %.4f m/s\n', ...
    metrics.background.mean);

fprintf('  Median:          %.4f m/s\n', ...
    metrics.background.median);

fprintf('  SD:              %.4f m/s\n', ...
    metrics.background.std);

fprintf('  Bias:            %.4f m/s\n', ...
    metrics.background.bias);

fprintf('  Relative bias:   %.2f %%\n', ...
    metrics.background.relative_bias_percent);

fprintf('  MAE:             %.4f m/s\n', ...
    metrics.background.MAE);

fprintf('  RMSE:            %.4f m/s\n', ...
    metrics.background.RMSE);
fprintf('\nGlobal\n');

fprintf('  Bias:            %.4f m/s\n', metrics.global.bias);

fprintf('  MAE:             %.4f m/s\n', metrics.global.MAE);

fprintf('  RMSE:            %.4f m/s\n', metrics.global.RMSE);

fprintf('\nImage quality\n');

fprintf('  CNR:                 %.4f\n', ...
    metrics.CNR);

fprintf('  Direct-KDE gCNR:     %.4f\n', ...
    metrics.gCNR_directKDE);

fprintf('  True contrast:       %.4f m/s\n', ...
    metrics.true_contrast_mps);

fprintf('  Estimated contrast:  %.4f m/s\n', ...
    metrics.estimated_contrast_mps);

fprintf('  Contrast recovery:   %.4f\n', ...
    metrics.contrast_recovery);

fprintf('  Edge accuracy:       %.4f mm\n', ...
    metrics.edge_accuracy_mm);

fprintf('  Edge sharpness:      %.6f m/s/pixel\n', ...
    metrics.edge_sharpness_mps_per_pixel);

fprintf('  Edge sharpness:      %.4f m/s/mm\n', ...
    metrics.edge_sharpness_mps_per_mm);

fprintf('  Edge pixels:         %d\n', ...
    metrics.detected_edge_pixels);

fprintf('========================================\n\n');
%% 11. Baseline table

baseline_table = table( ...
    shear_params.c_shear_incl, ...
    shear_params.c_shear_bkg, ...
    metrics.inclusion.mean, ...
    metrics.inclusion.median, ...
    metrics.inclusion.std, ...
    metrics.inclusion.bias, ...
    metrics.inclusion.relative_bias_percent, ...
    metrics.inclusion.MAE, ...
    metrics.inclusion.RMSE, ...
    metrics.background.mean, ...
    metrics.background.median, ...
    metrics.background.std, ...
    metrics.background.bias, ...
    metrics.background.relative_bias_percent, ...
    metrics.background.MAE, ...
    metrics.background.RMSE, ...
    metrics.global.bias, ...
    metrics.global.MAE, ...
    metrics.global.RMSE, ...
    metrics.true_contrast_mps, ...
    metrics.estimated_contrast_mps, ...
    metrics.contrast_recovery, ...
    metrics.CNR, ...
    metrics.gCNR_directKDE, ...
    metrics.edge_accuracy_mm, ...
    metrics.edge_sharpness_mps_per_pixel, ...
    metrics.edge_sharpness_mps_per_mm, ...
    metrics.detected_edge_pixels, ...
    'VariableNames', { ...
        'GT_Inclusion_SWS_mps', ...
        'GT_Background_SWS_mps', ...
        'InclusionMean_mps', ...
        'InclusionMedian_mps', ...
        'InclusionSD_mps', ...
        'InclusionBias_mps', ...
        'InclusionRelativeBias_percent', ...
        'InclusionMAE_mps', ...
        'InclusionRMSE_mps', ...
        'BackgroundMean_mps', ...
        'BackgroundMedian_mps', ...
        'BackgroundSD_mps', ...
        'BackgroundBias_mps', ...
        'BackgroundRelativeBias_percent', ...
        'BackgroundMAE_mps', ...
        'BackgroundRMSE_mps', ...
        'GlobalBias_mps', ...
        'GlobalMAE_mps', ...
        'GlobalRMSE_mps', ...
        'TrueContrast_mps', ...
        'EstimatedContrast_mps', ...
        'ContrastRecovery', ...
        'CNR', ...
        'DirectKDE_gCNR', ...
        'EdgeAccuracy_mm', ...
        'EdgeSharpness_mps_per_pixel', ...
        'EdgeSharpness_mps_per_mm', ...
        'DetectedEdgePixels'});

baseline_csv = fullfile( ...
    config.results_dir, ...
    [config.phantom_case '_SWS_baseline.csv']);

writetable(baseline_table, baseline_csv);
%% 12. Store baseline results

results = struct();
results.phantom_case = config.phantom_case;
results.metrics = metrics;
results.edge_metrics = edge_metrics;
results.gt_info = gt_info;
results.roi_info = roi_info;
results.roi = roi;
results.PRF = PRF;
results.target_frequency = target_frequency;
results.filter_speed = config.sws_filter_speed;
results.x_axis = x_axis;
results.z_axis = z_axis;
%% 13. Optional repeated SWS noise test

if config.run_sws_noise_test

    if ~isfield(data, 'channel_data')

        error(['SWS noise test requires channel_data in the dataset because ' ...
            'the autocorrelation displacement estimator uses it.']);

    end

    fprintf('\n');

    fprintf('========================================\n');

    fprintf('REPEATED SWS NOISE TEST\n');

    fprintf('SNR = %.1f dB\n', config.sws_noise_dB);

    fprintf('========================================\n');

    NR = config.n_sws_noise_realizations;

    snr_dB = config.sws_noise_dB;

    clean_b = b_data.data;

    signal_power = mean(abs(clean_b(:)).^2);

    noise_power = signal_power / (10^(snr_dB / 10));

    incl_rmse = nan(NR, 1);

    bkg_rmse = nan(NR, 1);

    global_rmse = nan(NR, 1);

    cnr_values = nan(NR, 1);

    incl_bias = nan(NR, 1);

    bkg_bias = nan(NR, 1);

    global_bias = nan(NR, 1);

    edge_accuracy = nan(NR, 1);

    edge_sharpness_pixel = nan(NR, 1);

    edge_sharpness_mm = nan(NR, 1);

    for r = 1:NR

        fprintf('Realization %d/%d\n', r, NR);

        rng(20000 + r, 'twister');

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

        u_noisy_raw = double(displacement_noisy.data);

        Nt_noisy = numel(u_noisy_raw) / (Nz * Nx);

        if mod(Nt_noisy, 1) ~= 0

            error('Unexpected noisy displacement dimensions.');

        end

        u_noisy = reshape(u_noisy_raw, Nz, Nx, Nt_noisy);

        u_noisy = -u_noisy;

        sws_noisy = run_phase_gradient_pipeline( ...
            u_noisy, ...
            scan, ...
            PRF, ...
            target_frequency, ...
            config.sws_filter_speed);

        m = compute_sws_metrics(sws_noisy, sws_gt, roi);

        em = compute_edge_metrics_dong( ...
            sws_noisy, ...
            sws_gt, ...
            x_axis, ...
            z_axis, ...
            shear_params.c_shear_bkg, ...
            shear_params.c_shear_incl);

        incl_rmse(r) = m.inclusion.RMSE;

        bkg_rmse(r) = m.background.RMSE;

        global_rmse(r) = m.global.RMSE;

        cnr_values(r) = m.CNR;

        incl_bias(r) = m.inclusion.bias;

        bkg_bias(r) = m.background.bias;

        global_bias(r) = m.global.bias;

        edge_accuracy(r) = em.edge_accuracy_mm;

        edge_sharpness_pixel(r) = em.edge_sharpness_mps_per_pixel;

        edge_sharpness_mm(r) = em.edge_sharpness_mps_per_mm;

    end

    % Restore clean data.

    b_data.data = clean_b;

    noise_table = table( ...
        snr_dB, ...
        NR, ...
        mean(incl_rmse, 'omitnan'), ...
        std(incl_rmse, 'omitnan'), ...
        mean(bkg_rmse, 'omitnan'), ...
        std(bkg_rmse, 'omitnan'), ...
        mean(global_rmse, 'omitnan'), ...
        std(global_rmse, 'omitnan'), ...
        mean(cnr_values, 'omitnan'), ...
        std(cnr_values, 'omitnan'), ...
        mean(edge_accuracy, 'omitnan'), ...
        std(edge_accuracy, 'omitnan'), ...
        mean(edge_sharpness_pixel, 'omitnan'), ...
        std(edge_sharpness_pixel, 'omitnan'), ...
        mean(edge_sharpness_mm, 'omitnan'), ...
        std(edge_sharpness_mm, 'omitnan'), ...
        mean(incl_bias, 'omitnan'), ...
        std(incl_bias, 'omitnan'), ...
        mean(bkg_bias, 'omitnan'), ...
        std(bkg_bias, 'omitnan'), ...
        mean(global_bias, 'omitnan'), ...
        std(global_bias, 'omitnan'), ...
        'VariableNames', { ...
            'SNR_dB', ...
            'N_realizations', ...
            'InclusionRMSE_mean_mps', ...
            'InclusionRMSE_std_mps', ...
            'BackgroundRMSE_mean_mps', ...
            'BackgroundRMSE_std_mps', ...
            'GlobalRMSE_mean_mps', ...
            'GlobalRMSE_std_mps', ...
            'CNR_mean', ...
            'CNR_std', ...
            'EdgeAccuracy_mean_mm', ...
            'EdgeAccuracy_std_mm', ...
            'EdgeSharpness_mean_mps_per_pixel', ...
            'EdgeSharpness_std_mps_per_pixel', ...
            'EdgeSharpness_mean_mps_per_mm', ...
            'EdgeSharpness_std_mps_per_mm', ...
            'InclusionBias_mean_mps', ...
            'InclusionBias_std_mps', ...
            'BackgroundBias_mean_mps', ...
            'BackgroundBias_std_mps', ...
            'GlobalBias_mean_mps', ...
            'GlobalBias_std_mps'});

    fprintf('\n');

    disp(noise_table);

    noise_csv = fullfile( ...
        config.results_dir, ...
        sprintf('%s_SWS_%gdB_repeated.csv', ...
        config.phantom_case, snr_dB));

    writetable(noise_table, noise_csv);

    results.noise_test.summary = noise_table;

    results.noise_test.inclusion_RMSE = incl_rmse;

    results.noise_test.background_RMSE = bkg_rmse;

    results.noise_test.global_RMSE = global_rmse;

    results.noise_test.CNR = cnr_values;

    results.noise_test.edge_accuracy_mm = edge_accuracy;

    results.noise_test.edge_sharpness_mps_per_pixel = edge_sharpness_pixel;

    results.noise_test.edge_sharpness_mps_per_mm = edge_sharpness_mm;

    results.noise_test.inclusion_bias = incl_bias;

    results.noise_test.background_bias = bkg_bias;

    results.noise_test.global_bias = global_bias;

end

%% 14. Save MAT results

save_file = fullfile( ...
    config.results_dir, ...
    [config.phantom_case '_SWS_quality_results.mat']);

save(save_file, 'results', 'sws_map', 'sws_gt');

fprintf('Saved SWS results:\n%s\n', save_file);

fprintf('Saved baseline table:\n%s\n', baseline_csv);

if config.run_sws_noise_test

    fprintf('Saved SWS noise summary:\n%s\n', noise_csv);

end

%% 15. Figures

[X, Z] = meshgrid(x_axis, z_axis);

figure;

imagesc(x_axis * 1e3, z_axis * 1e3, sws_gt);

axis image;

set(gca, 'YDir', 'normal');

colorbar;

xlabel('Lateral position [mm]');

ylabel('Depth [mm]');

title(sprintf('%s - Ground-truth SWS', ...
    strrep(config.phantom_case, '_', ' ')));

figure;

imagesc(x_axis * 1e3, z_axis * 1e3, sws_map);

axis image;

set(gca, 'YDir', 'normal');

colorbar;

xlabel('Lateral position [mm]');

ylabel('Depth [mm]');

title(sprintf('%s - Phase-gradient SWS', ...
    strrep(config.phantom_case, '_', ' ')));

figure;

imagesc(x_axis * 1e3, z_axis * 1e3, sws_map - sws_gt);

axis image;

set(gca, 'YDir', 'normal');

colorbar;

xlabel('Lateral position [mm]');

ylabel('Depth [mm]');

title('SWS error: estimate - ground truth');

figure;

imagesc(x_axis * 1e3, z_axis * 1e3, sws_map);

axis image;

set(gca, 'YDir', 'normal');

colorbar;

hold on;

contour(X * 1e3, Z * 1e3, roi.inclusion_core, [1 1], ...
    'k--', 'LineWidth', 1.2);

contour(X * 1e3, Z * 1e3, roi.background_core, [1 1], ...
    'w--', 'LineWidth', 1.2);

hold off;

xlabel('Lateral position [mm]');

ylabel('Depth [mm]');

title('SWS evaluation ROIs');

% Visualize GT and detected inclusion edges.

if edge_metrics.has_inclusion

    figure;

    imagesc( ...
        edge_metrics.x_iso_mm, ...
        edge_metrics.z_iso_mm, ...
        edge_metrics.sws_iso);

    axis image;

    set(gca, 'YDir', 'normal');

    colorbar;

    hold on;

    contour( ...
        edge_metrics.x_iso_mm, ...
        edge_metrics.z_iso_mm, ...
        double(edge_metrics.gt_edge), ...
        [0.5 0.5], ...
        'k--', ...
        'LineWidth', 1.5);

    contour( ...
        edge_metrics.x_iso_mm, ...
        edge_metrics.z_iso_mm, ...
        double(edge_metrics.estimated_edge), ...
        [0.5 0.5], ...
        'w-', ...
        'LineWidth', 1.0);

    hold off;

    xlabel('Lateral position [mm]');

    ylabel('Depth [mm]');

    title(sprintf( ...
        'Edge metrics: E_a = %.3f mm, E_s = %.3f m/s/mm', ...
        edge_metrics.edge_accuracy_mm, ...
        edge_metrics.edge_sharpness_mps_per_mm));

end

fprintf('\nSWS QUALITY ANALYSIS COMPLETE\n\n');

end

%% ========================================================================

% Local function: phase-gradient pipeline

% ========================================================================

function sws = run_phase_gradient_pipeline( ...
    displacement, ...
    scan, ...
    PRF, ...
    target_frequency, ...
    filter_speed)

elasto = uff.elastography();

elasto.PRF = PRF;

elasto.target_frequency = target_frequency;

elasto.scan = scan;

elasto.autocorrelation_field = uff.data_field();

elasto.autocorrelation_field.set_component( ...
    'u', ...
    reshape( ...
        displacement, ...
        size(displacement, 1), ...
        size(displacement, 2), ...
        1, ...
        size(displacement, 3)));

num_angles = 12;

elasto.directional_filtering_angles = ...
    2*pi*linspace(0, 1 - 1/num_angles, num_angles);

extractor = postprocess.elastography.extract_freq_component(elasto);

extractor.target_frequency = target_frequency;

extractor.component = 'u';

extractor.run();

omega = 2*pi*target_frequency;

elasto.radial_filter = uff.spatial_filter( ...
    'gaussian_charite', ...
    'radial');

elasto.radial_filter.sigma = filter_speed / omega;

elasto.angular_filter = uff.spatial_filter( ...
    'gaussian_charite', ...
    'angular');

elasto.angular_filter.deltaTheta = ...
    median(diff(elasto.directional_filtering_angles));

filter_operator = postprocess.elastography.apply_filter(elasto);

filter_operator.run();

pgm = postprocess.elastography.phase_gradient_method(elasto);

pgm.run();

sws = double(elasto.sws_map.u);

end

%% ========================================================================

% Local function: create ground-truth SWS map

% ========================================================================

function [sws_gt, info] = build_ground_truth_sws_map( ...
    shear_params, ...
    prb_center, ...
    prb_theta, ...
    x_axis, ...
    z_axis)

[kgrid, medium] = create_phantom(shear_params);

[TF, TF_rev] = coor_transformation(prb_center, prb_theta);

origin = [0; 0; 0; 1];

p1 = TF * origin;

p2 = TF_rev * origin;

err1 = norm(p1(1:3).' - double(prb_center(:).'));

err2 = norm(p2(1:3).' - double(prb_center(:).'));

if err1 <= err2

    T = TF;

    transform_mode = 1;

    origin_error = err1;

else

    T = TF_rev;

    transform_mode = 2;

    origin_error = err2;

end

fprintf('\nGROUND-TRUTH GEOMETRY\n');

fprintf('----------------------------------------\n');

fprintf('Probe center [mm]: [%.2f %.2f %.2f]\n', ...
    prb_center(1)*1e3, prb_center(2)*1e3, prb_center(3)*1e3);

fprintf('Probe theta [deg]: [%.2f %.2f]\n', ...
    rad2deg(prb_theta(1)), rad2deg(prb_theta(2)));

fprintf('Transform mode: %d\n', transform_mode);

fprintf('Origin error:   %.6g m\n', origin_error);

F = griddedInterpolant( ...
    { ...
        double(kgrid.x_vec(:)), ...
        double(kgrid.y_vec(:)), ...
        double(kgrid.z_vec(:)) ...
    }, ...
    double(medium.sound_speed_shear), ...
    'nearest', ...
    'none');

[X, Z] = meshgrid(x_axis, z_axis);

N = numel(X);

plane_maps = cell(3, 1);

plane_scores = zeros(3, 1);

plane_coverage = zeros(3, 1);

for mode = 1:3

    switch mode

        case 1

            local_points = [ ...
                X(:).'; ...
                zeros(1, N); ...
                Z(:).'; ...
                ones(1, N)];

        case 2

            local_points = [ ...
                X(:).'; ...
                Z(:).'; ...
                zeros(1, N); ...
                ones(1, N)];

        case 3

            local_points = [ ...
                zeros(1, N); ...
                X(:).'; ...
                Z(:).'; ...
                ones(1, N)];

    end

    global_points = T * local_points;

    values = F( ...
        global_points(1, :).', ...
        global_points(2, :).', ...
        global_points(3, :).');

    map_here = reshape(values, size(X));

    plane_maps{mode} = map_here;

    plane_coverage(mode) = nnz(isfinite(map_here));

    bkg = double(shear_params.c_shear_bkg);

    incl = double(shear_params.c_shear_incl);

    if abs(incl - bkg) > eps

        midpoint = 0.5 * (incl + bkg);

        if incl > bkg

            incl_pixels = isfinite(map_here) & map_here > midpoint;

        else

            incl_pixels = isfinite(map_here) & map_here < midpoint;

        end

        plane_scores(mode) = nnz(incl_pixels);

    else

        plane_scores(mode) = plane_coverage(mode);

    end

    fprintf('Plane mode %d: finite %.2f%%, inclusion pixels %d\n', ...
        mode, ...
        100 * plane_coverage(mode) / numel(map_here), ...
        plane_scores(mode));

end

if abs(shear_params.c_shear_incl - shear_params.c_shear_bkg) > eps

    [~, plane_mode] = max(plane_scores);

else

    [~, plane_mode] = max(plane_coverage);

end

sws_gt = plane_maps{plane_mode};

fprintf('Selected plane mode: %d\n', plane_mode);

fprintf('GT SWS range: %.4f to %.4f m/s\n', ...
    min(sws_gt(:), [], 'omitnan'), ...
    max(sws_gt(:), [], 'omitnan'));

fprintf('----------------------------------------\n\n');

info = struct();

info.transform_mode = transform_mode;

info.origin_error_m = origin_error;

info.plane_mode = plane_mode;

info.plane_scores = plane_scores;

info.plane_coverage = plane_coverage;

info.medium_size = size(medium.sound_speed_shear);

end

%% ========================================================================

% Local function: build inclusion/background ROIs

% ========================================================================

function [roi, info] = build_sws_rois( ...
    sws_gt, ...
    x_axis, ...
    z_axis, ...
    c_bkg, ...
    c_incl, ...
    roi_margin, ...
    image_margin)

[X, Z] = meshgrid(x_axis, z_axis);

finite_gt = isfinite(sws_gt);

if abs(c_incl - c_bkg) <= eps

    inclusion = false(size(sws_gt));

    background = finite_gt;

else

    midpoint = 0.5 * (double(c_bkg) + double(c_incl));

    if c_incl > c_bkg

        inclusion = finite_gt & sws_gt > midpoint;

    else

        inclusion = finite_gt & sws_gt < midpoint;

    end

    background = finite_gt & ~inclusion;

end

image_mask = ...
    X >= min(x_axis) + image_margin & ...
    X <= max(x_axis) - image_margin & ...
    Z >= min(z_axis) + image_margin & ...
    Z <= max(z_axis) - image_margin;

dx = median(abs(diff(x_axis)));

dz = median(abs(diff(z_axis)));

rx = max(0, ceil(roi_margin / dx));

rz = max(0, ceil(roi_margin / dz));

kernel = ones(2*rz + 1, 2*rx + 1);

kernel_n = numel(kernel);

if any(inclusion(:))

    inclusion_core = ...
        conv2(double(inclusion), kernel, 'same') >= kernel_n;

else

    inclusion_core = false(size(inclusion));

end

background_core = ...
    conv2(double(background), kernel, 'same') >= kernel_n;

inclusion_core = inclusion_core & image_mask;

background_core = background_core & image_mask;

roi = struct();

roi.inclusion = inclusion;

roi.background = background;

roi.inclusion_core = inclusion_core;

roi.background_core = background_core;

roi.image_mask = image_mask;

info = struct();

info.inclusion_pixels = nnz(inclusion);

info.background_pixels = nnz(background);

info.inclusion_core_pixels = nnz(inclusion_core);

info.background_core_pixels = nnz(background_core);

info.roi_margin_m = roi_margin;

info.image_margin_m = image_margin;

info.rx_pixels = rx;

info.rz_pixels = rz;

fprintf('SWS ROI SUMMARY\n');

fprintf('----------------------------------------\n');

fprintf('Inclusion pixels:       %d\n', info.inclusion_pixels);

fprintf('Background pixels:      %d\n', info.background_pixels);

fprintf('Inclusion core pixels:  %d\n', info.inclusion_core_pixels);

fprintf('Background core pixels: %d\n', info.background_core_pixels);

fprintf('----------------------------------------\n\n');

end

%% ========================================================================

% Local function: regional/global SWS metrics

% ========================================================================

function metrics = compute_sws_metrics(sws, sws_gt, roi)

valid_speed = ...
    isfinite(sws) & ...
    isfinite(sws_gt) & ...
    sws > 0 & ...
    sws < 10;

incl_mask = valid_speed & roi.inclusion_core;

bkg_mask = valid_speed & roi.background_core;

global_mask = valid_speed & roi.image_mask;

if nnz(incl_mask) == 0

    warning('No valid inclusion-core samples are available.');

end

if nnz(bkg_mask) == 0

    error('No valid background-core samples are available.');

end

metrics = struct();

metrics.inclusion = region_metrics(sws, sws_gt, incl_mask);

metrics.background = region_metrics(sws, sws_gt, bkg_mask);

metrics.global = region_metrics(sws, sws_gt, global_mask);

mu_i = metrics.inclusion.mean;

mu_b = metrics.background.mean;

sd_i = metrics.inclusion.std;

sd_b = metrics.background.std;

denom = sqrt(sd_i^2 + sd_b^2);

if isfinite(mu_i) && isfinite(mu_b) && denom > 0

    metrics.CNR = abs(mu_i - mu_b) / denom;

else

    metrics.CNR = NaN;

end

%% Contrast recovery

gt_i = metrics.inclusion.gt_mean;
gt_b = metrics.background.gt_mean;

true_contrast = gt_i - gt_b;
estimated_contrast = mu_i - mu_b;

metrics.true_contrast_mps = ...
    true_contrast;

metrics.estimated_contrast_mps = ...
    estimated_contrast;

if isfinite(true_contrast) && abs(true_contrast) > eps

    metrics.contrast_recovery = ...
        estimated_contrast / true_contrast;

else

    metrics.contrast_recovery = NaN;

end


%% Direct-KDE gCNR
% Exploratory application of gCNR to reconstructed SWS distributions.

pixels_i = double(sws(incl_mask));
pixels_b = double(sws(bkg_mask));

if ~isempty(pixels_i) && ~isempty(pixels_b)

    x_min = min([pixels_i; pixels_b]);
    x_max = max([pixels_i; pixels_b]);

    if x_max > x_min

        x_eval = ...
            linspace(x_min, x_max, 1000);

        f_i = ...
            ksdensity(pixels_i, x_eval);

        f_b = ...
            ksdensity(pixels_b, x_eval);

        % Normalize numerically so both represent PDFs.
        f_i = ...
            f_i ./ trapz(x_eval, f_i);

        f_b = ...
            f_b ./ trapz(x_eval, f_b);

        overlap = ...
            trapz( ...
                x_eval, ...
                min(f_i, f_b));

        metrics.gCNR_directKDE = ...
            1 - overlap;

    else

        metrics.gCNR_directKDE = 0;

    end

else

    metrics.gCNR_directKDE = NaN;

end


metrics.N_inclusion = nnz(incl_mask);
metrics.N_background = nnz(bkg_mask);
metrics.N_global = nnz(global_mask);
%% Generalized contrast-to-noise ratio (gCNR)
% Transform-KDE implementation from Arnestad et al.
% Default gamma = 0.7 is used for gCNR estimation.

pixels_i = double(sws(incl_mask));
pixels_b = double(sws(bkg_mask));


fprintf('\nSWS DISTRIBUTION CHECK\n');

fprintf('Inclusion: min %.4f, max %.4f m/s\n', ...
    min(pixels_i), max(pixels_i));

fprintf('Background: min %.4f, max %.4f m/s\n', ...
    min(pixels_b), max(pixels_b));

fprintf('Inclusion percentiles [1 5 50 95 99]:\n');
disp(prctile(pixels_i, [1 5 50 95 99]));

fprintf('Background percentiles [1 5 50 95 99]:\n');
disp(prctile(pixels_b, [1 5 50 95 99]));

if ~isempty(pixels_i) && ~isempty(pixels_b)

    %% Direct KDE gCNR sanity check

    x_min = min([pixels_i; pixels_b]);
    x_max = max([pixels_i; pixels_b]);

    x_eval = linspace(x_min, x_max, 1000);

    f_i_direct = ksdensity(pixels_i, x_eval);
    f_b_direct = ksdensity(pixels_b, x_eval);

    % Normalize numerically
    f_i_direct = f_i_direct ./ trapz(x_eval, f_i_direct);
    f_b_direct = f_b_direct ./ trapz(x_eval, f_b_direct);

    overlap_direct = trapz( ...
        x_eval, ...
        min(f_i_direct, f_b_direct));

    metrics.gCNR_directKDE = 1 - overlap_direct;

    fprintf('Direct-KDE gCNR sanity check: %.4f\n', ...
        metrics.gCNR_directKDE);

    %% Transform-KDE gCNR

    [metrics.gCNR, metrics.gCNR_info] = ...
        gCNR_tKDE( ...
            pixels_i, ...
            pixels_b, ...
            false);

else

    metrics.gCNR = NaN;
    metrics.gCNR_directKDE = NaN;
    metrics.gCNR_info = struct();

end

metrics.N_inclusion = nnz(incl_mask);

metrics.N_background = nnz(bkg_mask);

metrics.N_global = nnz(global_mask);

end

%% ========================================================================

% Local function: Dong edge accuracy / edge sharpness

% ========================================================================

function edge_metrics = compute_edge_metrics_dong( ...
    sws, ...
    sws_gt, ...
    x_axis, ...
    z_axis, ...
    c_bkg, ...
    c_incl)

%

% Dong et al. definitions:

%

%   Ea = r0/|dA| * sum_{P in dA} min_{P0 in dA0} ||P-P0||_2

%

%   Es = 1/|dA| * sum_{edge} |grad(SWS)|

%

% The article uses MATLAB edge() with default settings.

edge_metrics = struct();

edge_metrics.edge_accuracy_mm = NaN;

edge_metrics.edge_sharpness_mps_per_pixel = NaN;

edge_metrics.edge_sharpness_mps_per_mm = NaN;

edge_metrics.detected_edge_pixels = 0;

edge_metrics.has_inclusion = false;

edge_metrics.x_iso_mm = [];

edge_metrics.z_iso_mm = [];

edge_metrics.sws_iso = [];

edge_metrics.gt_edge = [];

edge_metrics.estimated_edge = [];

% No inclusion boundary exists for a homogeneous phantom.

if abs(c_incl - c_bkg) <= eps

    return;

end

dx = median(abs(diff(x_axis)));

dz = median(abs(diff(z_axis)));

% Finer native spacing used as isotropic evaluation spacing.

r0 = min(dx, dz);

x_iso = (min(x_axis):r0:max(x_axis)).';

z_iso = (min(z_axis):r0:max(z_axis)).';

[X_iso, Z_iso] = meshgrid(x_iso, z_iso);

% Estimated SWS: linear interpolation.

sws_iso = interp2( ...
    x_axis, ...
    z_axis, ...
    double(sws), ...
    X_iso, ...
    Z_iso, ...
    'linear', ...
    NaN);

% GT map: nearest interpolation preserves the discrete material classes.

gt_iso = interp2( ...
    x_axis, ...
    z_axis, ...
    double(sws_gt), ...
    X_iso, ...
    Z_iso, ...
    'nearest', ...
    NaN);

midpoint = 0.5 * (double(c_bkg) + double(c_incl));

if c_incl > c_bkg

    gt_inclusion = isfinite(gt_iso) & gt_iso > midpoint;

else

    gt_inclusion = isfinite(gt_iso) & gt_iso < midpoint;

end

if ~any(gt_inclusion(:))

    warning('No ground-truth inclusion was found for edge evaluation.');

    return;

end

edge_metrics.has_inclusion = true;

% Ground-truth boundary from the known phantom geometry.

gt_edge = bwperim(gt_inclusion, 8);

% Prepare SWS for edge(). No smoothing is introduced. Invalid values are

% replaced with nearest valid samples to avoid artificial NaN boundaries.

sws_edge_input = sws_iso;

invalid = ...
    ~isfinite(sws_edge_input) | ...
    sws_edge_input <= 0 | ...
    sws_edge_input >= 10;

sws_edge_input(invalid) = NaN;

sws_edge_input = fillmissing(sws_edge_input, 'nearest', 1);

sws_edge_input = fillmissing(sws_edge_input, 'nearest', 2);

sws_edge_input(~isfinite(sws_edge_input)) = double(c_bkg);

% -------------------------------------------------------------------------

% Normalize SWS image before edge detection

% -------------------------------------------------------------------------

%

% edge() operates more robustly when the image intensity range is

% normalized. This does not change the spatial structure of the SWS map;

% it only rescales the values for boundary detection.

%

valid_edge_values = ...
    sws_edge_input( ...
        isfinite(sws_edge_input));

if isempty(valid_edge_values)

    warning( ...
        'No valid SWS values available for edge detection.');

    return;

end

sws_min = ...
    min(valid_edge_values);

sws_max = ...
    max(valid_edge_values);

if sws_max <= sws_min

    warning( ...
        'SWS map has no dynamic range for edge detection.');

    return;

end

sws_edge_normalized = ...
    (sws_edge_input - sws_min) ./ ...
    (sws_max - sws_min);

sws_edge_normalized = ...
    max( ...
        0, ...
        min( ...
            1, ...
            sws_edge_normalized));

% Dong et al.: MATLAB edge() with default settings

estimated_edge_all = ...
    edge( ...
        sws_edge_normalized);

fprintf('\nEDGE DETECTION CHECK\n');

fprintf('All detected edge pixels: %d\n', nnz(estimated_edge_all));

% The paper refers specifically to the boundary of the inclusion. To avoid

% counting unrelated edges elsewhere in the SWS image, restrict the edge

% map to a reproducible bounding box around the known in-silico inclusion.

% The box is enlarged by 50% of the inclusion's largest dimension on each

% side. No distance-to-GT threshold is used, so spurious edges inside this

% region still penalize edge accuracy.

[gt_rows, gt_cols] = find(gt_inclusion);

row_min = min(gt_rows);

row_max = max(gt_rows);

col_min = min(gt_cols);

col_max = max(gt_cols);

height_px = row_max - row_min + 1;

width_px = col_max - col_min + 1;

margin_px = max(3, round(0.5 * max(height_px, width_px)));

row1 = max(1, row_min - margin_px);

row2 = min(size(gt_inclusion, 1), row_max + margin_px);

col1 = max(1, col_min - margin_px);

col2 = min(size(gt_inclusion, 2), col_max + margin_px);

inclusion_search_mask = false(size(gt_inclusion));

inclusion_search_mask(row1:row2, col1:col2) = true;

estimated_edge = estimated_edge_all & inclusion_search_mask;

figure;

imagesc( ...
    x_iso * 1e3, ...
    z_iso * 1e3, ...
    sws_edge_normalized);

axis image;

set(gca, 'YDir', 'normal');

colorbar;

hold on;

% All edges detected by MATLAB

contour( ...
    x_iso * 1e3, ...
    z_iso * 1e3, ...
    double(estimated_edge_all), ...
    [0.5 0.5], ...
    'w-', ...
    'LineWidth', ...
    1);

% Ground-truth inclusion boundary

contour( ...
    x_iso * 1e3, ...
    z_iso * 1e3, ...
    double(gt_edge), ...
    [0.5 0.5], ...
    'r--', ...
    'LineWidth', ...
    2);

% Search region

contour( ...
    x_iso * 1e3, ...
    z_iso * 1e3, ...
    double(inclusion_search_mask), ...
    [0.5 0.5], ...
    'g--', ...
    'LineWidth', ...
    1.5);

hold off;

xlabel('Lateral position [mm]');

ylabel('Depth [mm]');

title('Edge detection diagnostic');

legend( ...
    'Detected edges', ...
    'GT inclusion boundary', ...
    'Search region');

N_edge = nnz(estimated_edge);

edge_metrics.detected_edge_pixels = N_edge;

fprintf('Edges inside inclusion search region: %d\n', nnz(estimated_edge));

fprintf('----------------------------------------\n');

edge_metrics.x_iso_mm = x_iso * 1e3;

edge_metrics.z_iso_mm = z_iso * 1e3;

edge_metrics.sws_iso = sws_iso;

edge_metrics.gt_edge = gt_edge;

edge_metrics.estimated_edge = estimated_edge;

edge_metrics.r0_m = r0;

edge_metrics.r0_mm = r0 * 1e3;

edge_metrics.search_margin_pixels = margin_px;

if N_edge == 0

    warning(['No estimated inclusion edge was detected. ' ...
        'Edge accuracy and edge sharpness are NaN.']);

    return;

end

% Edge accuracy: bwdist gives Euclidean distance in pixels to nearest GT

% boundary pixel. Isotropic r0 converts directly to physical distance.

distance_to_gt_pixels = bwdist(gt_edge);

edge_accuracy_m = ...
    mean(distance_to_gt_pixels(estimated_edge), 'omitnan') * r0;

edge_metrics.edge_accuracy_mm = edge_accuracy_m * 1e3;

% Edge sharpness: source-like discrete spatial gradient magnitude.

[dSWS_dx_pixel, dSWS_dz_pixel] = gradient(sws_edge_input);

gradient_magnitude_pixel = sqrt( ...
    dSWS_dx_pixel.^2 + ...
    dSWS_dz_pixel.^2);

edge_sharpness_pixel = ...
    mean(gradient_magnitude_pixel(estimated_edge), 'omitnan');

edge_metrics.edge_sharpness_mps_per_pixel = edge_sharpness_pixel;

edge_metrics.edge_sharpness_mps_per_mm = ...
    edge_sharpness_pixel / edge_metrics.r0_mm;

end

%% ========================================================================

% Local function: regional error metrics

% ========================================================================

function m = region_metrics(sws, gt, mask)

v = double(sws(mask));
g = double(gt(mask));

valid = isfinite(v) & isfinite(g);

v = v(valid);
g = g(valid);

m = struct();

if isempty(v)

    m.N = 0;

    m.mean = NaN;
    m.median = NaN;
    m.std = NaN;

    m.gt_mean = NaN;

    m.bias = NaN;
    m.relative_bias_percent = NaN;

    m.MAE = NaN;
    m.RMSE = NaN;

    return;

end

e = v - g;

m.N = numel(v);

m.mean = mean(v);
m.median = median(v);
m.std = std(v);

m.gt_mean = mean(g);

m.bias = mean(e);
m.MAE = mean(abs(e));
m.RMSE = sqrt(mean(e.^2));

% Relative bias [%]
if abs(m.gt_mean) > eps

    m.relative_bias_percent = ...
        100 * m.bias / m.gt_mean;

else

    m.relative_bias_percent = NaN;

end

end
