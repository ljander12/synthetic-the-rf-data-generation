function summary_table = compare_SWS_methods()
%COMPARE_SWS_METHODS  Compare two SWS estimators on the same digital phantom.
%
% Reads the saved phase-gradient and spatial-autocorrelation SWS maps and
% evaluates both against the same ground truth, ROI masks and valid pixels.
% The existing phantom-analysis code and results are never modified.
%
% Place this file in:
%   synthetic-the-rf-data-generation/Analysis/compare_SWS_methods.m
% Then run from MATLAB:
%   compare_SWS_methods
%
% Required inputs (generated previously):
%   Results/method_comparison/high_contrast_method_comparison.mat
%   Results/high_contrast/high_contrast_SWS_quality_results.mat
%
% Output:
%   Results/method_comparison/high_contrast_method_metrics.csv
%   Results/method_comparison/high_contrast_method_maps.png
%   Results/method_comparison/high_contrast_method_errors.png
%   Results/method_comparison/high_contrast_method_metrics.mat

    clc;

    analysis_dir = fileparts(mfilename('fullpath'));
    repo_dir = fileparts(analysis_dir);
    results_dir = fullfile(repo_dir, 'Results');
    output_dir = fullfile(results_dir, 'method_comparison');

    comparison_file = fullfile(output_dir, ...
        'high_contrast_method_comparison.mat');
    baseline_file = fullfile(results_dir, 'high_contrast', ...
        'high_contrast_SWS_quality_results.mat');

    assert(isfile(comparison_file), ...
        ['Missing comparison MAT-file:\n%s\n' ...
         'Save sws_pg/sws_ac from the MATLAB debugger first.'], ...
        comparison_file);
    assert(isfile(baseline_file), ...
        ['Missing baseline SWS results:\n%s\n' ...
         'Run the standard high_contrast SWS analysis first.'], ...
        baseline_file);

    A = load(comparison_file, 'comparison');
    B = load(baseline_file, 'results', 'sws_gt');

    assert(isfield(A, 'comparison'), ...
        'Comparison MAT-file must contain a variable named comparison.');
    assert(isfield(B, 'results') && isfield(B.results, 'roi'), ...
        'Baseline result does not contain results.roi.');
    assert(isfield(B, 'sws_gt'), ...
        'Baseline SWS result is missing sws_gt.');

    C = A.comparison;
    required_fields = {'sws_gt', 'sws_phase_gradient', ...
        'sws_spatial_autocorr', 'x_axis', 'z_axis'};
    for k = 1:numel(required_fields)
        assert(isfield(C, required_fields{k}), ...
            'Comparison file lacks comparison.%s.', required_fields{k});
    end

    GT = double(C.sws_gt);
    PG = squeeze(double(C.sws_phase_gradient));
    AC = squeeze(double(C.sws_spatial_autocorr));
    roi = B.results.roi;
    x_axis = double(C.x_axis(:));
    z_axis = double(C.z_axis(:));

    assert(ismatrix(GT) && isequal(size(GT), size(PG), size(AC)), ...
        'Ground truth and reconstructed maps must have identical 2-D sizes.');
    assert(isequal(size(GT), [numel(z_axis), numel(x_axis)]), ...
        'SWS map size does not agree with x/z axes.');
    assert(isequaln(GT, double(B.sws_gt)), ...
        'Saved comparison ground truth differs from baseline ground truth.');
    assert(all(diff(x_axis) > 0) && all(diff(z_axis) > 0), ...
        'Scan axes must be strictly increasing.');

    mask_names = {'inclusion_core','background_core','image_mask'};
    for k = 1:numel(mask_names)
        name = mask_names{k};
        assert(isfield(roi, name), 'Missing ROI: %s.', name);
        assert(isequal(size(roi.(name)), size(GT)), ...
            'ROI %s has incorrect dimensions.', name);
    end

    % IMPORTANT: Fair comparison on exactly the SAME valid pixels.
    % The 0..10 m/s physical admissibility range matches the baseline code.
    common_valid = isfinite(GT) & ...
        isfinite(PG) & PG > 0 & PG < 10 & ...
        isfinite(AC) & AC > 0 & AC < 10;

    mask_global = common_valid & logical(roi.image_mask);
    mask_incl = common_valid & logical(roi.inclusion_core);
    mask_bkg = common_valid & logical(roi.background_core);

    assert(nnz(mask_global) > 0, 'No common valid global samples.');
    assert(nnz(mask_incl) > 1, 'Too few common valid inclusion samples.');
    assert(nnz(mask_bkg) > 1, 'Too few common valid background samples.');

    fprintf('\n========================================\n');
    fprintf('SWS METHOD COMPARISON - HIGH CONTRAST\n');
    fprintf('========================================\n');
    fprintf('Common global ROI pixels: %d\n', nnz(mask_global));
    fprintf('Common inclusion ROI pixels: %d\n', nnz(mask_incl));
    fprintf('Common background ROI pixels: %d\n', nnz(mask_bkg));

    method_names = ["Phase gradient"; "Spatial autocorrelation"];
    maps = {PG, AC};
    for k = 1:2
        metric = quality_metrics_common( ...
            maps{k}, GT, mask_global, mask_incl, mask_bkg);
        % Use identical valid-domain masking for each edge comparison.
        map_for_edges = maps{k};
        map_for_edges(~common_valid) = NaN;
        metric.edge = edge_metrics_common(map_for_edges, GT, x_axis, z_axis);
        if k == 1
            data = metric;
        else
            data(k, 1) = metric;
        end
    end

   
edge_accuracy = arrayfun( ...
    @(s) double(s.edge.edge_accuracy_mm), data);

edge_sharpness = arrayfun( ...
    @(s) double(s.edge.edge_sharpness_mps_per_mm), data);

edge_pixel_count = arrayfun( ...
    @(s) double(s.edge.detected_edge_pixels), data);


    summary_table = table( ...
        method_names, ...
        [data.global_rmse].', ...
        [data.inclusion_rmse].', ...
        [data.background_rmse].', ...
        [data.inclusion_bias].', ...
        [data.background_bias].', ...
        [data.CNR].', ...
        [data.gCNR_directKDE].', ...
        [data.contrast_recovery].', ...
        edge_accuracy(:), ...
        edge_sharpness(:), ...
        edge_pixel_count(:), ...
        'VariableNames', { ...
            'Method', 'GlobalRMSE_mps', 'InclusionRMSE_mps', ...
            'BackgroundRMSE_mps', 'InclusionBias_mps', ...
            'BackgroundBias_mps', 'CNR', 'DirectKDE_gCNR', ...
            'ContrastRecovery', 'EdgeAccuracy_mm', ...
            'EdgeSharpness_mps_per_mm', 'DetectedEdgePixels'});

    disp(summary_table);
    writetable(summary_table, fullfile(output_dir, ...
        'high_contrast_method_metrics.csv'));

    % Comparable figures: physical dimensions + identical color limits.
    x_mm = x_axis * 1e3;
    z_mm = z_axis * 1e3;
    labels = {'Ground truth', 'Phase gradient', 'Spatial autocorrelation'};
    plot_maps = {GT, PG, AC};

    f1 = figure('Name', 'High contrast - SWS comparison');
    tiledlayout(1, 3, 'TileSpacing', 'compact');
    for k = 1:3
        nexttile;
        imagesc(x_mm, z_mm, plot_maps{k});
        axis image;
        set(gca, 'YDir', 'normal');
        clim([0 8]);
        colorbar;
        title(labels{k});
        xlabel('Lateral [mm]');
        ylabel('Depth [mm]');
    end
    exportgraphics(f1, fullfile(output_dir, ...
        'high_contrast_method_maps.png'), 'Resolution', 300);

    f2 = figure('Name', 'High contrast - SWS estimation errors');
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    errors = {PG-GT, AC-GT};
    max_abs_error = max(abs([errors{1}(mask_global); errors{2}(mask_global)]));
    if ~isfinite(max_abs_error) || max_abs_error == 0
        max_abs_error = 1;
    end
    for k = 1:2
        nexttile;
        imagesc(x_mm, z_mm, errors{k});
        axis image;
        set(gca, 'YDir', 'normal');
        clim([-max_abs_error, max_abs_error]);
        colorbar;
        title(sprintf('%s - GT', char(method_names(k))));
        xlabel('Lateral [mm]');
        ylabel('Depth [mm]');
    end
    exportgraphics(f2, fullfile(output_dir, ...
        'high_contrast_method_errors.png'), 'Resolution', 300);

    save(fullfile(output_dir, 'high_contrast_method_metrics.mat'), ...
        'summary_table', 'data', 'common_valid', ...
        'mask_global', 'mask_incl', 'mask_bkg');

    fprintf('\nComparison completed. Outputs:\n%s\n', output_dir);
    fprintf(['NOTE: This is a comparison on one saved high-contrast phantom, ' ...
        'not a general performance ranking.\n']);
end

function m = quality_metrics_common(map, GT, global_mask, incl_mask, bkg_mask)
    m = struct();
    e_global = map(global_mask) - GT(global_mask);
    e_incl = map(incl_mask) - GT(incl_mask);
    e_bkg = map(bkg_mask) - GT(bkg_mask);

    m.global_rmse = sqrt(mean(e_global.^2));
    m.inclusion_rmse = sqrt(mean(e_incl.^2));
    m.background_rmse = sqrt(mean(e_bkg.^2));
    m.inclusion_bias = mean(e_incl);
    m.background_bias = mean(e_bkg);

    incl = map(incl_mask);
    bkg = map(bkg_mask);
    s_i = std(incl);
    s_b = std(bkg);
    denom = sqrt(s_i^2 + s_b^2);
    m.CNR = NaN;
    if denom > 0 && isfinite(denom)
        m.CNR = abs(mean(incl) - mean(bkg)) / denom;
    end

    true_contrast = mean(GT(incl_mask)) - mean(GT(bkg_mask));
    m.contrast_recovery = NaN;
    if isfinite(true_contrast) && abs(true_contrast) > eps
        m.contrast_recovery = (mean(incl) - mean(bkg)) / true_contrast;
    end

    % Direct-KDE gCNR (same definition used in baseline SWS analysis).
    m.gCNR_directKDE = NaN;
    x_lo = min([incl; bkg]);
    x_hi = max([incl; bkg]);
    if isfinite(x_lo) && isfinite(x_hi) && x_hi > x_lo
        x_eval = linspace(x_lo, x_hi, 1000);
        p_i = ksdensity(incl, x_eval);
        p_b = ksdensity(bkg, x_eval);
        norm_i = trapz(x_eval, p_i);
        norm_b = trapz(x_eval, p_b);
        if norm_i > 0 && norm_b > 0
            p_i = p_i / norm_i;
            p_b = p_b / norm_b;
            m.gCNR_directKDE = 1 - trapz(x_eval, min(p_i, p_b));
        end
    end
end

function e = edge_metrics_common(sws, GT, x_axis, z_axis)
% Mirror the existing Dong-style edge metric workflow, with the same
% isotropic grid, GT-relative search region and MATLAB edge() defaults.
    e = struct('edge_accuracy_mm', NaN, ...
               'edge_sharpness_mps_per_mm', NaN, ...
               'detected_edge_pixels', 0);

    vals = unique(GT(isfinite(GT)));
    if numel(vals) < 2
        return;
    end
    c_bkg = min(vals);
    c_incl = max(vals);
    if abs(c_incl - c_bkg) <= eps
        return;
    end

    dx = median(diff(x_axis));
    dz = median(diff(z_axis));
    r0 = min(dx, dz);
    x_iso = (min(x_axis):r0:max(x_axis)).';
    z_iso = (min(z_axis):r0:max(z_axis)).';
    [X_iso, Z_iso] = meshgrid(x_iso, z_iso);

    sws_iso = interp2(x_axis, z_axis, sws, ...
        X_iso, Z_iso, 'linear', NaN);
    gt_iso = interp2(x_axis, z_axis, GT, ...
        X_iso, Z_iso, 'nearest', NaN);

    gt_inclusion = isfinite(gt_iso) & ...
        gt_iso > 0.5*(c_bkg+c_incl);
    if ~any(gt_inclusion(:))
        return;
    end
    gt_edge = bwperim(gt_inclusion, 8);

    edge_input = sws_iso;
    invalid = ~isfinite(edge_input) | edge_input <= 0 | edge_input >= 10;
    edge_input(invalid) = NaN;
    edge_input = fillmissing(edge_input, 'nearest', 1);
    edge_input = fillmissing(edge_input, 'nearest', 2);
    edge_input(~isfinite(edge_input)) = c_bkg;

    m_min = min(edge_input(:));
    m_max = max(edge_input(:));
    if m_max <= m_min
        return;
    end
    edge_normalized = (edge_input - m_min) / (m_max - m_min);
    edge_normalized = max(0, min(1, edge_normalized));
    estimated_edge_all = edge(edge_normalized);

    [rr, cc] = find(gt_inclusion);
    row_min = min(rr); row_max = max(rr);
    col_min = min(cc); col_max = max(cc);
    height_px = row_max - row_min + 1;
    width_px = col_max - col_min + 1;
    margin_px = max(3, round(0.5*max(height_px,width_px)));
    r1 = max(1,row_min-margin_px);
    r2 = min(size(gt_inclusion,1),row_max+margin_px);
    c1 = max(1,col_min-margin_px);
    c2 = min(size(gt_inclusion,2),col_max+margin_px);
    search_mask = false(size(gt_inclusion));
    search_mask(r1:r2,c1:c2) = true;
    estimated_edge = estimated_edge_all & search_mask;

    e.detected_edge_pixels = nnz(estimated_edge);
    if e.detected_edge_pixels == 0
        return;
    end

    dist_px = bwdist(gt_edge);
    e.edge_accuracy_mm = mean(dist_px(estimated_edge)) * r0 * 1e3;
    [grad_x,grad_z] = gradient(edge_input);
    grad_magnitude_px = hypot(grad_x,grad_z);
    e.edge_sharpness_mps_per_mm = ...
        mean(grad_magnitude_px(estimated_edge)) / (r0 * 1e3);
end
