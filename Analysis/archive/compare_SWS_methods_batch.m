function batch_summary = compare_SWS_methods_batch(phantom_cases)
%COMPARE_SWS_METHODS_BATCH  Automatic comparison of two SWS estimators.
%
% Uses the SAME saved displacement data, frequency extraction and directional
% filtering for phase-gradient and spatial-autocorrelation SWS estimation.
% Evaluates both maps with the SAME ground truth, baseline ROIs, and shared
% valid-pixel mask. The existing analysis scripts and datasets are READ-ONLY.
%
% Place this function in:
%   synthetic-the-rf-data-generation/Analysis/compare_SWS_methods_batch.m
%
% Example (first verify two cases):
%   compare_SWS_methods_batch
%
% Example (all currently configured nonhomogeneous phantoms):
%   compare_SWS_methods_batch(["chaoran_baseline", "low_contrast", ...
%       "high_contrast", "small_inclusions", "large_inclusions"])
%
% PREREQUISITE for each phantom:
%   Results/<case>/<case>_SWS_quality_results.mat, including results.roi,
%   sws_gt, and the original SWS baseline; and the displacement MAT dataset.
%
% Outputs go ONLY to:
%   Results/method_comparison/batch/<case>/
%   Results/method_comparison/batch/batch_summary.csv
%
% This is METHOD comparison at each dataset's original excitation frequency.
% It is NOT a frequency sweep or a validation of k-MDEV/full-waveform inversion.

    if nargin < 1 || isempty(phantom_cases)
        phantom_cases = ["high_contrast", "large_inclusions"];
    end
    phantom_cases = string(phantom_cases(:));
    assert(~isempty(phantom_cases), 'At least one phantom must be specified.');

    analysis_dir = fileparts(mfilename('fullpath'));
    repo_dir = fileparts(analysis_dir);
    master_dir = fileparts(repo_dir);
    data_root = fullfile(repo_dir, 'Data');
    results_root = fullfile(repo_dir, 'Results');
    batch_dir = fullfile(results_root, 'method_comparison', 'batch');

    % Match the working run_phantom_analysis.m setup.
    addpath(genpath(repo_dir));
    toolbox_dir = fullfile(master_dir, 'elastography', 'ElastographyToolbox');
    assert(isfolder(toolbox_dir), 'Missing ElastographyToolbox: %s', toolbox_dir);
    addpath(genpath(toolbox_dir));

    assert(~isempty(which('uff.elastography')), 'USTB uff.elastography not found.');
    assert(~isempty(which('postprocess.elastography.phase_gradient_method')), ...
        'USTB phase_gradient_method not found.');
    assert(~isempty(which('postprocess.elastography.autocorrelation_SWS_method')), ...
        'USTB autocorrelation_SWS_method not found.');
    assert(~isempty(which('ksdensity')) && ~isempty(which('edge')) && ...
        ~isempty(which('bwperim')) && ~isempty(which('bwdist')), ...
        'Required Statistics/Image Processing Toolbox functions are missing.');

    if ~isfolder(batch_dir), mkdir(batch_dir); end
    batch_summary = table();

    for case_idx = 1:numel(phantom_cases)
        case_name = char(phantom_cases(case_idx));
        fprintf('\n==============================================\n');
        fprintf('METHOD COMPARISON %d/%d: %s\n', ...
            case_idx, numel(phantom_cases), case_name);
        fprintf('==============================================\n');

        data_file = locate_phantom_data(data_root, repo_dir, analysis_dir, case_name);
        baseline_file = fullfile(results_root, case_name, ...
            [case_name '_SWS_quality_results.mat']);
        assert(isfile(baseline_file), ...
            ['Missing baseline result: %s\n' ...
             'First run run_phantom_analysis with RUN_SWS_ANALYSIS=true ' ...
             'for this phantom.'], baseline_file);

        % Read-only loads; NO saving into existing dataset or baseline result.
        D = load(data_file, 'b_data', 'displacement_estimation_autocorr', ...
            'shear_params', 'prf', 'channel_data', 'sws_gt');
        B = load(baseline_file, 'results', 'sws_gt', 'sws_map');
        assert(isfield(B, 'results') && isfield(B.results, 'roi') && ...
            isfield(B, 'sws_gt'), ...
            'Baseline file lacks required results.roi or sws_gt: %s', baseline_file);
        assert(isfield(D, 'b_data') && ...
            isfield(D, 'displacement_estimation_autocorr') && ...
            isfield(D, 'shear_params'), ...
            'Displacement data missing required fields: %s', data_file);

        GT = double(B.sws_gt);
        roi = B.results.roi;
        if isfield(D, 'sws_gt')
            assert(isequaln(double(D.sws_gt), GT), ...
                'Ground truth in dataset differs from baseline results for %s.', case_name);
        end

        scan = D.b_data.scan;
        x_axis = double(scan.x_axis(:));
        z_axis = double(scan.z_axis(:));
        Nz = numel(z_axis);
        Nx = numel(x_axis);
        assert(isequal(size(GT), [Nz, Nx]), ...
            'Ground truth size does not match scan axes for %s.', case_name);
        assert(all(diff(x_axis) > 0) && all(diff(z_axis) > 0), ...
            'Scan axes must be strictly increasing.');
        assert(all(isfield(roi, {'image_mask','inclusion_core','background_core'})), ...
            'Missing expected ROI masks in %s.', baseline_file);
        assert(isequal(size(roi.image_mask), size(GT)) && ...
            isequal(size(roi.inclusion_core), size(GT)) && ...
            isequal(size(roi.background_core), size(GT)), ...
            'ROI dimensions do not match ground truth.');

        if isfield(D, 'prf') && ~isempty(D.prf)
            PRF = double(D.prf);
        elseif isfield(D, 'channel_data') && ...
                isprop(D.channel_data, 'PRF')
            PRF = double(D.channel_data.PRF);
        else
            PRF = 3000;  % Same documented legacy fallback as main analysis.
            warning('No PRF stored for %s: using legacy PRF = 3000 Hz.', case_name);
        end
        target_frequency = double(D.shear_params.source_freq);
        filter_speed = 3.2;
        if isfield(B.results, 'filter_speed')
            filter_speed = double(B.results.filter_speed);
        end
        assert(isscalar(PRF) && isfinite(PRF) && PRF > 0, 'Invalid PRF.');
        assert(isscalar(target_frequency) && isfinite(target_frequency) && ...
            target_frequency > 0, 'Invalid excitation frequency.');

        u_raw = double(D.displacement_estimation_autocorr.data);
        assert(mod(numel(u_raw), Nz*Nx) == 0, ...
            'Displacement does not match scan dimensions.');
        Nt = numel(u_raw)/(Nz*Nx);
        assert(Nt >= 2, 'Displacement time series is too short.');
        u = -reshape(u_raw, Nz, Nx, Nt); % Match original sign convention.
        clear u_raw D; % Free the loaded ultrasound objects before reconstruction.

        % Construct one common preprocessed elastography object.
        [PG, AC] = estimate_both_sws(u, scan, PRF, ...
            target_frequency, filter_speed);
        clear u;
        assert(isequal(size(GT), size(PG), size(AC)), ...
            'Reconstructed SWS maps do not match the ground-truth grid.');

        if isfield(B, 'sws_map')
            baseline_PG = double(B.sws_map);
            if isequal(size(baseline_PG), size(PG))
                both_finite = isfinite(PG) & isfinite(baseline_PG);
                if any(both_finite(:))
                    max_difference = max(abs(PG(both_finite) - ...
                        baseline_PG(both_finite)));
                    fprintf('Regenerated PG vs stored baseline max difference: %.6g m/s\n', ...
                        max_difference);
                    if max_difference > 1e-3
                        warning(['Recomputed phase-gradient map differs from saved ' ...
                            'baseline. Verify toolbox version/settings before interpretation.']);
                    end
                end
            end
        end

        % Exact same valid samples for both methods (0 < SWS < 10 m/s).
        common_valid = isfinite(GT) & ...
            isfinite(PG) & PG > 0 & PG < 10 & ...
            isfinite(AC) & AC > 0 & AC < 10;
        mask_global = common_valid & logical(roi.image_mask);
        mask_incl = common_valid & logical(roi.inclusion_core);
        mask_bkg = common_valid & logical(roi.background_core);
        assert(nnz(mask_global) > 1 && nnz(mask_bkg) > 1, ...
            'Too few common valid global/background pixels in %s.', case_name);
        assert(nnz(mask_incl) > 1, ...
            'Too few common valid inclusion pixels in %s.', case_name);

        fprintf('Common ROI pixels: global=%d, inclusion=%d, background=%d\n', ...
            nnz(mask_global), nnz(mask_incl), nnz(mask_bkg));

        maps = {PG, AC};
        method_names = ["Phase gradient"; "Spatial autocorrelation"];
        metrics = struct([]);
        for k = 1:2
            m = quality_metrics_common(maps{k}, GT, mask_global, mask_incl, mask_bkg);
            edge_map = maps{k};
            edge_map(~common_valid) = NaN;
            m.edge = edge_metrics_common(edge_map, GT, x_axis, z_axis);
            if k == 1
                metrics = m;
            else
                metrics(k,1) = m;
            end
        end

        cases = repmat(string(case_name), 2, 1);
        this_table = table( ...
            cases, method_names, ...
            [metrics.global_rmse].', [metrics.inclusion_rmse].', ...
            [metrics.background_rmse].', [metrics.inclusion_bias].', ...
            [metrics.background_bias].', [metrics.CNR].', ...
            [metrics.gCNR_directKDE].', [metrics.contrast_recovery].', ...
            arrayfun(@(m) double(m.edge.edge_accuracy_mm), metrics), ...
            arrayfun(@(m) double(m.edge.edge_sharpness_mps_per_mm), metrics), ...
            arrayfun(@(m) double(m.edge.detected_edge_pixels), metrics), ...
            repmat(nnz(mask_global),2,1), ...
            repmat(nnz(mask_incl),2,1), ...
            repmat(nnz(mask_bkg),2,1), ...
            repmat(target_frequency,2,1), ...
            'VariableNames', { ...
                'Phantom', 'Method', 'GlobalRMSE_mps', 'InclusionRMSE_mps', ...
                'BackgroundRMSE_mps', 'InclusionBias_mps', 'BackgroundBias_mps', ...
                'CNR', 'DirectKDE_gCNR', 'ContrastRecovery', ...
                'EdgeAccuracy_mm', 'EdgeSharpness_mps_per_mm', ...
                'DetectedEdgePixels', 'GlobalPixels', 'InclusionPixels', ...
                'BackgroundPixels', 'Frequency_Hz'});
        % Avoid single/double MATLAB arrayfun issues by converting all
        % metrics to double before table assembly.
        disp(this_table);
        if isempty(batch_summary)
            batch_summary = this_table;
        else
            batch_summary = [batch_summary; this_table]; %#ok<AGROW>
        end

        out_dir = fullfile(batch_dir, case_name);
        if ~isfolder(out_dir), mkdir(out_dir); end
        writetable(this_table, fullfile(out_dir, 'method_metrics.csv'));
        save(fullfile(out_dir, 'method_comparison.mat'), ...
            'PG', 'AC', 'GT', 'roi', 'common_valid', ...
            'x_axis', 'z_axis', 'target_frequency', 'PRF', ...
            'filter_speed', 'metrics', 'this_table', '-v7.3');
        write_method_figures(out_dir, case_name, GT, PG, AC, ...
            x_axis, z_axis, mask_global);
        fprintf('Saved %s\n', out_dir);
    end

    writetable(batch_summary, fullfile(batch_dir, 'batch_summary.csv'));
    save(fullfile(batch_dir, 'batch_summary.mat'), 'batch_summary');
    fprintf('\nBATCH COMPLETE: %d phantoms.\nResults: %s\n', ...
        numel(phantom_cases), batch_dir);
end

function path = locate_phantom_data(data_root, repo_dir, analysis_dir, case_name)
    switch lower(case_name)
        case 'low_contrast'
            data_name = 'low_contrast_displacement_with_quality.mat';
        case 'large_inclusions'
            data_name = 'large_inclusions_displacement_with_quality.mat';
        case 'small_inclusions'
            data_name = 'small_inclusions_displacement_with_quality.mat';
        case 'high_contrast'
            data_name = 'high_contrast_displacement_with_quality.mat';
        case 'chaoran_baseline'
            data_name = 'chaoran_baseline_displacement_with_quality.mat';
        otherwise
            error('Unknown phantom case: %s', case_name);
    end
    candidates = {fullfile(data_root,case_name,data_name), ...
        fullfile(data_root,data_name), fullfile(analysis_dir,data_name), ...
        fullfile(repo_dir,data_name)};
    path = '';
    for k = 1:numel(candidates)
        if isfile(candidates{k})
            path = candidates{k};
            break;
        end
    end
    assert(~isempty(path), 'No displacement dataset found for %s (%s).', ...
        case_name, data_name);
end

function [PG, AC] = estimate_both_sws(u, scan, PRF, target_frequency, filter_speed)
    elasto = uff.elastography();
    elasto.PRF = PRF;
    elasto.target_frequency = target_frequency;
    elasto.scan = scan;
    elasto.autocorrelation_field = uff.data_field();
    elasto.autocorrelation_field.set_component('u', ...
        reshape(u, size(u,1), size(u,2), 1, size(u,3)));

    num_angles = 12; % Same as the existing phase-gradient pipeline.
    elasto.directional_filtering_angles = ...
        2*pi*linspace(0, 1-1/num_angles, num_angles);

    extractor = postprocess.elastography.extract_freq_component(elasto);
    extractor.target_frequency = target_frequency;
    extractor.component = 'u';
    extractor.run();

    omega = 2*pi*target_frequency;
    elasto.radial_filter = uff.spatial_filter('gaussian_charite', 'radial');
    elasto.radial_filter.sigma = filter_speed/omega;
    elasto.angular_filter = uff.spatial_filter('gaussian_charite', 'angular');
    elasto.angular_filter.deltaTheta = ...
        median(diff(elasto.directional_filtering_angles));

    filter_operator = postprocess.elastography.apply_filter(elasto);
    filter_operator.run();

    pgm = postprocess.elastography.phase_gradient_method(elasto);
    pgm.run();
    PG = squeeze(double(elasto.sws_map.u));

    acm = postprocess.elastography.autocorrelation_SWS_method(elasto);
    acm.run();
    assert(~isempty(elasto.sws_map), ...
        'Spatial autocorrelation returned no SWS map.');
    AC = squeeze(double(elasto.sws_map.u));
end

function write_method_figures(out_dir, case_name, GT, PG, AC, x_axis, z_axis, mask_global)
    x_mm = x_axis*1e3;
    z_mm = z_axis*1e3;
    maps = {GT, PG, AC};
    labels = {'Ground truth', 'Phase gradient', 'Spatial autocorrelation'};
    max_gt = max(GT(isfinite(GT)));
    color_max = max(8, ceil(max_gt));
    fig = figure('Visible','off', 'Name',[case_name ' - SWS methods']);
    tiledlayout(1,3,'TileSpacing','compact');
    for k = 1:3
        nexttile;
        imagesc(x_mm,z_mm,maps{k});
        axis image;
        set(gca,'YDir','normal');
        clim([0,color_max]); colorbar;
        title(labels{k}); xlabel('Lateral [mm]'); ylabel('Depth [mm]');
    end
    exportgraphics(fig,fullfile(out_dir,'method_maps.png'),'Resolution',300);
    close(fig);

    errors = {PG-GT, AC-GT};
    err_values = [errors{1}(mask_global); errors{2}(mask_global)];
    max_abs_error = max(abs(err_values));
    if ~isfinite(max_abs_error) || max_abs_error == 0
        max_abs_error = 1;
    end
    fig = figure('Visible','off', 'Name',[case_name ' - SWS errors']);
    tiledlayout(1,2,'TileSpacing','compact');
    for k = 1:2
        nexttile;
        imagesc(x_mm,z_mm,errors{k});
        axis image;
        set(gca,'YDir','normal');
        clim([-max_abs_error,max_abs_error]); colorbar;
        title([labels{k+1} ' - GT']);
        xlabel('Lateral [mm]'); ylabel('Depth [mm]');
    end
    exportgraphics(fig,fullfile(out_dir,'method_errors.png'),'Resolution',300);
    close(fig);
end

function m = quality_metrics_common(map, GT, global_mask, incl_mask, bkg_mask)
    m = struct();
    e_global = map(global_mask) - GT(global_mask);
    e_incl = map(incl_mask) - GT(incl_mask);
    e_bkg = map(bkg_mask) - GT(bkg_mask);
    m.global_rmse = double(sqrt(mean(e_global.^2)));
    m.inclusion_rmse = double(sqrt(mean(e_incl.^2)));
    m.background_rmse = double(sqrt(mean(e_bkg.^2)));
    m.inclusion_bias = double(mean(e_incl));
    m.background_bias = double(mean(e_bkg));

    incl = map(incl_mask);
    bkg = map(bkg_mask);
    denom = hypot(std(incl),std(bkg));
    m.CNR = NaN;
    if isfinite(denom) && denom > 0
        m.CNR = double(abs(mean(incl)-mean(bkg))/denom);
    end
    true_contrast = mean(GT(incl_mask))-mean(GT(bkg_mask));
    m.contrast_recovery = NaN;
    if isfinite(true_contrast) && abs(true_contrast)>eps
        m.contrast_recovery = double((mean(incl)-mean(bkg))/true_contrast);
    end

    % Same direct-KDE gCNR definition as the original comparison script.
    m.gCNR_directKDE = NaN;
    x_lo = min([incl; bkg]);
    x_hi = max([incl; bkg]);
    if isfinite(x_lo) && isfinite(x_hi) && x_hi>x_lo
        x_eval = linspace(x_lo,x_hi,1000);
        p_i = ksdensity(incl,x_eval);
        p_b = ksdensity(bkg,x_eval);
        norm_i = trapz(x_eval,p_i);
        norm_b = trapz(x_eval,p_b);
        if norm_i>0 && norm_b>0
            p_i = p_i/norm_i;
            p_b = p_b/norm_b;
            m.gCNR_directKDE = double(1-trapz(x_eval,min(p_i,p_b)));
        end
    end
end

function e = edge_metrics_common(sws, GT, x_axis, z_axis)
    % Identical methodology to compare_SWS_methods.m.
    e = struct('edge_accuracy_mm', NaN, ...
               'edge_sharpness_mps_per_mm', NaN, ...
               'detected_edge_pixels', 0);
    vals = unique(GT(isfinite(GT)));
    if numel(vals)<2, return; end
    c_bkg = min(vals); c_incl = max(vals);
    if abs(c_incl-c_bkg)<=eps, return; end

    r0 = min(median(diff(x_axis)),median(diff(z_axis)));
    x_iso = (min(x_axis):r0:max(x_axis)).';
    z_iso = (min(z_axis):r0:max(z_axis)).';
    [X_iso,Z_iso] = meshgrid(x_iso,z_iso);
    sws_iso = interp2(x_axis,z_axis,sws,X_iso,Z_iso,'linear',NaN);
    gt_iso  = interp2(x_axis,z_axis,GT,X_iso,Z_iso,'nearest',NaN);
    gt_inclusion = isfinite(gt_iso) & gt_iso>0.5*(c_bkg+c_incl);
    if ~any(gt_inclusion(:)), return; end
    gt_edge = bwperim(gt_inclusion,8);

    edge_input = sws_iso;
    invalid = ~isfinite(edge_input) | edge_input<=0 | edge_input>=10;
    edge_input(invalid) = NaN;
    edge_input = fillmissing(edge_input,'nearest',1);
    edge_input = fillmissing(edge_input,'nearest',2);
    edge_input(~isfinite(edge_input)) = c_bkg;
    m_min = min(edge_input(:));
    m_max = max(edge_input(:));
    if m_max<=m_min, return; end
    edge_normalized = (edge_input-m_min)/(m_max-m_min);
    edge_normalized = max(0,min(1,edge_normalized));
    estimated_edge_all = edge(edge_normalized);

    [rr,cc] = find(gt_inclusion);
    row_min=min(rr); row_max=max(rr);
    col_min=min(cc); col_max=max(cc);
    margin_px=max(3,round(0.5*max(row_max-row_min+1,col_max-col_min+1)));
    r1=max(1,row_min-margin_px);
    r2=min(size(gt_inclusion,1),row_max+margin_px);
    c1=max(1,col_min-margin_px);
    c2=min(size(gt_inclusion,2),col_max+margin_px);
    search_mask=false(size(gt_inclusion));
    search_mask(r1:r2,c1:c2)=true;
    estimated_edge=estimated_edge_all & search_mask;
    e.detected_edge_pixels=double(nnz(estimated_edge));
    if e.detected_edge_pixels==0, return; end

    dist_px=bwdist(gt_edge);
    e.edge_accuracy_mm=double(mean(dist_px(estimated_edge))*r0*1e3);
    [grad_x,grad_z]=gradient(edge_input);
    grad_magnitude_px=hypot(grad_x,grad_z);
    e.edge_sharpness_mps_per_mm=double( ...
        mean(grad_magnitude_px(estimated_edge))/(r0*1e3));
end
