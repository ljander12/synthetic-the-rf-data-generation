function summary_table = compare_gt_displacement_sws()
%COMPARE_GT_DISPLACEMENT_SWS  Diagnostic: ideal vs ultrasound-estimated displacement.
%
% Runs the same phase-gradient SWS reconstruction on PRF-sampled k-Wave
% ground-truth displacement and compares it to the previously saved SWS
% reconstruction from ultrasound-estimated displacement.
%
% Read-only inputs (for the existing large_inclusions case):
%   large_inclusions_displacement_with_quality.mat
%   Results/large_inclusions/large_inclusions_SWS_quality_results.mat
%
% Outputs ONLY to Results/gt_displacement_diagnostic/large_inclusions/.
%
% IMPORTANT limitations:
%  - Ground truth has only 46x39 spatial points, so interpolation to the
%    dense ultrasound grid does NOT restore missing spatial information.
%  - The first N packet-averaged ground-truth frames are assumed to align
%    with the N estimated displacement frames, matching the current
%    analyze_displacement_quality.m methodology. Confirm temporal alignment
%    against acquisition indexing before using results as final evidence.
%  - This separates candidate error contributions; it does not uniquely
%    attribute the source of bias.

    phantom_case = 'large_inclusions';
    analysis_dir = fileparts(mfilename('fullpath'));
    repo_dir = fileparts(analysis_dir);
    toolbox_dir = fullfile(fileparts(repo_dir), ...
        'elastography', 'ElastographyToolbox');

    addpath(genpath(repo_dir));
    assert(isfolder(toolbox_dir), 'ElastographyToolbox not found: %s', toolbox_dir);
    addpath(genpath(toolbox_dir));
    assert(~isempty(which('create_phantom')), ...
        'create_phantom.m not on the MATLAB path.');
    assert(~isempty(which('kWaveGrid')), ...
        'k-Wave not on the MATLAB path; add the installed k-Wave folder first.');

    filename = [phantom_case '_displacement_with_quality.mat'];
    candidates = {fullfile(repo_dir, 'Data', phantom_case, filename), ...
        fullfile(repo_dir, 'Data', filename), ...
        fullfile(analysis_dir, filename), fullfile(repo_dir, filename)};
    data_file = '';
    for k = 1:numel(candidates)
        if isfile(candidates{k})
            data_file = candidates{k};
            break;
        end
    end
    assert(~isempty(data_file), 'Cannot locate %s.', filename);

    baseline_file = fullfile(repo_dir, 'Results', phantom_case, ...
        [phantom_case '_SWS_quality_results.mat']);
    assert(isfile(baseline_file), 'Missing baseline file: %s', baseline_file);

    fprintf('\n=== SWS diagnostic: GT vs estimated displacement ===\n');
    fprintf('Phantom: %s\n', phantom_case);
    D = load(data_file, 'disp_gt', 'b_data', 'shear_params', ...
        'autocorr_params', 'displacement_estimation_autocorr');
    B = load(baseline_file, 'sws_map', 'sws_gt', 'results');

    assert(isfield(D, 'disp_gt') && isfield(D, 'b_data') && ...
        isfield(D, 'displacement_estimation_autocorr') && ...
        isfield(D, 'shear_params'), 'Missing required displacement data.');
    assert(isfield(B, 'results') && isfield(B.results, 'roi') && ...
        isfield(B, 'sws_gt') && isfield(B, 'sws_map'), ...
        'Baseline file must contain results.roi, sws_gt, and sws_map.');
    assert(isfield(B.results, 'PRF') && ...
        isfield(B.results, 'target_frequency') && ...
        isfield(B.results, 'filter_speed'), ...
        'Baseline missing PRF, target frequency, or filter speed.');

    scan_original = D.b_data.scan;
    x_est_full = double(scan_original.x_axis(:));
    z_est_full = double(scan_original.z_axis(:));
    Nx = numel(x_est_full);
    Nz = numel(z_est_full);
    GT_SWS_full = double(B.sws_gt);
    PG_EST_full = double(B.sws_map);
    assert(isequal(size(GT_SWS_full), [Nz,Nx]) && ...
        isequal(size(PG_EST_full), [Nz,Nx]), ...
        'Baseline maps do not match ultrasound scan geometry.');

    if isfield(D, 'autocorr_params') && ...
            isfield(D.autocorr_params, 'packet_size')
        packet_size = double(D.autocorr_params.packet_size);
    else
        packet_size = 4;
        warning('Using legacy packet_size = 4; verify against simulation.');
    end
    packet_intervals = packet_size - 1;
    assert(packet_intervals >= 1 && packet_intervals == fix(packet_intervals), ...
        'Invalid autocorrelation packet size.');

    % The same estimated frame count as in analyze_displacement_quality.m.
    est_data = D.displacement_estimation_autocorr.data;
    assert(mod(numel(est_data), Nz*Nx) == 0, ...
        'Estimated displacement data do not match ultrasound grid.');
    N_est = numel(est_data)/(Nz*Nx);
    clear est_data;
    N_GT = size(D.disp_gt,3);
    assert(N_GT - packet_intervals + 1 >= N_est, ...
        'Not enough GT frames to align with estimated displacement.');

    % These constants are taken from the existing SIM_PHANTOM_MAIN settings:
    %   prb.N = 128; prb.pitch = 0.300e-3; scene_depth = 50e-3.
    % The coordinate vectors match get_GT_kwave.m exactly.
    probe_N = 128;
    probe_pitch_m = 0.300e-3;
    array_length_m = (probe_N - 1)*probe_pitch_m;
    scene_depth_m = 50e-3;
    [kgrid, ~] = create_phantom(D.shear_params); % Grid only; NO wave simulation.
    x_gt = (-array_length_m/2:kgrid.dx:array_length_m/2).';
    z_gt = (5e-3:kgrid.dz:scene_depth_m).';

    fprintf('GT field: %s, estimate: %dx%dx%d\n', ...
        mat2str(size(D.disp_gt)), Nz, Nx, N_est);
    fprintf('GT x: %.3f to %.3f mm; z: %.3f to %.3f mm\n', ...
        1e3*x_gt(1), 1e3*x_gt(end), 1e3*z_gt(1), 1e3*z_gt(end));
    fprintf('Scan x: %.3f to %.3f mm; z: %.3f to %.3f mm\n', ...
        1e3*x_est_full(1), 1e3*x_est_full(end), ...
        1e3*z_est_full(1), 1e3*z_est_full(end));

    assert(isequal(size(D.disp_gt, 1), numel(z_gt)) && ...
           isequal(size(D.disp_gt, 2), numel(x_gt)), ...
        ['Ground-truth array dimensions do not match reconstructed ' ...
         'coordinates. Confirm original probe pitch and k-Wave grid.']);
    assert(all(diff(x_gt)>0) && all(diff(z_gt)>0), ...
        'Ground-truth coordinates must be strictly increasing.');

    % Restrict BOTH pipelines to the spatial overlap before filtering.
    % The GT x-axis ends at +18.95 mm, whereas the ultrasound scan ends
    % at +19.05 mm. Using the common scan avoids extrapolation and makes
    % the PG reconstruction comparable for GT and estimated displacement.
    tol = 1e-9;  % numerical tolerance only (1 nanometre)
    keep_x = x_est_full >= x_gt(1)-tol & ...
             x_est_full <= x_gt(end)+tol;
    keep_z = z_est_full >= z_gt(1)-tol & ...
             z_est_full <= z_gt(end)+tol;
    assert(nnz(keep_x)>3 && nnz(keep_z)>3, ...
        'Ground-truth and ultrasound scan have insufficient spatial overlap.');

    x_est = x_est_full(keep_x);
    z_est = z_est_full(keep_z);
    % Clip only floating-point tolerance at the GT endpoints (NOT data
    % extrapolation). The physical comparison grid remains x_est/z_est.
    x_q = min(max(x_est, x_gt(1)), x_gt(end));
    z_q = min(max(z_est, z_gt(1)), z_gt(end));
    [Xq,Zq] = meshgrid(x_q,z_q);
    scan_common = uff.linear_scan('x_axis', x_est, 'z_axis', z_est);

    fprintf(['Using common scan: %d x %d (removed %d lateral and ' ...
        '%d axial pixels outside GT coverage).\n'], ...
        numel(z_est),numel(x_est),Nx-nnz(keep_x),Nz-nnz(keep_z));
    GT_SWS = GT_SWS_full(keep_z,keep_x);
    PG_EST_reference = PG_EST_full(keep_z,keep_x);
    clear GT_SWS_full PG_EST_full;

    % Recreate the estimated-displacement PG map on the same CROPPED scan.
    % Reusing the saved full-scan PG map would produce unequal boundary
    % conditions during directional filtering, biasing the comparison.
    u_EST_full = reshape(double(D.displacement_estimation_autocorr.data), ...
                         Nz, Nx, N_est);
    u_EST = -u_EST_full(keep_z,keep_x,:);  % existing sign convention
    clear u_EST_full;

    % Packet averaging matches the existing displacement quality analysis.
    % GT and estimated displacement already represent increments per PRF.
    gt_packet = zeros(size(D.disp_gt,1), size(D.disp_gt,2), N_est);
    for t = 1:N_est
        gt_packet(:,:,t) = mean( ...
            double(D.disp_gt(:,:,t:(t+packet_intervals-1))), 3);
    end
    clear D;

    u_GT = zeros(numel(z_est),numel(x_est),N_est);
    for t = 1:N_est
        u_GT(:,:,t) = interp2(x_gt,z_gt,gt_packet(:,:,t), ...
            Xq,Zq,'linear',NaN);
    end
    clear gt_packet;
    coverage = nnz(isfinite(u_GT))/numel(u_GT);
    fprintf('GT displacement interpolation coverage: %.2f %%\n', coverage*100);
    assert(coverage == 1, ...
        'GT displacement interpolation has NaNs; inspect spatial axes.');

    PRF = double(B.results.PRF);
    target_frequency = double(B.results.target_frequency);
    filter_speed = double(B.results.filter_speed);
    assert(isfinite(PRF) && PRF>0 && isfinite(target_frequency) && ...
        target_frequency>0 && isfinite(filter_speed) && filter_speed>0, ...
        'Invalid stored reconstruction parameters.');

    fprintf('Reconstructing GT displacement with phase gradient...\n');
    PG_GT = phase_gradient_reconstruct(u_GT, scan_common, ...
        PRF, target_frequency, filter_speed);
    clear u_GT;
    fprintf('Reconstructing estimated displacement on the SAME scan...\n');
    PG_EST = phase_gradient_reconstruct(u_EST, scan_common, ...
        PRF, target_frequency, filter_speed);
    clear u_EST;
    assert(isequal(size(PG_GT),size(PG_EST),size(GT_SWS)), ...
        'Reconstructed SWS maps do not agree with common scan size.');

    % Diagnostic only: cropping changes the filtering boundary conditions.
    reference_valid = isfinite(PG_EST) & isfinite(PG_EST_reference);
    if any(reference_valid(:))
        fprintf('Cropped-vs-original PG difference RMSE: %.5f m/s\n', ...
            sqrt(mean((PG_EST(reference_valid) - ...
                       PG_EST_reference(reference_valid)).^2)));
    end
    clear PG_EST_reference;

    % Fair comparison: same ROI and same set of valid SWS pixels.
    common = isfinite(GT_SWS) & ...
        isfinite(PG_GT) & PG_GT>0 & PG_GT<10 & ...
        isfinite(PG_EST) & PG_EST>0 & PG_EST<10;
    roi_full = B.results.roi;
    roi = struct();
    roi.image_mask = roi_full.image_mask(keep_z,keep_x);
    roi.inclusion_core = roi_full.inclusion_core(keep_z,keep_x);
    roi.background_core = roi_full.background_core(keep_z,keep_x);
    mask_names = {'image_mask', 'inclusion_core', 'background_core'};
    for k=1:numel(mask_names)
        assert(isfield(roi,mask_names{k}) && ...
            isequal(size(roi.(mask_names{k})), size(GT_SWS)), ...
            'Missing or mismatched ROI: %s', mask_names{k});
    end
    masks = {common & roi.image_mask, ...
             common & roi.inclusion_core, ...
             common & roi.background_core};
    assert(all(cellfun(@nnz,masks)>1), ...
        'Not enough common valid pixels in all ROIs.');

    maps = {PG_GT,PG_EST};
    method = ["GT displacement";"Estimated displacement"];
    bias = nan(2,3);
    rmse = nan(2,3);
    for j=1:2
        for k=1:3
            e=maps{j}(masks{k})-GT_SWS(masks{k});
            bias(j,k)=mean(e);
            rmse(j,k)=sqrt(mean(e.^2));
        end
    end
    summary_table = table(method, ...
        rmse(:,1),rmse(:,2),rmse(:,3), ...
        bias(:,1),bias(:,2),bias(:,3), ...
        'VariableNames',{'Method','GlobalRMSE_mps', ...
        'InclusionRMSE_mps','BackgroundRMSE_mps', ...
        'GlobalBias_mps','InclusionBias_mps','BackgroundBias_mps'});
    fprintf('Common valid pixels: global=%d, inclusion=%d, background=%d\n', ...
        nnz(masks{1}),nnz(masks{2}),nnz(masks{3}));
    disp(summary_table);

    out_dir=fullfile(repo_dir,'Results','gt_displacement_diagnostic',phantom_case);
    if ~isfolder(out_dir),mkdir(out_dir);end
    writetable(summary_table,fullfile(out_dir,'gt_vs_estimated_metrics.csv'));
    save(fullfile(out_dir,'gt_vs_estimated_maps.mat'), ...
        'PG_GT','PG_EST','GT_SWS','roi','x_est','z_est', ...
        'summary_table','PRF','target_frequency','filter_speed', ...
        'keep_x','keep_z','-v7.3');
    fig=figure('Name','GT vs estimated displacement - phase gradient');
    map_titles={'SWS ground truth', 'From GT displacement', ...
        'From estimated displacement'};
    figures={GT_SWS,PG_GT,PG_EST};
    tiledlayout(1,3,'TileSpacing','compact');
    for k=1:3
        nexttile;
        imagesc(1e3*x_est,1e3*z_est,figures{k});
        axis image;set(gca,'YDir','normal');clim([0 8]);colorbar;
        title(map_titles{k});xlabel('Lateral [mm]');ylabel('Depth [mm]');
    end
    exportgraphics(fig,fullfile(out_dir,'gt_vs_estimated_maps.png'), ...
        'Resolution',300);
    fprintf('Results saved to:\n%s\n',out_dir);
end

function sws = phase_gradient_reconstruct(displacement,scan,PRF,target_frequency,filter_speed)
% Exact reconstruction sequence used by analyze_SWS_quality.m.
    elasto=uff.elastography();
    elasto.PRF=PRF;
    elasto.target_frequency=target_frequency;
    elasto.scan=scan;
    elasto.autocorrelation_field=uff.data_field();
    elasto.autocorrelation_field.set_component('u', ...
        reshape(displacement,size(displacement,1), ...
        size(displacement,2),1,size(displacement,3)));
    num_angles=12;
    elasto.directional_filtering_angles= ...
        2*pi*linspace(0,1-1/num_angles,num_angles);
    extractor=postprocess.elastography.extract_freq_component(elasto);
    extractor.target_frequency=target_frequency;
    extractor.component='u';
    extractor.run();
    omega=2*pi*target_frequency;
    elasto.radial_filter=uff.spatial_filter('gaussian_charite','radial');
    elasto.radial_filter.sigma=filter_speed/omega;
    elasto.angular_filter=uff.spatial_filter('gaussian_charite','angular');
    elasto.angular_filter.deltaTheta= ...
        median(diff(elasto.directional_filtering_angles));
    filter_operator=postprocess.elastography.apply_filter(elasto);
    filter_operator.run();
    pgm=postprocess.elastography.phase_gradient_method(elasto);
    pgm.run();
    sws=squeeze(double(elasto.sws_map.u));
end
