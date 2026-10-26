function demo_elastography_master(mode)
%DEMO_ELASTOGRAPHY_MASTER  One-click demonstration of the MSc THE workflow.
%
% Put this file in synthetic-the-rf-data-generation/Analysis/.
%
% FAST MEETING DEMO (reads existing results; no re-simulation):
%   demo_elastography_master
%
% OPTIONAL LIVE DEMO (runs two SWS estimators on high_contrast again):
%   demo_elastography_master('live')
%
% WHAT IS SHOWN:
%   1. Ground truth / phase gradient / spatial autocorrelation (two phantoms)
%   2. Method comparison metrics for the four simulated phantoms
%   3. Stress tests: global bias, map-domain noise, blur, sharpness
%   4. GT displacement vs estimated displacement (diagnostic)
%   5. Stored displacement quality metrics, if available
%
% IMPORTANT:
% - Does not run k-Wave / Field II or alter original datasets.
% - Standard mode only reads previously saved analysis output.
% - 'live' mode performs USTB SWS reconstruction from saved displacement
%   data and writes results to Results/meeting_demo ONLY.
% - This is NOT a frequency sweep, k-MDEV/full-waveform comparison, or
%   conclusive proof of the source of SWS bias.
%
% Tested against the file/variable conventions of the user's existing
% compare_SWS_methods_batch.m and compare_gt_displacement_sws.m.
% The script itself has not been executed in MATLAB in this environment.

    if nargin < 1 || isempty(mode)
        mode = 'show';
    end
    mode = lower(string(mode));
    assert(isscalar(mode) && any(mode == ["show","live"]), ...
        'Use demo_elastography_master or demo_elastography_master(''live'').');

    analysis_dir = fileparts(mfilename('fullpath'));
    repo_dir = fileparts(analysis_dir);
    results_dir = fullfile(repo_dir, 'Results');
    batch_dir = fullfile(results_dir, 'method_comparison', 'batch');

    assert(isfolder(results_dir), 'Missing Results folder: %s', results_dir);
    addpath(analysis_dir);

    fprintf('\n=====================================================\n');
    fprintf('  MASTER THESIS - TIME-HARMONIC ELASTOGRAPHY\n');
    fprintf('  LIVE MATLAB DEMONSTRATION\n');
    fprintf('=====================================================\n');
    fprintf('Repository: %s\n', repo_dir);
    fprintf('Mode: %s\n\n', upper(char(mode)));
    fprintf('Pipeline: simulated phantom -> RF/beamforming -> displacement\n');
    fprintf('          -> harmonic extraction/directional filtering -> SWS\n');
    fprintf('          -> quality metrics and controlled stress testing\n\n');

    % The optional live step is intentionally independent of batch_summary.csv
    % to avoid replacing the existing four-phantom evaluation with one case.
    if mode == "live"
        fprintf('[LIVE] Re-running two estimators on stored high-contrast displacement...\n');
        try
            run_live_reconstruction(repo_dir, results_dir, batch_dir);
        catch ME
            warning('Live reconstruction failed, but saved-result demo continues: %s', ME.message);
        end
    end

    % Section 1: visual maps (saved batch output from the working pipeline).
    fprintf('\n[1/5] SWS ground truth + two reconstruction methods\n');
    try
        show_method_maps(batch_dir);
    catch ME
        warning('Method map overview was skipped: %s', ME.message);
    end

    % Section 2: consolidated quantitative evaluation.
    fprintf('\n[2/5] Method comparison across existing phantom cases\n');
    try
        show_method_metrics(batch_dir);
    catch ME
        warning('Method metrics were skipped: %s', ME.message);
    end

    % Section 3: saved independent quality-metric stress tests.
    fprintf('\n[3/5] Controlled metric stress tests\n');
    try
        show_stress_tests(results_dir);
    catch ME
        warning('Stress-test overview was skipped: %s', ME.message);
    end

    % Section 4: GT displacement vs estimated displacement diagnostic.
    fprintf('\n[4/5] GT displacement vs estimated displacement diagnostic\n');
    try
        show_displacement_diagnostic(results_dir);
    catch ME
        warning('GT displacement diagnostic was skipped: %s', ME.message);
    end

    % Section 5: existing displacement quality table (no computation).
    fprintf('\n[5/5] Stored displacement quality metrics\n');
    disp_csv = fullfile(results_dir,'large_inclusions', ...
        'large_inclusions_displacement_baseline.csv');
    if isfile(disp_csv)
        disp(readtable(disp_csv));
    else
        fprintf('No displacement baseline CSV found at: %s\n', disp_csv);
    end

    fprintf('\n=====================================================\n');
    fprintf('  DEMO COMPLETE\n');
    fprintf('=====================================================\n');
    fprintf('Completed: 4 simulated phantoms + 2 SWS methods.\n');
    fprintf('Metric stress tests: bias / map-domain noise / spatial blur.\n');
    fprintf('GT-displacement SWS diagnostic: exploratory, not final attribution.\n');
    fprintf('Not included yet: small_inclusions simulation or frequency sweep.\n');
    
end

% ========================================================================
% 1. Maps and residuals for two representative phantom configurations
% ========================================================================
function show_method_maps(batch_dir)
    cases = {'high_contrast','large_inclusions'};
    titles = {'High contrast','Large inclusions'};
    arrays = cell(2,1);
    for k = 1:2
        f = fullfile(batch_dir,cases{k},'method_comparison.mat');
        assert(isfile(f), 'Cannot find %s. Run compare_SWS_methods_batch first.', f);
        arrays{k} = load(f,'GT','PG','AC','x_axis','z_axis');
    end

    fig = figure('Name','THE | 1. SWS map comparison', ...
        'NumberTitle','off','Color','w');
    layout = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    for k = 1:2
        S = arrays{k};
        A = {double(S.GT),double(S.PG),double(S.AC)};
        subtitles = {'Ground truth','Phase gradient','Spatial autocorrelation'};
        for j = 1:3
            ax = nexttile(layout,(k-1)*3+j);
            imagesc(ax,double(S.x_axis(:))*1e3, ...
                double(S.z_axis(:))*1e3,A{j});
            axis(ax,'image'); set(ax,'YDir','normal');
            clim(ax,[0 8]); colorbar(ax);
            xlabel(ax,'Lateral [mm]'); ylabel(ax,'Depth [mm]');
            title(ax,sprintf('%s - %s',titles{k},subtitles{j}),'Interpreter','none');
        end
    end
    title(layout,'Same spatial grid and colour scale | 200 Hz');

    % Extra figure: a difference image makes under/overestimation clear.
    S = arrays{2};
    fig2 = figure('Name','THE | 2. Large-inclusion estimation errors', ...
        'NumberTitle','off','Color','w');
    L = tiledlayout(fig2,1,2,'TileSpacing','compact');
    errors = {double(S.PG)-double(S.GT),double(S.AC)-double(S.GT)};
    maxerr = max(abs([errors{1}(:);errors{2}(:)]),[],'omitnan');
    if ~isfinite(maxerr) || maxerr == 0, maxerr=1;end
    for j=1:2
        ax=nexttile(L);
        imagesc(ax,double(S.x_axis(:))*1e3,double(S.z_axis(:))*1e3,errors{j});
        axis(ax,'image');set(ax,'YDir','normal');
        clim(ax,[-maxerr maxerr]);colorbar(ax);
        xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
        if j==1, title(ax,'Phase gradient - ground truth');
        else, title(ax,'Spatial autocorrelation - ground truth');end
    end
    title(L,'Systematic SWS errors | large inclusions');
end

% ========================================================================
% 2. Four-phantom quantitative comparison, using existing batch_summary
% ========================================================================
function show_method_metrics(batch_dir)
    csv_path = fullfile(batch_dir,'batch_summary.csv');
    assert(isfile(csv_path), 'Missing %s. Run four-phantom batch comparison.',csv_path);
    T = readtable(csv_path);
    desired = { ...
        'Phantom','Method','GlobalRMSE_mps','InclusionRMSE_mps', ...
        'BackgroundBias_mps','CNR','DirectKDE_gCNR','ContrastRecovery'};
    assert(all(ismember(desired,T.Properties.VariableNames)), ...
        'batch_summary.csv has unexpected columns.');
    fprintf('FOUR-PHANTOM SWS RESULTS\n');
    disp(T(:,desired));

    phantoms = ["chaoran_baseline","low_contrast", ...
        "high_contrast","large_inclusions"];
    method_list = ["Phase gradient","Spatial autocorrelation"];
    field_names = {'InclusionRMSE_mps','CNR', ...
        'ContrastRecovery','BackgroundBias_mps'};
    field_titles = {'Inclusion RMSE [m/s] (lower is better)', ...
        'CNR (higher is better)', ...
        'Contrast recovery [%] (higher is better)', ...
        'Background bias [m/s] (closer to zero is better)'};
    fig = figure('Name','THE | 3. Quantitative results - four phantoms', ...
        'NumberTitle','off','Color','w');
    L = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
    for m=1:numel(field_names)
        vals = nan(numel(phantoms),2);
        for i=1:numel(phantoms)
            for j=1:2
                ix = find(string(T.Phantom)==phantoms(i) & ...
                    string(T.Method)==method_list(j),1);
                if ~isempty(ix)
                    vals(i,j)=double(T.(field_names{m})(ix));
                end
            end
        end
        if strcmp(field_names{m},'ContrastRecovery')
            vals=100*vals;
        end
        ax=nexttile(L);
        bar(ax,vals,'grouped');
        xticks(ax,1:4);
        xticklabels(ax,{'Baseline','Low contrast','High contrast','Large incl.'});
        title(ax,field_titles{m});ylabel(ax,'Metric value');
        grid(ax,'on');
        if m==1
            legend(ax,{'Phase gradient','Spatial autocorrelation'},'Location','best');
        end
    end
    title(L,'Same ground truth, ROI definitions and valid-pixel masks');
end

% ========================================================================
% 3. Stress tests from previously saved CSV outputs
% ========================================================================
function show_stress_tests(results_dir)
    root = fullfile(results_dir,'metric_stress_test');
    if ~isfolder(root)
        warning('No metric_stress_test folder; run your existing stress-test script first.');
        return;
    end
    files = {fullfile(root,'global_bias_stress_test.csv'), ...
             fullfile(root,'map_noise_stress_test.csv'), ...
             fullfile(root,'spatial_blur_stress_test.csv'), ...
             fullfile(root,'sharpness_noise_sensitivity.csv')};
    for k=1:numel(files)
        if ~isfile(files{k})
            warning('Stress-test output missing: %s',files{k});
        end
    end
    fig=figure('Name','THE | 4. Quality-metric stress tests', ...
        'NumberTitle','off','Color','w');
    L=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

    ax=nexttile(L);
    if isfile(files{1})
        T=readtable(files{1});
        plot(ax,T.AppliedBias_mps,T.GlobalRMSE_mps,'-o','LineWidth',1.4);
        xlabel(ax,'Added global SWS bias [m/s]');ylabel(ax,'RMSE [m/s]');
        title(ax,'Global bias stress test');grid(ax,'on');
    else, show_missing(ax,'Global bias CSV unavailable');end

    ax=nexttile(L);
    if isfile(files{2})
        T=readtable(files{2});
        errorbar(ax,T.NoiseSigma_mps,T.GlobalRMSE_mean_mps, ...
            T.GlobalRMSE_std_mps,'-o','LineWidth',1.3);
        xlabel(ax,'Added SWS map noise SD [m/s]');ylabel(ax,'RMSE [m/s]');
        title(ax,'Map-domain Gaussian noise (10 realizations)');grid(ax,'on');
    else, show_missing(ax,'Noise CSV unavailable');end

    ax=nexttile(L);
    if isfile(files{3})
        T=readtable(files{3});
        plot(ax,T.BlurSigma_mm,T.EdgeSharpness_mps_per_mm, ...
            '-o','LineWidth',1.4);hold(ax,'on');
        plot(ax,T.BlurSigma_mm,T.GTBoundarySharpness_mps_per_mm, ...
            '-s','LineWidth',1.4);hold(ax,'off');
        xlabel(ax,'Gaussian blur SD [mm]');ylabel(ax,'Sharpness [m/s/mm]');
        title(ax,'Blur: detected vs GT-referenced boundary');
        legend(ax,{'Detected edge','Fixed GT boundary'},'Location','best');
        grid(ax,'on');
    else, show_missing(ax,'Blur CSV unavailable');end

    ax=nexttile(L);
    if isfile(files{4})
        T=readtable(files{4});
        plot(ax,T.NoiseSigma_mps,T.DongSharpness_normalized, ...
            '-o','LineWidth',1.4);hold(ax,'on');
        plot(ax,T.NoiseSigma_mps,T.GTBoundarySharpness_normalized, ...
            '-s','LineWidth',1.4);hold(ax,'off');
        xlabel(ax,'Added SWS map noise SD [m/s]');
        ylabel(ax,'Normalised sharpness');
        title(ax,'Sharpness-metric noise sensitivity');
        legend(ax,{'Detected edge','Fixed GT boundary'},'Location','best');
        grid(ax,'on');
    else, show_missing(ax,'Sharpness-noise CSV unavailable');end

    title(L,'Metric validation: synthetic perturbations applied after SWS reconstruction');
end

function show_missing(ax,msg)
    axis(ax,'off');text(ax,0.05,0.5,msg,'Units','normalized');
end

% ========================================================================
% 4. Displacement-source diagnostic, read-only
% ========================================================================
function show_displacement_diagnostic(results_dir)
    root=fullfile(results_dir,'gt_displacement_diagnostic','large_inclusions');
    mat_path=fullfile(root,'gt_vs_estimated_maps.mat');
    if ~isfile(mat_path)
        warning('Missing %s; run compare_gt_displacement_sws first.',mat_path);
        return;
    end
    D=load(mat_path,'GT_SWS','PG_GT','PG_EST','x_est','z_est','summary_table');
    fprintf('GT-DISPLACEMENT DIAGNOSTIC (EXPLORATORY)\n');
    disp(D.summary_table);
    fig=figure('Name','THE | 5. GT vs estimated displacement diagnostic', ...
        'NumberTitle','off','Color','w');
    L=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
    maps={D.GT_SWS,D.PG_GT,D.PG_EST};
    labels={'SWS ground truth','From GT displacement','From estimated displacement'};
    for k=1:3
        ax=nexttile(L);
        imagesc(ax,1e3*double(D.x_est(:)),1e3*double(D.z_est(:)),double(maps{k}));
        axis(ax,'image');set(ax,'YDir','normal');clim(ax,[0 8]);colorbar(ax);
        title(ax,labels{k});xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
    end
    title(L,'Exploratory: temporal alignment and cropped-grid sensitivity require verification');
end

% ========================================================================
% OPTIONAL LIVE reconstruction from *already simulated* displacement data
% ========================================================================
function run_live_reconstruction(repo_dir,results_dir,batch_dir)
    phantom='high_contrast';
    data_name=[phantom '_displacement_with_quality.mat'];
    candidates={fullfile(repo_dir,'Data',phantom,data_name), ...
        fullfile(repo_dir,'Data',data_name), ...
        fullfile(repo_dir,'Analysis',data_name),fullfile(repo_dir,data_name)};
    data_file='';
    for k=1:numel(candidates)
        if isfile(candidates{k}),data_file=candidates{k};break;end
    end
    assert(~isempty(data_file),'Missing saved high-contrast displacement dataset.');

    toolbox_dir=fullfile(fileparts(repo_dir),'elastography','ElastographyToolbox');
    assert(isfolder(toolbox_dir),'Missing ElastographyToolbox.');
    addpath(genpath(toolbox_dir));
    assert(~isempty(which('postprocess.elastography.autocorrelation_SWS_method')), ...
        'USTB autocorrelation SWS method not available.');

    baseline_path=fullfile(results_dir,phantom,[phantom '_SWS_quality_results.mat']);
    assert(isfile(baseline_path),'Missing baseline SWS result.');
    B=load(baseline_path,'sws_gt','results');
    D=load(data_file,'b_data','displacement_estimation_autocorr','shear_params');
    scan=D.b_data.scan;
    Nx=numel(scan.x_axis);Nz=numel(scan.z_axis);
    raw=double(D.displacement_estimation_autocorr.data);
    assert(mod(numel(raw),Nx*Nz)==0,'Unexpected displacement dimensions.');
    Nt=numel(raw)/(Nx*Nz);
    u=-reshape(raw,Nz,Nx,Nt);
    clear raw D;

    PRF=double(B.results.PRF);
    freq=double(B.results.target_frequency);
    filter_speed=double(B.results.filter_speed);
    elasto=uff.elastography();
    elasto.PRF=PRF;elasto.target_frequency=freq;elasto.scan=scan;
    elasto.autocorrelation_field=uff.data_field();
    elasto.autocorrelation_field.set_component('u', ...
        reshape(u,Nz,Nx,1,Nt));
    num_angles=12;
    elasto.directional_filtering_angles= ...
        2*pi*linspace(0,1-1/num_angles,num_angles);
    ext=postprocess.elastography.extract_freq_component(elasto);
    ext.target_frequency=freq;ext.component='u';ext.run();
    omega=2*pi*freq;
    elasto.radial_filter=uff.spatial_filter('gaussian_charite','radial');
    elasto.radial_filter.sigma=filter_speed/omega;
    elasto.angular_filter=uff.spatial_filter('gaussian_charite','angular');
    elasto.angular_filter.deltaTheta= ...
        median(diff(elasto.directional_filtering_angles));
    filter_op=postprocess.elastography.apply_filter(elasto);
    filter_op.run();
    pg=postprocess.elastography.phase_gradient_method(elasto);pg.run();
    PG=squeeze(double(elasto.sws_map.u));
    ac=postprocess.elastography.autocorrelation_SWS_method(elasto);ac.run();
    AC=squeeze(double(elasto.sws_map.u));
    GT=double(B.sws_gt);
    assert(isequal(size(PG),size(AC),size(GT)), ...
        'SWS map sizes are inconsistent.');

    % Comparison to last saved batch (when available) proves reproducibility.
    batch_path=fullfile(batch_dir,phantom,'method_comparison.mat');
    if isfile(batch_path)
        Ref=load(batch_path,'PG','AC');
        common=isfinite(PG) & isfinite(Ref.PG);
        if any(common(:))
            fprintf('LIVE vs saved phase-gradient maximum difference: %.6g m/s\n', ...
                max(abs(PG(common)-Ref.PG(common))));
        end
    end
    output_dir=fullfile(results_dir,'meeting_demo');
    if ~isfolder(output_dir),mkdir(output_dir);end
    save(fullfile(output_dir,'live_high_contrast_maps.mat'), ...
        'PG','AC','GT','freq','PRF','filter_speed','-v7.3');
    fprintf('[LIVE] Two SWS reconstructions succeeded. Saved only to: %s\n', output_dir);

    fig=figure('Name','THE | LIVE: Both SWS estimators', ...
        'NumberTitle','off','Color','w');
    L=tiledlayout(fig,1,3,'TileSpacing','compact');
    maps={GT,PG,AC}; labels={'Ground truth','LIVE phase gradient', ...
        'LIVE spatial autocorrelation'};
    for k=1:3
        ax=nexttile(L);
        imagesc(ax,1e3*double(scan.x_axis(:)), ...
            1e3*double(scan.z_axis(:)),maps{k});
        axis(ax,'image');set(ax,'YDir','normal');clim(ax,[0 8]);colorbar(ax);
        xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');title(ax,labels{k});
    end
    title(L,sprintf('LIVE reconstruction from previously simulated displacement | %.0f Hz',freq));
end
