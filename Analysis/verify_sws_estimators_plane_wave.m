function summary = verify_sws_estimators_plane_wave(direction, grid_mode)
%VERIFY_SWS_ESTIMATORS_PLANE_WAVE Check PG and AC against a known plane wave.
%
% PURPOSE
%   Isolate the SWS reconstruction stage from k-Wave, Field II, beamforming,
%   displacement tracking, and inclusions. This is a DIAGNOSTIC, not a new
%   quality metric and not a validation of the physical phantom simulation.
%
% RUN (put this .m file in the repo's Analysis/ directory):
%   T = verify_sws_estimators_plane_wave();                 % right, full grid
%   T = verify_sws_estimators_plane_wave('both','full');    % two directions
%   T = verify_sws_estimators_plane_wave('right','quick');  % smaller preview
%
% INPUTS
%   direction : 'right' (default), 'left', or 'both'.
%   grid_mode : 'full' (default) uses EXACT saved ultrasound x/z axes;
%               'quick' downsamples those axes and is for debugging only.
%
% Uses metadata ONLY from:
%   Results/method_comparison/batch/large_inclusions/method_comparison.mat
% (or high_contrast if large_inclusions is unavailable).
%
% Uses precisely the USTB operations/settings of compare_SWS_methods_batch:
%   uff.elastography -> extract_freq_component -> apply_filter ->
%   phase_gradient_method -> autocorrelation_SWS_method.
%
% SYNTHETIC INPUT: time increments from an ideal plane wave
%   d(x,t) = A*[cos(w*(t+dt) - s*k*x) - cos(w*t - s*k*x)]
%   k = w/c,  c = 2.32 m/s, f = 200 Hz, PRF = 3000 Hz.
% Both phase directions have the same magnitude of SWS.
%
% No data files or existing results are overwritten. Output goes under:
%   Results/sws_estimator_diagnostic/analytic_plane_wave/<grid_mode>/
%
% NOTE: USTB is not available in this environment, so the *USTB* part must
% be run on your MATLAB machine. Any numerical thresholds printed here are
% screening criteria, not published acceptance limits.

    if nargin < 1 || isempty(direction), direction = 'right'; end
    if nargin < 2 || isempty(grid_mode), grid_mode = 'full'; end
    direction = lower(char(string(direction)));
    grid_mode = lower(char(string(grid_mode)));
    assert(any(strcmp(direction, {'right','left','both'})), ...
        'direction must be right, left, or both.');
    assert(any(strcmp(grid_mode, {'full','quick'})), ...
        'grid_mode must be full or quick.');

    analysis_dir = fileparts(mfilename('fullpath'));
    repo_dir = fileparts(analysis_dir);
    master_dir = fileparts(repo_dir);

    toolbox_dir = fullfile(master_dir, 'elastography', 'ElastographyToolbox');
    assert(isfolder(toolbox_dir), 'Missing ElastographyToolbox: %s', toolbox_dir);
    addpath(genpath(toolbox_dir));

    classes = {'uff.elastography', 'uff.linear_scan', 'uff.data_field', ...
        'postprocess.elastography.extract_freq_component', ...
        'postprocess.elastography.apply_filter', ...
        'postprocess.elastography.phase_gradient_method', ...
        'postprocess.elastography.autocorrelation_SWS_method'};
    for i = 1:numel(classes)
        assert(~isempty(which(classes{i})), ...
            'Missing USTB class/function %s on MATLAB path.', classes{i});
    end

    % Read the EXACT source grid and acquisition settings used in the
    % previously executed estimator comparison, without loading big maps.
    result_dir = fullfile(repo_dir, 'Results', 'method_comparison', 'batch');
    source_file = fullfile(result_dir, 'large_inclusions', 'method_comparison.mat');
    if ~isfile(source_file)
        source_file = fullfile(result_dir, 'high_contrast', 'method_comparison.mat');
    end
    assert(isfile(source_file), ...
        'Could not find method_comparison.mat in Results/method_comparison/batch.');
    M = load(source_file, 'x_axis', 'z_axis', 'PRF', ...
        'target_frequency', 'filter_speed');
    required = {'x_axis','z_axis','PRF','target_frequency','filter_speed'};
    assert(all(isfield(M, required)), 'Saved comparison lacks scan/settings metadata.');
    x_full = double(M.x_axis(:));
    z_full = double(M.z_axis(:));
    assert(all(diff(x_full)>0) && all(diff(z_full)>0), ...
        'Scan coordinates must be increasing.');

    if strcmp(grid_mode, 'quick')
        % Non-production spacing: quick sanity check only.
        x = x_full(1:2:end);
        z = z_full(1:4:end);
    else
        x = x_full;
        z = z_full;
    end
    Nx = numel(x); Nz = numel(z);
    assert(Nx>20 && Nz>20, 'Too few scan points for a meaningful plane wave.');

    cfg = struct();
    cfg.c_true_mps = 2.32;
    cfg.PRF_Hz = double(M.PRF);
    cfg.f0_Hz = double(M.target_frequency);
    cfg.filter_speed_mps = double(M.filter_speed);
    cfg.Nt = 59; % Same as previously analyzed displacement datasets.
    cfg.absolute_amplitude_m = 20e-6;
    cfg.num_angles = 12;
    cfg.grid_mode = grid_mode;
    cfg.source_metadata = source_file;
    assert(all(isfinite([cfg.PRF_Hz,cfg.f0_Hz,cfg.filter_speed_mps])) ...
        && cfg.PRF_Hz>2*cfg.f0_Hz && cfg.f0_Hz>0 ...
        && cfg.filter_speed_mps>0, 'Invalid sampling/filter settings.');
    assert(cfg.Nt>=8, 'Too few time samples.');

    if abs(cfg.f0_Hz - 200)>1e-8 || abs(cfg.PRF_Hz - 3000)>1e-8
        warning(['Source settings are PRF=%.3f, f0=%.3f rather than ' ...
            'previously observed PRF=3000, f0=200. ' ...
            'Using the SAVED metadata, not hardcoded values.'], ...
            cfg.PRF_Hz, cfg.f0_Hz);
    end

    omega = 2*pi*cfg.f0_Hz;
    k = omega/cfg.c_true_mps;
    dx = median(diff(x));
    dz = median(diff(z));
    assert(abs(k*dx)<pi, ...
        'Lateral phase sampling aliases: k*dx=%.2f rad >= pi.', abs(k*dx));
    scan = uff.linear_scan('x_axis', x, 'z_axis', z);
    out_dir = fullfile(repo_dir, 'Results', ...
        'sws_estimator_diagnostic', 'analytic_plane_wave', grid_mode);
    if ~isfolder(out_dir), mkdir(out_dir); end

    if strcmp(direction,'both')
        directions = {'right','left'};
    else
        directions = {direction};
    end

    summary = table();
    fprintf('\n=== ANALYTIC PLANE-WAVE SWS CHECK ===\n');
    fprintf('Source: %s\n', source_file);
    fprintf('Grid: %d axial x %d lateral; dx=%.5f mm, dz=%.5f mm\n', ...
        Nz,Nx,dx*1e3,dz*1e3);
    fprintf('f0=%.1f Hz; PRF=%.1f Hz; c=%.4f m/s; wavelength=%.3f mm\n', ...
        cfg.f0_Hz,cfg.PRF_Hz,cfg.c_true_mps,1e3*cfg.c_true_mps/cfg.f0_Hz);
    fprintf('Temporal frames: %d; estimated raw input memory: %.2f GB (single)\n', ...
        cfg.Nt, single(Nz)*single(Nx)*single(cfg.Nt)*4/1e9);
    if strcmp(grid_mode,'quick')
        fprintf('NOTE: QUICK grid has different spacing/boundaries; use FULL for conclusions.\n');
    end

    % Fixed interior ROI, independent of reconstructed values. Helps
    % separate finite-window/filter boundary effects from interior behavior.
    margin_x = max(2, round(0.15*Nx));
    margin_z = max(2, round(0.15*Nz));
    interior = false(Nz,Nx);
    interior(1+margin_z:Nz-margin_z,1+margin_x:Nx-margin_x) = true;

    for j = 1:numel(directions)
        label = directions{j};
        sgn = 1;
        if strcmp(label,'left'), sgn = -1; end
        fprintf('\n--- Plane wave propagating %s (+/- x) ---\n', label);

        % Generate displacement increments on the ACTUAL USTB scan grid.
        % Identical across depth; only lateral phase varies.
        dt = 1/cfg.PRF_Hz;
        xx = reshape(x,1,Nx);
        u = zeros(Nz,Nx,cfg.Nt,'single');
        for n = 1:cfg.Nt
            t = (n-1)*dt;
            sample = cfg.absolute_amplitude_m * ( ...
                cos(omega*(t+dt)-sgn*k*xx) - ...
                cos(omega*t-sgn*k*xx));
            u(:,:,n) = repmat(single(sample),Nz,1);
        end

        % INDEPENDENT numerical sanity test: recover complex phasor from
        % the sampled time series by least-squares fitting cos(wt),sin(wt).
        % This does not invoke any USTB SWS estimator.
        t = ((0:cfg.Nt-1).'/cfg.PRF_Hz);
        basis = [cos(omega*t), sin(omega*t)];
        iz = round(Nz/2);
        lateral_series = double(squeeze(u(iz,:,:))).'; % [Nt, Nx]
        coef = basis \ lateral_series;
        phasor = coef(1,:) - 1i*coef(2,:);
        dphi = angle(conj(phasor(1:end-1)).*phasor(2:end));
        local_dx = diff(x).';
        c_oracle = omega ./ abs(dphi ./ local_dx);
        oracle_error = abs(median(c_oracle(isfinite(c_oracle))) - cfg.c_true_mps);
        fprintf('Independent temporal-fit/phase-difference oracle: %.6f m/s\n', ...
            median(c_oracle(isfinite(c_oracle))));
        assert(oracle_error < 0.005*cfg.c_true_mps, ...
            'Synthetic signal or numerical phase sampling is inconsistent.');
        clear lateral_series coef phasor dphi c_oracle;

        fprintf('Running USTB extract_freq_component + spatial filters + PG/AC...\n');
        [PG,AC] = run_original_sws_pipeline(u,scan,cfg);
        clear u;
        assert(isequal(size(PG),[Nz,Nx]) && isequal(size(AC),[Nz,Nx]), ...
            'SWS output dimensions differ from scan.');

        % Do not hide invalid pixels: coverage is reported explicitly.
        values = {PG,AC};
        names = {'PG','AC'};
        for m = 1:2
            map = values{m};
            good = interior & isfinite(map) & map>0 & map<10;
            N = nnz(good);
            coverage = N/nnz(interior);
            if N > 0
                e = map(good)-cfg.c_true_mps;
                med = median(map(good));
                bias = mean(e);
                rmse = sqrt(mean(e.^2));
                mae = mean(abs(e));
            else
                med=NaN; bias=NaN; rmse=NaN; mae=NaN;
            end
            % Informal initial screening only, NOT a pass/fail guarantee.
            screening_ok = coverage>=0.8 && ...
                abs(med-cfg.c_true_mps)/cfg.c_true_mps<0.05 && ...
                rmse/cfg.c_true_mps<0.10;
            r = table(string(label),string(grid_mode),string(names{m}), ...
                cfg.c_true_mps, med, bias, rmse, mae, coverage, ...
                100*rmse/cfg.c_true_mps, screening_ok, ...
                'VariableNames',{'Direction','Grid','Estimator', ...
                'TrueSWS_mps','MedianSWS_mps','Bias_mps','RMSE_mps', ...
                'MAE_mps','ValidFraction','RelativeRMSE_percent','ScreeningOK'});
            summary = [summary; r]; %#ok<AGROW>
            fprintf('%s: valid=%.1f%%, median=%.4f, bias=%+.4f, RMSE=%.4f m/s', ...
                names{m},100*coverage,med,bias,rmse);
            if screening_ok
                fprintf(' [screening OK]\n');
            else
                fprintf(' [INVESTIGATE]\n');
            end
        end

        % Figure: constant known SWS, two reconstructions and their errors.
        GT = cfg.c_true_mps * ones(Nz,Nx);
        fig = figure('Visible','off','Color','w','Position',[50 50 1470 800]);
        tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
        maps = {GT,PG,AC};
        titles = {'Known SWS (analytic)','PG estimated SWS','AC estimated SWS'};
        for m=1:3
            nexttile(tl,m);
            imagesc(x*1e3,z*1e3,maps{m});
            axis image; set(gca,'YDir','normal');
            clim([0,8]); colorbar;
            xlabel('Lateral [mm]'); ylabel('Depth [mm]'); title(titles{m});
        end
        nexttile(tl,4);
        imagesc(x*1e3,z*1e3,double(interior));
        axis image; set(gca,'YDir','normal');
        clim([0,1]); colorbar;
        xlabel('Lateral [mm]'); ylabel('Depth [mm]');
        title('Fixed interior evaluation ROI');
        for m=1:2
            nexttile(tl,4+m);
            imagesc(x*1e3,z*1e3,values{m}-GT);
            axis image; set(gca,'YDir','normal');
            clim([-3,3]); colorbar;
            xlabel('Lateral [mm]'); ylabel('Depth [mm]');
            title([names{m} ' error [m/s]']);
        end
        title(tl,sprintf('Analytic plane wave (%s) | true SWS %.2f m/s | %s grid', ...
            label,cfg.c_true_mps,grid_mode));
        exportgraphics(fig,fullfile(out_dir, ...
            ['plane_wave_' label '_maps.png']),'Resolution',180);
        close(fig);

        % Lateral profile at the center depth; watch for edge artifacts.
        fig=figure('Visible','off','Color','w');
        plot(x*1e3,PG(iz,:),'LineWidth',1.4); hold on;
        plot(x*1e3,AC(iz,:),'LineWidth',1.4);
        yline(cfg.c_true_mps,'k--','LineWidth',1.4);
        xline(x(1+margin_x)*1e3,':');
        xline(x(Nx-margin_x)*1e3,':');
        xlabel('Lateral [mm]'); ylabel('SWS [m/s]');
        title(sprintf('Central-depth SWS profile: %s',label));
        legend('PG','AC','Analytic truth','Location','best');
        grid on; ylim([0,8]);
        exportgraphics(fig,fullfile(out_dir, ...
            ['plane_wave_' label '_profile.png']),'Resolution',180);
        close(fig);

        save(fullfile(out_dir,['plane_wave_' label '_results.mat']), ...
            'GT','PG','AC','x','z','interior','cfg','oracle_error','-v7.3');
        clear PG AC GT values;
    end

    disp(summary);
    writetable(summary,fullfile(out_dir,'plane_wave_summary.csv'));
    fprintf('\nSaved diagnostics: %s\n',out_dir);
    fprintf(['Interpretation: if the numerical oracle succeeds but PG/AC ' ...
        'deviate inside the ROI, investigate extraction/filter/estimator ' ...
        'settings. If both succeed, move next to the homogeneous ' ...
        'k-Wave case; this does NOT yet prove inclusion simulations correct.\n']);
end

function [PG,AC] = run_original_sws_pipeline(u,scan,cfg)
    % Duplicates the operative sequence in compare_SWS_methods_batch.m.
    elasto = uff.elastography();
    elasto.PRF = cfg.PRF_Hz;
    elasto.target_frequency = cfg.f0_Hz;
    elasto.scan = scan;
    elasto.autocorrelation_field = uff.data_field();
    elasto.autocorrelation_field.set_component('u', ...
        reshape(u,size(u,1),size(u,2),1,size(u,3)));

    elasto.directional_filtering_angles = ...
        2*pi*linspace(0,1-1/cfg.num_angles,cfg.num_angles);
    extractor = postprocess.elastography.extract_freq_component(elasto);
    extractor.target_frequency = cfg.f0_Hz;
    extractor.component = 'u';
    extractor.run();

    omega = 2*pi*cfg.f0_Hz;
    elasto.radial_filter = uff.spatial_filter('gaussian_charite','radial');
    elasto.radial_filter.sigma = cfg.filter_speed_mps/omega;
    elasto.angular_filter = uff.spatial_filter('gaussian_charite','angular');
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
        'AC did not produce a SWS map.');
    AC = squeeze(double(elasto.sws_map.u));
end
