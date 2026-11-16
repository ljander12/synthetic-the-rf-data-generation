function T = diagnose_high_contrast_phase_filter()
%DIAGNOSE_HIGH_CONTRAST_PHASE_FILTER  Is high-contrast phase information lost by filtering?
%
% Run from synthetic-the-rf-data-generation/Analysis/:
%   T = diagnose_high_contrast_phase_filter();
%
% ONE focused diagnostic, no new k-Wave/Field II simulations:
%   1. Native k-Wave displacement -> direct 200-Hz complex phasor.
%   2. SAME displacement, packet average and cropped US scan as in
%      diagnose_high_contrast_gt_vs_us.m -> USTB 200-Hz extraction.
%   3. Inspect USTB's ACTUAL directional-filter output for ALL angles.
%   4. Compare local phase-derived apparent SWS before/after filtering,
%      inside and outside the known inclusion, with saved PG/AC SWS.
%
% IMPORTANT: Phase-derived apparent speed is DIAGNOSTIC, not a new SWS
% estimator or quality metric. Single-frequency local phase gradients do
% not necessarily equal true local material SWS under interference.
% A 2-D slice cannot characterize out-of-plane wave propagation.
% k-Wave GT is native 1-mm grid; dense interpolation restores NO resolution.
%
% Requires prior result:
%   Results/method_comparison/batch/high_contrast/method_comparison.mat
% Optional previous controlled comparison:
%   Results/sws_estimator_diagnostic/gt_vs_us_high_contrast/comparison_results.mat
%
% Output ONLY to:
%   Results/sws_estimator_diagnostic/high_contrast_phase_filter/
%   phase_filter_summary.csv, phase_filter_angles.csv,
%   phase_filter_maps.png, phase_filter_profiles.png,
%   phase_filter_results.mat
%
% Written 2026-11-11; MATLAB/USTB runtime NOT available for execution here.

    here = fileparts(mfilename('fullpath'));
    if isfolder(fullfile(here,'Data')) && isfolder(fullfile(here,'Results'))
        repo = here;
    else
        repo = fileparts(here);
    end
    assert(isfolder(fullfile(repo,'Data')) && isfolder(fullfile(repo,'Results')), ...
        'Place diagnose_high_contrast_phase_filter.m in repository/Analysis/.');
    addpath(genpath(repo));
    toolbox_dir = fullfile(fileparts(repo),'elastography','ElastographyToolbox');
    assert(isfolder(toolbox_dir),'Cannot find ElastographyToolbox: %s', toolbox_dir);
    addpath(genpath(toolbox_dir));
    assert(~isempty(which('create_phantom')), 'create_phantom is missing from MATLAB path.');
    assert(~isempty(which('uff.elastography')), 'uff.elastography is missing from MATLAB path.');
    assert(~isempty(which('postprocess.elastography.apply_filter')), ...
        'USTB apply_filter is missing from MATLAB path.');

    name = 'high_contrast_displacement_with_quality.mat';
    paths = {fullfile(repo,'Data','high_contrast',name), ...
        fullfile(repo,'Data',name),fullfile(repo,name),fullfile(here,name)};
    datafile = '';
    for i=1:numel(paths)
        if isfile(paths{i}), datafile=paths{i}; break; end
    end
    assert(~isempty(datafile),'Cannot find original %s.',name);
    mfile = fullfile(repo,'Results','method_comparison','batch', ...
        'high_contrast','method_comparison.mat');
    assert(isfile(mfile),'Missing saved PG/AC method comparison: %s',mfile);
    variables = whos('-file',datafile);
    varnames = {variables.name};
    must = {'disp_gt','b_data','shear_params','displacement_estimation_autocorr'};
    assert(all(ismember(must,varnames)), ...
        'Original MAT lacks a saved k-Wave displacement or required metadata.');
    opts = {'disp_gt','b_data','shear_params','autocorr_params', ...
        'array_length','scene_depth','prf','prb_theta','displacement_estimation_autocorr'};
    opts = opts(ismember(opts,varnames));
    D = load(datafile,opts{:});
    M = load(mfile,'GT','PG','AC','roi','x_axis','z_axis', ...
        'PRF','target_frequency','filter_speed');
    req = {'GT','PG','AC','roi','x_axis','z_axis', ...
        'PRF','target_frequency','filter_speed'};
    assert(all(isfield(M,req)),'Saved method comparison has missing fields.');

    x0 = double(M.x_axis(:)); z0 = double(M.z_axis(:));
    [Nz,Nx] = size(M.GT);
    assert(isequal(size(M.PG),[Nz Nx]) && isequal(size(M.AC),[Nz Nx]) && ...
        Nx==numel(x0) && Nz==numel(z0),'Saved SWS axes/maps differ.');
    assert(isequal(size(D.b_data.scan.x_axis(:)),size(x0)) && ...
        max(abs(double(D.b_data.scan.x_axis(:))-x0))<1e-10 && ...
        max(abs(double(D.b_data.scan.z_axis(:))-z0))<1e-10, ...
        'Displacement US scan and saved SWS scan differ.');
    assert(all(diff(x0)>0) && all(diff(z0)>0),'Axes must increase.');
    assert(all(isfield(M.roi,{'image_mask','inclusion_core','background_core'})), ...
        'ROI masks missing from saved comparison.');

    f0 = double(M.target_frequency);
    PRF = double(M.PRF);
    vfilter = double(M.filter_speed);
    assert(f0>0 && PRF>2*f0 && vfilter>0 && ...
        abs(double(D.shear_params.source_freq)-f0)<1e-6, ...
        'Target frequency/PRF/filter metadata inconsistent.');
    if isfield(D,'prf') && ~isempty(D.prf)
        assert(abs(double(D.prf)-PRF)<1e-6,'Saved PRFs disagree.');
    end
    c_bg = double(D.shear_params.c_shear_bkg);
    c_inc = double(D.shear_params.c_shear_incl);
    assert(c_inc>c_bg,'Expected positive high-contrast inclusion.');
    if isfield(D,'prb_theta')
        fprintf('Probe orientation: %s degrees\n', ...
            mat2str(rad2deg(double(D.prb_theta))));
    end
    if isfield(D,'autocorr_params') && isfield(D.autocorr_params,'packet_size')
        packet = double(D.autocorr_params.packet_size);
    else
        packet = 4;
        warning('Using packet=4 from previous GT-vs-US diagnostic.');
    end
    assert(isscalar(packet) && packet>=2 && fix(packet)==packet, ...
        'Invalid packet size.');
    n_avg = packet-1; % Match EXACT previous GT-vs-US diagnostic convention.
    [ngz,ngx,ngt] = size(D.disp_gt);
    n_est = numel(D.displacement_estimation_autocorr.data)/(Nz*Nx);
    assert(n_est==fix(n_est) && ngt>=n_est+n_avg-1, ...
        'GT time frames cannot be matched to the saved US-displacement frames.');

    % Recover native k-Wave spatial coordinates from original setup.
    [kg,~] = create_phantom(D.shear_params); % NO simulation
    dx = double(kg.dx); dz = double(kg.dz);
    clear kg;
    if isfield(D,'array_length')
        array_length = double(D.array_length);
    else
        array_length = (128-1)*0.300e-3;
        warning('array_length missing: using documented 128-element 0.3 mm probe.');
    end
    if isfield(D,'scene_depth')
        scene_depth = double(D.scene_depth);
    else
        scene_depth = 50e-3;
        warning('scene_depth missing: using documented 50 mm.');
    end
    xn = (-array_length/2:dx:array_length/2).';
    zn = (5e-3:dz:scene_depth).';
    assert(numel(xn)==ngx && numel(zn)==ngz, ...
        'k-Wave native-grid axes do not match saved disp_gt.');

    keep_x = find(x0>=xn(1)-1e-9 & x0<=xn(end)+1e-9);
    keep_z = find(z0>=zn(1)-1e-9 & z0<=zn(end)+1e-9);
    assert(numel(keep_x)>25 && numel(keep_z)>25,'Insufficient GT/US grid overlap.');
    x = x0(keep_x); z = z0(keep_z);
    [Xq,Zq] = meshgrid(min(max(x,xn(1)),xn(end)), ...
        min(max(z,zn(1)),zn(end)));
    ctrue = double(M.GT(keep_z,keep_x));
    gt_inc = ctrue>(c_bg+c_inc)/2;
    inc = logical(M.roi.inclusion_core(keep_z,keep_x));
    bg = logical(M.roi.background_core(keep_z,keep_x));
    valid_roi = logical(M.roi.image_mask(keep_z,keep_x));
    assert(nnz(inc)>100 && nnz(bg)>100,'Too few GT inclusion/background pixels.');
    % Fixed diagnostic masks; do not optimize after seeing results.
    inc = inc & valid_roi;
    bg = bg & valid_roi;

    fprintf('\n=== HIGH CONTRAST: RAW vs FILTERED k-Wave PHASE ===\n');
    fprintf('Source: %s\n',datafile);
    fprintf('Material SWS: background %.3f, inclusion %.3f m/s\n',c_bg,c_inc);
    fprintf('Native GT %dx%dx%d, dx=%.3f mm dz=%.3f mm\n', ...
        ngz,ngx,ngt,1e3*dx,1e3*dz);
    fprintf('Common scan %dx%d; GT packet=%d frames, used=%d\n', ...
        numel(z),numel(x),n_avg,n_est);
    fprintf('f0=%.1f Hz PRF=%.1f Hz filter_speed=%.3f m/s\n',f0,PRF,vfilter);
    fprintf('ROI: inclusion %d, background %d pixels\n',nnz(inc),nnz(bg));

    % Form SAME k-Wave GT input as previous diagnostic: mean of packet-1
    % successive PRF-sampled displacement increments, then spatial interpolation.
    native_u = zeros(ngz,ngx,n_est,'single');
    dense_u = zeros(numel(z),numel(x),n_est,'single');
    for t=1:n_est
        frame = single(mean(double(D.disp_gt(:,:,t:t+n_avg-1)),3));
        native_u(:,:,t) = frame;
        interp_frame = interp2(xn,zn,double(frame),Xq,Zq,'linear',NaN);
        assert(all(isfinite(interp_frame(:))), ...
            'GT interpolation failed for time frame %d.',t);
        dense_u(:,:,t)=single(interp_frame);
    end
    clear D Xq Zq interp_frame frame;

    % Independent direct harmonic fit on native 1-mm grid, a control for
    % interpolation-related changes (not a substitute for the USTB extractor).
    U_native = harmonic_phasor(native_u,PRF,f0);
    clear native_u;
    [knx,knz] = phase_gradient(U_native,xn,zn,1e-3);
    c_native = phase_speed(knx,f0);
    c2_native = phase_speed(hypot(knx,knz),f0);

    % Reproduce the SWS pipeline and read EXACT USTB intermediate fields.
    elasto = uff.elastography();
    elasto.PRF = PRF;
    elasto.target_frequency = f0;
    elasto.scan = uff.linear_scan('x_axis',x,'z_axis',z);
    elasto.autocorrelation_field = uff.data_field();
    elasto.autocorrelation_field.set_component('u', ...
        reshape(dense_u,numel(z),numel(x),1,n_est));
    clear dense_u;
    angles = 2*pi*linspace(0,1-1/12,12);
    elasto.directional_filtering_angles = angles;
    extractor = postprocess.elastography.extract_freq_component(elasto);
    extractor.target_frequency = f0;
    extractor.component = 'u';
    extractor.run();
    assert(isprop(elasto,'unfiltered_phase_image_xz'), ...
        'Elastography object lacks unfiltered_phase_image_xz property.');
    U_raw = fetch_phasor(elasto.unfiltered_phase_image_xz, ...
        numel(z),numel(x),'unfiltered_phase_image_xz');
    assert(size(U_raw,3)==1,'Expected one unfiltered complex phase field.');
    U_raw = U_raw(:,:,1);

    omega = 2*pi*f0;
    elasto.radial_filter = uff.spatial_filter('gaussian_charite','radial');
    elasto.radial_filter.sigma = vfilter/omega;
    elasto.angular_filter = uff.spatial_filter('gaussian_charite','angular');
    elasto.angular_filter.deltaTheta = median(diff(angles));
    filter_op = postprocess.elastography.apply_filter(elasto);
    filter_op.run();
    assert(isprop(elasto,'directional_filtered_field_xz'), ...
        'Elastography object lacks directional_filtered_field_xz property.');
    U_filt = fetch_phasor(elasto.directional_filtered_field_xz, ...
        numel(z),numel(x),'directional_filtered_field_xz');
    nang = size(U_filt,3);
    assert(nang>=1,'USTB filtering produced no direction channels.');
    if nang~=numel(angles)
        warning('Filtered field has %d channels, configured for %d directions.', ...
            nang,numel(angles));
    end

    % Validate that the re-created SWS pipeline matches saved/cropped GT.
    pg = postprocess.elastography.phase_gradient_method(elasto);
    pg.run();
    PG = squeeze(double(elasto.sws_map.u));
    ac = postprocess.elastography.autocorrelation_SWS_method(elasto);
    ac.run();
    AC = squeeze(double(elasto.sws_map.u));
    assert(isequal(size(PG),size(ctrue),size(AC)), ...
        'Reconstructed PG/AC map dimensions differ from comparison scan.');

    % Fixed 1-mm half-width in phase derivative: equivalent to a 2-mm
    % centered spatial difference on both native and dense grids. This
    % reduces tiny-dx numerical errors in interpolated data.
    [krx,krz] = phase_gradient(U_raw,x,z,1e-3);
    c_raw = phase_speed(krx,f0);
    c2_raw = phase_speed(hypot(krx,krz),f0);

    scores = nan(nang,1);
    cmaps = cell(nang,1);
    c2maps = cell(nang,1);
    field_amp = nan(nang,1);
    for j=1:nang
        A = U_filt(:,:,j);
        m = abs(A(bg));
        m = m(isfinite(m));
        if ~isempty(m), scores(j)=mean(m.^2);end
        field_amp(j)=mean(abs(A(inc)),'omitnan');
        [kx,kz] = phase_gradient(A,x,z,1e-3);
        cmaps{j} = phase_speed(kx,f0);
        c2maps{j} = phase_speed(hypot(kx,kz),f0);
    end
    [~, chosen] = max(scores);
    assert(isfinite(scores(chosen)),'No finite filtered-field amplitude in background ROI.');
    c_filt = cmaps{chosen};
    c2_filt = c2maps{chosen};
    fprintf('USTB directional channels: %d; chosen channel %d by BG energy.\n', ...
        nang,chosen);
    if nang==numel(angles)
        fprintf('Selected filter angle: %.1f degrees.\n',rad2deg(angles(chosen)));
    end

    % Compare on the same valid inclusion/background pixels. A local
    % phase-speed estimate can become unbounded if phase slope is ~0;
    % keep coverage as a separate, explicit statistic.
    common = inc | bg;
    same = common & isfinite(c_raw) & isfinite(c_filt) & ...
        c_raw>0 & c_raw<12 & c_filt>0 & c_filt<12;
    native_inc = interp2(x,z,double(inc), ...
        repmat(xn.',numel(zn),1),repmat(zn,1,numel(xn)), ...
        'nearest',0)>0.5;
    native_bg = interp2(x,z,double(bg), ...
        repmat(xn.',numel(zn),1),repmat(zn,1,numel(xn)), ...
        'nearest',0)>0.5;
    native_inc = native_inc & isfinite(c_native);
    native_bg = native_bg & isfinite(c_native);

    stages = ["Native k-Wave raw";"USTB before filter"; ...
              "USTB filtered / BG-dominant";"USTB PG";"USTB AC"];
    vals = {c_native,c_raw,c_filt,PG,AC};
    maski = {native_inc,inc,inc,inc,inc};
    maskb = {native_bg,bg,bg,bg,bg};
    nrows=numel(vals);
    medI = nan(nrows,1); medB = medI; meanI=medI; meanB=medI;
    fracI = medI; fracB=medI;
    rec = medI;
    for j=1:nrows
        [medI(j),meanI(j),fracI(j)] = stage_stat(vals{j},maski{j});
        [medB(j),meanB(j),fracB(j)] = stage_stat(vals{j},maskb{j});
        rec(j)=(meanI(j)-meanB(j))/(c_inc-c_bg);
    end
    T=table(stages,medI,medB,meanI,meanB,rec,fracI,fracB, ...
        'VariableNames',{'Stage','InclMedian_mps','BgMedian_mps', ...
            'InclMean_mps','BgMean_mps','ContrastRecovery', ...
            'InclValidFraction','BgValidFraction'});

    id = (1:nang).';
    angledeg = nan(nang,1);
    if nang==numel(angles),angledeg=rad2deg(angles(:));end
    powerBG = scores;
    inclMean = nan(nang,1); bgMean = inclMean;
    inclMedian = inclMean; bgMedian = inclMean;
    fracIncl = inclMean;
    for j=1:nang
        [inclMedian(j),inclMean(j),fracIncl(j)] = stage_stat(cmaps{j},inc);
        [bgMedian(j),bgMean(j)] = stage_stat(cmaps{j},bg);
    end
    contrast = (inclMean-bgMean)/(c_inc-c_bg);
    Aangles=table(id,angledeg,powerBG,inclMedian,bgMedian, ...
        inclMean,bgMean,contrast,fracIncl, ...
        'VariableNames',{'Channel','Angle_deg','BackgroundEnergy', ...
          'InclMedian_mps','BgMedian_mps','InclMean_mps', ...
          'BgMean_mps','ContrastRecovery','InclValidFraction'});
    fprintf('\n--- RAW vs FILTERED PHASE SPEED (lateral) ---\n');
    disp(T);
    fprintf('\n--- ALL FILTER DIRECTIONS ---\n');
    disp(Aangles);
    fprintf('True inclusion speed %.3f m/s, BG %.3f m/s.\n',c_inc,c_bg);
    fprintf('Raw-filter common valid ROI: %d/%d pixels.\n', ...
        nnz(same),nnz(common));

    % Cross-check with the previous GT-vs-US controlled run (same cropped
    % ultrasound scan and packet convention). A discrepancy signals a
    % pipeline mismatch before any physical interpretation.
    prior = fullfile(repo,'Results','sws_estimator_diagnostic', ...
        'gt_vs_us_high_contrast','comparison_results.mat');
    if isfile(prior)
        P=load(prior,'pg_gt','ac_gt','x','z');
        if isfield(P,'pg_gt') && isfield(P,'ac_gt') && ...
                isequal(size(P.pg_gt),size(PG)) && ...
                max(abs(double(P.x(:))-x))<1e-10 && ...
                max(abs(double(P.z(:))-z))<1e-10
            ix=inc & isfinite(P.pg_gt) & isfinite(PG);
            ia=inc & isfinite(P.ac_gt) & isfinite(AC);
            dpg=max(abs(PG(ix)-double(P.pg_gt(ix))));
            dac=max(abs(AC(ia)-double(P.ac_gt(ia))));
            fprintf('Reproduced GT-input PG max difference: %.5g m/s\n',dpg);
            fprintf('Reproduced GT-input AC max difference: %.5g m/s\n',dac);
            if dpg>0.01 || dac>0.01
                warning('Current filter/SWS setup differs from previous GT-input maps.');
            end
        else
            warning('Previous GT-vs-US scan differs; skipped numerical cross-check.');
        end
    end

    outdir=fullfile(repo,'Results','sws_estimator_diagnostic', ...
        'high_contrast_phase_filter');
    if ~isfolder(outdir),mkdir(outdir);end
    writetable(T,fullfile(outdir,'phase_filter_summary.csv'));
    writetable(Aangles,fullfile(outdir,'phase_filter_angles.csv'));

    % Maps in physical units; all speeds share [0,8] m/s, no auto-contrast.
    fig=figure('Color','w','Position',[50 70 1850 870],'Visible','off');
    tl=tiledlayout(fig,2,5,'TileSpacing','compact','Padding','compact');
    plot_sws(nexttile(tl,1),x,z,ctrue,gt_inc, ...
        'True material SWS',[0 8]);
    plot_sws(nexttile(tl,2),xn,zn,c_native,native_inc, ...
        'Native raw phase c_x (1-mm grid)',[0 8]);
    plot_sws(nexttile(tl,3),x,z,c_raw,gt_inc, ...
        'Raw 200-Hz phase c_x',[0 8]);
    plot_sws(nexttile(tl,4),x,z,c_filt,gt_inc, ...
        sprintf('Filtered phase c_x (dir %d)',chosen),[0 8]);
    plot_sws(nexttile(tl,5),x,z,PG,gt_inc, ...
        'USTB PG (after filter)',[0 8]);
    plot_sws(nexttile(tl,6),x,z,AC,gt_inc, ...
        'USTB AC (after filter)',[0 8]);
    plot_sws(nexttile(tl,7),x,z,c2_raw,gt_inc, ...
        'Raw phase c_{xz}',[0 8]);
    plot_sws(nexttile(tl,8),x,z,c2_filt,gt_inc, ...
        'Filtered phase c_{xz}',[0 8]);
    plot_sws(nexttile(tl,9),x,z,abs(U_raw),gt_inc, ...
        'Raw phase amplitude',[]);
    plot_sws(nexttile(tl,10),x,z,abs(U_filt(:,:,chosen)),gt_inc, ...
        'Filtered phase amplitude',[]);
    title(tl,'high contrast | k-Wave GT: raw vs filtered phase | same scan', ...
        'Interpreter','none');
    exportgraphics(fig,fullfile(outdir,'phase_filter_maps.png'),'Resolution',170);
    close(fig);

    % Lateral profile at geometrical center of known inclusion.
    zs=find(any(gt_inc,2));
    assert(~isempty(zs),'Inclusion contour is missing.');
    zmid=zs(round(numel(zs)/2));
    fig2=figure('Color','w','Position',[110 110 1280 520],'Visible','off');
    ax=axes('Parent',fig2);
    plot(ax,x*1e3,ctrue(zmid,:),'k-','LineWidth',2.1);hold(ax,'on');
    plot(ax,x*1e3,c_raw(zmid,:),'b-','LineWidth',1.6);
    plot(ax,x*1e3,c_filt(zmid,:),'r-','LineWidth',1.6);
    plot(ax,x*1e3,PG(zmid,:),'Color',[0.15 0.60 0.37],'LineWidth',1.5);
    plot(ax,x*1e3,AC(zmid,:),'Color',[0.52 0.28 0.62],'LineWidth',1.4);
    grid(ax,'on');ylim(ax,[0 10]);
    xlabel(ax,'Lateral [mm]');ylabel(ax,'Apparent speed [m/s]');
    title(ax,sprintf('Through inclusion at z = %.2f mm | filtered dir %d', ...
        1e3*z(zmid),chosen));
    legend(ax,{'True material','Raw phase','Filtered phase', ...
        'PG after filter','AC after filter'},'Location','best');
    exportgraphics(fig2,fullfile(outdir,'phase_filter_profiles.png'), ...
        'Resolution',180);
    close(fig2);

    metadata=struct('datafile',datafile,'comparisonFile',mfile, ...
        'PRF',PRF,'targetFrequency',f0,'filterSpeed',vfilter, ...
        'packetAveragedIncrements',n_avg,'GTNativeSpacing',[dx dz], ...
        'spatialPhaseHalfLag_m',1e-3,'selectedDirectionIndex',chosen, ...
        'directionSelection','maximum background phasor energy', ...
        'c_true_inclusion',c_inc,'c_true_background',c_bg, ...
        'limitations',['Phase slope measures apparent local propagation ' ...
        'not material speed in mixed wavefields; GT native resolution ' ...
        '1mm, interpolation supplies no new spatial detail; packet ' ...
        'alignment follows previous diagnostic assumption.']);
    save(fullfile(outdir,'phase_filter_results.mat'), ...
        'T','Aangles','c_raw','c2_raw','c_filt','c2_filt', ...
        'c_native','c2_native','PG','AC','ctrue','gt_inc', ...
        'inc','bg','same','x','z','xn','zn','metadata','-v7.3');
    fprintf('\nSaved to: %s\n',outdir);
    fprintf(['Interpretation: native -> raw dense tests sampling/interpolation; ' ...
        'raw dense -> filtered tests directional processing; ' ...
        'filtered -> PG/AC tests phase-to-SWS estimation.\n']);
    fprintf('Do NOT conclude that filtering is causal solely from c_x differences.\n');
end

function U = harmonic_phasor(u,PRF,f0)
    [nz,nx,nt]=size(u);
    t=(0:nt-1)'/PRF;
    B=[ones(nt,1),cos(2*pi*f0*t),sin(2*pi*f0*t)];
    Y=reshape(double(u),[],nt).';
    coeff=B\Y;
    U=reshape(coeff(2,:)-1i*coeff(3,:),nz,nx);
end

function [kx,kz] = phase_gradient(U,x,z,half_lag_m)
% Central phase difference with ~2-mm total spacing, using phasor products.
    [nz,nx]=size(U);
    sx=max(1,round(half_lag_m/median(diff(x))));
    sz=max(1,round(half_lag_m/median(diff(z))));
    assert(nx>2*sx && nz>2*sz, ...
        'Grid too small for chosen phase-difference lag.');
    kx=nan(nz,nx);
    kz=nan(nz,nx);
    left=1:nx-2*sx; ctr=left+sx; right=left+2*sx;
    dx=(x(right)-x(left)).';
    kx(:,ctr)=angle(conj(U(:,left)).*U(:,right))./dx;
    top=1:nz-2*sz; mid=top+sz; bottom=top+2*sz;
    dz=z(bottom)-z(top);
    kz(mid,:)=angle(conj(U(top,:)).*U(bottom,:))./dz;
end

function c = phase_speed(k,f0)
    c=2*pi*f0./abs(k);
    c(~isfinite(c) | c<=0)=NaN;
end

function [med,avg,coverage] = stage_stat(c,mask)
    total=nnz(mask);
    v=c(mask);
    ok=isfinite(v) & v>0 & v<12;
    coverage=nnz(ok)/max(total,1);
    if any(ok)
        med=median(v(ok));avg=mean(v(ok));
    else
        med=NaN;avg=NaN;
    end
end

function A = fetch_phasor(v,nz,nx,propname)
% Read a complex numeric component from USTB data_field or numeric field.
% USTB custom branches vary in accessors; handle common representations.
    a=[];
    if isnumeric(v)
        a=v;
    elseif isstruct(v)
        if isfield(v,'u'),a=v.u;
        elseif isfield(v,'data'),a=v.data;
        end
    elseif isobject(v)
        if ismethod(v,'get_component')
            try, a=v.get_component('u'); catch, end
        end
        if isempty(a) && isprop(v,'u')
            try, a=v.u; catch, end
        end
        if isempty(a) && isprop(v,'data')
            try, a=v.data; catch, end
        end
    end
    if isstruct(a) && isfield(a,'u'),a=a.u;end
    if iscell(a) && numel(a)==1,a=a{1};end
    if isempty(a) || ~isnumeric(a)
        error(['Could not read numeric complex field "%s" (class %s). ' ...
            'Inspect properties(elasto.%s) and contact with that output.'], ...
            propname,class(v),propname);
    end
    if isvector(a)
        assert(mod(numel(a),nz*nx)==0, ...
            'Vectorized %s field does not match scan.',propname);
        a=reshape(a,nz,nx,[]);
    else
        dims=size(a);
        assert(numel(dims)>=2 && dims(1)==nz && dims(2)==nx, ...
            '%s size %s does not match scan [%d %d].', ...
            propname,mat2str(dims),nz,nx);
        a=reshape(a,nz,nx,[]);
    end
    A=double(a);
    if isreal(A)
        % Some USTB versions save phase ANGLE rather than complex phasor.
        av=A(isfinite(A));
        if ~isempty(av) && min(av)>=-pi-0.05 && max(av)<=2*pi+0.05
            warning('%s appears to store wrapped phase angles; using exp(i*phase).', ...
                propname);
            A=exp(1i*A);
        else
            error(['Field %s is real-valued, not an identifiable phase ' ...
                'image or complex phasor. Stop: cannot assess phase gradients.'], ...
                propname);
        end
    end
    assert(all(isfinite(real(A(:))) | isnan(real(A(:)))) && ...
           all(isfinite(imag(A(:))) | isnan(imag(A(:)))), ...
           '%s contains Inf values.',propname);
end

function plot_sws(ax,x,z,m,outline,ttl,limits)
    imagesc(ax,x*1e3,z*1e3,m);axis(ax,'image');
    set(ax,'YDir','normal');
    colorbar(ax);
    if ~isempty(limits),clim(ax,limits);end
    xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
    title(ax,ttl,'Interpreter','none');
    if ~isempty(outline)
        hold(ax,'on');
        contour(ax,x*1e3,z*1e3,double(outline),[0.5 0.5], ...
            'w','LineWidth',1.2);
    end
end
