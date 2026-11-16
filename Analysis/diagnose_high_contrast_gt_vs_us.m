function T = diagnose_high_contrast_gt_vs_us()
%DIAGNOSE_HIGH_CONTRAST_GT_VS_US  Where is high-contrast SWS information lost?
% Revision 2026-11-11: Fix two invalid multiline fprintf format strings.
%
% Run from the synthetic-the-rf-data-generation repository:
%   T = diagnose_high_contrast_gt_vs_us();
%
% One controlled experiment, FOUR maps:
%       k-Wave displacement -> PG / AC
%       US-estimated displacement -> PG / AC
% SAME scan, temporal frequency, direction filters and evaluation pixels.
%
% No new simulation, no RF generation and NO changes to previous outputs.
% Earlier GT-vs-US implementation: compare_gt_displacement_sws.m
% Earlier PG/AC implementation: compare_SWS_methods_batch.m
%
% Cautions:
%  * k-Wave GT is typically a 46x39 grid, which is spatially INTERPOLATED
%    to the existing ultrasound grid. Interpolation does not recover
%    missing wavefield detail. If GT also fails, this test alone cannot
%    distinguish wave physics from limited original spatial sampling.
%  * Temporal packet alignment follows the prior project's assumption:
%    average packet_size-1 consecutive PRF displacement increments starting
%    with frame 1. This is NOT a verified acquisition-index alignment.
%  * Comparing 4 reconstructions on a cropped scan alters filtering edges
%    compared to saved full-scan results. The script prints this difference.
%
% Outputs:
%   Results/sws_estimator_diagnostic/gt_vs_us_high_contrast/
%      comparison_summary.csv
%      comparison_maps.png
%      inclusion_profile.png
%      comparison_results.mat
%
% Only reads previously saved data. Requires MATLAB, USTB, k-Wave/create_phantom.

    case_name = 'high_contrast';
    here = fileparts(mfilename('fullpath'));
    if isfolder(fullfile(here,'Data')) && isfolder(fullfile(here,'Results'))
        repo = here;
    else
        repo = fileparts(here);
    end
    assert(isfolder(fullfile(repo,'Data')) && isfolder(fullfile(repo,'Results')), ...
        'Put this file in <repository>/Analysis/.');

    addpath(genpath(repo));
    ext = fullfile(fileparts(repo), 'elastography', 'ElastographyToolbox');
    assert(isfolder(ext), 'Cannot locate ElastographyToolbox: %s',ext);
    addpath(genpath(ext));
    for s = {'create_phantom','kWaveGrid','uff.linear_scan', ...
             'uff.elastography','postprocess.elastography.phase_gradient_method', ...
             'postprocess.elastography.autocorrelation_SWS_method'}
        assert(~isempty(which(s{1})), 'Required class/function not found: %s',s{1});
    end

    data_name = [case_name '_displacement_with_quality.mat'];
    data_candidates = {fullfile(repo,'Data',case_name,data_name), ...
        fullfile(repo,'Data',data_name),fullfile(repo,data_name), ...
        fullfile(here,data_name)};
    idx = find(cellfun(@isfile,data_candidates),1,'first');
    assert(~isempty(idx),'Cannot find %s in repository.',data_name);
    data_file = data_candidates{idx};
    comparison_file = fullfile(repo,'Results','method_comparison', ...
        'batch',case_name,'method_comparison.mat');
    assert(isfile(comparison_file), ...
        'Missing prior method comparison: %s',comparison_file);

    fprintf('\n=== HIGH CONTRAST: k-Wave GT vs ultrasound displacement ===\n');
    fprintf('Input: %s\n',data_file);
    fprintf('Saved PG/AC: %s\n',comparison_file);

    % Check metadata before loading large USTB objects.
    listing = whos('-file',data_file);
    names = {listing.name};
    required = {'disp_gt','b_data','displacement_estimation_autocorr', ...
        'shear_params'};
    missing = required(~ismember(required,names));
    if ~isempty(missing)
        error(['Original displacement MAT lacks required fields: %s. ' ...
               'Do NOT substitute sws_gt for a time-dependent wavefield. ' ...
               'Find the saved k-Wave disp_gt sequence first.'], ...
               strjoin(missing,', '));
    end
    wanted = {'disp_gt','b_data','displacement_estimation_autocorr', ...
        'shear_params','autocorr_params','array_length', ...
        'scene_depth','prb_theta','prf'};
    wanted = wanted(ismember(wanted,names));
    D = load(data_file,wanted{:});
    M = load(comparison_file,'PG','AC','GT','roi','x_axis','z_axis', ...
        'PRF','target_frequency','filter_speed');
    need = {'PG','AC','GT','roi','x_axis','z_axis', ...
            'PRF','target_frequency','filter_speed'};
    assert(all(isfield(M,need)), ...
        'Saved method comparison lacks expected variables.');
    assert(all(isfield(M.roi,{'image_mask','inclusion_core', ...
                            'background_core'})), ...
        'Method comparison lacks required ROI masks.');

    x_full = double(M.x_axis(:));
    z_full = double(M.z_axis(:));
    Nx = numel(x_full); Nz = numel(z_full);
    assert(isequal(size(M.GT),[Nz,Nx],size(M.PG),size(M.AC)), ...
        'Saved SWS maps/axes have inconsistent sizes.');
    assert(all(diff(x_full)>0) && all(diff(z_full)>0), ...
        'US scan coordinate axes must be strictly increasing.');
    scan_in = D.b_data.scan;
    ax_x = double(scan_in.x_axis(:));
    ax_z = double(scan_in.z_axis(:));
    assert(isequal(size(ax_x),size(x_full)) && ...
           isequal(size(ax_z),size(z_full)) && ...
           max(abs(ax_x-x_full))<1e-10 && ...
           max(abs(ax_z-z_full))<1e-10, ...
        'Saved SWS and estimated displacement use different ultrasound grids.');
    clear scan_in;

    PRF=double(M.PRF);
    f0=double(M.target_frequency);
    filter_speed=double(M.filter_speed);
    assert(all(isfinite([PRF,f0,filter_speed])) && ...
        PRF>2*f0 && f0>0 && filter_speed>0, ...
        'Invalid saved PRF, frequency or filter-speed parameters.');
    assert(abs(f0-double(D.shear_params.source_freq))<1e-6, ...
        'Saved SWS target frequency differs from original k-Wave source.');
    if isfield(D,'prf') && ~isempty(D.prf)
        assert(abs(double(D.prf)-PRF)<1e-6, ...
            'Original displacement PRF differs from prior method comparison.');
    end
    if isfield(D,'prb_theta') && ~isempty(D.prb_theta)
        fprintf('Stored probe orientation: %s degrees\n', ...
            mat2str(rad2deg(double(D.prb_theta))));
        if abs(rad2deg(double(D.prb_theta(1)))-90)>0.1
            warning('Probe angle is not 90 degrees in original MAT.');
        end
    else
        fprintf('Probe angle not saved directly; do not infer from filename.\n');
    end

    if isfield(D,'autocorr_params') && ...
            isfield(D.autocorr_params,'packet_size')
        packet_size=double(D.autocorr_params.packet_size);
    else
        packet_size=4;  % identical to prior GT-vs-US script
        warning('packet_size missing; using previous project value 4.');
    end
    assert(isscalar(packet_size) && packet_size>=2 && ...
        packet_size==round(packet_size),'Invalid packet_size.');
    n_avg = packet_size-1;

    u_vec = D.displacement_estimation_autocorr.data;
    assert(mod(numel(u_vec),Nz*Nx)==0, ...
        'Estimated displacement cube does not match US scan.');
    N_est=numel(u_vec)/(Nz*Nx);
    [Ngz,Ngx,Nraw]=size(D.disp_gt);
    assert(Nraw>=N_est+n_avg-1, ...
        ['Not enough GT frames for existing packet convention: ' ...
         'GT=%d, estimated=%d, packet intervals=%d.'],Nraw,N_est,n_avg);
    if Nraw~=N_est+n_avg-1
        warning(['GT and estimated frame counts do not exactly satisfy ' ...
                 'packet alignment. Only the first %d GT frames are used.'], ...
                N_est+n_avg-1);
    end

    % Recover the exact native k-Wave spacing from the existing phantom
    % configuration; creates grid/material arrays ONLY, no simulation.
    [kg,~] = create_phantom(D.shear_params);
    dx=double(kg.dx); dz=double(kg.dz);
    clear kg;
    % These are the same documented defaults as get_GT_kwave and the
    % earlier compare_gt_displacement_sws.m. Demand consistency with dims.
    if isfield(D,'array_length') && ~isempty(D.array_length)
        array_length=double(D.array_length);
    else
        array_length=(128-1)*0.300e-3;
        warning('array_length absent: using documented probe N=128,pitch=0.3 mm.');
    end
    if isfield(D,'scene_depth') && ~isempty(D.scene_depth)
        scene_depth=double(D.scene_depth);
    else
        scene_depth=50e-3;
        warning('scene_depth absent: using documented 50 mm.');
    end
    x_gt=(-array_length/2:dx:array_length/2).';
    z_gt=(5e-3:dz:scene_depth).';
    assert(numel(x_gt)==Ngx && numel(z_gt)==Ngz, ...
        ['Inferred GT axes (%dx%d) differ from saved disp_gt (%dx%d). ' ...
         'STOP: check array_length, scene_depth, k-Wave spacing and ' ...
         'get_GT_kwave before interpreting speeds.'], ...
         numel(z_gt),numel(x_gt),Ngz,Ngx);
    assert(all(diff(x_gt)>0) && all(diff(z_gt)>0), ...
        'GT axes must increase.');

    % Same cropped scan FOR BOTH wavefield sources. No extrapolation.
    tol=1e-9;
    ix=find(x_full>=x_gt(1)-tol & x_full<=x_gt(end)+tol);
    iz=find(z_full>=z_gt(1)-tol & z_full<=z_gt(end)+tol);
    assert(numel(ix)>20 && numel(iz)>20, 'No sufficient grid overlap.');
    x=x_full(ix); z=z_full(iz);
    xq=min(max(x,x_gt(1)),x_gt(end));
    zq=min(max(z,z_gt(1)),z_gt(end));
    [Xq,Zq]=meshgrid(xq,zq);
    sc=uff.linear_scan('x_axis',x,'z_axis',z);

    c_true=double(M.GT(iz,ix));
    r=M.roi;
    roi_global=logical(r.image_mask(iz,ix));
    roi_incl=logical(r.inclusion_core(iz,ix));
    roi_bg=logical(r.background_core(iz,ix));
    gt_speed_core=median(c_true(roi_incl & isfinite(c_true)));
    bg_speed_core=median(c_true(roi_bg & isfinite(c_true)));
    fprintf('Material GT: inclusion %.4f, background %.4f m/s\n', ...
        gt_speed_core,bg_speed_core);
    fprintf('GT wavefield %dx%dx%d; US grid %dx%dx%d\n', ...
        Ngz,Ngx,Nraw,Nz,Nx,N_est);
    fprintf('Comparison scan %dx%d (cropped before BOTH filters).\n', ...
        numel(iz),numel(ix));
    fprintf('PRF=%.2f Hz; f0=%.2f Hz; filter speed=%.3f m/s; packet=%d.\n', ...
        PRF,f0,filter_speed,packet_size);
    fprintf(['NOTE: GT wavefield is coarsely sampled (dx=%.3f mm, ' ...
             'dz=%.3f mm) and interpolated spatially.\n'],1000*dx,1000*dz);

    % Build the GT displacement on the common scan; packet averages over
    % the same PRF increments used by the previous project diagnostic.
    u_gt=zeros(numel(iz),numel(ix),N_est,'single');
    for t=1:N_est
        g=mean(double(D.disp_gt(:,:,t:t+n_avg-1)),3);
        q=interp2(x_gt,z_gt,g,Xq,Zq,'linear',NaN);
        assert(all(isfinite(q(:))), ...
            'GT interpolation returned nonfinite values (frame %d).',t);
        u_gt(:,:,t)=single(q);
    end
    clear g q Xq Zq;

    % Match original estimated displacement sign; a global sign does not
    % change phase-gradient magnitudes, but preserves established convention.
    fprintf('Reconstructing PG + AC from k-Wave displacement ...\n');
    [pg_gt,ac_gt]=sws_pair(u_gt,sc,PRF,f0,filter_speed);
    clear u_gt;

    % Reshape raw US estimated displacement and crop to EXACT same scan.
    % Keep D only until the cropped cube is assembled, to limit memory.
    u_est_full=reshape(u_vec,Nz,Nx,N_est);
    u_est=single(-u_est_full(iz,ix,:));
    clear u_est_full u_vec D;
    fprintf('Reconstructing PG + AC from US-estimated displacement ...\n');
    [pg_us,ac_us]=sws_pair(u_est,sc,PRF,f0,filter_speed);
    clear u_est;

    sz=size(c_true);
    assert(isequal(size(pg_gt),sz,size(ac_gt),sz, ...
                  size(pg_us),sz,size(ac_us),sz), ...
        'Reconstructed SWS dimensions are not the common scan dimensions.');
    all_valid=isfinite(c_true) & ...
        isfinite(pg_gt) & pg_gt>0 & pg_gt<10 & ...
        isfinite(ac_gt) & ac_gt>0 & ac_gt<10 & ...
        isfinite(pg_us) & pg_us>0 & pg_us<10 & ...
        isfinite(ac_us) & ac_us>0 & ac_us<10;
    global_mask=all_valid & roi_global;
    incl_mask=all_valid & roi_incl;
    bg_mask=all_valid & roi_bg;
    assert(nnz(global_mask)>100 && nnz(incl_mask)>100 && ...
           nnz(bg_mask)>100, ...
           'Too few common valid global/inclusion/background pixels.');
    fprintf('IDENTICAL valid pixels: global %d, inclusion %d, background %d.\n', ...
        nnz(global_mask),nnz(incl_mask),nnz(bg_mask));

    % One mask for all four combinations (no unfair excluded pixels).
    entries={pg_gt,ac_gt,pg_us,ac_us};
    estimator=["PG";"AC";"PG";"AC"];
    source=["k-Wave GT";"k-Wave GT";"Ultrasound";"Ultrasound"];
    n=4;
    global_rmse=nan(n,1); inclusion_rmse=nan(n,1);
    background_rmse=nan(n,1); inclusion_bias=nan(n,1);
    background_bias=nan(n,1); inclusion_mean=nan(n,1);
    background_mean=nan(n,1); contrast_recovery=nan(n,1);
    true_contrast=gt_speed_core-bg_speed_core;
    for j=1:n
        e=entries{j}-c_true;
        global_rmse(j)=sqrt(mean(e(global_mask).^2));
        inclusion_rmse(j)=sqrt(mean(e(incl_mask).^2));
        background_rmse(j)=sqrt(mean(e(bg_mask).^2));
        inclusion_bias(j)=mean(e(incl_mask));
        background_bias(j)=mean(e(bg_mask));
        inclusion_mean(j)=mean(entries{j}(incl_mask));
        background_mean(j)=mean(entries{j}(bg_mask));
        contrast_recovery(j)=(inclusion_mean(j)-background_mean(j)) ...
            /true_contrast;
    end
    T=table(source,estimator,global_rmse,inclusion_rmse,background_rmse, ...
        inclusion_bias,background_bias,inclusion_mean,background_mean, ...
        contrast_recovery, ...
        'VariableNames',{'Source','Estimator','GlobalRMSE_mps', ...
            'InclusionRMSE_mps','BackgroundRMSE_mps', ...
            'InclusionBias_mps','BackgroundBias_mps', ...
            'InclusionMean_mps','BackgroundMean_mps','ContrastRecovery'});
    fprintf('\n--- SAME-PIXEL COMPARISON ---\n');
    disp(T);

    % Audit difference vs original full-scan saved method comparison.
    origpg=double(M.PG(iz,ix));
    origac=double(M.AC(iz,ix));
    report_map_diff('PG',origpg,pg_us,global_mask);
    report_map_diff('AC',origac,ac_us,global_mask);
    clear M;

    outdir=fullfile(repo,'Results','sws_estimator_diagnostic', ...
        'gt_vs_us_high_contrast');
    if ~isfolder(outdir),mkdir(outdir);end
    writetable(T,fullfile(outdir,'comparison_summary.csv'));

    % Map overview with fixed units and GT inclusion outline.
    gt_inc=isfinite(c_true) & abs(c_true-gt_speed_core) < ...
        abs(c_true-bg_speed_core);
    fig=figure('Color','w','Position',[80 60 1620 930], ...
               'Visible','off');
    tiled=tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
    pmap(nexttile(tiled,1),x,z,c_true,[0 8],gt_inc,'Material GT SWS');
    pmap(nexttile(tiled,2),x,z,pg_gt,[0 8],gt_inc,'PG from k-Wave GT');
    pmap(nexttile(tiled,3),x,z,pg_us,[0 8],gt_inc,'PG from ultrasound');
    pmap(nexttile(tiled,4),x,z,pg_us-pg_gt,[-4 4],gt_inc, ...
        'PG ultrasound - GT input');
    pmap(nexttile(tiled,5),x,z,double(global_mask),[0 1],[], ...
        'Common valid ROI');
    pmap(nexttile(tiled,6),x,z,ac_gt,[0 8],gt_inc,'AC from k-Wave GT');
    pmap(nexttile(tiled,7),x,z,ac_us,[0 8],gt_inc,'AC from ultrasound');
    pmap(nexttile(tiled,8),x,z,ac_us-ac_gt,[-4 4],gt_inc, ...
        'AC ultrasound - GT input');
    title(tiled,'high contrast | original k-Wave GT vs US displacement | same SWS settings', ...
        'Interpreter','none');
    exportgraphics(fig,fullfile(outdir,'comparison_maps.png'),'Resolution',155);
    close(fig);

    % One lateral slice through the TRUE inclusion: interpretable overlap.
    idx_z=find(any(gt_inc,2));
    assert(~isempty(idx_z),'No GT inclusion found in common image.');
    iz0=idx_z(round(numel(idx_z)/2));
    fig2=figure('Color','w','Position',[100 100 1240 510], ...
                'Visible','off');
    tl=tiledlayout(fig2,1,2,'TileSpacing','compact');
    for j=1:2
        ax=nexttile(tl,j);
        if j==1
            gt_recon=pg_gt; us_recon=pg_us;
        else
            gt_recon=ac_gt; us_recon=ac_us;
        end
        plot(ax,x*1e3,c_true(iz0,:),'k-','LineWidth',2);hold(ax,'on');
        plot(ax,x*1e3,gt_recon(iz0,:),'b-','LineWidth',1.5);
        plot(ax,x*1e3,us_recon(iz0,:),'r-','LineWidth',1.5);
        xlabel(ax,'Lateral [mm]');ylabel(ax,'SWS [m/s]');
        title(ax,sprintf('%s at depth %.2f mm',char(estimator(j)), ...
              z(iz0)*1e3));
        ylim(ax,[0 9]); grid(ax,'on');
        legend(ax,{'Material GT','k-Wave GT input','US input'}, ...
            'Location','best');
    end
    exportgraphics(fig2,fullfile(outdir,'inclusion_profile.png'), ...
        'Resolution',170);
    close(fig2);

    meta=struct('data_file',data_file,'method_file',comparison_file, ...
        'PRF',PRF,'target_frequency',f0,'filter_speed',filter_speed, ...
        'packet_size',packet_size,'GT_native_spacing_dx_m',dx, ...
        'GT_native_spacing_dz_m',dz, ...
        'native_size',[Ngz Ngx Nraw], ...
        'US_source_size',[Nz Nx N_est], ...
        'comparison_size',size(c_true), ...
        'notes',['GT wavefield is interpolated from native coarse grid; ' ...
          'packet alignment assumes packet_size-1 consecutive increments. ' ...
          'Do not attribute errors uniquely to physics, RF or estimator ' ...
          'without checking this assumption.']);
    save(fullfile(outdir,'comparison_results.mat'), ...
        'pg_gt','ac_gt','pg_us','ac_us','c_true','x','z', ...
        'global_mask','incl_mask','bg_mask','T','meta','-v7.3');
    fprintf('Outputs: %s\n',outdir);
    fprintf(['\nDECISION: If GT input also has poor inclusion recovery, ' ...
        'the loss occurs upstream of RF or in GT-grid/interpolation/SWS ' ...
        'model. If GT input succeeds but US fails, focus on RF/displacement.\n']);
end

function [PG,AC]=sws_pair(u,scan,PRF,f0,filter_speed)
% EXACT sequence/parameters from compare_SWS_methods_batch.m.
    elasto=uff.elastography();
    elasto.PRF=PRF;
    elasto.target_frequency=f0;
    elasto.scan=scan;
    elasto.autocorrelation_field=uff.data_field();
    elasto.autocorrelation_field.set_component('u', ...
        reshape(u,size(u,1),size(u,2),1,size(u,3)));
    num_angles=12;
    elasto.directional_filtering_angles= ...
        2*pi*linspace(0,1-1/num_angles,num_angles);
    extractor=postprocess.elastography.extract_freq_component(elasto);
    extractor.target_frequency=f0;
    extractor.component='u';
    extractor.run();
    omega=2*pi*f0;
    elasto.radial_filter=uff.spatial_filter('gaussian_charite','radial');
    elasto.radial_filter.sigma=filter_speed/omega;
    elasto.angular_filter=uff.spatial_filter('gaussian_charite','angular');
    elasto.angular_filter.deltaTheta= ...
        median(diff(elasto.directional_filtering_angles));
    filt=postprocess.elastography.apply_filter(elasto);
    filt.run();
    pgm=postprocess.elastography.phase_gradient_method(elasto);
    pgm.run();
    PG=squeeze(double(elasto.sws_map.u));
    acm=postprocess.elastography.autocorrelation_SWS_method(elasto);
    acm.run();
    AC=squeeze(double(elasto.sws_map.u));
end

function pmap(ax,x,z,map,lim,outline,ttl)
    imagesc(ax,x*1e3,z*1e3,map);
    set(ax,'YDir','normal');axis(ax,'image');
    clim(ax,lim);colorbar(ax);
    xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
    title(ax,ttl,'Interpreter','none');
    if ~isempty(outline)
        hold(ax,'on');
        contour(ax,x*1e3,z*1e3,double(outline),[.5 .5], ...
            'w','LineWidth',1.2);
    end
end

function report_map_diff(label,original,new,mask)
    ok=mask & isfinite(original) & isfinite(new);
    if ~any(ok(:))
        fprintf('Original %s vs cropped %s: no common valid pixels.\n', ...
            label,label);
        return;
    end
    d=new(ok)-original(ok);
    fprintf(['Original full-scan vs reprocessed/cropped %s: RMSE %.5f m/s, ' ...
        'median difference %.5f m/s (filter boundary change).\n'], ...
        label,sqrt(mean(d.^2)),median(d));
    if sqrt(mean(d.^2))>.3
        warning('Cropped %s results differ materially from stored full-scan map.',label);
    end
end
