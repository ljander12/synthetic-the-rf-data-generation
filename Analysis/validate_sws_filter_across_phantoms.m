function T = validate_sws_filter_across_phantoms()
%VALIDATE_SWS_FILTER_ACROSS_PHANTOMS  Fixed-filter holdout validation.
%
% T = validate_sws_filter_across_phantoms();
%
% PURPOSE
%   Validate the filter candidate found on high_contrast without tuning it.
%   Holdout phantoms: large_inclusions and low_contrast.
%   Compare ORIGINAL radial parameter 3.2 m/s and CANDIDATE 1.6 m/s,
%   with 12 directions, unchanged 30-degree angular-filter width and
%   EXACTLY the same PG / spatial-autocorrelation USTB implementations.
%
%   Two displacement inputs are evaluated:
%      GT: packet-averaged k-Wave ground-truth displacement interpolated to
%          the overlapping ultrasound scan, identical to the development
%          experiment used to identify the high_contrast candidate.
%      US: ultrasound-estimated displacement on the same cropped scan;
%          this is the relevant check for actual analysis-pipeline use.
%
%   All 8 combinations within each phantom (2 inputs x 2 filters x
%   2 estimators) are evaluated on an identical, fixed intersection of
%   valid pixels, based on *saved* ROI masks. No parameter retuning.
%
% INPUTS (read-only)
%   <case>_displacement_with_quality.mat
%   Results/method_comparison/batch/<case>/method_comparison.mat
%
% OUTPUTS (new folder; never modifies prior results)
%   Results/sws_estimator_diagnostic/filter_holdout_validation/
%       filter_holdout_summary.csv, filter_holdout_decision.csv,
%       filter_holdout_comparison.png, filter_holdout_results.mat
%       <case>/sws_filter_maps.png
%       <case>/filter_results.csv
%
% NOTE
%   MATLAB/USTB cannot be executed here; run this in the project's
%   configured MATLAB installation. Stop on metadata/grid mismatch.

here = fileparts(mfilename('fullpath'));
if isfolder(fullfile(here,'Data')) && isfolder(fullfile(here,'Results'))
    repo = here;
else
    repo = fileparts(here);
end
assert(isfolder(fullfile(repo,'Data')) && isfolder(fullfile(repo,'Results')), ...
    'Place this file in the repository Analysis/ folder.');
addpath(genpath(repo));
toolbox_dir = fullfile(fileparts(repo),'elastography','ElastographyToolbox');
assert(isfolder(toolbox_dir),'ElastographyToolbox not found: %s',toolbox_dir);
addpath(genpath(toolbox_dir));

required = {'create_phantom','uff.elastography','uff.data_field', ...
    'uff.linear_scan','postprocess.elastography.extract_freq_component', ...
    'postprocess.elastography.apply_filter', ...
    'postprocess.elastography.phase_gradient_method', ...
    'postprocess.elastography.autocorrelation_SWS_method'};
for k=1:numel(required)
    assert(~isempty(which(required{k})), ...
        'Missing MATLAB/USTB function: %s',required{k});
end

% LOCKED settings, selected BEFORE inspecting the two holdout phantoms.
cases = {'large_inclusions','low_contrast'};
filter_values = [3.2, 1.6];  % m/s as radial FILTER parameters, not tissue c
nangles = 12;
angles = 2*pi*(0:nangles-1)/nangles;
angular_width = median(diff(angles)); % 30 degrees, unchanged
sources = {'GT','US'};
methods = {'PG','AC'};

outroot=fullfile(repo,'Results','sws_estimator_diagnostic', ...
    'filter_holdout_validation');
if ~isfolder(outroot), mkdir(outroot); end
T=table();
coverage_rows=table();

for caseid=1:numel(cases)
    name=cases{caseid};
    fprintf('\n========================================================\n');
    fprintf('FILTER HOLDOUT %d/%d: %s\n',caseid,numel(cases),name);
    fprintf('========================================================\n');
    [D,M,datafile] = load_case(repo,name);
    [uGT,uUS,GT,x,z,rois,meta] = build_inputs(D,M,name);
    clear D M;

    fprintf('GT material: BG %.3f m/s, INCL %.3f m/s | f0 %.1f Hz\n', ...
        meta.c_bg,meta.c_inc,meta.f0);
    fprintf('Probe orientation %.3f deg | scan %d x %d x %d\n', ...
        meta.probe_deg,numel(z),numel(x),size(uGT,3));
    fprintf('Filters: original %.1f, candidate %.1f m/s; angle width %.1f deg\n', ...
        filter_values(1),filter_values(2),rad2deg(angular_width));

    maps=cell(2,2,2); % source, filter, method
    for s=1:2
        if s==1, u=uGT; else, u=uUS; end
        for f=1:2
            fprintf('  %s | radial %.1f m/s | PG + AC ...\n', ...
                sources{s},filter_values(f));
            [maps{s,f,1},maps{s,f,2}]=run_estimators( ...
                u,x,z,meta.PRF,meta.f0,angles, ...
                filter_values(f)/(2*pi*meta.f0),angular_width);
        end
    end
    clear uGT uUS u;

    % ONE fixed pixel mask, independent of candidate, data source and method.
    valid=isfinite(GT);
    for j=1:numel(maps)
        v=maps{j};
        valid=valid & isfinite(v) & v>0 & v<10;
    end
    mask_g=rois.global & valid;
    mask_i=rois.incl & valid;
    mask_b=rois.bg & valid;
    cover_g=100*nnz(mask_g)/nnz(rois.global);
    cover_i=100*nnz(mask_i)/nnz(rois.incl);
    cover_b=100*nnz(mask_b)/nnz(rois.bg);
    fprintf('Fixed-mask coverage: global %.2f%%, inclusion %.2f%%, BG %.2f%%\n', ...
        cover_g,cover_i,cover_b);
    assert(cover_g>=80 && cover_i>=80 && cover_b>=80, ...
        ['Too few common pixels. Do NOT compare. Check which branch ' ...
         'produces out-of-range or invalid SWS maps.']);
    assert(nnz(mask_g)>100 && nnz(mask_i)>50 && nnz(mask_b)>100, ...
        'Too few valid pixels for meaningful regional metrics.');

    % Metrics over the exact same pixels. Positive and negative contrast
    % handled by a signed denominator. No per-filter ROI redefinition.
    den=mean(GT(mask_i))-mean(GT(mask_b));
    assert(abs(den)>0.01,'Material contrast is too small to score.');
    Tc=table();
    for s=1:2
        for f=1:2
            for m=1:2
                v=maps{s,f,m};
                Phantom=string(name);
                Source=string(sources{s});
                Estimator=string(methods{m});
                FilterSetting_mps=filter_values(f);
                GlobalRMSE_mps=sqrt(mean((v(mask_g)-GT(mask_g)).^2));
                InclusionRMSE_mps=sqrt(mean((v(mask_i)-GT(mask_i)).^2));
                BackgroundRMSE_mps=sqrt(mean((v(mask_b)-GT(mask_b)).^2));
                InclusionMean_mps=mean(v(mask_i));
                BackgroundMean_mps=mean(v(mask_b));
                ContrastRecovery_pct=100*(InclusionMean_mps-BackgroundMean_mps)/den;
                CommonGlobal_pct=cover_g;
                CommonInclusion_pct=cover_i;
                CommonBackground_pct=cover_b;
                row=table(Phantom,Source,Estimator,FilterSetting_mps, ...
                    GlobalRMSE_mps,InclusionRMSE_mps,BackgroundRMSE_mps, ...
                    InclusionMean_mps,BackgroundMean_mps,ContrastRecovery_pct, ...
                    CommonGlobal_pct,CommonInclusion_pct,CommonBackground_pct);
                Tc=[Tc;row]; %#ok<AGROW>
            end
        end
    end
    fprintf('\nSAME-PIXEL RESULTS: %s\n',name);
    disp(Tc(:,{'Source','Estimator','FilterSetting_mps', ...
        'GlobalRMSE_mps','InclusionRMSE_mps', ...
        'BackgroundRMSE_mps','ContrastRecovery_pct'}));
    T=[T;Tc]; %#ok<AGROW>

    dest=fullfile(outroot,name);
    if ~isfolder(dest),mkdir(dest);end
    writetable(Tc,fullfile(dest,'filter_results.csv'));
    fig=figure('Color','w','Visible','off','Position',[60 60 1500 980]);
    tl=tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
    outline=GT>0.5*(meta.c_bg+meta.c_inc);
    plotmap(nexttile(tl,1),x,z,GT,outline, ...
        sprintf('%s: material GT',strrep(name,'_',' ')),[0 8]);
    plotmap(nexttile(tl,2),x,z,maps{2,1,1},outline, ...
        'US: PG original 3.2',[0 8]);
    plotmap(nexttile(tl,3),x,z,maps{2,2,1},outline, ...
        'US: PG candidate 1.6',[0 8]);
    plotmap(nexttile(tl,4),x,z,double(mask_g),[], ...
        sprintf('Fixed common ROI (%.1f%%)',cover_g),[0 1]);
    plotmap(nexttile(tl,5),x,z,maps{2,1,2},outline, ...
        'US: AC original 3.2',[0 8]);
    plotmap(nexttile(tl,6),x,z,maps{2,2,2},outline, ...
        'US: AC candidate 1.6',[0 8]);
    plotmap(nexttile(tl,7),x,z,maps{1,1,1},outline, ...
        'GT: PG original 3.2',[0 8]);
    plotmap(nexttile(tl,8),x,z,maps{1,2,1},outline, ...
        'GT: PG candidate 1.6',[0 8]);
    plotmap(nexttile(tl,9),x,z,maps{1,2,2},outline, ...
        'GT: AC candidate 1.6',[0 8]);
    title(tl,sprintf('%s | Filter holdout, 200 Hz, probe 90 deg', ...
        strrep(name,'_',' ')));
    exportgraphics(fig,fullfile(dest,'sws_filter_maps.png'),'Resolution',160);
    close(fig);
    save(fullfile(dest,'filter_maps.mat'),'maps','GT','x','z','mask_g', ...
        'mask_i','mask_b','meta','-v7.3');
end

writetable(T,fullfile(outroot,'filter_holdout_summary.csv'));

% Pre-registered, transparent rule applied to ULTRASOUND source only:
% For BOTH holdout phantoms AND both estimators, candidate must not worsen
% global or inclusion RMSE, and BG RMSE must not increase > 0.10 m/s.
% This is a screening rule, not statistical significance or proof of
% generalization to different anatomy/excitation frequencies.
D=table();
for caseid=1:numel(cases)
    for m=1:numel(methods)
        id0=T.Phantom==string(cases{caseid}) & ...
            T.Source=="US" & T.Estimator==string(methods{m}) & ...
            abs(T.FilterSetting_mps-filter_values(1))<1e-8;
        id1=T.Phantom==string(cases{caseid}) & ...
            T.Source=="US" & T.Estimator==string(methods{m}) & ...
            abs(T.FilterSetting_mps-filter_values(2))<1e-8;
        assert(nnz(id0)==1 && nnz(id1)==1,'Missing comparison table row.');
        row0=T(id0,:);row1=T(id1,:);
        Phantom=string(cases{caseid}); Estimator=string(methods{m});
        DeltaGlobalRMSE_mps=row1.GlobalRMSE_mps-row0.GlobalRMSE_mps;
        DeltaInclusionRMSE_mps=row1.InclusionRMSE_mps-row0.InclusionRMSE_mps;
        DeltaBackgroundRMSE_mps=row1.BackgroundRMSE_mps-row0.BackgroundRMSE_mps;
        DeltaContrastRecovery_pct=row1.ContrastRecovery_pct-row0.ContrastRecovery_pct;
        Pass=DeltaGlobalRMSE_mps<=0 && ...
            DeltaInclusionRMSE_mps<=0 && DeltaBackgroundRMSE_mps<=0.10;
        D=[D;table(Phantom,Estimator,DeltaGlobalRMSE_mps, ...
            DeltaInclusionRMSE_mps,DeltaBackgroundRMSE_mps, ...
            DeltaContrastRecovery_pct,Pass)]; %#ok<AGROW>
    end
end
writetable(D,fullfile(outroot,'filter_holdout_decision.csv'));

fprintf('\n=== HOLDOUT VALIDATION DECISION (ULTRASOUND INPUT) ===\n');
disp(D);
if all(D.Pass)
    fprintf(['Both holdout phantoms support 1.6 m/s for PG and AC under ' ...
        'the predefined RMSE criterion. CANDIDATE for pipeline update; ' ...
        'still verify full-scan behavior and preserve baseline.\n']);
else
    fprintf(['1.6 m/s FAILS the predefined criterion for at least one ' ...
        'phantom/estimator. Do NOT change the production filter globally.\n']);
end

fig=figure('Color','w','Visible','off','Position',[70 70 1300 740]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
metric={'GlobalRMSE_mps','InclusionRMSE_mps', ...
    'BackgroundRMSE_mps','ContrastRecovery_pct'};
labels={'Global RMSE [m/s]','Inclusion RMSE [m/s]', ...
    'Background RMSE [m/s]','Contrast recovery [%]'};
for q=1:4
    ax=nexttile(tl,q);
    vals=nan(4,2);
    names=cell(4,1);
    p=0;
    for caseid=1:2
        for m=1:2
            p=p+1;
            names{p}=sprintf('%s %s',strrep(cases{caseid},'_',' '),methods{m});
            for f=1:2
                ix=T.Phantom==string(cases{caseid}) & T.Source=="US" & ...
                    T.Estimator==string(methods{m}) & ...
                    abs(T.FilterSetting_mps-filter_values(f))<1e-8;
                vals(p,f)=T.(metric{q})(ix);
            end
        end
    end
    bar(ax,vals);grid(ax,'on');
    set(ax,'XTick',1:4,'XTickLabel',names);
    xtickangle(ax,20);ylabel(ax,labels{q});
    if q==1, legend(ax,{'3.2 original','1.6 candidate'},'Location','best');end
end
title(tl,'Filter holdout validation | Ultrasound displacement | fixed-mask');
exportgraphics(fig,fullfile(outroot,'filter_holdout_comparison.png'), ...
    'Resolution',160);
close(fig);
save(fullfile(outroot,'filter_holdout_results.mat'),'T','D', ...
    'cases','filter_values','angles','angular_width');
fprintf('Saved results: %s\n',outroot);
fprintf('No simulations or production MATLAB files were modified.\n');
end

function [D,M,f]=load_case(repo,name)
    fn=[name '_displacement_with_quality.mat'];
    choices={fullfile(repo,'Data',name,fn),fullfile(repo,'Data',fn), ...
        fullfile(repo,fn)};
    idx=find(cellfun(@isfile,choices),1);
    assert(~isempty(idx),'Missing data file for %s: %s',name,fn);
    f=choices{idx};
    result_file=fullfile(repo,'Results','method_comparison','batch', ...
        name,'method_comparison.mat');
    assert(isfile(result_file),'Missing method comparison: %s',result_file);
    v=whos('-file',f); names={v.name};
    must={'disp_gt','b_data','displacement_estimation_autocorr','shear_params'};
    assert(all(ismember(must,names)), ...
        'Dataset for %s is missing one of disp_gt, b_data, displacement_estimation_autocorr, shear_params.',name);
    wanted=[must,{'autocorr_params','array_length','scene_depth','prb_theta','prf'}];
    wanted=wanted(ismember(wanted,names));
    D=load(f,wanted{:});
    M=load(result_file,'GT','PG','AC','roi','x_axis','z_axis', ...
        'PRF','target_frequency','filter_speed');
    need={'GT','PG','AC','roi','x_axis','z_axis', ...
        'PRF','target_frequency','filter_speed'};
    assert(all(isfield(M,need)),'Missing method-comparison fields for %s.',name);
    assert(all(isfield(M.roi,{'image_mask','inclusion_core','background_core'})), ...
        'Missing original ROI masks for %s.',name);
end

function [uGT,uUS,GT,x,z,rois,meta]=build_inputs(D,M,name)
    xfull=double(M.x_axis(:));zfull=double(M.z_axis(:));
    [Nz,Nx]=size(M.GT);
    assert(isequal(size(M.PG),[Nz Nx],size(M.AC),[Nz Nx]) && ...
        Nx==numel(xfull) && Nz==numel(zfull), ...
        'Comparison grid mismatch for %s.',name);
    assert(max(abs(double(D.b_data.scan.x_axis(:))-xfull))<1e-10 && ...
        max(abs(double(D.b_data.scan.z_axis(:))-zfull))<1e-10, ...
        'B-mode scan and saved SWS axes differ for %s.',name);
    PRF=double(M.PRF);f0=double(M.target_frequency);
    assert(PRF>2*f0 && f0>0,'Invalid timing/frequency settings for %s.',name);
    assert(abs(f0-double(D.shear_params.source_freq))<1e-6, ...
        'Source frequency mismatch for %s.',name);
    assert(abs(double(M.filter_speed)-3.2)<1e-6, ...
        'Original reconstruction did not use a 3.2 m/s radial parameter.');
    if isfield(D,'prf') && ~isempty(D.prf)
        assert(abs(double(D.prf)-PRF)<1e-6,'PRF mismatch for %s.',name);
    end
    if isfield(D,'prb_theta') && ~isempty(D.prb_theta)
        probe_deg=rad2deg(double(D.prb_theta(1)));
        assert(abs(probe_deg-90)<0.2, ...
            'Probe orientation is not 90 degrees for %s.',name);
    else
        warning('No probe angle metadata for %s; verify simulation configuration.',name);
        probe_deg=NaN;
    end
    c_bg=double(D.shear_params.c_shear_bkg);
    c_inc=double(D.shear_params.c_shear_incl);
    assert(abs(c_inc-c_bg)>0.05,'No material contrast in %s.',name);
    if isfield(D,'autocorr_params') && ...
            isfield(D.autocorr_params,'packet_size')
        packet=double(D.autocorr_params.packet_size);
    else
        packet=4;
        warning('%s missing packet size; using legacy 4.',name);
    end
    assert(isscalar(packet) && packet>=2 && packet==fix(packet), ...
        'Invalid packet size.');
    navg=packet-1;
    [ngz,ngx,ngt]=size(D.disp_gt);
    nEst=numel(D.displacement_estimation_autocorr.data)/(Nz*Nx);
    assert(nEst>=2 && nEst==fix(nEst) && ngt>=nEst+navg-1, ...
        'Temporal dimension mismatch for %s.',name);
    [kg,~]=create_phantom(D.shear_params);  % grid only, no simulation
    dx=double(kg.dx);dz=double(kg.dz);clear kg;
    if isfield(D,'array_length') && ~isempty(D.array_length)
        L=double(D.array_length);
    else
        L=(128-1)*0.3e-3;
        warning('%s: legacy array_length fallback used.',name);
    end
    if isfield(D,'scene_depth') && ~isempty(D.scene_depth)
        depth=double(D.scene_depth);
    else
        depth=50e-3;
        warning('%s: legacy scene_depth fallback used.',name);
    end
    xgt=(-L/2:dx:L/2).';
    zgt=(5e-3:dz:depth).';
    assert(ngz==numel(zgt) && ngx==numel(xgt), ...
        'Original GT displacement axes mismatch for %s.',name);
    ix=find(xfull>=xgt(1)-1e-9 & xfull<=xgt(end)+1e-9);
    iz=find(zfull>=zgt(1)-1e-9 & zfull<=zgt(end)+1e-9);
    assert(numel(ix)>25 && numel(iz)>25,'Insufficient GT/US scan overlap.');
    x=xfull(ix);z=zfull(iz);
    [XX,ZZ]=meshgrid(min(max(x,xgt(1)),xgt(end)), ...
        min(max(z,zgt(1)),zgt(end)));
    GT=double(M.GT(iz,ix));
    rois.global=logical(M.roi.image_mask(iz,ix));
    rois.incl=rois.global & logical(M.roi.inclusion_core(iz,ix));
    rois.bg=rois.global & logical(M.roi.background_core(iz,ix));
    assert(nnz(rois.incl)>50 && nnz(rois.bg)>100, ...
        'Insufficient inclusion/background ROIs for %s.',name);
    uUS=-single(reshape(D.displacement_estimation_autocorr.data, ...
        Nz,Nx,nEst));
    uUS=uUS(iz,ix,:);
    uGT=zeros(numel(z),numel(x),nEst,'single');
    for t=1:nEst
        frame=mean(double(D.disp_gt(:,:,t:t+navg-1)),3);
        interp=interp2(xgt,zgt,frame,XX,ZZ,'linear',NaN);
        assert(all(isfinite(interp(:))), ...
            'Invalid GT interpolation at t=%d for %s.',t,name);
        uGT(:,:,t)=single(interp);
    end
    assert(all(isfinite(uUS(:))),'Nonfinite US displacements for %s.',name);
    meta=struct('f0',f0,'PRF',PRF,'c_bg',c_bg,'c_inc',c_inc, ...
        'probe_deg',probe_deg,'packet_size',packet,'source',name);
end

function [pg,ac]=run_estimators(u,x,z,PRF,f0,angles,sigma,dtheta)
    nz=numel(z);nx=numel(x);nt=size(u,3);
    elasto=uff.elastography();
    elasto.PRF=PRF;
    elasto.target_frequency=f0;
    elasto.scan=uff.linear_scan('x_axis',x,'z_axis',z);
    elasto.autocorrelation_field=uff.data_field();
    elasto.autocorrelation_field.set_component('u', ...
        reshape(u,nz,nx,1,nt));
    elasto.directional_filtering_angles=angles;
    extract=postprocess.elastography.extract_freq_component(elasto);
    extract.target_frequency=f0;
    extract.component='u';
    extract.run();
    elasto.radial_filter=uff.spatial_filter('gaussian_charite','radial');
    elasto.radial_filter.sigma=sigma;
    elasto.angular_filter=uff.spatial_filter('gaussian_charite','angular');
    elasto.angular_filter.deltaTheta=dtheta;
    filt=postprocess.elastography.apply_filter(elasto);
    filt.run();
    p=postprocess.elastography.phase_gradient_method(elasto);
    p.run();
    pg=squeeze(double(elasto.sws_map.u));
    a=postprocess.elastography.autocorrelation_SWS_method(elasto);
    a.run();
    ac=squeeze(double(elasto.sws_map.u));
    assert(isequal(size(pg),[nz nx],size(ac)), ...
        'Estimator output dimensions do not match ultrasound scan.');
end

function plotmap(ax,x,z,im,outline,ttl,lims)
    imagesc(ax,x*1e3,z*1e3,im);
    axis(ax,'image');set(ax,'YDir','normal');
    colormap(ax,'parula');colorbar(ax);clim(ax,lims);
    xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
    title(ax,ttl,'Interpreter','none');
    if ~isempty(outline)
        hold(ax,'on');
        contour(ax,x*1e3,z*1e3,double(outline),[.5 .5], ...
            'w','LineWidth',1.0);
    end
end
