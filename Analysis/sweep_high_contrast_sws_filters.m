function T = sweep_high_contrast_sws_filters()
%SWEEP_HIGH_CONTRAST_SWS_FILTERS  One-factor-at-a-time USTB filter study.
%
%   T = sweep_high_contrast_sws_filters();
%
% Purpose: Find whether adjusting the *existing* radial and angular
% directional-filter parameters improves SWS recovery for high_contrast.
% This is a controlled DEVELOPMENT screening, NOT a new SWS estimator.
%
% Input: saved k-Wave GT displacement from high_contrast, 90-deg probe.
% Exactly the same GT packet averaging, 200-Hz extraction, scan crop,
% 12 angles, PG/AC classes, masks, and SWS evaluation as the validated
% high_contrast_filter_on_off diagnostic. No new k-Wave or RF simulation.
%
% Test one factor at a time:
%  - Radial parameter: elasto.radial_filter.sigma = speed_parameter/(2*pi*f0)
%  - Angular parameter: elasto.angular_filter.deltaTheta
% Parameter values are USTB filter settings, NOT measured tissue speeds.
%
% Strict controls:
%  1. First reconstruction at original settings must agree with the prior
%     controlled FILTER-ON maps (if that result file exists).
%  2. All variants are evaluated using the SAME intersection of valid
%     pixels. In addition, valid fractions are reported per variant.
%  3. No existing data/results/scripts are overwritten by this diagnostic.
%  4. A candidate is only a candidate if BOTH background and global error
%     constraints are met and inclusion error improves. Validate on other
%     phantoms BEFORE changing a production pipeline.
%
% MATLAB and user's local USTB version not executable in this environment.

here = fileparts(mfilename('fullpath'));
if isfolder(fullfile(here,'Data')) && isfolder(fullfile(here,'Results'))
    repo = here;
else
    repo = fileparts(here);
end
assert(isfolder(fullfile(repo,'Data')) && isfolder(fullfile(repo,'Results')), ...
    'Place this script inside the repository Analysis/ directory.');
addpath(genpath(repo));
toolbox_dir = fullfile(fileparts(repo),'elastography','ElastographyToolbox');
assert(isfolder(toolbox_dir),'Cannot find ElastographyToolbox: %s',toolbox_dir);
addpath(genpath(toolbox_dir));

required_functions = {'create_phantom','uff.elastography','uff.data_field', ...
    'uff.linear_scan','postprocess.elastography.extract_freq_component', ...
    'postprocess.elastography.apply_filter', ...
    'postprocess.elastography.phase_gradient_method', ...
    'postprocess.elastography.autocorrelation_SWS_method'};
for j=1:numel(required_functions)
    assert(~isempty(which(required_functions{j})), ...
        'Missing MATLAB/USTB component: %s',required_functions{j});
end

filename='high_contrast_displacement_with_quality.mat';
locations={fullfile(repo,'Data','high_contrast',filename), ...
    fullfile(repo,'Data',filename),fullfile(repo,filename)};
found=find(cellfun(@isfile,locations),1);
assert(~isempty(found),'Cannot locate high_contrast input MAT.');
datafile=locations{found};
priorfile=fullfile(repo,'Results','method_comparison','batch', ...
    'high_contrast','method_comparison.mat');
assert(isfile(priorfile),'Missing original method_comparison.mat: %s',priorfile);
info=whos('-file',datafile); available={info.name};
needed={'disp_gt','b_data','displacement_estimation_autocorr','shear_params'};
assert(all(ismember(needed,available)),'Dataset lacks GT displacement or metadata.');
toload={'disp_gt','b_data','displacement_estimation_autocorr', ...
    'shear_params','autocorr_params','array_length','scene_depth','prb_theta','prf'};
toload=toload(ismember(toload,available));
D=load(datafile,toload{:});
M=load(priorfile,'GT','PG','AC','roi','x_axis','z_axis', ...
    'PRF','target_frequency','filter_speed');
req={'GT','PG','AC','roi','x_axis','z_axis', ...
    'PRF','target_frequency','filter_speed'};
assert(all(isfield(M,req)),'Method comparison metadata incomplete.');
assert(all(isfield(M.roi,{'image_mask','inclusion_core','background_core'})), ...
    'Saved GT ROI masks missing.');

xfull=double(M.x_axis(:)); zfull=double(M.z_axis(:));
[Nz,Nx]=size(M.GT);
assert(isequal(size(M.PG),[Nz Nx],size(M.AC),[Nz Nx]) && ...
    numel(xfull)==Nx && numel(zfull)==Nz, ...
    'Saved SWS dimensions disagree.');
assert(max(abs(double(D.b_data.scan.x_axis(:))-xfull))<1e-10 && ...
    max(abs(double(D.b_data.scan.z_axis(:))-zfull))<1e-10, ...
    'Ultrasound scan coordinates differ from method comparison.');
assert(all(diff(xfull)>0) && all(diff(zfull)>0), ...
    'Scan axes must be strictly increasing.');

PRF=double(M.PRF);
f0=double(M.target_frequency);
speed_base=double(M.filter_speed);
assert(PRF>2*f0 && f0>0 && speed_base>0, 'Invalid timing/filter metadata.');
assert(abs(f0-double(D.shear_params.source_freq))<1e-6, ...
    'Metadata mismatch: source frequency.');
if isfield(D,'prf') && ~isempty(D.prf)
    assert(abs(double(D.prf)-PRF)<1e-6,'Saved PRF does not match.');
end
if isfield(D,'prb_theta') && ~isempty(D.prb_theta)
    fprintf('Probe angle (NOT filter angle): %s deg\n', ...
        mat2str(rad2deg(double(D.prb_theta))));
    assert(abs(rad2deg(double(D.prb_theta(1)))-90)<0.2, ...
        'Unexpected probe orientation. Stop before interpreting this test.');
end
c_bg=double(D.shear_params.c_shear_bkg);
c_inc=double(D.shear_params.c_shear_incl);
assert(c_inc>c_bg,'Expected fast inclusion.');
if isfield(D,'autocorr_params') && isfield(D.autocorr_params,'packet_size')
    packet=double(D.autocorr_params.packet_size);
else
    packet=4;
    warning('packet_size missing; using legacy 4.');
end
assert(isscalar(packet) && packet>=2 && packet==fix(packet), ...
    'Unexpected packet_size.');
navg=packet-1;
[ngz,ngx,ngt]=size(D.disp_gt);
nt=numel(D.displacement_estimation_autocorr.data)/(Nz*Nx);
assert(nt==fix(nt) && ngt>=nt+navg-1, 'Displacement frame count mismatch.');

[kg,~]=create_phantom(D.shear_params); % Only grid/material; NO simulation
DX=double(kg.dx); DZ=double(kg.dz);
clear kg;
if isfield(D,'array_length') && ~isempty(D.array_length)
    array_length=double(D.array_length);
else
    array_length=(128-1)*0.3e-3;
    warning('Missing array_length: assuming 128 x 0.3 mm legacy probe.');
end
if isfield(D,'scene_depth') && ~isempty(D.scene_depth)
    scene_depth=double(D.scene_depth);
else
    scene_depth=50e-3;
    warning('Missing scene_depth: assuming 50 mm.');
end
xgt=(-array_length/2:DX:array_length/2).';
zgt=(5e-3:DZ:scene_depth).';
assert(numel(xgt)==ngx && numel(zgt)==ngz, ...
    'Native GT axes do not agree with disp_gt dimensions.');
ix=find(xfull>=xgt(1)-1e-9 & xfull<=xgt(end)+1e-9);
iz=find(zfull>=zgt(1)-1e-9 & zfull<=zgt(end)+1e-9);
assert(numel(ix)>25 && numel(iz)>25,'Not enough overlap between scans.');
x=xfull(ix); z=zfull(iz);
nz=numel(z); nx=numel(x);
[XX,ZZ]=meshgrid(min(max(x,xgt(1)),xgt(end)), ...
                 min(max(z,zgt(1)),zgt(end)));
GT=double(M.GT(iz,ix));
roi=logical(M.roi.image_mask(iz,ix));
inc=roi & logical(M.roi.inclusion_core(iz,ix));
bg=roi & logical(M.roi.background_core(iz,ix));
assert(nnz(inc)>100 && nnz(bg)>100,'Missing GT inclusion/background pixels.');
outline=(GT>(c_bg+c_inc)/2);

% Build SAME input as validated GT-vs-US/on-off control.
fprintf('\n=== HIGH CONTRAST: RADIAL / ANGULAR FILTER SWEEP ===\n');
fprintf('Input: %s\n',datafile);
fprintf('Material: background %.3f, inclusion %.3f m/s\n',c_bg,c_inc);
fprintf('Native %dx%dx%d (%.3f mm); scan %dx%d; %d frames\n', ...
    ngz,ngx,ngt,DX*1e3,nz,nx,nt);
fprintf('f0 %.1f Hz, PRF %.1f Hz, radial parameter base %.3f m/s\n', ...
    f0,PRF,speed_base);

u=zeros(nz,nx,nt,'single');
for t=1:nt
    fr=mean(double(D.disp_gt(:,:,t:t+navg-1)),3);
    q=interp2(xgt,zgt,fr,XX,ZZ,'linear',NaN);
    assert(all(isfinite(q(:))), 'Interpolation failed at frame %d.',t);
    u(:,:,t)=single(q);
end
clear D XX ZZ fr q;

nang=12;
angles=2*pi*linspace(0,1-1/nang,nang);
base_theta=median(diff(angles));
omega=2*pi*f0;

% Fixed, pre-specified ONE-PARAMETER sweeps. Baseline is run FIRST.
% RadialParameter_mps is an input to sigma=parameter/omega, NOT tissue c.
Mode=["Baseline";repmat("Radial",4,1);repmat("Angular",4,1)];
RadialParameter_mps=[speed_base; 1.6; 2.4; 4.8; 6.4; ...
    repmat(speed_base,4,1)];
AngularWidth_deg=[rad2deg(base_theta)*ones(5,1);15;45;60;90];
nrun=numel(Mode);
PG=cell(nrun,1);AC=cell(nrun,1);

fprintf('Nine configurations: baseline + 4 radial + 4 angular.\n');
% Reference was already validated against the original filtered pipeline.
prioron=fullfile(repo,'Results','sws_estimator_diagnostic', ...
    'high_contrast_filter_on_off','filter_on_off_results.mat');
assert(isfile(prioron), ...
    'Run compare_high_contrast_filter_on_off() first (reference missing).');
R=load(prioron,'PG_on','AC_on','x','z');
dpg=NaN; dac=NaN;
for j=1:nrun
    fprintf('(%d/%d) %-8s  radial %.2f m/s, angular %.1f deg ...\n', ...
        j,nrun,char(Mode(j)),RadialParameter_mps(j),AngularWidth_deg(j));
    timer=tic;
    [PG{j},AC{j}]=run_one( ...
        u,x,z,PRF,f0,angles, ...
        RadialParameter_mps(j)/omega,deg2rad(AngularWidth_deg(j)));
    fprintf('    finished (%.1f s)\n',toc(timer));
    if j==1
        assert(isequal(size(PG{1}),size(R.PG_on),size(R.AC_on)) && ...
            numel(R.x)==nx && numel(R.z)==nz && ...
            max(abs(double(R.x(:))-x))<1e-10 && ...
            max(abs(double(R.z(:))-z))<1e-10, ...
            'Reference baseline grid mismatch.');
        dpg=maxdiff(PG{1},double(R.PG_on));
        dac=maxdiff(AC{1},double(R.AC_on));
        fprintf('Baseline vs prior filter ON: PG %.4g, AC %.4g m/s\n',dpg,dac);
        assert(dpg<0.02 && dac<0.02, ...
            'Baseline output changed. STOP: not comparable with prior validated test.');
        clear R;
    end
end
clear u;

% Common pixels throughout sweep and both estimators. Separate per-run
% validity is reported so invalid pixels cannot silently look better.
validall=isfinite(GT);
for j=1:nrun
    validall=validall & isfinite(PG{j}) & PG{j}>0 & PG{j}<10 & ...
        isfinite(AC{j}) & AC{j}>0 & AC{j}<10;
end
common=roi & validall;
corei=inc & common;
coreb=bg & common;
fraction_global=nnz(common)/nnz(roi);
fraction_i=nnz(corei)/nnz(inc);
fraction_b=nnz(coreb)/nnz(bg);
fprintf('Common coverage: global %.2f%%, incl %.2f%%, BG %.2f%%\n', ...
    100*fraction_global,100*fraction_i,100*fraction_b);
assert(fraction_global>0.90 && fraction_i>0.90 && fraction_b>0.90, ...
    ['Too many invalid pixels in one of the filter variants. ' ...
     'Reduce parameter range before comparing results.']);

nrow=2*nrun;
ModeAll=strings(nrow,1);
Estimator=strings(nrow,1);
RadialSetting_mps=zeros(nrow,1);
AngularSetting_deg=zeros(nrow,1);
GlobalRMSE_mps=zeros(nrow,1);
InclusionRMSE_mps=zeros(nrow,1);
BackgroundRMSE_mps=zeros(nrow,1);
InclusionMean_mps=zeros(nrow,1);
BackgroundMean_mps=zeros(nrow,1);
ContrastRecovery_pct=zeros(nrow,1);
GlobalValid_pct=zeros(nrow,1);
InclusionValid_pct=zeros(nrow,1);
BackgroundValid_pct=zeros(nrow,1);
pass=zeros(nrow,1);
for j=1:nrun
    for k=1:2
        p=2*j-2+k;
        if k==1, v=PG{j}; name="PG"; else, v=AC{j}; name="AC"; end
        ModeAll(p)=Mode(j);
        Estimator(p)=name;
        RadialSetting_mps(p)=RadialParameter_mps(j);
        AngularSetting_deg(p)=AngularWidth_deg(j);
        GlobalRMSE_mps(p)=sqrt(mean((v(common)-GT(common)).^2));
        InclusionRMSE_mps(p)=sqrt(mean((v(corei)-GT(corei)).^2));
        BackgroundRMSE_mps(p)=sqrt(mean((v(coreb)-GT(coreb)).^2));
        InclusionMean_mps(p)=mean(v(corei));
        BackgroundMean_mps(p)=mean(v(coreb));
        ContrastRecovery_pct(p)=100*(InclusionMean_mps(p)-BackgroundMean_mps(p)) / ...
            (mean(GT(corei))-mean(GT(coreb)));
        finitev=isfinite(v) & v>0 & v<10;
        GlobalValid_pct(p)=100*nnz(finitev & roi)/nnz(roi);
        InclusionValid_pct(p)=100*nnz(finitev & inc)/nnz(inc);
        BackgroundValid_pct(p)=100*nnz(finitev & bg)/nnz(bg);
    end
end

% Predefined screening criteria. Corrected parameter settings must improve
% inclusion RMSE by >= 0.20 m/s while permitting at most +0.10 m/s BG RMSE
% and +0.02 m/s global RMSE relative to baseline of the SAME estimator.
% Screening on this ONE phantom is not evidence of generalization.
for p=3:nrow
    base=1;
    if Estimator(p)=="AC",base=2;end
    pass(p)= ...
        InclusionRMSE_mps(p)<=InclusionRMSE_mps(base)-0.20 && ...
        BackgroundRMSE_mps(p)<=BackgroundRMSE_mps(base)+0.10 && ...
        GlobalRMSE_mps(p)<=GlobalRMSE_mps(base)+0.02 && ...
        InclusionValid_pct(p)>=95 && BackgroundValid_pct(p)>=95;
end
ScreenPass=logical(pass);
T=table(ModeAll,Estimator,RadialSetting_mps,AngularSetting_deg, ...
    GlobalRMSE_mps,InclusionRMSE_mps,BackgroundRMSE_mps, ...
    InclusionMean_mps,BackgroundMean_mps,ContrastRecovery_pct, ...
    GlobalValid_pct,InclusionValid_pct,BackgroundValid_pct,ScreenPass);

fprintf('\n=== SAME-PIXEL FILTER SWEEP RESULTS ===\n');
disp(T(:,{'ModeAll','Estimator','RadialSetting_mps','AngularSetting_deg', ...
    'GlobalRMSE_mps','InclusionRMSE_mps', ...
    'BackgroundRMSE_mps','ContrastRecovery_pct','ScreenPass'}));

% Best provisional per estimator ONLY among candidates meeting pre-set
% criteria. Otherwise baseline remains the default.
selected=[1 2];
for k=1:2
    if k==1, kind="PG"; else, kind="AC"; end
    ii=find(ScreenPass & Estimator==kind);
    if ~isempty(ii)
        [~,best]=min(GlobalRMSE_mps(ii));
        selected(k)=ii(best);
        fprintf('%s: provisional candidate row %d, %s (radial %.2f, angular %.1f).\n', ...
            char(Estimator(selected(k))),selected(k),char(ModeAll(selected(k))), ...
            RadialSetting_mps(selected(k)),AngularSetting_deg(selected(k)));
    else
        fprintf('%s: NO candidate meets all criteria; KEEP original for now.\n', ...
            char(Estimator(k)));
    end
end

outdir=fullfile(repo,'Results','sws_estimator_diagnostic', ...
    'high_contrast_filter_sweep');
if ~isfolder(outdir),mkdir(outdir);end
writetable(T,fullfile(outdir,'filter_sweep_summary.csv'));

% Figures: one line per estimator for radial/angular one-factor sweeps.
fig=figure('Color','w','Position',[80 80 1550 900],'Visible','off');
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
for k=1:2
    inds=k:2:nrow;
    ix=inds(:).';
    makechart(nexttile(tl,(k-1)*3+1),T,ix, ...
        'InclusionRMSE_mps','Inclusion RMSE [m/s]');
    makechart(nexttile(tl,(k-1)*3+2),T,ix, ...
        'BackgroundRMSE_mps','Background RMSE [m/s]');
    makechart(nexttile(tl,(k-1)*3+3),T,ix, ...
        'GlobalRMSE_mps','Global RMSE [m/s]');
end
title(tl,'High contrast | radial and angular USTB filter parameter sweep');
exportgraphics(fig,fullfile(outdir,'filter_sweep_metrics.png'),'Resolution',170);
close(fig);

fig2=figure('Color','w','Position',[40 80 1550 840],'Visible','off');
tl=tiledlayout(fig2,2,3,'TileSpacing','compact','Padding','compact');
showmap(nexttile(tl,1),x,z,GT,outline,[0 8],'Material GT');
showmap(nexttile(tl,2),x,z,PG{1},outline,[0 8],'PG baseline');
showmap(nexttile(tl,3),x,z,PG{ceil(selected(1)/2)},outline,[0 8], ...
    sprintf('PG screened choice: %s',char(ModeAll(selected(1)))));
showmap(nexttile(tl,4),x,z,double(common),[],[0 1],'Fixed common ROI');
showmap(nexttile(tl,5),x,z,AC{1},outline,[0 8],'AC baseline');
showmap(nexttile(tl,6),x,z,AC{ceil(selected(2)/2)},outline,[0 8], ...
    sprintf('AC screened choice: %s',char(ModeAll(selected(2)))));
title(tl,'High contrast | screened setting only if RMSE constraints pass');
exportgraphics(fig2,fullfile(outdir,'filter_sweep_maps.png'),'Resolution',170);
close(fig2);

meta=struct('datafile',datafile,'referencefile',prioron, ...
    'Nangles',nang,'f0_Hz',f0,'PRF_Hz',PRF,'GT_grid_m',[DX DZ], ...
    'baselineRadialSpeed_mps',speed_base, ...
    'baselineAngularWidth_deg',rad2deg(base_theta), ...
    'baselineErrorPG_mps',dpg,'baselineErrorAC_mps',dac, ...
    'commonGlobalFraction',fraction_global, ...
    'commonInclusionFraction',fraction_i, ...
    'commonBackgroundFraction',fraction_b, ...
    'screening',['Inclusion RMSE improvement >= 0.20 m/s; ' ...
    'background RMSE increase <= 0.10 m/s; ' ...
    'global RMSE increase <= 0.02 m/s; valid ROI >=95%.'], ...
    'warning',['One-phantom screening is exploratory. Do not change ' ...
    'production pipeline before validation on other phantoms.']);
save(fullfile(outdir,'filter_sweep_results.mat'), ...
    'T','PG','AC','GT','x','z','common','corei','coreb', ...
    'selected','meta','-v7.3');
fprintf('Saved in: %s\n',outdir);
fprintf(['No original data or production code changed. The sweep is ' ...
    'DEVELOPMENT ONLY, NOT validation on independent phantoms.\n']);
end

function [pg,ac]=run_one(u,x,z,PRF,f0,angles,sigma,deltaTheta)
    nz=numel(z);nx=numel(x);nt=size(u,3);
    elasto=uff.elastography();
    elasto.PRF=PRF;
    elasto.target_frequency=f0;
    elasto.scan=uff.linear_scan('x_axis',x,'z_axis',z);
    elasto.autocorrelation_field=uff.data_field();
    elasto.autocorrelation_field.set_component('u',reshape(u,nz,nx,1,nt));
    elasto.directional_filtering_angles=angles;
    extract=postprocess.elastography.extract_freq_component(elasto);
    extract.target_frequency=f0;
    extract.component='u';
    extract.run();
    elasto.radial_filter=uff.spatial_filter('gaussian_charite','radial');
    elasto.radial_filter.sigma=sigma;
    elasto.angular_filter=uff.spatial_filter('gaussian_charite','angular');
    elasto.angular_filter.deltaTheta=deltaTheta;
    filt=postprocess.elastography.apply_filter(elasto);
    filt.run();
    pg_obj=postprocess.elastography.phase_gradient_method(elasto);
    pg_obj.run();
    pg=squeeze(double(elasto.sws_map.u));
    ac_obj=postprocess.elastography.autocorrelation_SWS_method(elasto);
    ac_obj.run();
    ac=squeeze(double(elasto.sws_map.u));
    assert(isequal(size(pg),[nz nx],size(ac)), ...
        'USTB estimator returned wrong scan dimensions.');
end

function mx=maxdiff(A,B)
    if ~isequal(size(A),size(B)),mx=Inf;return;end
    valid=isfinite(A); second=isfinite(B);
    if ~isequal(valid,second),mx=Inf;return;end
    if ~any(valid(:)),mx=Inf;return;end
    mx=max(abs(A(valid)-B(valid)));
end

function makechart(ax,T,indices,columnName,ylab)
    t=T(indices,:);
    y=t.(columnName);
    labs=cell(height(t),1);
    for j=1:height(t)
        if t.ModeAll(j)=="Baseline"
            labs{j}='Original';
        elseif t.ModeAll(j)=="Radial"
            labs{j}=sprintf('R %.1f',t.RadialSetting_mps(j));
        else
            labs{j}=sprintf('A %g deg',t.AngularSetting_deg(j));
        end
    end
    bar(ax,1:numel(y),y,'FaceColor',[.28 .53 .74]);
    hold(ax,'on');
    plot(ax,[.5 numel(y)+.5],[y(1) y(1)],'k--','LineWidth',1);
    passed=find(t.ScreenPass);
    if ~isempty(passed)
        plot(ax,passed,y(passed),'kp','MarkerSize',13, ...
            'MarkerFaceColor',[.92 .75 .22]);
    end
    set(ax,'XTick',1:numel(y),'XTickLabel',labs);
    xtickangle(ax,40);
    ylabel(ax,ylab);
    grid(ax,'on');
    title(ax,sprintf('%s | %s',char(t.Estimator(1)),ylab),'Interpreter','none');
end

function showmap(ax,x,z,M,outline,lims,ttl)
    imagesc(ax,x*1e3,z*1e3,M);
    axis(ax,'image');set(ax,'YDir','normal');
    colormap(ax,'parula'); colorbar(ax);clim(ax,lims);
    xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
    title(ax,ttl,'Interpreter','none');
    if ~isempty(outline)
        hold(ax,'on');
        contour(ax,x*1e3,z*1e3,double(outline),[.5 .5], ...
            'w','LineWidth',1.2);
    end
end
