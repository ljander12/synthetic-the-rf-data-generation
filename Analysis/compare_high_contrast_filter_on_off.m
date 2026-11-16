function T = compare_high_contrast_filter_on_off()
%COMPARE_HIGH_CONTRAST_FILTER_ON_OFF  Controlled filter bypass (high_contrast).
%
%  T = compare_high_contrast_filter_on_off();
%
% Uses the SAVED k-Wave GT displacement. Does NOT run k-Wave or Field II.
% Both branches use the exact same spatial scan, time samples, harmonic
% extraction, USTB PG/AC classes, and shared evaluation ROI.
%
% FILTER ON: the established 12-angle directional-filter pipeline.
% FILTER OFF (diagnostic bypass): the SAME unfiltered 200-Hz complex field
% is placed in each of the 12 directional-input slots, retaining the same
% estimator classes and direction-angle metadata. This bypasses directional
% separation; it is NOT a physically equivalent directional decomposition.
% A favorable OFF result therefore proves sensitivity to filtering, but
% DOES NOT by itself justify adopting this bypass as a production method.
%
% IMPORTANT: We explicitly check that replacing the directional field with
% its OWN original values reproduces the original PG and AC maps. If that
% check fails, we STOP rather than silently misinterpret the comparison.
%
% Requires an existing dataset and saved method comparison for high_contrast.
% Saves ONLY to Results/sws_estimator_diagnostic/high_contrast_filter_on_off/
%
% Diagnostic script, untested against the user's local MATLAB USTB branch.

here = fileparts(mfilename('fullpath'));
if isfolder(fullfile(here,'Data')) && isfolder(fullfile(here,'Results'))
    repo = here;
else
    repo = fileparts(here);
end
assert(isfolder(fullfile(repo,'Data')) && isfolder(fullfile(repo,'Results')), ...
    'Put the script in the repository Analysis/ folder.');
addpath(genpath(repo));
toolbox_dir = fullfile(fileparts(repo),'elastography','ElastographyToolbox');
assert(isfolder(toolbox_dir),'ElastographyToolbox not found: %s',toolbox_dir);
addpath(genpath(toolbox_dir));
needed = {'create_phantom','uff.elastography','uff.data_field', ...
    'uff.linear_scan','postprocess.elastography.extract_freq_component', ...
    'postprocess.elastography.apply_filter', ...
    'postprocess.elastography.phase_gradient_method', ...
    'postprocess.elastography.autocorrelation_SWS_method'};
for j=1:numel(needed)
    assert(~isempty(which(needed{j})),'Missing required function: %s',needed{j});
end

name='high_contrast_displacement_with_quality.mat';
opts={fullfile(repo,'Data','high_contrast',name), ...
      fullfile(repo,'Data',name),fullfile(repo,name),fullfile(here,name)};
pos=find(cellfun(@isfile,opts),1,'first');
assert(~isempty(pos),'No %s found.',name);
datafile=opts{pos};
priorfile=fullfile(repo,'Results','method_comparison','batch', ...
    'high_contrast','method_comparison.mat');
assert(isfile(priorfile),'Missing original PG/AC result: %s',priorfile);
info=whos('-file',datafile); names={info.name};
required={'disp_gt','b_data','displacement_estimation_autocorr','shear_params'};
assert(all(ismember(required,names)), ...
    'Data file must include disp_gt, b_data, displacement_estimation_autocorr, shear_params.');
wanted={'disp_gt','b_data','displacement_estimation_autocorr', ...
    'shear_params','autocorr_params','array_length','scene_depth','prb_theta','prf'};
wanted=wanted(ismember(wanted,names));
D=load(datafile,wanted{:});
M=load(priorfile,'GT','PG','AC','roi','x_axis','z_axis', ...
    'PRF','target_frequency','filter_speed');
required={'GT','PG','AC','roi','x_axis','z_axis', ...
    'PRF','target_frequency','filter_speed'};
assert(all(isfield(M,required)),'Prior method_comparison.mat lacks required fields.');
assert(all(isfield(M.roi,{'image_mask','inclusion_core','background_core'})), ...
    'Saved ROI does not contain expected masks.');

xfull=double(M.x_axis(:)); zfull=double(M.z_axis(:));
[Nz,Nx]=size(M.GT);
assert(numel(xfull)==Nx && numel(zfull)==Nz && ...
    isequal(size(M.PG),[Nz Nx],size(M.AC),[Nz Nx]), ...
    'Saved SWS maps differ in size.');
assert(all(diff(xfull)>0) && all(diff(zfull)>0), ...
    'SWS axes not increasing.');
assert(max(abs(double(D.b_data.scan.x_axis(:))-xfull))<1e-10 && ...
    max(abs(double(D.b_data.scan.z_axis(:))-zfull))<1e-10, ...
    'Saved US displacement grid differs from saved SWS grid.');
PRF=double(M.PRF); f0=double(M.target_frequency);
filterSpeed=double(M.filter_speed);
assert(PRF>2*f0 && f0>0 && filterSpeed>0,'Invalid saved PRF/f0/filter speed.');
assert(abs(f0-double(D.shear_params.source_freq))<1e-6, ...
    'Frequency metadata mismatch.');
if isfield(D,'prf') && ~isempty(D.prf)
    assert(abs(double(D.prf)-PRF)<1e-6,'Original PRF differs.');
end
c_bg=double(D.shear_params.c_shear_bkg);
c_inc=double(D.shear_params.c_shear_incl);
assert(c_inc>c_bg,'This diagnostic expects a fast inclusion.');
if isfield(D,'prb_theta')
    fprintf('Probe angle: %s degrees\n',mat2str(rad2deg(double(D.prb_theta))));
end
if isfield(D,'autocorr_params') && isfield(D.autocorr_params,'packet_size')
    packet=double(D.autocorr_params.packet_size);
else
    packet=4;
    warning('autocorr_params.packet_size absent; using previous diagnostic value 4.');
end
assert(packet>=2 && packet==round(packet),'Invalid packet_size.');
navg=packet-1;
[ngz,ngx,ngt]=size(D.disp_gt);
N_est=numel(D.displacement_estimation_autocorr.data)/(Nz*Nx);
assert(N_est==round(N_est) && ngt>=N_est+navg-1, ...
    'Number of GT wavefield frames does not match US-displacement frames.');

% Recover axes of the original 1-mm k-Wave grid. No new simulation.
[kg,~]=create_phantom(D.shear_params);
dx=double(kg.dx); dz=double(kg.dz); clear kg;
if isfield(D,'array_length') && ~isempty(D.array_length)
    array_length=double(D.array_length);
else
    array_length=(128-1)*0.3e-3;
    warning('Using legacy array_length fallback; verify original grid.');
end
if isfield(D,'scene_depth') && ~isempty(D.scene_depth)
    scene_depth=double(D.scene_depth);
else
    scene_depth=50e-3;
    warning('Using legacy scene_depth fallback; verify original grid.');
end
xgt=(-array_length/2:dx:array_length/2).';
zgt=(5e-3:dz:scene_depth).';
assert(numel(xgt)==ngx && numel(zgt)==ngz, ...
    'Recovered native k-Wave axes do not match disp_gt.');
ix=find(xfull>=xgt(1)-1e-9 & xfull<=xgt(end)+1e-9);
iz=find(zfull>=zgt(1)-1e-9 & zfull<=zgt(end)+1e-9);
assert(numel(ix)>25 && numel(iz)>25,'Insufficient grid overlap.');
x=xfull(ix);z=zfull(iz);
[XX,ZZ]=meshgrid(min(max(x,xgt(1)),xgt(end)), ...
                 min(max(z,zgt(1)),zgt(end)));
nz=numel(z);nx=numel(x);
GT=double(M.GT(iz,ix));
roi=logical(M.roi.image_mask(iz,ix));
inc=roi & logical(M.roi.inclusion_core(iz,ix));
bg=roi & logical(M.roi.background_core(iz,ix));
assert(nnz(inc)>100 && nnz(bg)>100, ...
    'Insufficient pixels in inclusion/background ROI.');
outline = GT > 0.5*(c_bg+c_inc);

fprintf('\n=== HIGH CONTRAST: FILTER ON/OFF — SAME k-Wave INPUT ===\n');
fprintf('Saved source: %s\n',datafile);
fprintf('GT: inclusion %.3f m/s, background %.3f m/s\n',c_inc,c_bg);
fprintf('Native GT %dx%dx%d (dx %.3f mm); cropped %dx%d, Nframes=%d\n', ...
    ngz,ngx,ngt,dx*1e3,nz,nx,N_est);
fprintf('f0 %.1f Hz, PRF %.1f Hz, filter speed %.3f m/s, packet=%d\n', ...
    f0,PRF,filterSpeed,packet);

u=zeros(nz,nx,N_est,'single');
for t=1:N_est
    native_frame=mean(double(D.disp_gt(:,:,t:t+navg-1)),3);
    sampled=interp2(xgt,zgt,native_frame,XX,ZZ,'linear',NaN);
    assert(all(isfinite(sampled(:))), ...
        'GT spatial interpolation failed at frame %d.',t);
    u(:,:,t)=single(sampled);
end
clear D XX ZZ sampled native_frame;

elasto=uff.elastography();
elasto.PRF=PRF;
elasto.target_frequency=f0;
elasto.scan=uff.linear_scan('x_axis',x,'z_axis',z);
elasto.autocorrelation_field=uff.data_field();
elasto.autocorrelation_field.set_component('u',reshape(u,nz,nx,1,N_est));
clear u;
nang=12;
angles=2*pi*linspace(0,1-1/nang,nang);
elasto.directional_filtering_angles=angles;
extractor=postprocess.elastography.extract_freq_component(elasto);
extractor.target_frequency=f0;
extractor.component='u';
extractor.run();
raw=fetch_component(elasto.unfiltered_phase_image_xz,nz,nx, ...
    'unfiltered_phase_image_xz');
assert(size(raw,3)==1, 'Expected one unfiltered 200-Hz phasor.');
raw=raw(:,:,1);

omega=2*pi*f0;
elasto.radial_filter=uff.spatial_filter('gaussian_charite','radial');
elasto.radial_filter.sigma=filterSpeed/omega;
elasto.angular_filter=uff.spatial_filter('gaussian_charite','angular');
elasto.angular_filter.deltaTheta=median(diff(angles));
filt=postprocess.elastography.apply_filter(elasto);
filt.run();
filtered_shape=field_shape(elasto.directional_filtered_field_xz);
filtered=fetch_component(elasto.directional_filtered_field_xz,nz,nx, ...
    'directional_filtered_field_xz');
assert(size(filtered,3)==nang, ...
    'Expected exactly %d directional filter channels, got %d.', ...
    nang,size(filtered,3));

fprintf('Running ORIGINAL PG/AC with directional filter ON...\n');
[PG_on,AC_on]=run_both(elasto,nz,nx);

% Verify the injection mechanism BEFORE bypassing the filter. Different
% USTB branches may use a read-only field, hidden cache, or alternate size.
fprintf('Control: re-injecting original filtered phasors...\n');
put_component(elasto,filtered,nz,nx,filtered_shape);
[PG_recheck,AC_recheck]=run_both(elasto,nz,nx);
check1=diff_finite(PG_on,PG_recheck);
check2=diff_finite(AC_on,AC_recheck);
fprintf('Filtered field re-injection max difference: PG %.4g, AC %.4g m/s\n', ...
    check1,check2);
assert(check1<1e-6 && check2<1e-6, ...
    ['USTB changed after re-injecting the unchanged filtered field. ' ...
     'STOP: bypass is not validated for this toolbox version.']);

% OFF: no directional filter. The same raw field is provided in each of
% the 12 direction slots because the current USTB PG/AC classes consume
% directional fields. The ANGLE METADATA and estimators remain the same.
% This is a controlled bypass, NOT a new physically decomposed wavefield.
fprintf('Running PG/AC with directional filter BYPASSED...\n');
raw12=repmat(raw,1,1,nang);
put_component(elasto,raw12,nz,nx,filtered_shape);
check_raw=fetch_component(elasto.directional_filtered_field_xz,nz,nx, ...
    'directional_filtered_field_xz');
field_err=max(abs(check_raw(:)-raw12(:)));
assert(isfinite(field_err) && field_err<1e-8*max(1,max(abs(raw(:)))), ...
    'Raw phasor injection failed to round-trip.');
clear raw12 check_raw filtered;
[PG_off,AC_off]=run_both(elasto,nz,nx);
clear elasto;

% Identical masks for all FOUR maps, independently of the branch ranking.
valid=isfinite(GT) & isfinite(PG_on) & isfinite(AC_on) & ...
    isfinite(PG_off) & isfinite(AC_off) & ...
    PG_on>0 & PG_on<10 & AC_on>0 & AC_on<10 & ...
    PG_off>0 & PG_off<10 & AC_off>0 & AC_off<10;
mask_g=roi & valid;
mask_i=inc & valid;
mask_b=bg & valid;
fprintf('Common valid pixels: all %d / %d; inclusion %d / %d; BG %d / %d\n', ...
    nnz(mask_g),nnz(roi),nnz(mask_i),nnz(inc),nnz(mask_b),nnz(bg));
assert(nnz(mask_i)>100 && nnz(mask_b)>100 && nnz(mask_g)>100, ...
    'Insufficient same-pixel coverage. Do NOT interpret branch comparison.');

maps={PG_on,AC_on,PG_off,AC_off};
Mode=["Filter ON";"Filter ON";"Filter OFF (bypass)";"Filter OFF (bypass)"];
Estimator=["PG";"AC";"PG";"AC"];
GlobalRMSE_mps=zeros(4,1); InclusionRMSE_mps=zeros(4,1);
BackgroundRMSE_mps=zeros(4,1); InclusionMean_mps=zeros(4,1);
BackgroundMean_mps=zeros(4,1); ContrastRecovery_pct=zeros(4,1);
InclusionBias_mps=zeros(4,1);
for j=1:4
    v=maps{j};
    GlobalRMSE_mps(j)=sqrt(mean((v(mask_g)-GT(mask_g)).^2));
    InclusionRMSE_mps(j)=sqrt(mean((v(mask_i)-GT(mask_i)).^2));
    BackgroundRMSE_mps(j)=sqrt(mean((v(mask_b)-GT(mask_b)).^2));
    InclusionMean_mps(j)=mean(v(mask_i));
    BackgroundMean_mps(j)=mean(v(mask_b));
    InclusionBias_mps(j)=mean(v(mask_i)-GT(mask_i));
    contrastGT=mean(GT(mask_i))-mean(GT(mask_b));
    ContrastRecovery_pct(j)=100*(InclusionMean_mps(j)-BackgroundMean_mps(j))/contrastGT;
end
T=table(Mode,Estimator,GlobalRMSE_mps,InclusionRMSE_mps, ...
    BackgroundRMSE_mps,InclusionBias_mps,InclusionMean_mps, ...
    BackgroundMean_mps,ContrastRecovery_pct);

fprintf('\n--- SAME-PIXEL FILTER ON/OFF RESULTS ---\n');
disp(T);
for j=1:2
    fprintf(['%s OFF minus ON: inclusion RMSE %+.4f m/s, background ' ...
        'RMSE %+.4f m/s, contrast recovery %+.2f percentage points\n'], ...
        char(Estimator(j)), ...
        InclusionRMSE_mps(j+2)-InclusionRMSE_mps(j), ...
        BackgroundRMSE_mps(j+2)-BackgroundRMSE_mps(j), ...
        ContrastRecovery_pct(j+2)-ContrastRecovery_pct(j));
end

% Cross-check filter ON with existing previous GT-input reconstruction.
prior=fullfile(repo,'Results','sws_estimator_diagnostic', ...
    'gt_vs_us_high_contrast','comparison_results.mat');
if isfile(prior)
    q=load(prior,'pg_gt','ac_gt','x','z');
    if all(isfield(q,{'pg_gt','ac_gt','x','z'})) && ...
            isequal(size(q.pg_gt),[nz nx]) && ...
            isequal(size(q.ac_gt),[nz nx]) && ...
            max(abs(double(q.x(:))-x))<1e-10 && ...
            max(abs(double(q.z(:))-z))<1e-10
        chk_pg=diff_finite(PG_on,double(q.pg_gt));
        chk_ac=diff_finite(AC_on,double(q.ac_gt));
        fprintf('ON vs prior GT-input PG/AC max difference: %.4g / %.4g m/s\n', ...
            chk_pg,chk_ac);
        assert(chk_pg<0.02 && chk_ac<0.02, ...
            'Filter-ON output differs from validated previous GT-input run.');
    else
        warning('Prior GT-vs-US diagnostic scan differs; baseline cross-check skipped.');
    end
else
    warning('Prior GT-vs-US result absent; only local reinjection check performed.');
end

outdir=fullfile(repo,'Results','sws_estimator_diagnostic', ...
    'high_contrast_filter_on_off');
if ~isfolder(outdir),mkdir(outdir);end
writetable(T,fullfile(outdir,'filter_on_off_summary.csv'));

fig=figure('Color','w','Position',[50 60 1700 900],'Visible','off');
tile=tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
showmap(nexttile(tile,1),x,z,GT,outline,[0 8], ...
    'Ground-truth SWS');
showmap(nexttile(tile,2),x,z,PG_on,outline,[0 8], ...
    'PG / filter ON');
showmap(nexttile(tile,3),x,z,PG_off,outline,[0 8], ...
    'PG / filter OFF (bypass)');
showmap(nexttile(tile,4),x,z,PG_off-PG_on,outline,[-4 4], ...
    'PG / OFF - ON');
showmap(nexttile(tile,5),x,z,double(mask_g),[],[0 1], ...
    'Identical valid ROI');
showmap(nexttile(tile,6),x,z,AC_on,outline,[0 8], ...
    'AC / filter ON');
showmap(nexttile(tile,7),x,z,AC_off,outline,[0 8], ...
    'AC / filter OFF (bypass)');
showmap(nexttile(tile,8),x,z,AC_off-AC_on,outline,[-4 4], ...
    'AC / OFF - ON');
title(tile,'high contrast | SAME k-Wave input | filter ON vs bypass', ...
    'Interpreter','none');
exportgraphics(fig,fullfile(outdir,'filter_on_off_maps.png'),'Resolution',165);
close(fig);

zs=find(any(outline,2));
zi=zs(round(numel(zs)/2));
fig2=figure('Color','w','Position',[100 110 1300 490],'Visible','off');
tl=tiledlayout(fig2,1,2,'TileSpacing','compact');
for j=1:2
    ax=nexttile(tl,j);
    plot(ax,x*1e3,GT(zi,:),'k-','LineWidth',2);hold(ax,'on');
    plot(ax,x*1e3,maps{j}(zi,:),'b-','LineWidth',1.5);
    plot(ax,x*1e3,maps{j+2}(zi,:),'r-','LineWidth',1.5);
    xlabel(ax,'Lateral [mm]');ylabel(ax,'SWS [m/s]');
    title(ax,sprintf('%s / z=%.2f mm',char(Estimator(j)),z(zi)*1e3));
    legend(ax,{'Material GT','Filter ON','Filter OFF (bypass)'}, ...
        'Location','best');
    grid(ax,'on');ylim(ax,[0 10]);
end
exportgraphics(fig2,fullfile(outdir,'filter_on_off_profile.png'), ...
    'Resolution',170);
close(fig2);

metadata=struct('datafile',datafile,'methodFile',priorfile, ...
    'PRF',PRF,'f0',f0,'filterSpeed',filterSpeed,'Nangles',nang, ...
    'comparisonGrid',[nz nx],'nativeGridSpacing',[dx dz], ...
    'GTpacketAveragingFrames',navg, ...
    'reinjectMaxErrorPG',check1,'reinjectMaxErrorAC',check2, ...
    'limitations',['OFF bypass repeats one unfiltered complex field over ' ...
    '12 directional-input channels: it does not reproduce a physically ' ...
    'direction-separated field. It is a method-sensitivity CONTROL only. ' ...
    'Native k-Wave 1 mm grid interpolated to denser ultrasound scan: ' ...
    'no new spatial information is added.']);
save(fullfile(outdir,'filter_on_off_results.mat'), ...
    'T','PG_on','PG_off','AC_on','AC_off','GT','x','z', ...
    'inc','bg','mask_g','mask_i','mask_b','metadata','-v7.3');

fprintf('\nSaved outputs to: %s\n',outdir);
fprintf(['CAUTION: OFF = diagnostic filter bypass with replicated raw fields; ' ...
    'not a production-ready estimator. Better contrast alone does not ' ...
    'prove that the original directional filter is wrong.\n']);
end

function [pg_map,ac_map]=run_both(elasto,nz,nx)
    pgop=postprocess.elastography.phase_gradient_method(elasto);
    pgop.run();
    pg_map=squeeze(double(elasto.sws_map.u));
    assert(isequal(size(pg_map),[nz nx]),'Unexpected PG map dimensions.');
    acop=postprocess.elastography.autocorrelation_SWS_method(elasto);
    acop.run();
    ac_map=squeeze(double(elasto.sws_map.u));
    assert(isequal(size(ac_map),[nz nx]),'Unexpected AC map dimensions.');
end

function put_component(elasto,A,nz,nx,original_shape)
    f=elasto.directional_filtered_field_xz;
    if ~(isobject(f) && ismethod(f,'set_component'))
        error(['USTB directional_filtered_field_xz is class %s and has no ' ...
            'set_component method. Cannot safely bypass filter. ' ...
            'Send output of class(elasto.directional_filtered_field_xz).'], ...
            class(f));
    end
    assert(size(A,1)==nz && size(A,2)==nx, ...
        'Incompatible bypass field dimensions.');
    assert(prod(original_shape)==numel(A), ...
        'Cannot preserve original USTB directional-field size.');
    f.set_component('u',reshape(A,original_shape));
end

function shp=field_shape(v)
    a=raw_component(v);
    assert(isnumeric(a) && ~isempty(a), ...
        'Cannot inspect USTB directional field dimensions.');
    shp=size(a);
end

function a=raw_component(v)
    a=[];
    if isnumeric(v)
        a=v;
    elseif isstruct(v)
        if isfield(v,'u'),a=v.u;elseif isfield(v,'data'),a=v.data;end
    elseif isobject(v)
        if ismethod(v,'get_component')
            try,a=v.get_component('u');catch,end
        end
        if isempty(a) && isprop(v,'u')
            try,a=v.u;catch,end
        end
        if isempty(a) && isprop(v,'data')
            try,a=v.data;catch,end
        end
    end
    if isstruct(a) && isfield(a,'u'),a=a.u;end
    if iscell(a) && numel(a)==1,a=a{1};end
end

function A=fetch_component(v,nz,nx,name)
    a=raw_component(v);
    assert(isnumeric(a) && ~isempty(a), ...
        'Cannot read %s from USTB class %s.',name,class(v));
    assert(size(a,1)==nz && size(a,2)==nx, ...
        'USTB field %s size mismatch.',name);
    A=reshape(double(a),nz,nx,[]);
    if isreal(A)
        finite=A(isfinite(A));
        if ~isempty(finite) && min(finite)>=-pi-.05 && ...
                max(finite)<=2*pi+.05
            A=exp(1i*A);
            warning('%s appears to contain phase angles; converted to phasors.',name);
        else
            error('USTB %s is real and not identifiable as a phasor.',name);
        end
    end
end

function md=diff_finite(A,B)
    if ~isequal(size(A),size(B)),md=Inf;return;end
    fa=isfinite(A);fb=isfinite(B);
    if ~isequal(fa,fb),md=Inf;return;end
    if any(fa(:)),md=max(abs(A(fa)-B(fa)));else,md=Inf;end
end

function showmap(ax,x,z,M,outline,lim,ttl)
    imagesc(ax,x*1e3,z*1e3,M);
    axis(ax,'image');set(ax,'YDir','normal');
    colorbar(ax);clim(ax,lim);
    xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
    title(ax,ttl,'Interpreter','none');
    if ~isempty(outline)
        hold(ax,'on');
        contour(ax,x*1e3,z*1e3,double(outline),[.5 .5], ...
            'w','LineWidth',1.2);
    end
end
