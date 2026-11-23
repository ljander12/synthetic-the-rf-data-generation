%% RUN_ARTIFICIAL_SWS_VALIDATION.m
% One-click, controlled artificial SWS and displacement quality validation.
% Put this file in repository/Analysis/ and press RUN in MATLAB.
%
% NO k-Wave, Field II, RF synthesis, USTB or precomputed MAT files required.
% This is an analytical / kinematic laboratory, NOT a physical wave solver.
%
% Metrics evaluated:
%   GT-based: RMSE, MAE, bias, CNR, contrast recovery, gCNR,
%             edge accuracy / edge sharpness (synthetic-grid approximation)
%   Controlled displacement SNR: measured against known analytic CLEAN u
%             ONLY in cases where stochastic displacement noise is injected.
%             NaN for other cases (NOT an SWS-map SNR or clinical SNR).
%   GT-free:  Correlation Q (Schoen et al. 2025 mapping, zero-lag adaptation)
%             RMI-1D (Choi et al. 2015 Eq. 6 adaptation; fixed calibration)
%   GT-based validation OF Q/RMI: ROC AUC, false confidence, acceptance.
%
% Phase Planarity Q / Log-RMI intentionally EXCLUDED.
% IMPORTANT: Artificially modifying c_est without modifying displacement
% cannot influence Correlation Q. This is a deliberate negative control.
% RMI is sensitive to BOTH displacement and estimated speed.
%
% Existing project production functions / results are NOT overwritten.

%% USER SETTINGS
rng(42, 'twister');
showFigures = false;                % avoid opening 24 panels; summary opens
saveMAT = true;
rho_kg_m3 = 1000;                   % for approximate Young modulus
thresholdCorr = 0.63;               % illustrative pixel-wise only
thresholdRMI = 0.85;                % illustrative; NOT validated for these data
errorThreshold_mps = 0.5;
frequency_Hz = 200;
PRF_Hz = 3000;
Nt = 75;

% Adjust case list to run a smaller experiment.
caseNames = { ...
    'homogeneous_ideal', ...        % correct speed, ideal plane wave
    'inclusion_ideal', ...          % correct speed, kinematic inclusion wave
    'high_contrast_ideal', ...      % correct speed, 6 m/s inclusion
    'global_speed_bias', ...       % change ONLY reconstructed SWS
    'localized_speed_bias', ...    % change ONLY SWS in a patch
    'random_speed_error', ...      % change ONLY SWS pixelwise
    'blurred_inclusion', ...       % change ONLY SWS at boundary
    'lateral_speed_drift', ...     % change ONLY SWS smoothly
    'displacement_noise', ...      % corrupt BOTH displacement and SWS
    'decorrelated_patch', ...      % corrupt BOTH in a localized patch
    'amplitude_attenuation', ...   % lower waveform amplitude, correct SWS
    'coherent_wrong_speed' ...     % wrong SWS AND a coherent matching wave
    };

%% PATHS AND GRID
here = fileparts(mfilename('fullpath'));
if isfolder(fullfile(here,'Results'))
    repo = here;
elseif isfolder(fullfile(fileparts(here),'Results'))
    repo = fileparts(here);
else
    repo = here;  % standalone use outside a Git repo
end
if isfolder(fullfile(repo,'Analysis'))
    addpath(fullfile(repo,'Analysis')); % find optional project helpers
end
outDir = fullfile(repo,'Results','artificial_sws_validation');
if ~isfolder(outDir), mkdir(outDir); end

% Nearly isotropic spatial grid; x is shear-wave propagation direction.
x = linspace(-0.0165,0.0165,111);
z = linspace(0.008,0.046,128).';
t = (0:Nt-1)/PRF_Hz;
[X,Z] = meshgrid(x,z);
[Nz,Nx] = size(X);
dx = mean(diff(x));

GT_hom = 2.32*ones(Nz,Nx);
inclusion = (X/0.0065).^2 + ((Z-0.023)/0.0065).^2 <= 1;
GT_inc = GT_hom;
GT_inc(inclusion) = 4.67;
GT_high = GT_hom;
GT_high(inclusion) = 6.00;

% Kinematic (travel-time) wavefields: these are NOT k-Wave solutions.
uHom = make_kinematic_wave(GT_hom,x,z,t,frequency_Hz);
uInc = make_kinematic_wave(GT_inc,x,z,t,frequency_Hz);
uHigh = make_kinematic_wave(GT_high,x,z,t,frequency_Hz);

% Fixed RMI normalization calibrated once WITHOUT using SWS-error labels.
% The calibration wavefield is noisy on purpose so numerical limits span a
% useful range. Unlike the main study, this calibration is SYNTHETIC.
ucal = uInc + 0.09*randn(size(uInc));
[rescal, ampcal] = wave_eq_residual_1d(ucal,GT_inc,x,PRF_Hz);
rv = rescal(isfinite(rescal));
av = ampcal(isfinite(ampcal));
assert(numel(rv)>100 && numel(av)>100,'Synthetic calibration failed.');
limits = struct('res_min',percentile_simple(rv,5), ...
                'res_max',percentile_simple(rv,95), ...
                'u_min',percentile_simple(av,5), ...
                'u_max',percentile_simple(av,95));
assert(limits.res_max>limits.res_min && limits.u_max>limits.u_min, ...
    'RMI calibration bounds collapsed.');
save(fullfile(outDir,'synthetic_rmi_calibration.mat'),'limits', ...
    'thresholdCorr','thresholdRMI','frequency_Hz','PRF_Hz');

fprintf('\n=== CONTROLLED ARTIFICIAL SWS VALIDATION ===\n');
fprintf('Cases: %d | Grid: %d x %d x %d | f=%.0f Hz\n', ...
    numel(caseNames),Nz,Nx,Nt,frequency_Hz);
fprintf('Q: zero-lag adapted Schoen (2025). RMI: 1D Choi Eq. 6.\n');
fprintf('RMI calibration frozen across all cases; NOT main-study calibration.\n');
fprintf('GT used only for scoring, never to compute Q/RMI.\n\n');

%% MAIN EXPERIMENTS
summary = table();
for ic = 1:numel(caseNames)
    name = caseNames{ic};
    [gt,c_est,u,category,description] = make_case(name, ...
        GT_hom,GT_inc,GT_high,uHom,uInc,uHigh, ...
        X,Z,x,z,t,frequency_Hz,inclusion);
    assert(isequal(size(gt),size(c_est),[Nz Nx]),'SWS map size mismatch.');
    assert(isequal(size(u),[Nz Nx Nt]),'Wavefield size mismatch.');
    assert(all(isfinite(c_est(:))) && all(c_est(:)>0), ...
        'Artificial SWS must be positive and finite.');

    % === GT-INDEPENDENT quality calculations ===
    Qcorr = correlation_q_zero_lag(u);
    [residual,urep] = wave_eq_residual_1d(u,c_est,x,PRF_Hz);
    RMI = make_rmi_1d(residual,urep,limits);

    % === Ground truth used ONLY in evaluation below ===
    spatialMask = abs(X)<0.014 & Z>0.011 & Z<0.043;
    valid = spatialMask & isfinite(Qcorr) & isfinite(RMI) & ...
        isfinite(c_est) & isfinite(gt);
    assert(nnz(valid)>100, 'Insufficient common pixels in %s.',name);
    absErr = abs(c_est-gt);

    % Noise SNR uses the known clean analytical displacement field, NOT GT
    % SWS. It is applicable only when we explicitly injected displacement
    % noise into ui in make_case(). SWS-only corruption is NOT signal noise.
    injectedDisplacementSNR_dB = NaN;
    if strcmp(name,'displacement_noise') || strcmp(name,'decorrelated_patch')
        injectedDisplacementSNR_dB = displacement_snr_from_clean( ...
            u, uInc, spatialMask);
    end

    isBad = absErr(valid)>errorThreshold_mps;
    qc = Qcorr(valid); qr = RMI(valid);
    fracBad=mean(isBad);

    % Metrics used elsewhere in the project. Any optional unavailable
    % dependency is returned as NaN, never silently replaced by a new metric.
    m = sws_reference_metrics(c_est,gt,inclusion,valid,x,z);
    corrAUC = quality_auc(qc,~isBad);
    rmiAUC = quality_auc(qr,~isBad);
    corrFC = false_confidence_top_half(qc,isBad);
    rmiFC = false_confidence_top_half(qr,isBad);
    corrAccepted = qc>=thresholdCorr;
    rmiAccepted = qr>=thresholdRMI;
    corrBadAccepted = accepted_bad_rate(isBad,corrAccepted);
    rmiBadAccepted = accepted_bad_rate(isBad,rmiAccepted);

    row = table(string(name),string(category),nnz(valid),fracBad, ...
        m.RMSE,m.MAE,m.Bias,m.InclRMSE,m.BkgRMSE, ...
        m.CNR,injectedDisplacementSNR_dB, ...
        m.gCNR,m.gCNR_hist,m.ContrastRecoveryPct, ...
        m.EdgeAccuracy_mm,m.EdgeSharpness_mpsPerMm, ...
        mean(qc),mean(qr),mean(qc>=1-1e-10), ...
        mean(qr>=1-1e-10),corrAUC,rmiAUC,corrFC,rmiFC, ...
        mean(corrAccepted),mean(rmiAccepted),corrBadAccepted,rmiBadAccepted, ...
        'VariableNames',{'Case','Category','N','BadFraction', ...
        'GlobalRMSE_mps','GlobalMAE_mps','GlobalBias_mps', ...
        'InclusionRMSE_mps','BackgroundRMSE_mps', ...
        'CNR','InjectedDisplacementSNR_dB','gCNR_tKDE','gCNR_hist_exploratory', ...
        'ContrastRecovery_percent','EdgeAccuracy_mm', ...
        'EdgeSharpness_mps_per_mm','CorrelationQmean','RMImean', ...
        'CorrelationQsaturation','RMIsaturation', ...
        'CorrelationAUC','RMIAUC','CorrelationFalseConfidenceTop50', ...
        'RMIFalseConfidenceTop50','CorrelationCoverage', 'RMICoverage', ...
        'CorrelationWrongAmongAccepted','RMIWrongAmongAccepted'});
    summary = [summary;row]; %#ok<AGROW>

    make_case_figures(outDir,name,description,x,z,gt,c_est,absErr, ...
        Qcorr,RMI,valid,rho_kg_m3,thresholdCorr,thresholdRMI,showFigures);
    if saveMAT
        save(fullfile(outDir,[name '_data.mat']), ...
            'gt','c_est','Qcorr','RMI','residual','urep','valid', ...
            'x','z','category','description','-v7.3');
    end
    fprintf(['%-25s RMSE %.3f | CNR %.3f | u-SNR %.2f dB | ' ...
        'Edge accuracy %.2f mm | Edge sharpness %.2f m/s/mm | ' ...
        'Q AUC %.3f | RMI AUC %.3f\n'], ...
        name,m.RMSE,m.CNR,injectedDisplacementSNR_dB, ...
        m.EdgeAccuracy_mm,m.EdgeSharpness_mpsPerMm,corrAUC,rmiAUC);
end

%% EXPORT SUMMARY AND COMPARISON CHARTS
writetable(summary,fullfile(outDir,'artificial_sws_metrics.csv'));
save(fullfile(outDir,'artificial_sws_summary.mat'),'summary','caseNames','limits');
make_summary_figure(summary,outDir,true); % always open ONE summary window
make_standard_metrics_figure(summary,outDir,false);
fprintf('\nResults saved to:\n  %s\n',outDir);
fprintf(['NOTE: AUC=NaN for a case containing only good or only bad SWS.\n' ...
    'Correlation Q cannot detect SWS-only corruptions by design.\n' ...
    'InjectedDisplacementSNR_dB: NaN except for explicit waveform-noise cases.\n' ...
    'Edge accuracy and edge sharpness apply only when an inclusion exists.\n' ...
    'gCNR_tKDE is NaN if its source is unavailable; histogram gCNR is exploratory.\n' ...
    'These are CONTROLLED synthetic tests, NOT physical simulation proof.\n']);

%% ------------------------------------------------------------------------
function snr_dB = displacement_snr_from_clean(u_noisy,u_clean,spatialMask)
% Controlled injected-displacement SNR, not reconstruction SNR.
% Uses KNOWN analytical signal (available only in synthetic experiments).
% Sum over the same central spatial ROI and over all time frames.
    assert(isequal(size(u_noisy),size(u_clean)), ...
        'Clean/noisy displacement dimensions differ.');
    assert(isequal(size(spatialMask),size(u_noisy(:,:,1))), ...
        'Displacement SNR ROI dimensions differ.');
    roi3 = repmat(logical(spatialMask),1,1,size(u_noisy,3));
    reference = double(u_clean(roi3));
    noise = double(u_noisy(roi3)) - reference;
    ok = isfinite(reference) & isfinite(noise);
    reference = reference(ok); noise = noise(ok);
    if isempty(reference),snr_dB=NaN;return;end
    ps = mean(reference.^2); pn = mean(noise.^2);
    if ~(ps>0) || ~(pn>0)
        snr_dB=NaN;return;
    end
    snr_dB=10*log10(ps/pn);
end

%% ------------------------------------------------------------------------
function [gt,s,u,category,description] = make_case(name,hom,inc,high,uh,ui,uhigh,X,Z,x,z,t,f,incl)
    gt=inc; s=inc; u=ui;
    category='SWS_only';
    description='Synthetic SWS-map error with fixed displacement';
    switch name
        case 'homogeneous_ideal'
            gt=hom;s=hom;u=uh;category='reference';
            description='Ideal homogeneous traveling wave; no imposed SWS error';
        case 'inclusion_ideal'
            category='reference';
            description='Ideal inclusion SWS; approximate kinematic travel-time wavefield';
        case 'high_contrast_ideal'
            gt=high;s=high;u=uhigh;category='reference';
            description='6 m/s inclusion with kinematic wavefield (not wave scattering)';
        case 'global_speed_bias'
            s=inc+0.75;
            description='Uniform +0.75 m/s SWS error; unchanged displacement';
        case 'localized_speed_bias'
            region = X>0.004 & X<0.012 & Z>0.028 & Z<0.039;
            s(region)=s(region)+1.3;
            description='Localized +1.3 m/s error, displacement unchanged';
        case 'random_speed_error'
            s=max(0.8,inc+0.6*randn(size(inc)));
            description='Pixelwise artificial SWS noise; displacement unchanged';
        case 'blurred_inclusion'
            s=gaussian_blur_replicate(inc,7);
            description='Blurred SWS boundary; displacement unchanged';
        case 'lateral_speed_drift'
            s=max(0.8,inc+0.85*X/max(abs(x)));
            description='Smooth lateral SWS bias; displacement unchanged';
        case 'displacement_noise'
            category='coupled_perturbation';
            A = repmat(max(abs(ui),[],3),1,1,size(ui,3));
            u = ui + 0.90*A.*randn(size(ui));
            s=max(0.8,inc+0.65*randn(size(inc)));
            description='Independent noise imposed on displacement and SWS';
        case 'decorrelated_patch'
            category='coupled_perturbation';
            region = X>0.002 & X<0.010 & Z>0.029 & Z<0.038;
            s(region)=s(region)+1.4;
            A = repmat(max(abs(ui),[],3),1,1,size(ui,3));
            mask=repmat(region,1,1,numel(t));
            jitter=2.3*A.*randn(size(ui));
            u(mask)=ui(mask)+jitter(mask);
            description='Spatially localized speed error + heavily decorrelated waveforms';
        case 'amplitude_attenuation'
            category='wavefield_only';
            attenuation=0.05+0.95*(1-(Z-min(z))/(max(z)-min(z)));
            u=bsxfun(@times,ui,attenuation);
            description='Depth attenuation of wave amplitude; correct SWS';
        case 'coherent_wrong_speed'
            category='coherent_false_confidence';
            gt=hom;s=3.20*ones(size(hom));
            u=make_kinematic_wave(s,x,z,t,f);
            description='Coherent plane wave fits WRONG estimated SWS, not material GT';
        otherwise
            error('Unknown artificial case: %s',name);
    end
    assert(any(incl(:)),'Inclusion mask missing.');
end

function u = make_kinematic_wave(c,x,z,t,f)
% The wave phase is constructed by lateral line integration of slowness.
% This is an analytical test field, NOT a physical heterogeneous solver.
    [Nz,Nx]=size(c);
    dx=diff(x);
    tau=zeros(Nz,Nx);
    for j=2:Nx
        tau(:,j)=tau(:,j-1)+0.5*dx(j-1)* ...
            (1./c(:,j-1)+1./c(:,j));
    end
    [X,Z]=meshgrid(x,z);
    amplitude=(1+0.10*cos(2*pi*X/(max(x)-min(x)))) .* ...
        exp(-0.15*(Z-min(z))/(max(z)-min(z)));
    u=zeros(Nz,Nx,numel(t),'single');
    omega=2*pi*f;
    for it=1:numel(t)
        u(:,:,it)=single(amplitude.*cos(omega*(t(it)-tau)));
    end
end

function Q = correlation_q_zero_lag(u)
% Match the project's test_sws_measurement_quality.m adaptation:
% zero-lag temporally centered Pearson correlation between adjacent x.
% Map pairwise dB correlation through Schoen Jr. et al. (2025) Q function.
    [Nz,Nx,Nt]=size(u);
    sa=zeros(Nz,Nx-1);sb=sa;saa=sa;sbb=sa;sab=sa;
    for it=1:Nt
        a=double(u(:,1:end-1,it));b=double(u(:,2:end,it));
        sa=sa+a;sb=sb+b;saa=saa+a.^2;sbb=sbb+b.^2;sab=sab+a.*b;
    end
    num=sab-sa.*sb/Nt;
    va=max(0,saa-sa.^2/Nt);vb=max(0,sbb-sb.^2/Nt);
    den=sqrt(va.*vb);
    rho=num./den;
    rho(den<=eps(max(den(:))))=NaN;
    rho=max(-1,min(1,rho));
    rmsPair=sqrt((saa+sbb)/(2*Nt));
    v=rmsPair(isfinite(rmsPair)&rmsPair>0);
    if ~isempty(v),rho(rmsPair<median(v)*1e-4)=NaN;end
    C=20*log10(max(abs(rho),realmin));
    pair=nan(size(C));
    pair(C>=-10)=1;
    mid=C>-50 & C<-10;
    pair(mid)=cos((pi/2)*(C(mid)+10)/40);
    pair(C<=-50)=0;
    Q=nan(Nz,Nx);
    Q(:,1)=pair(:,1);Q(:,end)=pair(:,end);
    Q(:,2:end-1)=min(pair(:,1:end-1),pair(:,2:end));
end

function [res,amp] = wave_eq_residual_1d(u,c,x,prf)
% Same Choi (2015) Eq. 6 finite differences as validate_sws_rmi.m.
    [Nz,Nx,Nt]=size(u);
    num=zeros(Nz,Nx);den=zeros(Nz,Nx);
    amp=max(abs(double(u(:,:,1))),abs(double(u(:,:,end))));
    dr=x(3:end)-x(2:end-1);
    dl=x(2:end-1)-x(1:end-2);
    span=x(3:end)-x(1:end-2);
    c2=c.^2;
    for it=2:Nt-1
        u0=double(u(:,:,it));
        up=double(u(:,:,it+1));um=double(u(:,:,it-1));
        utt=(up-2*u0+um)*prf^2;
        ux2=nan(Nz,Nx);
        forward=bsxfun(@rdivide,u0(:,3:end)-u0(:,2:end-1),dr);
        backward=bsxfun(@rdivide,u0(:,2:end-1)-u0(:,1:end-2),dl);
        ux2(:,2:end-1)=bsxfun(@rdivide,2*(forward-backward),span);
        rr=utt(:,2:end-1)-c2(:,2:end-1).*ux2(:,2:end-1);
        num(:,2:end-1)=num(:,2:end-1)+rr.^2;
        den(:,2:end-1)=den(:,2:end-1)+ux2(:,2:end-1).^2;
        amp=max(amp,abs(u0));
    end
    good=den(isfinite(den)&den>0);
    assert(~isempty(good),'Zero spatial curvature in all waveforms.');
    floorDen=median(good)*1e-12;
    res=num./den;
    res(:,[1,end])=NaN;
    res(~isfinite(den)|den<=floorDen)=NaN;
end

function Q = make_rmi_1d(res,amp,L)
    qr=(L.res_max-res)/(L.res_max-L.res_min);
    qa=(amp-L.u_min)/(L.u_max-L.u_min);
    qr=max(0,min(1,qr));qa=max(0,min(1,qa));
    Q=0.8*qr+0.2*qa;
    Q(~isfinite(res)|~isfinite(amp))=NaN;
end

function m = sws_reference_metrics(est,gt,incl,valid,x,z)
    e=double(est(valid)-gt(valid));
    m.RMSE=sqrt(mean(e.^2));
    m.MAE=mean(abs(e));m.Bias=mean(e);
    validIncl=valid & incl;validBg=valid & ~incl;
    m.InclRMSE=NaN;m.BkgRMSE=NaN;
    m.CNR=NaN;m.gCNR=NaN;m.gCNR_hist=NaN;m.ContrastRecoveryPct=NaN;
    m.EdgeAccuracy_mm=NaN;m.EdgeSharpness_mpsPerMm=NaN;
    if ~any(incl(:)) || max(gt(:))-min(gt(:))<1e-10
        % Homogeneous GT has NO real inclusion, even if a geometric
        % placeholder mask exists for the other synthetic cases.
        m.BkgRMSE=m.RMSE;
        return;
    end
    a=double(est(validIncl));b=double(est(validBg));
    ga=double(gt(validIncl));gb=double(gt(validBg));
    if isempty(a)||isempty(b),return;end
    m.InclRMSE=sqrt(mean((a-ga).^2));
    m.BkgRMSE=sqrt(mean((b-gb).^2));
    den=sqrt(var(a)+var(b));
    if den>0,m.CNR=abs(mean(a)-mean(b))/den;end
    trueContrast=mean(ga)-mean(gb);
    if abs(trueContrast)>eps
        m.ContrastRecoveryPct=100*(mean(a)-mean(b))/trueContrast;
    end
    % Always compute a SEPARATELY LABELLED exploratory histogram gCNR.
    edges=linspace(min([a(:);b(:)]),max([a(:);b(:)]),61);
    if edges(end)>edges(1)
        ca=histcounts(a,edges,'Normalization','probability');
        cb=histcounts(b,edges,'Normalization','probability');
        m.gCNR_hist=1-sum(min(ca,cb));
    end
    % Use project's actual gCNR_tKDE estimator if available.
    if exist('gCNR_tKDE','file')==2
        try
            m.gCNR=gCNR_tKDE(a,b,false);
        catch ME
            warning('gCNR_tKDE failed (%s); gCNR left NaN.',ME.message);
        end
    end
    % Approx. Dong-style edge metrics on the nearly isotropic synthetic grid.
    if exist('edge','file')~=2 || exist('bwperim','file')~=2 || ...
            exist('bwdist','file')~=2
        return;
    end
    gtEdge=bwperim(incl,8);
    if ~any(gtEdge(:)),return;end
    emin=min(est(:));emax=max(est(:));
    if emax<=emin,return;end
    estNorm=(est-emin)/(emax-emin);
    detected=edge(estNorm); % MATLAB default detector
    [rr,cc]=find(incl);
    margin=max(3,round(0.5*max(max(rr)-min(rr)+1,max(cc)-min(cc)+1)));
    box=false(size(incl));
    box(max(1,min(rr)-margin):min(size(incl,1),max(rr)+margin), ...
        max(1,min(cc)-margin):min(size(incl,2),max(cc)+margin))=true;
    detected=detected & box;
    if ~any(detected(:)),return;end
    r0=min(median(abs(diff(x))),median(abs(diff(z))));
    dist=bwdist(gtEdge);
    m.EdgeAccuracy_mm=mean(dist(detected))*r0*1e3;
    [gx,gz]=gradient(est);
    edgeGrad=hypot(gx,gz);
    m.EdgeSharpness_mpsPerMm=mean(edgeGrad(detected))/(r0*1e3);
end

function B = gaussian_blur_replicate(A,sigmaPixels)
% Spatial blur, without periodic/zero boundary effects.
    half=ceil(3*sigmaPixels);
    j=-half:half;
    g=exp(-(j.^2)/(2*sigmaPixels^2));g=g/sum(g);
    [nrow,ncol]=size(A);
    tmp=[repmat(A(:,1),1,half),A,repmat(A(:,end),1,half)];
    tmp=conv2(tmp,g,'valid');
    tmp=[repmat(tmp(1,:),half,1);tmp;repmat(tmp(end,:),half,1)];
    B=conv2(tmp,g','valid');
    assert(isequal(size(B),[nrow ncol]),'Blurred map has wrong size.');
end

function auc=quality_auc(q,isGood)
    q=double(q(:));good=logical(isGood(:));
    ok=isfinite(q);q=q(ok);good=good(ok);
    ng=nnz(good);nb=nnz(~good);
    if ng==0 || nb==0,auc=NaN;return;end
    [v,idx]=sort(q,'ascend');
    y=double(good(idx));n=numel(v);
    starts=[1;find(diff(v)~=0)+1];
    ends=[starts(2:end)-1;n];
    counts=cumsum([0;y]);
    groupGood=counts(ends+1)-counts(starts);
    avgRank=(starts+ends)/2;
    auc=(sum(groupGood.*avgRank)-ng*(ng+1)/2)/(ng*nb);
    auc=max(0,min(1,auc));
end

function fc=false_confidence_top_half(q,isBad)
    q=double(q(:));isBad=double(isBad(:));
    if isempty(q)||~all(isfinite(q)),fc=NaN;return;end
    n=numel(q);v=sort(q,'descend');cut=v(ceil(n/2));
    above=q>cut;tie=q==cut;
    weight=(n/2-nnz(above))/nnz(tie);
    weights=double(above)+max(0,min(1,weight))*double(tie);
    fc=sum(weights.*isBad)/sum(weights);
end

function x=accepted_bad_rate(isBad,accepted)
    if any(accepted),x=mean(isBad(accepted));else,x=NaN;end
end

function q=percentile_simple(values,p)
    v=sort(double(values(isfinite(values))));
    assert(numel(v)>1,'Not enough values for percentile.');
    t=1+(numel(v)-1)*p/100;
    lo=floor(t);hi=ceil(t);
    q=v(lo)+(t-lo)*(v(hi)-v(lo));
end

function make_case_figures(outDir,name,description,x,z,gt,est,err,Q,RMI,valid,rho,corrTh,rmiTh,showFigures)
    visibility='off';if showFigures,visibility='on';end
    f=figure('Visible',visibility,'Color','w', ...
        'Position',[50 80 1650 930],'Name',['Artificial SWS: ' name]);
    tl=tiledlayout(f,2,3,'Padding','compact','TileSpacing','compact');
    vmax=max([6.5;gt(:);est(:)]);
    plot_panel(nexttile(tl),x,z,gt,[0 vmax],'Known material SWS [m/s]');
    plot_panel(nexttile(tl),x,z,est,[0 vmax],'Artificial estimated SWS [m/s]');
    plot_panel(nexttile(tl),x,z,err,[0 max(2,max(err(:)))], ...
        'Absolute SWS error [m/s] (GT validation)');
    qc=Q;qc(~valid)=NaN;
    qr=RMI;qr(~valid)=NaN;
    plot_panel(nexttile(tl),x,z,qc,[0 1],'Correlation Q (no GT)');
    plot_panel(nexttile(tl),x,z,qr,[0 1],'RMI-1D (no GT)');
    highWrong=nan(size(err));
    highWrong(valid) = double(err(valid)>0.5 & RMI(valid)>=rmiTh);
    plot_panel(nexttile(tl),x,z,highWrong,[0 1], ...
        'Wrong SWS and RMI accepted (GT validation)');
    title(tl,[strrep(name,'_',' ') ' | ' description], 'Interpreter','none');
    exportgraphics(f,fullfile(outDir,[name '_quality.png']),'Resolution',150);
    if ~showFigures,close(f);end

    % All three stiffness panels use the SAME estimated SWS.
    stiffness=3*rho*est.^2/1000;
    stiffness(~valid)=NaN;
    maskedCorr=stiffness;maskedCorr(Q<corrTh)=NaN;
    maskedRMI=stiffness;maskedRMI(RMI<rmiTh)=NaN;
    f=figure('Visible',visibility,'Color','w', ...
        'Position',[70 90 1570 560],'Name',['Masked stiffness: ' name]);
    tl=tiledlayout(f,1,3,'Padding','compact','TileSpacing','compact');
    cmax=max(20,ceil(percentile_simple(stiffness,99)/10)*10);
    plot_panel(nexttile(tl),x,z,stiffness,[0 cmax],'Original Young modulus [kPa]');
    plot_panel(nexttile(tl),x,z,maskedCorr,[0 cmax], ...
        sprintf('Correlation Q >= %.2f',corrTh));
    plot_panel(nexttile(tl),x,z,maskedRMI,[0 cmax], ...
        sprintf('RMI-1D >= %.2f',rmiTh));
    title(tl,[strrep(name,'_',' ') ' | Same SWS estimate; quality masking only'], ...
        'Interpreter','none');
    exportgraphics(f,fullfile(outDir,[name '_stiffness.png']),'Resolution',150);
    if ~showFigures,close(f);end
end

function plot_panel(ax,x,z,map,limits,label)
    h=imagesc(ax,x*1e3,z*1e3,map);
    set(h,'AlphaData',double(isfinite(map)));
    set(ax,'Color',[0.92 0.92 0.92],'YDir','normal');
    axis(ax,'image');
    clim(ax,limits);colorbar(ax);
    xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
    title(ax,label,'Interpreter','none');
end

function make_summary_figure(T,outDir,showFigures)
    visibility='off';if showFigures,visibility='on';end
    n=height(T);
    f=figure('Visible',visibility,'Color','w','Position',[50 60 1700 850]);
    tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
    ax=nexttile(tl);
    bar(ax,T.GlobalRMSE_mps);ylabel(ax,'Global RMSE [m/s]');grid(ax,'on');
    title(ax,'GT-based reconstruction error');
    ax=nexttile(tl);
    bar(ax,[T.CorrelationAUC T.RMIAUC]);ylim(ax,[0 1]);grid(ax,'on');
    yline(ax,0.5,'k--');ylabel(ax,'AUC (high Q predicts correct SWS)');
    title(ax,'Quality-score discrimination');
    legend(ax,{'Correlation Q','RMI-1D'},'Location','best');
    ax=nexttile(tl);
    bar(ax,100*[T.CorrelationCoverage T.RMICoverage]);ylim(ax,[0 100]);
    ylabel(ax,'Accepted pixels [%]');grid(ax,'on');
    title(ax,'Illustrative acceptance thresholds');
    ax=nexttile(tl);
    bar(ax,100*[T.CorrelationFalseConfidenceTop50 T.RMIFalseConfidenceTop50]);
    ylim(ax,[0 100]);grid(ax,'on');
    ylabel(ax,'Wrong SWS among top 50% Q [%]');
    title(ax,'False confidence (GT validation)');
    for k=1:4
        ax=nexttile(tl,k);
        set(ax,'XTick',1:n,'XTickLabel',strrep(cellstr(T.Case),'_',' '), ...
            'XTickLabelRotation',35);
    end
    exportgraphics(f,fullfile(outDir,'all_cases_summary.png'),'Resolution',170);
    if ~showFigures,close(f);end
end


function make_standard_metrics_figure(T,outDir,showFigures)
% Additional four-metric comparison requested for the thesis.  Each
% series has its own units and applicability; never normalize scores into
% a misleading shared scale.  NaN values are intentionally left missing.
    visibility='off';if showFigures,visibility='on';end
    f=figure('Visible',visibility,'Color','w', ...
        'Position',[60 70 1720 880], ...
        'Name','Artificial validation: CNR, displacement SNR and edges');
    tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
    vals = {T.CNR, T.InjectedDisplacementSNR_dB, ...
            T.EdgeAccuracy_mm, T.EdgeSharpness_mps_per_mm};
    headings = {'CNR: inclusion versus background', ...
                'Injected displacement SNR (known clean reference)', ...
                'Edge accuracy (lower is better)', ...
                'Edge sharpness (larger = sharper detected boundary)'};
    ylabs = {'CNR [dimensionless]', 'SNR [dB]', ...
             'Edge accuracy [mm]', 'Edge sharpness [m/s/mm]'};
    for k=1:4
        ax=nexttile(tl);
        v=double(vals{k});
        bar(ax,v);grid(ax,'on');
        title(ax,headings{k}); ylabel(ax,ylabs{k});
        set(ax,'XTick',1:height(T), ...
            'XTickLabel',strrep(cellstr(T.Case),'_',' '), ...
            'XTickLabelRotation',40);
        xlim(ax,[0.4 height(T)+0.6]);
        if k==2
            % NaN outside injected-displacement-noise cases is intentional.
            text(ax,0.02,0.96,'Only cases with injected displacement noise', ...
                'Units','normalized','VerticalAlignment','top','FontSize',9);
        end
    end
    title(tl,'Controlled synthetic tests | GT-based spatial metrics and injected-displacement SNR');
    exportgraphics(f,fullfile(outDir,'cnr_snr_edge_metrics_comparison.png'), ...
        'Resolution',170);
    if ~showFigures,close(f);end
end
