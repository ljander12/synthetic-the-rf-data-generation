function T = analyze_rmi_speed_dependence(phantom_case)
%ANALYZE_RMI_SPEED_DEPENDENCE Does c_est inside RMI mask SWS errors?
%
% CONTROL ANALYSIS OF EXISTING METRIC 3 (RMI-1D): NOT A NEW METRIC.
%
% Compare on identical phase-ramped displacement u_beta:
%   Adaptive: RMI(u_beta, c_est_beta)  [existing calculation]
%   Frozen:   RMI(u_beta, c_est_0)     [diagnostic control]
%
% At beta=0 both variants must agree pixel for pixel. Therefore:
%   Delta RMI(adaptive vs baseline)
%      = [RMI(frozen,beta) - RMI(baseline,0)]
%      + [RMI(adaptive,beta) - RMI(frozen,beta)]
% The first bracket isolates displacement changes at fixed c; the second
% isolates the effect of changing c for the same displacement field.
%
% RMI implementation: Choi et al., IEEE IUS (2015), Eq. 4-6 and 8.
% This is the existing 1D Eq. 6 ADAPTATION, not the full radial Eq. 7.
% The frozen-speed variant is a CONTROL ONLY, not a deployable metric.
% Ground truth is used ONLY to calculate validation RMSE/risk, never Q.
%
% INSTALL: put this file into repo/Analysis/ (do not modify prior scripts).
% PREREQUISITES: previously run validate_sws_rmi and phase-ramp scripts.
%
% RUN:
%   analyze_rmi_speed_dependence('large_inclusions')
%   analyze_rmi_speed_dependence('high_contrast')
%
% OUTPUT (Results/sws_measurement_quality/<case>/rmi_speed_diagnostic):
%   rmi_speed_summary.csv
%   rmi_speed_trends.png
%   rmi_speed_decomposition.png
%   rmi_speed_maps_pg.png, rmi_speed_maps_ac.png
%   rmi_speed_risk_coverage.png
%   rmi_speed_results.mat
%
    if nargin < 1 || isempty(phantom_case)
        phantom_case = 'large_inclusions';
    end
    phantom_case = char(string(phantom_case));
    root = fileparts(mfilename('fullpath'));
    if ~isfolder(fullfile(root,'Results')), root=fileparts(root); end
    top = fullfile(root,'Results','sws_measurement_quality');
    pdir = fullfile(top,phantom_case);
    source_file = fullfile(pdir,'phase_ramp_sws_validation','phase_ramp_metrics.mat');
    baseline_file = fullfile(pdir,'rmi_1d_diagnostic','rmi_1d_results.mat');
    calib_file = fullfile(top,'rmi_1d_calibration.mat');
    assert(isfile(source_file),'Missing phase-ramp reconstruction: %s',source_file);
    assert(isfile(baseline_file),'Missing original RMI: %s',baseline_file);
    assert(isfile(calib_file),'Missing fixed calibration: %s',calib_file);

    S = load(source_file,'maps','GT','common','betas','prf','f0','x','z');
    R = load(baseline_file,'RmiPG','RmiAC','basic','masks','regions','GT','x','z');
    K = load(calib_file,'limits');
    lim = K.limits;
    assert(strcmp(lim.algorithm,'Choi2015_Eq6_RMI_1D_adaptation'), ...
        'Calibration belongs to a different RMI algorithm.');
    assert(lim.res_max>lim.res_min && lim.u_max>lim.u_min, ...
        'Invalid fixed RMI calibration.');
    betas = double(S.betas(:));
    assert(numel(S.maps)==numel(betas) && nnz(betas==0)==1, ...
        'Expected phase-ramp maps and one beta=0 baseline.');
    k0 = find(betas==0,1);
    [Nz,Nx] = size(S.GT);
    x = double(S.x(:)); z=double(S.z(:));
    assert(numel(x)==Nx && numel(z)==Nz && ...
        all(diff(x)>0) && Nx>=5, 'Invalid scan axes.');
    assert(isequal(size(S.common),[Nz Nx],size(R.GT),size(R.basic)), ...
        'Saved results do not have matching dimensions.');
    assert(max(abs(double(S.GT(:))-double(R.GT(:))))<1e-6, ...
        'Ground-truth arrays differ between saved analyses.');
    assert(isequal(size(R.RmiPG),[Nz Nx],size(R.RmiAC)), ...
        'Baseline RMI dimensions mismatch.');
    assert(numel(R.regions)==numel(R.masks), 'Region mask metadata mismatch.');
    assert(max(abs(x-double(R.x(:))))<1e-9 && ...
        max(abs(z-double(R.z(:))))<1e-9, 'Saved x/z axes differ.');
    for k=1:numel(betas)
        assert(abs(double(S.maps(k).beta)-betas(k))<1e-12, ...
            'Stored beta and phase-ramp map indices differ.');
        assert(isequal(size(S.maps(k).PG),[Nz Nx], ...
            size(S.maps(k).AC),size(S.maps(k).Q)), ...
            'Stored PG/AC/phase-Q sizes differ.');
    end

    % Read the SAME displacement data and USTB objects as prior scripts.
    if isempty(which('uff.beamformed_data'))
        ustb_root=fullfile(fileparts(root),'elastography', ...
            'ElastographyToolbox','USTB');
        if isfolder(ustb_root), addpath(genpath(ustb_root)); rehash; end
    end
    assert(~isempty(which('uff.beamformed_data')), ...
        'USTB not on MATLAB path. Add ElastographyToolbox/USTB first.');
    data_name = [phantom_case '_displacement_with_quality.mat'];
    candidates = {fullfile(root,'Data',data_name), ...
        fullfile(root,'Data',phantom_case,data_name), ...
        fullfile(root,data_name),fullfile(root,'Analysis',data_name)};
    data_file='';
    for i=1:numel(candidates)
        if isfile(candidates{i}), data_file=candidates{i};break;end
    end
    assert(~isempty(data_file),'Cannot find %s',data_name);

    fprintf('\n=== RMI SPEED-COUPLING CONTROL: %s ===\n',phantom_case);
    fprintf('Using FIXED calibration from %s\n',lim.reference_phantom);
    fprintf('Displacement source: %s\n',data_file);
    D=load(data_file,'channel_data','displacement_estimation_autocorr', ...
        'shear_params');
    prf=double(D.channel_data.PRF);
    f0=double(D.shear_params.source_freq);
    assert(isscalar(prf) && prf>0 && isfinite(prf) && ...
        abs(prf-double(S.prf))<1e-8 && abs(f0-double(S.f0))<1e-8, ...
        'PRF or excitation frequency differs from prior phase-ramp test.');
    raw=D.displacement_estimation_autocorr.data;
    assert(isnumeric(raw) && isreal(raw) && ...
        mod(numel(raw),Nz*Nx)==0, 'Unexpected displacement dimensions.');
    Nt=numel(raw)/(Nz*Nx);
    assert(Nt>=5,'Too few displacement frames.');
    u=-reshape(single(raw),Nz,Nx,Nt);
    clear D raw;
    fprintf('Displacement %dx%dx%d; PRF %.1f Hz; f0 %.1f Hz\n', ...
        Nz,Nx,Nt,prf,f0);

    % Use exactly the harmonic least-squares perturbation convention of
    % validate_rmi_phase_ramp.m and analyze_rmi_components.m.
    tt=(0:Nt-1)'/prf;
    cos_t=cos(2*pi*f0*tt);
    sin_t=sin(2*pi*f0*tt);
    design=[cos_t sin_t ones(Nt,1)];
    assert(rank(design)==3,'Harmonic LS fit is rank deficient.');
    weights=pinv(design);
    a=zeros(Nz,Nx);b=zeros(Nz,Nx);
    for it=1:Nt
        frame=double(u(:,:,it));
        a=a+weights(1,it)*frame;
        b=b+weights(2,it)*frame;
    end
    [Xmm,~]=meshgrid(x*1e3,z*1e3);
    methods={'PG','AC'};
    c0=cell(1,2);
    for mi=1:2
        v=double(S.maps(k0).(methods{mi}));
        v(~isfinite(v)|v<=0|v>=10)=NaN;
        c0{mi}=v;
    end

    nBeta=numel(betas);
    adaptive=cell(nBeta,2); frozen=cell(nBeta,2);
    qresA=cell(nBeta,2);qresF=cell(nBeta,2);
    qamp=cell(nBeta,1);
    for k=1:nBeta
        beta=betas(k);
        theta=beta*(Xmm-mean(x)*1e3);
        da=a.*(cos(theta)-1)+b.*sin(theta);
        db=b.*(cos(theta)-1)-a.*sin(theta);
        cBeta=cell(1,2);
        for mi=1:2
            v=double(S.maps(k).(methods{mi}));
            v(~isfinite(v)|v<=0|v>=10)=NaN;
            cBeta{mi}=v;
        end
        fprintf('Computing adaptive/frozen RMI: beta=%+.2f rad/mm ...\n',beta);
        [resA,resF,amp]=wave_residual_both(u,da,db,cos_t,sin_t, ...
            x,prf,cBeta,c0);
        qamp{k}=clip01((amp-lim.u_min)/(lim.u_max-lim.u_min));
        for mi=1:2
            qresA{k,mi}=clip01((lim.res_max-resA{mi}) / ...
                (lim.res_max-lim.res_min));
            qresF{k,mi}=clip01((lim.res_max-resF{mi}) / ...
                (lim.res_max-lim.res_min));
            adaptive{k,mi}=0.8*qresA{k,mi}+0.2*qamp{k};
            frozen{k,mi}=0.8*qresF{k,mi}+0.2*qamp{k};
            adaptive{k,mi}(~isfinite(resA{mi})|~isfinite(amp))=NaN;
            frozen{k,mi}(~isfinite(resF{mi})|~isfinite(amp))=NaN;
        end
        clear da db resA resF amp;
    end
    clear u a b;

    % BASELINE GUARANTEES: adaptive/frozen must coincide exactly at beta=0.
    for mi=1:2
        equalmask=isfinite(adaptive{k0,mi}) & isfinite(frozen{k0,mi});
        assert(nnz(equalmask)>100,'Too few valid baseline RMI pixels.');
        diff0=max(abs(adaptive{k0,mi}(equalmask)-frozen{k0,mi}(equalmask)));
        assert(diff0<1e-10, 'Frozen c does not reproduce adaptive beta=0.');
        original=double(R.(['Rmi' methods{mi}]));
        check=equalmask & logical(R.basic) & isfinite(original);
        assert(nnz(check)>100,'No common pixels with saved baseline RMI.');
        delta=max(abs(adaptive{k0,mi}(check)-original(check)));
        fprintf('%s baseline vs saved RMI: max abs difference %.4g\n', ...
            methods{mi},delta);
        assert(delta<2e-3,'Cannot reproduce saved baseline RMI.');
    end

    % Verify each adaptive map against PREVIOUS phase-ramp validation,
    % when its MAT artifact exists. Fail loudly on methodological drift.
    previous_file=fullfile(pdir,'rmi_phase_ramp_validation', ...
        'rmi_phase_ramp_results.mat');
    if isfile(previous_file)
        P=load(previous_file,'rmaps');
        if isfield(P,'rmaps')
            assert(isequal(size(P.rmaps),size(adaptive)), ...
                'Prior phase-ramp maps have different dimensions.');
            for k=1:nBeta
                for mi=1:2
                    old=double(P.rmaps{k,mi});
                    valid_old=isfinite(old) & isfinite(adaptive{k,mi});
                    assert(nnz(valid_old)>100,'No overlap with prior RMI map.');
                    discrepancy=max(abs(old(valid_old)-adaptive{k,mi}(valid_old)));
                    assert(discrepancy<2e-3, ...
                        'Adaptive RMI mismatch with prior validation.');
                end
            end
            fprintf('Adaptive RMI agrees with saved phase-ramp experiment.\n');
        end
    end

    % SAME pixels for all betas, both estimators and both c variants.
    valid=logical(S.common)&logical(R.basic)&isfinite(double(S.GT));
    for k=1:nBeta
        M=S.maps(k);
        valid=valid & isfinite(M.PG)&M.PG>0&M.PG<10 & ...
            isfinite(M.AC)&M.AC>0&M.AC<10 & isfinite(M.Q);
        for mi=1:2
            valid=valid & isfinite(adaptive{k,mi}) & ...
                isfinite(frozen{k,mi});
        end
    end
    fprintf('Common comparison pixels: %d (%.2f %% of prior common)\n', ...
        nnz(valid),100*nnz(valid)/max(nnz(S.common & R.basic),1));
    assert(nnz(valid)>100,'Too few common pixels.');

    % Compute ground-truth error only AFTER all Q maps are constructed.
    % For each (region, estimator, beta), compare adaptive/frozen on SAME
    % pixels. Block-level analysis can be added later if required.
    T=table();
    cov=(0.05:0.05:1)';
    for ri=1:numel(R.regions)
        rm=logical(R.masks{ri})&valid;
        if nnz(rm)<30
            warning('Skipped region %s (fewer than 30 pixels).',R.regions{ri});
            continue;
        end
        for mi=1:2
            for k=1:nBeta
                c=double(S.maps(k).(methods{mi}));
                gt=double(S.GT);
                e=abs(c(rm)-gt(rm));
                qa=double(adaptive{k,mi}(rm));
                qf=double(frozen{k,mi}(rm));
                ra=double(qresA{k,mi}(rm));
                rf=double(qresF{k,mi}(rm));
                aa=double(qamp{k}(rm));
                rA=tie_aware_risk(qa,e,cov);
                rF=tie_aware_risk(qf,e,cov);
                riskRandom=mean(e);
                aurcA=trapz(cov,rA)/(cov(end)-cov(1))/riskRandom;
                aurcF=trapz(cov,rF)/(cov(end)-cov(1))/riskRandom;
                row=table(string(phantom_case),string(R.regions{ri}), ...
                    string(methods{mi}),betas(k),nnz(rm), ...
                    sqrt(mean((c(rm)-gt(rm)).^2)),riskRandom, ...
                    mean(qa),mean(qf),mean(ra),mean(rf),mean(aa), ...
                    mean(qa-qf),mean(ra-rf), ...
                    spearman_safe(qa,e),spearman_safe(qf,e), ...
                    aurcA,aurcF,rA(10),rF(10), ...
                    mean(ra<=1e-12),mean(ra>=1-1e-12), ...
                    mean(rf<=1e-12),mean(rf>=1-1e-12), ...
                    'VariableNames',{'Phantom','Region','Estimator', ...
                    'Beta_rad_per_mm','N_Pixels','RMSE_mps','MAE_mps', ...
                    'MeanRMI_Adaptive','MeanRMI_Frozen', ...
                    'MeanQres_Adaptive','MeanQres_Frozen','MeanQamp', ...
                    'SpeedEffect_MeanRMI','SpeedEffect_MeanQres', ...
                    'Spearman_Adaptive','Spearman_Frozen', ...
                    'AURCratio_Adaptive','AURCratio_Frozen', ...
                    'Risk50_Adaptive_mps','Risk50_Frozen_mps', ...
                    'FracResAdaptive0','FracResAdaptive1', ...
                    'FracResFrozen0','FracResFrozen1'});
                T=[T;row]; %#ok<AGROW>
            end
        end
    end

    % Exact mean decomposition relative to beta=0 on fixed pixel masks:
    % deltaAdaptive = deltaFrozen + (adaptive - frozen) at each beta.
    T.DeltaRMI_Adaptive=nan(height(T),1);
    T.DisplacementEffect=nan(height(T),1);
    T.DeltaRMSE_mps=nan(height(T),1);
    for j=1:height(T)
        sel=T.Region==T.Region(j)&T.Estimator==T.Estimator(j) & ...
            T.Beta_rad_per_mm==0;
        assert(nnz(sel)==1, 'Expected exactly one beta=0 per region/estimator.');
        T.DeltaRMI_Adaptive(j)=T.MeanRMI_Adaptive(j)-T.MeanRMI_Adaptive(sel);
        T.DisplacementEffect(j)=T.MeanRMI_Frozen(j)-T.MeanRMI_Frozen(sel);
        T.DeltaRMSE_mps(j)=T.RMSE_mps(j)-T.RMSE_mps(sel);
    end
    discrepancy=max(abs(T.DeltaRMI_Adaptive - ...
        (T.DisplacementEffect+T.SpeedEffect_MeanRMI)));
    assert(discrepancy<1e-9, ...
        'Displacement + speed effects do not sum to total change.');
    fprintf('Mean RMI decomposition identity: max residual %.3g\n',discrepancy);

    outdir=fullfile(pdir,'rmi_speed_diagnostic');
    if ~isfolder(outdir),mkdir(outdir);end
    writetable(T,fullfile(outdir,'rmi_speed_summary.csv'));
    save(fullfile(outdir,'rmi_speed_results.mat'), ...
        'T','betas','valid','lim','adaptive','frozen', ...
        'qresA','qresF','qamp','-v7.3');
    plot_trends(T,betas,outdir,phantom_case);
    plot_decomposition(T,betas,outdir,phantom_case);
    plot_maps(adaptive,frozen,betas,x,z,outdir,phantom_case,methods);
    plot_risk(S,adaptive,frozen,betas,valid,cov,outdir,phantom_case);

    fprintf('\nGLOBAL SPEED-COUPLING RESULTS (fixed common pixels):\n');
    cols={'Beta_rad_per_mm','Estimator','RMSE_mps', ...
        'MeanRMI_Adaptive','MeanRMI_Frozen','DisplacementEffect', ...
        'SpeedEffect_MeanRMI','DeltaRMI_Adaptive', ...
        'AURCratio_Adaptive','AURCratio_Frozen'};
    disp(T(T.Region=="Global",cols));
    fprintf('Saved results in: %s\n',outdir);
    fprintf(['Interpretation: Frozen-c is an analytical CONTROL, not a new ' ...
        'RMI metric. Phase ramps are synthetic, and RMI remains ' ...
        'the Eq. 6 1D adaptation.\n']);
end

function [rA,rF,amp]=wave_residual_both(u,da,db,cos_t,sin_t,x,prf,cBeta,c0)
% Same finite differences/denominator/temporal frames as earlier RMI code.
% Reuse each u_tt and u_xx for both c choices, preventing inconsistencies.
    [Nz,Nx,Nt]=size(u);
    na={zeros(Nz,Nx),zeros(Nz,Nx)};
    nf={zeros(Nz,Nx),zeros(Nz,Nx)};
    den=zeros(Nz,Nx);
    cb2=cellfun(@(c)c.^2,cBeta,'UniformOutput',false);
    c02=cellfun(@(c)c.^2,c0,'UniformOutput',false);
    dxr=x(3:end)-x(2:end-1);
    dxl=x(2:end-1)-x(1:end-2);
    dxspan=x(3:end)-x(1:end-2);
    frame=@(t)double(u(:,:,t))+da*cos_t(t)+db*sin_t(t);
    um=frame(1);uz=frame(2);
    amp=max(abs(um),abs(uz));
    for it=2:Nt-1
        up=frame(it+1);
        utt=(up-2*uz+um)*prf^2;
        fwd=bsxfun(@rdivide,uz(:,3:end)-uz(:,2:end-1),dxr');
        bwd=bsxfun(@rdivide,uz(:,2:end-1)-uz(:,1:end-2),dxl');
        uxx=bsxfun(@rdivide,2*(fwd-bwd),dxspan');
        den(:,2:end-1)=den(:,2:end-1)+uxx.^2;
        for mi=1:2
            eA=utt(:,2:end-1)-cb2{mi}(:,2:end-1).*uxx;
            eF=utt(:,2:end-1)-c02{mi}(:,2:end-1).*uxx;
            na{mi}(:,2:end-1)=na{mi}(:,2:end-1)+eA.^2;
            nf{mi}(:,2:end-1)=nf{mi}(:,2:end-1)+eF.^2;
        end
        amp=max(amp,abs(up));
        um=uz;uz=up;
    end
    positive=den(isfinite(den)&den>0);
    assert(~isempty(positive),'Zero spatial curvature in displacement.');
    floorDen=median(positive)*1e-12;
    invalid=~isfinite(den)|den<=floorDen;
    rA=cell(1,2);rF=cell(1,2);
    for mi=1:2
        rA{mi}=na{mi}./den;
        rF{mi}=nf{mi}./den;
        rA{mi}(:,[1 end])=NaN;rF{mi}(:,[1 end])=NaN;
        rA{mi}(invalid)=NaN;rF{mi}(invalid)=NaN;
    end
end

function y=clip01(x)
    y=max(0,min(1,x));
end

function risk=tie_aware_risk(q,e,coverage)
% Fractional inclusion of a tie group; robust when many Q scores are 1.
    q=double(q(:));e=double(e(:));
    assert(~isempty(q)&&numel(q)==numel(e) && ...
        all(isfinite(q))&&all(isfinite(e))&&all(e>=0), ...
        'Risk-coverage requires finite Q and nonnegative errors.');
    [sq,ord]=sort(q,'descend');e=e(ord);
    group=cumsum([true;diff(sq)~=0]);
    counts=accumarray(group,1);errsum=accumarray(group,e);
    cs=cumsum(counts);es=cumsum(errsum);
    risk=zeros(size(coverage));
    for i=1:numel(coverage)
        nWanted=coverage(i)*numel(q);
        g=find(cs>=nWanted-1e-8,1,'first');
        if isempty(g),g=numel(counts);end
        if g==1,nPrevious=0;ePrevious=0;
        else,nPrevious=cs(g-1);ePrevious=es(g-1);end
        nPartial=max(0,min(counts(g),nWanted-nPrevious));
        risk(i)=(ePrevious+nPartial*errsum(g)/counts(g))/nWanted;
    end
    assert(abs(risk(end)-mean(e))<1e-8*max(1,mean(e)), ...
        'Risk at 100 pct must equal the unfiltered MAE.');
end

function rho=spearman_safe(q,e)
    if numel(q)<5||numel(unique(q))<3||numel(unique(e))<3
        rho=NaN;return;
    end
    rho=corr(q(:),e(:),'Type','Spearman');
end

function plot_trends(T,betas,outdir,p)
    f=figure('Color','w','Name',['Speed dependence: ' p], ...
        'Position',[30 40 1300 740]);
    tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
    for mi=1:2
        method=string({'PG','AC'});method=method(mi);
        G=T(T.Region=="Global" & T.Estimator==method,:);
        [~,order]=sort(G.Beta_rad_per_mm);G=G(order,:);
        ax=nexttile(tl,mi);hold(ax,'on');
        plot(ax,G.Beta_rad_per_mm,G.MeanRMI_Adaptive,'-o', ...
            'LineWidth',1.7,'DisplayName','Adaptive c_{beta}');
        plot(ax,G.Beta_rad_per_mm,G.MeanRMI_Frozen,'-s', ...
            'LineWidth',1.7,'DisplayName','Frozen c_0');
        ylabel(ax,'Mean RMI');xlabel(ax,'Phase ramp [rad/mm]');
        title(ax,[char(method) ' / same displacement, different c']);
        ylim(ax,[0 1]);grid(ax,'on');legend(ax,'Location','best');
        ax=nexttile(tl,mi+2);hold(ax,'on');
        plot(ax,G.Beta_rad_per_mm,G.RMSE_mps,'-o','LineWidth',1.7);
        ylabel(ax,'SWS RMSE [m/s] (validation only)');
        xlabel(ax,'Phase ramp [rad/mm]');title(ax,[char(method) ' / SWS error']);
        grid(ax,'on');
    end
    title(tl,['Does c_{est} affect RMI? ' strrep(p,'_',' ')]);
    exportgraphics(f,fullfile(outdir,'rmi_speed_trends.png'),'Resolution',170);
end

function plot_decomposition(T,betas,outdir,p)
    f=figure('Color','w','Name',['RMI delta decomposition: ' p], ...
        'Position',[40 60 1300 460]);
    tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
    for mi=1:2
        method=string({'PG','AC'});method=method(mi);
        G=T(T.Region=="Global" & T.Estimator==method,:);
        [~,order]=sort(G.Beta_rad_per_mm);G=G(order,:);
        ax=nexttile(tl,mi);hold(ax,'on');
        bar(ax,G.Beta_rad_per_mm, ...
            [G.DisplacementEffect G.SpeedEffect_MeanRMI],'stacked');
        plot(ax,G.Beta_rad_per_mm,G.DeltaRMI_Adaptive,'ko-', ...
            'LineWidth',1.5,'MarkerFaceColor','k');
        yline(ax,0,'k:');grid(ax,'on');
        xlabel(ax,'Phase ramp [rad/mm]');ylabel(ax,'Change in mean RMI from beta=0');
        title(ax,[char(method) ' / exact RMI change decomposition']);
        legend(ax,{'Displacement effect (frozen c)', ...
            'Speed-estimate effect','Total adaptive delta'},'Location','best');
    end
    title(tl,['Displacement vs reconstructed-speed contributions: ' ...
        strrep(p,'_',' ')]);
    exportgraphics(f,fullfile(outdir,'rmi_speed_decomposition.png'), ...
        'Resolution',170);
end

function plot_maps(adaptive,frozen,betas,x,z,outdir,p,methods)
    for mi=1:2
        f=figure('Color','w','Name',[methods{mi} ' speed-coupling maps'], ...
            'Position',[40 60 1480 1000]);
        tl=tiledlayout(f,3,numel(betas),'TileSpacing','compact', ...
            'Padding','compact');
        for k=1:numel(betas)
            ax=nexttile(tl,k);
            imagesc(ax,x*1e3,z*1e3,adaptive{k,mi});
            set(ax,'YDir','normal');axis(ax,'image');clim(ax,[0 1]);colorbar(ax);
            title(ax,sprintf('Adaptive: beta=%+.2f',betas(k)));
            xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
            ax=nexttile(tl,k+numel(betas));
            imagesc(ax,x*1e3,z*1e3,frozen{k,mi});
            set(ax,'YDir','normal');axis(ax,'image');clim(ax,[0 1]);colorbar(ax);
            title(ax,sprintf('Frozen c0: beta=%+.2f',betas(k)));
            xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
            ax=nexttile(tl,k+2*numel(betas));
            imagesc(ax,x*1e3,z*1e3,adaptive{k,mi}-frozen{k,mi});
            set(ax,'YDir','normal');axis(ax,'image');clim(ax,[-1 1]);colorbar(ax);
            title(ax,sprintf('Adaptive - frozen: beta=%+.2f',betas(k)));
            xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
        end
        title(tl,[methods{mi} ' / fixed-speed control: ' strrep(p,'_',' ')]);
        exportgraphics(f,fullfile(outdir,['rmi_speed_maps_' ...
            lower(methods{mi}) '.png']),'Resolution',160);
    end
end

function plot_risk(S,adaptive,frozen,betas,valid,cov,outdir,p)
    f=figure('Color','w','Name',['RMI risk curves: ' p], ...
        'Position',[40 55 1450 850]);
    tl=tiledlayout(f,2,numel(betas),'TileSpacing','compact', ...
        'Padding','compact');
    methods={'PG','AC'};
    for mi=1:2
        for k=1:numel(betas)
            m=S.maps(k);gt=double(S.GT(valid));
            est=double(m.(methods{mi}));
            e=abs(est(valid)-gt);
            rA=tie_aware_risk(adaptive{k,mi}(valid),e,cov);
            rF=tie_aware_risk(frozen{k,mi}(valid),e,cov);
            ax=nexttile(tl,(mi-1)*numel(betas)+k);hold(ax,'on');
            plot(ax,cov*100,rA,'LineWidth',1.5,'DisplayName','Adaptive c_{beta}');
            plot(ax,cov*100,rF,'LineWidth',1.5,'DisplayName','Frozen c_0');
            yline(ax,mean(e),'k--','DisplayName','Random');
            xlabel(ax,'Accepted coverage [%]');ylabel(ax,'Accepted MAE [m/s]');
            title(ax,sprintf('%s: beta=%+.2f rad/mm',methods{mi},betas(k)));
            grid(ax,'on');legend(ax,'Location','best');
        end
    end
    title(tl,['Same-pixel risk-coverage: ' strrep(p,'_',' ')]);
    exportgraphics(f,fullfile(outdir,'rmi_speed_risk_coverage.png'), ...
        'Resolution',170);
end
