function T = validate_rmi_phase_ramp(phantom_case)
%VALIDATE_RMI_PHASE_RAMP Test RMI-1D vs Phase Planarity Q under phase ramps.
%
% THIS IS A NEW VALIDATION TEST, NOT A NEW QUALITY METRIC.
% Compared quality metrics:
%   (1) Phase Planarity Q: our experimental candidate.
%   (2) RMI-1D: adapted Choi et al. 2015 Eq. (6), using fixed calibration.
%       DOI 10.1109/ULTSYM.2015.0382.
%
% Both quality metrics are ground-truth independent. Ground truth enters
% ONLY the evaluation (MAE, RMSE, Spearman, risk-coverage and AURC).
% RMI contains the reconstructed PG/AC velocity: its residual can share
% estimator biases. The perturbations are controlled, not physical waves.
%
% REQUIRED FIRST:
%   validate_phase_ramp_sws('large_inclusions')  % PG/AC maps for each beta
%   validate_sws_rmi('large_inclusions','calibrate') % only ONCE, already done
%   validate_sws_rmi('large_inclusions')            % saved RMI baseline
% RUN:
%   validate_rmi_phase_ramp('large_inclusions')
%   validate_rmi_phase_ramp('high_contrast')
%
% Reads:
%   Results/sws_measurement_quality/rmi_1d_calibration.mat
%   Results/sws_measurement_quality/<case>/rmi_1d_diagnostic/rmi_1d_results.mat
%   Results/sws_measurement_quality/<case>/phase_ramp_sws_validation/
%      phase_ramp_metrics.mat
%   Data/<case>_displacement_with_quality.mat
%
% Saves:
%   Results/sws_measurement_quality/<case>/rmi_phase_ramp_validation/
%      rmi_phase_ramp_summary.csv, rmi_phase_ramp_curves.csv,
%      rmi_phase_ramp_results.mat, rmi_phase_ramp_maps.png,
%      rmi_phase_ramp_trends.png, rmi_phase_ramp_risk_coverage.png
%
% Uses EXACT SAME finite-difference residual, empirical calibration, and
% displacement amplitude conventions as validate_sws_rmi.m. The harmonic
% phase ramp follows validate_phase_ramp_sws.m (LS harmonic fit and
% replacement of its cosine/sine component; remaining signal untouched).
%
% NOTE: acceptance at tied Q scores is fractional and unbiased. AURC
% covers 5-100%% (not 0-100%%). Pixels are spatially correlated.
%
    if nargin < 1 || isempty(phantom_case)
        phantom_case = 'large_inclusions';
    end
    phantom_case = char(string(phantom_case));
    root = fileparts(mfilename('fullpath'));
    if ~isfolder(fullfile(root,'Results')), root = fileparts(root); end
    top = fullfile(root,'Results','sws_measurement_quality');
    pdir = fullfile(top,phantom_case);
    source_file = fullfile(pdir,'phase_ramp_sws_validation','phase_ramp_metrics.mat');
    rmi_file = fullfile(pdir,'rmi_1d_diagnostic','rmi_1d_results.mat');
    calib_file = fullfile(top,'rmi_1d_calibration.mat');
    assert(isfile(source_file),'Missing phase-ramp maps: %s',source_file);
    assert(isfile(rmi_file),'Missing baseline RMI maps: %s',rmi_file);
    assert(isfile(calib_file),'Missing fixed RMI calibration: %s',calib_file);

    S=load(source_file,'maps','GT','common','betas','prf','f0','x','z');
    R=load(rmi_file,'RmiPG','RmiAC','basic','masks','regions','GT','x','z');
    K=load(calib_file,'limits');
    lim=K.limits;
    assert(strcmp(lim.algorithm,'Choi2015_Eq6_RMI_1D_adaptation'), ...
        'Wrong calibration algorithm. Do NOT substitute a different RMI calibration.');
    required={'res_min','res_max','u_min','u_max'};
    assert(all(isfield(lim,required)) && ...
        lim.res_max>lim.res_min && lim.u_max>lim.u_min, ...
        'Invalid fixed calibration bounds.');
    x=double(S.x(:)); z=double(S.z(:));
    [Nz,Nx]=size(S.GT);
    assert(isequal(size(R.GT),[Nz,Nx],size(R.basic),size(S.common)), ...
        'Saved maps have incompatible dimensions.');
    assert(max(abs(x-double(R.x(:))))<1e-9 && ...
        max(abs(z-double(R.z(:))))<1e-9, 'Lateral/depth axes mismatch.');
    assert(max(abs(double(S.GT(:))-double(R.GT(:))))<1e-6, ...
        'Ground truth differs between saved analyses.');
    betas=double(S.betas(:));
    assert(numel(betas)==numel(S.maps) && ...
        nnz(betas==0)==1 && numel(betas)>=3, ...
        'Need phase-ramp maps with exactly one beta=0 baseline.');
    assert(numel(R.masks)==numel(R.regions),'Invalid stored region masks.');
    for k=1:numel(S.maps)
        assert(abs(S.maps(k).beta-betas(k))<1e-12, ...
            'Phase-ramp beta does not match stored reconstruction.');
        assert(isequal(size(S.maps(k).PG),[Nz Nx], ...
            size(S.maps(k).AC),size(S.maps(k).Q)), ...
            'Inconsistent PG/AC/Phase Q map dimensions.');
    end

    % Ensure USTB objects can be loaded.
    if isempty(which('uff.beamformed_data'))
        ustb_root=fullfile(fileparts(root),'elastography','ElastographyToolbox','USTB');
        if isfolder(ustb_root), addpath(genpath(ustb_root)); rehash; end
    end
    assert(~isempty(which('uff.beamformed_data')), ...
        'USTB classes missing. Add ElastographyToolbox/USTB to MATLAB path.');
    data_name=[phantom_case '_displacement_with_quality.mat'];
    candidates={fullfile(root,'Data',data_name), ...
        fullfile(root,'Data',phantom_case,data_name),fullfile(root,data_name), ...
        fullfile(root,'Analysis',data_name)};
    datafile='';
    for i=1:numel(candidates)
        if isfile(candidates{i}),datafile=candidates{i};break;end
    end
    assert(~isempty(datafile),'Cannot find %s in project.',data_name);

    fprintf('\n=== RMI vs PHASE Q | PHASE-RAMP VALIDATION: %s ===\n',phantom_case);
    fprintf('Using already reconstructed PG/AC maps and fixed RMI calibration from %s\n', ...
        lim.reference_phantom);
    fprintf('Loading displacement from %s\n',datafile);
    D=load(datafile,'channel_data','displacement_estimation_autocorr', ...
        'shear_params');
    prf=double(D.channel_data.PRF);
    f0=double(D.shear_params.source_freq);
    assert(abs(prf-double(S.prf))<1e-8 && ...
           abs(f0-double(S.f0))<1e-8, ...
        'Dataset PRF or excitation frequency differs from phase-ramp source.');
    raw=D.displacement_estimation_autocorr.data;
    assert(isnumeric(raw) && isreal(raw) && ...
        mod(numel(raw),Nx*Nz)==0,'Invalid displacement array.');
    Nt=numel(raw)/(Nx*Nz);
    assert(Nt>=5,'Need at least 5 displacement frames.');
    % Negative sign matches validate_phase_ramp_sws. RMI residual is
    % invariant to flipping displacement globally (squared terms).
    u=-reshape(single(raw),Nz,Nx,Nt);
    clear D raw;
    fprintf('Displacement %dx%dx%d | PRF %.1f Hz | f0 %.1f Hz\n', ...
        Nz,Nx,Nt,prf,f0);

    % Reproduce harmonic-fit convention from validate_phase_ramp_sws.m.
    tt=(0:Nt-1)'/prf; omega=2*pi*f0;
    cos_t=cos(omega*tt); sin_t=sin(omega*tt);
    design=[cos_t sin_t ones(Nt,1)];
    assert(rank(design)==3,'Harmonic LS design rank deficient.');
    pinv_design=pinv(design);
    a=zeros(Nz,Nx); b=zeros(Nz,Nx);
    for it=1:Nt
        frame=double(u(:,:,it));
        a=a+pinv_design(1,it)*frame;
        b=b+pinv_design(2,it)*frame;
    end
    [Xmm,~]=meshgrid(x*1e3,z*1e3);
    rmaps=cell(numel(betas),2);
    % Run RMI on phase-ramped DISPLACEMENT, not on SWS maps alone.
    for k=1:numel(betas)
        beta=betas(k);
        theta=beta*(Xmm-mean(x)*1e3);
        da=a.*(cos(theta)-1)+b.*sin(theta);
        db=b.*(cos(theta)-1)-a.*sin(theta);
        fprintf('RMI finite-difference residual for beta=%+.2f rad/mm...\n',beta);
        cPG=double(S.maps(k).PG);
        cAC=double(S.maps(k).AC);
        % Match validate_sws_rmi's invalid-speed guard.
        cPG(~isfinite(cPG)|cPG<=0|cPG>=10)=NaN;
        cAC(~isfinite(cAC)|cAC<=0|cAC>=10)=NaN;
        [rp,ra,amp]=residual_ramped_1d(u,da,db,cos_t,sin_t,x,prf,cPG,cAC);
        rmaps{k,1}=rmi_from_residual(rp,amp,lim);
        rmaps{k,2}=rmi_from_residual(ra,amp,lim);
        clear da db rp ra amp;
    end
    clear u a b;

    k0=find(betas==0,1);
    baseline_mask=logical(R.basic) & ...
        isfinite(rmaps{k0,1}) & isfinite(rmaps{k0,2});
    assert(nnz(baseline_mask)>100,'Not enough baseline RMI pixels.');
    refPG=double(R.RmiPG); refAC=double(R.RmiAC);
    bdiffPG=max(abs(rmaps{k0,1}(baseline_mask)-refPG(baseline_mask)));
    bdiffAC=max(abs(rmaps{k0,2}(baseline_mask)-refAC(baseline_mask)));
    fprintf('RMI baseline max-map discrepancy: PG %.5g; AC %.5g\n', ...
        bdiffPG,bdiffAC);
    % A mismatch here indicates inconsistency between the code versions.
    assert(bdiffPG<2e-3 && bdiffAC<2e-3, ...
        ['Rebuilt RMI does not reproduce saved original RMI map. ' ...
         'Check saved code versions, frame ordering and calibration.']);

    % Same FIXED pixels for all beta values, estimators and Q methods.
    valid=logical(S.common) & logical(R.basic) & isfinite(double(S.GT));
    for k=1:numel(betas)
        M=S.maps(k);
        valid=valid & isfinite(M.Q) & ...
            isfinite(M.PG) & M.PG>0 & M.PG<10 & ...
            isfinite(M.AC) & M.AC>0 & M.AC<10 & ...
            isfinite(rmaps{k,1}) & isfinite(rmaps{k,2});
    end
    fraction=nnz(valid)/max(nnz(S.common & R.basic),1);
    fprintf('Joint mask: %d pixels (%.2f%% of previously common pixels)\n', ...
        nnz(valid),100*fraction);
    assert(nnz(valid)>=100,'Too few pixels remain valid across ramps.');
    if fraction<0.9
        warning('Under 90%% common coverage; report excluded-pixel fraction.');
    end

    [X,Z]=meshgrid(x,z);
    bx=floor((X-min(x))/(1e-3));
    bz=floor((Z-min(z))/(1e-3));
    blockid=bz*(max(bx(:))+1)+bx+1;
    regions=R.regions;
    methods={'PG','AC'};
    quality_names={'Phase Planarity Q','RMI-1D'};
    coverage=(0.05:0.05:1)';
    T=table(); Curves=table();
    for ri=1:numel(regions)
        mask=logical(R.masks{ri}) & valid;
        if nnz(mask)<30
            warning('Skipping region %s: fewer than 30 pixels.',regions{ri});
            continue;
        end
        ids=blockid(mask);
        for k=1:numel(betas)
            for mi=1:2
                est=double(S.maps(k).(methods{mi}));
                err=abs(est(mask)-double(S.GT(mask)));
                err_sign=est(mask)-double(S.GT(mask));
                Qs={double(S.maps(k).Q(mask)),double(rmaps{k,mi}(mask))};
                for qi=1:2
                    q=Qs{qi};
                    for unit=1:2
                        if unit==1
                            qq=q; ee=err; label='Pixel';
                        else
                            [qq,ee]=spatial_means(ids,q,err);
                            label='Block1mm';
                        end
                        risk=tie_aware_risk(qq,ee,coverage);
                        random_risk=mean(ee);
                        aurc=trapz(coverage,risk)/(coverage(end)-coverage(1));
                        ratio=aurc/random_risk;
                        spearman=spearman_safe(qq,ee);
                        if unit==1
                            RMSE=sqrt(mean(err_sign.^2));
                            MAE=mean(err);
                        else
                            RMSE=NaN; % pixel RMSE is not defined for block-mean absolute errors
                            MAE=mean(ee);
                        end
                        row=table(string(phantom_case),string(regions{ri}), ...
                            betas(k),string(methods{mi}),string(quality_names{qi}), ...
                            string(label),numel(ee),nnz(mask),mean(qq), ...
                            RMSE,MAE,spearman,ratio,risk(4),risk(10),risk(16), ...
                            'VariableNames',{'Phantom','Region','Beta_rad_per_mm', ...
                            'Estimator','QualityMetric','Unit','N_Units','N_Pixels', ...
                            'MeanQ','RMSE_mps','MAE_mps','Spearman_Q_vs_AbsError', ...
                            'AURC_over_Random','Risk20_mps','Risk50_mps','Risk80_mps'});
                        T=[T;row]; %#ok<AGROW>
                        n=numel(coverage);
                        part=table(repmat(string(phantom_case),n,1), ...
                            repmat(string(regions{ri}),n,1), ...
                            repmat(betas(k),n,1),repmat(string(methods{mi}),n,1), ...
                            repmat(string(quality_names{qi}),n,1), ...
                            repmat(string(label),n,1),coverage*100,risk, ...
                            repmat(random_risk,n,1), ...
                            'VariableNames',{'Phantom','Region','Beta_rad_per_mm', ...
                            'Estimator','QualityMetric','Unit','Coverage_pct', ...
                            'RiskMAE_mps','RandomRisk_mps'});
                        Curves=[Curves;part]; %#ok<AGROW>
                    end
                end
            end
        end
    end
    % Add delta columns to make Q-vs-error response explicit.
    T.DeltaMeanQ=nan(height(T),1);
    T.DeltaRMSE_mps=nan(height(T),1);
    T.DeltaAURC_ratio=nan(height(T),1);
    for i=1:height(T)
        same=T.Region==T.Region(i)&T.Estimator==T.Estimator(i)& ...
            T.QualityMetric==T.QualityMetric(i)&T.Unit==T.Unit(i)& ...
            T.Beta_rad_per_mm==0;
        assert(nnz(same)==1,'Expected exactly one baseline row.');
        T.DeltaMeanQ(i)=T.MeanQ(i)-T.MeanQ(same);
        T.DeltaRMSE_mps(i)=T.RMSE_mps(i)-T.RMSE_mps(same);
        T.DeltaAURC_ratio(i)=T.AURC_over_Random(i)-T.AURC_over_Random(same);
    end
    G=T(T.Region=="Global" & T.Unit=="Pixel",:);
    disp(G(:,{'Beta_rad_per_mm','Estimator','QualityMetric','RMSE_mps', ...
        'MeanQ','DeltaMeanQ','Spearman_Q_vs_AbsError','AURC_over_Random', ...
        'Risk50_mps'}));
    outdir=fullfile(pdir,'rmi_phase_ramp_validation');
    if ~isfolder(outdir),mkdir(outdir);end
    writetable(T,fullfile(outdir,'rmi_phase_ramp_summary.csv'));
    writetable(Curves,fullfile(outdir,'rmi_phase_ramp_curves.csv'));
    save(fullfile(outdir,'rmi_phase_ramp_results.mat'),'T','Curves', ...
        'rmaps','betas','valid','lim','bdiffPG','bdiffAC', ...
        'fraction','coverage','-v7.3');
    plot_maps(S,rmaps,betas,k0,x,z,outdir,phantom_case);
    plot_trends(T,methods,quality_names,betas,outdir,phantom_case);
    plot_risk(Curves,methods,quality_names,betas,outdir,phantom_case);
    fprintf('Saved validation output: %s\n',outdir);
    fprintf(['Interpretation: this is a synthetic phase-ramp perturbation; ' ...
        'RMI uses reconstructed speed; no GT enters Q calculation.\n']);
end

function [rp,ra,amp]=residual_ramped_1d(u,da,db,cos_t,sin_t,x,prf,cPG,cAC)
% The finite-difference operator matches validate_sws_rmi.m Eq. 6.
    [Nz,Nx,Nt]=size(u);
    numPG=zeros(Nz,Nx); numAC=zeros(Nz,Nx); den=zeros(Nz,Nx);
    cp2=cPG.^2; ca2=cAC.^2;
    dxr=x(3:end)-x(2:end-1);
    dxl=x(2:end-1)-x(1:end-2);
    dxspan=x(3:end)-x(1:end-2);
    getframe=@(k) double(u(:,:,k))+da*cos_t(k)+db*sin_t(k);
    um=getframe(1); u0=getframe(2);
    amp=max(abs(um),abs(u0));
    for it=2:Nt-1
        up=getframe(it+1);
        utt=(up-2*u0+um)*prf^2;
        fwd=bsxfun(@rdivide,u0(:,3:end)-u0(:,2:end-1),dxr');
        bwd=bsxfun(@rdivide,u0(:,2:end-1)-u0(:,1:end-2),dxl');
        uxx=bsxfun(@rdivide,2*(fwd-bwd),dxspan');
        p=utt(:,2:end-1)-cp2(:,2:end-1).*uxx;
        q=utt(:,2:end-1)-ca2(:,2:end-1).*uxx;
        numPG(:,2:end-1)=numPG(:,2:end-1)+p.^2;
        numAC(:,2:end-1)=numAC(:,2:end-1)+q.^2;
        den(:,2:end-1)=den(:,2:end-1)+uxx.^2;
        amp=max(amp,abs(up));
        um=u0;u0=up;
    end
    finiteDen=den(isfinite(den)&den>0);
    assert(~isempty(finiteDen),'Zero spatial second derivative throughout wavefield.');
    floorDen=median(finiteDen)*1e-12;
    rp=numPG./den; ra=numAC./den;
    rp(:,[1 end])=NaN; ra(:,[1 end])=NaN;
    bad=~isfinite(den)|den<=floorDen;
    rp(bad)=NaN;ra(bad)=NaN;
end

function q=rmi_from_residual(res,amp,lim)
    r=(lim.res_max-res)/(lim.res_max-lim.res_min);
    a=(amp-lim.u_min)/(lim.u_max-lim.u_min);
    r=max(0,min(1,r));a=max(0,min(1,a));
    q=0.8*r+0.2*a;
    q(~isfinite(res)|~isfinite(amp))=NaN;
end

function [qb,eb]=spatial_means(ids,q,e)
    [~,~,g]=unique(double(ids(:)));
    counts=accumarray(g,1);
    qb=accumarray(g,double(q(:)))./counts;
    eb=accumarray(g,double(e(:)))./counts;
end

function risk=tie_aware_risk(q,e,cov)
% Fair fractional inclusion when all / many quality values are tied.
    q=double(q(:));e=double(e(:));
    assert(~isempty(q)&&numel(q)==numel(e) && ...
        all(isfinite(q))&&all(isfinite(e))&&all(e>=0));
    [s,ix]=sort(q,'descend');e=e(ix);
    gid=cumsum([true;diff(s)~=0]);
    count=accumarray(gid,1); errsum=accumarray(gid,e);
    cn=cumsum(count);ce=cumsum(errsum);
    risk=zeros(size(cov));
    for k=1:numel(cov)
        wanted=double(cov(k))*numel(q);
        g=find(cn>=wanted-1e-8,1,'first');
        if isempty(g),g=numel(count);end
        if g==1
            nprev=0; eprev=0;
        else
            nprev=cn(g-1); eprev=ce(g-1);
        end
        part=max(0,min(count(g),wanted-nprev));
        risk(k)=(eprev+part*errsum(g)/count(g))/wanted;
    end
    assert(abs(risk(end)-mean(e))<1e-8*max(1,mean(e)), ...
        'Tie-aware risk at 100%% should be unfiltered MAE.');
end

function r=spearman_safe(a,b)
    a=double(a(:));b=double(b(:));
    if numel(a)<5||numel(unique(a))<3||numel(unique(b))<3
        r=NaN;return;
    end
    r=corr(a,b,'Type','Spearman');
end

function plot_maps(S,rmaps,betas,k0,x,z,outdir,p)
    order=unique([1 k0 numel(betas)],'stable');
    f=figure('Color','w','Name',['RMI phase ramp maps: ' p], ...
        'Position',[40 70 1450 730]);
    tl=tiledlayout(f,2,numel(order),'TileSpacing','compact','Padding','compact');
    for j=1:numel(order)
        k=order(j);
        ax=nexttile(tl,j);
        imagesc(ax,x*1000,z*1000,rmaps{k,1});
        set(ax,'YDir','normal');axis(ax,'image');clim(ax,[0 1]);colorbar(ax);
        title(ax,sprintf('RMI/PG: beta=%+.2f rad/mm',betas(k)));
        xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
        ax=nexttile(tl,numel(order)+j);
        imagesc(ax,x*1000,z*1000,rmaps{k,2});
        set(ax,'YDir','normal');axis(ax,'image');clim(ax,[0 1]);colorbar(ax);
        title(ax,sprintf('RMI/AC: beta=%+.2f rad/mm',betas(k)));
        xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
    end
    title(tl,['RMI-1D on phase-ramped displacement: ' strrep(p,'_',' ')]);
    exportgraphics(f,fullfile(outdir,'rmi_phase_ramp_maps.png'),'Resolution',160);
end

function plot_trends(T,methods,quality_names,betas,outdir,p)
    f=figure('Color','w','Name',['RMI vs Phase Q trends: ' p], ...
        'Position',[40 60 1250 820]);
    tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
    for mi=1:2
        g=T(T.Region=="Global"&T.Unit=="Pixel"& ...
            T.Estimator==string(methods{mi}),:);
        subplot_idx=mi;
        ax=nexttile(tl,subplot_idx);hold(ax,'on');
        for qi=1:2
            h=g(g.QualityMetric==string(quality_names{qi}),:);
            [~,idx]=sort(h.Beta_rad_per_mm);
            plot(ax,h.Beta_rad_per_mm(idx),h.MeanQ(idx),'-o', ...
                'DisplayName',quality_names{qi});
        end
        xlabel(ax,'Phase ramp [rad/mm]');ylabel(ax,'Mean quality Q');
        title(ax,[methods{mi} ' / Mean Q']);legend(ax,'Location','best');grid(ax,'on');
        ax=nexttile(tl,mi+2);hold(ax,'on');
        for qi=1:2
            h=g(g.QualityMetric==string(quality_names{qi}),:);
            [~,idx]=sort(h.Beta_rad_per_mm);
            plot(ax,h.Beta_rad_per_mm(idx),h.AURC_over_Random(idx),'-o', ...
                'DisplayName',quality_names{qi});
        end
        yline(ax,1,'k--','Random');grid(ax,'on');
        title(ax,[methods{mi} ' / Risk-coverage']);
        xlabel(ax,'Phase ramp [rad/mm]');ylabel(ax,'AURC / random');
        legend(ax,'Location','best');
    end
    title(tl,['Quality response to phase ramps: ' strrep(p,'_',' ')]);
    exportgraphics(f,fullfile(outdir,'rmi_phase_ramp_trends.png'),'Resolution',160);
end

function plot_risk(C,methods,quality_names,betas,outdir,p)
    f=figure('Color','w','Name',['RMI phase ramp risk-coverage: ' p], ...
        'Position',[40 60 1550 870]);
    tl=tiledlayout(f,2,numel(betas),'TileSpacing','compact','Padding','compact');
    for mi=1:2
        for k=1:numel(betas)
            ax=nexttile(tl,(mi-1)*numel(betas)+k);hold(ax,'on');
            subset=C(C.Region=="Global" & C.Unit=="Pixel" & ...
                C.Estimator==string(methods{mi}) & ...
                abs(C.Beta_rad_per_mm-betas(k))<1e-10,:);
            for qi=1:2
                h=subset(subset.QualityMetric==string(quality_names{qi}),:);
                plot(ax,h.Coverage_pct,h.RiskMAE_mps,'LineWidth',1.3, ...
                    'DisplayName',quality_names{qi});
            end
            yline(ax,subset.RandomRisk_mps(1),'k--','Random');
            title(ax,sprintf('%s: beta=%+.2f rad/mm',methods{mi},betas(k)));
            xlabel(ax,'Accepted coverage [%]');ylabel(ax,'Accepted MAE [m/s]');
            grid(ax,'on');if mi==1&&k==1,legend(ax,'Location','best');end
        end
    end
    title(tl,['Fixed-pixel risk-coverage: ' strrep(p,'_',' ')]);
    exportgraphics(f,fullfile(outdir,'rmi_phase_ramp_risk_coverage.png'), ...
        'Resolution',160);
end
