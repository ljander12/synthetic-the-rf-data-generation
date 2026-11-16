function T = diagnose_rmi_residual_saturation(phantom_case)
%DIAGNOSE_RMI_RESIDUAL_SATURATION Inspect clipping in existing RMI-1D.
%
% This is a DIAGNOSTIC of quality metric RMI-1D, NOT a new quality metric.
% Reuses the EXACT finite-difference / denominator and harmonic phase-ramp
% convention from analyze_rmi_speed_dependence.m (Choi et al. 2015 Eq. 6
% adaptation). Fixed calibration is loaded, not refitted.
%
% Run AFTER analyze_rmi_speed_dependence.m and its prerequisites:
%   T = diagnose_rmi_residual_saturation('high_contrast');
%   T = diagnose_rmi_residual_saturation('large_inclusions');
%
% Expected outputs in Results/sws_measurement_quality/<case>/
%                         rmi_residual_saturation_diagnostic/
%   residual_saturation_summary.csv
%   residual_saturation_trends.png
%   residual_saturation_maps_pg.png
%   residual_saturation_maps_ac.png
%   residual_saturation_results.mat
%
% IMPORTANT: The 'res' below is Choi's normalized PDE residual (Eq. 6),
% BEFORE its affine mapping/clipping into the Q_res score. It is NOT the
% raw, unnormalized PDE numerator. Two versions: adaptive c_beta and frozen c0.
% Ground truth is NOT used to calculate residual or quality scores.

    if nargin<1 || isempty(phantom_case)
        phantom_case='high_contrast';
    end
    phantom_case=char(string(phantom_case));
    here=fileparts(mfilename('fullpath'));
    if isfolder(fullfile(here,'Results'))
        repo=here;
    else
        repo=fileparts(here);
    end
    top=fullfile(repo,'Results','sws_measurement_quality');
    pdir=fullfile(top,phantom_case);
    src=fullfile(pdir,'phase_ramp_sws_validation','phase_ramp_metrics.mat');
    base=fullfile(pdir,'rmi_1d_diagnostic','rmi_1d_results.mat');
    cal=fullfile(top,'rmi_1d_calibration.mat');
    assert(isfile(src),'Missing prerequisite: %s',src);
    assert(isfile(base),'Missing prerequisite: %s',base);
    assert(isfile(cal),'Missing fixed calibration: %s',cal);
    S=load(src,'maps','GT','common','betas','prf','f0','x','z');
    R=load(base,'RmiPG','RmiAC','basic','GT','x','z');
    K=load(cal,'limits');
    lim=K.limits;
    assert(strcmp(lim.algorithm,'Choi2015_Eq6_RMI_1D_adaptation'), ...
        'Incorrect RMI calibration algorithm.');
    assert(lim.res_min<lim.res_max && lim.u_min<lim.u_max, ...
        'Invalid calibration bounds.');

    betas=double(S.betas(:));
    k0=find(abs(betas)<1e-12);
    assert(numel(k0)==1 && numel(S.maps)==numel(betas), ...
        'Expected exactly one beta=0 and matching SWS maps.');
    [Nz,Nx]=size(S.GT);
    x=double(S.x(:));z=double(S.z(:));
    assert(numel(x)==Nx && numel(z)==Nz && all(diff(x)>0) && Nx>=5, ...
        'Scan axes or dimensions are invalid.');
    assert(isequal(size(S.common),[Nz Nx],size(R.GT),size(R.basic)), ...
        'Different saved masks/shapes.');
    assert(max(abs(double(S.GT(:))-double(R.GT(:))))<1e-6, ...
        'Ground truth differs between prerequisite datasets.');
    assert(max(abs(x-double(R.x(:))))<1e-9 && ...
        max(abs(z-double(R.z(:))))<1e-9, 'Scan axes mismatch.');
    for k=1:numel(betas)
        assert(abs(double(S.maps(k).beta)-betas(k))<1e-12, ...
            'Beta index mismatch.');
        assert(isequal(size(S.maps(k).PG),[Nz Nx], ...
            size(S.maps(k).AC),size(S.maps(k).Q)), ...
            'Saved SWS map shape mismatch.');
    end

    if isempty(which('uff.beamformed_data'))
        ustb=fullfile(fileparts(repo),'elastography','ElastographyToolbox','USTB');
        if isfolder(ustb),addpath(genpath(ustb));rehash;end
    end
    assert(~isempty(which('uff.beamformed_data')), ...
        'USTB not on MATLAB path. Add the USTB folder before running.');
    fn=[phantom_case '_displacement_with_quality.mat'];
    candidates={fullfile(repo,'Data',fn), ...
        fullfile(repo,'Data',phantom_case,fn),fullfile(repo,fn), ...
        fullfile(repo,'Analysis',fn)};
    datafile='';
    for j=1:numel(candidates)
        if isfile(candidates{j}),datafile=candidates{j};break;end
    end
    assert(~isempty(datafile),'Cannot find displacement data: %s',fn);
    D=load(datafile,'channel_data','displacement_estimation_autocorr', ...
        'shear_params');
    prf=double(D.channel_data.PRF);
    f0=double(D.shear_params.source_freq);
    assert(isscalar(prf) && prf>0 && isfinite(prf) && ...
        abs(prf-double(S.prf))<1e-8 && abs(f0-double(S.f0))<1e-8, ...
        'Saved PRF/frequency mismatch.');
    raw=D.displacement_estimation_autocorr.data;
    assert(isnumeric(raw) && isreal(raw) && mod(numel(raw),Nz*Nx)==0, ...
        'Unexpected displacement data dimensions.');
    Nt=numel(raw)/(Nz*Nx);
    assert(Nt>=5,'Insufficient displacement frames.');
    u=-reshape(single(raw),Nz,Nx,Nt); % SAME sign as original script
    clear D raw;

    % Same harmonic least-squares fit as analyze_rmi_speed_dependence.m.
    tt=(0:Nt-1)'/prf;
    cos_t=cos(2*pi*f0*tt);
    sin_t=sin(2*pi*f0*tt);
    weights=pinv([cos_t sin_t ones(Nt,1)]);
    a=zeros(Nz,Nx);b=zeros(Nz,Nx);
    for it=1:Nt
        frame=double(u(:,:,it));
        a=a+weights(1,it)*frame;
        b=b+weights(2,it)*frame;
    end
    [Xmm,~]=meshgrid(x*1e3,z*1e3);
    methods={'PG','AC'};
    frozenC=cell(1,2);
    for mi=1:2
        v=double(S.maps(k0).(methods{mi}));
        v(~isfinite(v)|v<=0|v>=10)=NaN;
        frozenC{mi}=v;
    end

    nB=numel(betas);
    resAdaptive=cell(nB,2);
    resFrozen=cell(nB,2);
    qAdaptive=cell(nB,2);
    qFrozen=cell(nB,2);
    ampScores=cell(nB,1);
    fprintf('\n=== RMI RESIDUAL SATURATION: %s ===\n',phantom_case);
    fprintf('Source: %s\n',datafile);
    fprintf('PRF = %.1f Hz; f0 = %.1f Hz; size = %dx%dx%d\n', ...
        prf,f0,Nz,Nx,Nt);
    fprintf('Fixed res_min = %.8g, res_max = %.8g (%s)\n', ...
        lim.res_min,lim.res_max,lim.reference_phantom);
    for k=1:nB
        beta=betas(k);
        theta=beta*(Xmm-mean(x)*1e3);
        da=a.*(cos(theta)-1)+b.*sin(theta);
        db=b.*(cos(theta)-1)-a.*sin(theta);
        adaptiveC=cell(1,2);
        for mi=1:2
            v=double(S.maps(k).(methods{mi}));
            v(~isfinite(v)|v<=0|v>=10)=NaN;
            adaptiveC{mi}=v;
        end
        fprintf('Residual for beta=%+.2f rad/mm ...\n',beta);
        [ra,rf,amp]=same_wave_residual(u,da,db,cos_t,sin_t, ...
            x,prf,adaptiveC,frozenC);
        ampScores{k}=clip01((amp-lim.u_min)/(lim.u_max-lim.u_min));
        for mi=1:2
            resAdaptive{k,mi}=ra{mi};
            resFrozen{k,mi}=rf{mi};
            qAdaptive{k,mi}=clip01((lim.res_max-ra{mi})/ ...
                (lim.res_max-lim.res_min));
            qFrozen{k,mi}=clip01((lim.res_max-rf{mi})/ ...
                (lim.res_max-lim.res_min));
        end
        clear da db ra rf amp;
    end
    clear u a b;

    % SAME pixels for all betas, speed variants, and both estimators.
    valid=logical(S.common)&logical(R.basic);
    for k=1:nB
        valid=valid & isfinite(S.maps(k).PG) & ...
            (S.maps(k).PG>0) & (S.maps(k).PG<10) & ...
            isfinite(S.maps(k).AC) & (S.maps(k).AC>0) & ...
            (S.maps(k).AC<10) & isfinite(S.maps(k).Q) & ...
            isfinite(ampScores{k});
        for mi=1:2
            valid=valid & isfinite(resAdaptive{k,mi}) & ...
                isfinite(resFrozen{k,mi});
        end
    end
    assert(nnz(valid)>100,'Too few common valid pixels.');
    fprintf('Identical mask for ALL comparisons: %d pixels\n',nnz(valid));

    % Reproduce existing saved baseline RMI (checks exact implementation).
    for mi=1:2
        zeroDifference=max(abs(resAdaptive{k0,mi}(valid)- ...
            resFrozen{k0,mi}(valid)));
        assert(zeroDifference<1e-10,'Frozen/adaptive baseline mismatch.');
        original=double(R.(['Rmi' methods{mi}]));
        chk=valid & isfinite(original);
        assert(nnz(chk)>100,'Saved RMI has no valid overlap.');
        local=0.8*qAdaptive{k0,mi}+0.2*ampScores{k0};
        diff0=max(abs(local(chk)-original(chk)));
        fprintf('%s check against original baseline RMI: %.6g max difference\n', ...
            methods{mi},diff0);
        assert(diff0<2e-3, ...
            'Reconstruction differs from saved RMI. Do not interpret!');
    end

    T=table();
    for mi=1:2
        for k=1:nB
            for mode=1:2
                if mode==1
                    type='Adaptive';
                    r=double(resAdaptive{k,mi}(valid));
                    q=double(qAdaptive{k,mi}(valid));
                    baseR=double(resAdaptive{k0,mi}(valid));
                    baseQ=double(qAdaptive{k0,mi}(valid));
                else
                    type='Frozen';
                    r=double(resFrozen{k,mi}(valid));
                    q=double(qFrozen{k,mi}(valid));
                    baseR=double(resFrozen{k0,mi}(valid));
                    baseQ=double(qFrozen{k0,mi}(valid));
                end
                % The 'res' in Eq. 6 is already a denominator-normalized
                % wave-equation residual; these are PRE-clipping values.
                row=table(string(phantom_case),string(methods{mi}), ...
                    string(type),betas(k),nnz(valid), ...
                    median(r),prctile(r,10),prctile(r,90), ...
                    mean(r),mean(abs(r-baseR)), ...
                    mean(q),median(q),mean(abs(q-baseQ)), ...
                    mean(r<=lim.res_min),mean(r>=lim.res_max), ...
                    mean(q>=1-1e-12),mean(q<=1e-12), ...
                    'VariableNames',{'Phantom','Estimator','SpeedMode', ...
                    'Beta_rad_per_mm','CommonPixels','MedianResidual', ...
                    'P10Residual','P90Residual','MeanResidual', ...
                    'MeanAbsDeltaResidual','MeanQres','MedianQres', ...
                    'MeanAbsDeltaQres','FracResidualBelowMin', ...
                    'FracResidualAboveMax','FracQresAt1','FracQresAt0'});
                T=[T;row]; %#ok<AGROW>
            end
        end
    end

    outdir=fullfile(pdir,'rmi_residual_saturation_diagnostic');
    if ~isfolder(outdir),mkdir(outdir);end
    writetable(T,fullfile(outdir,'residual_saturation_summary.csv'));
    save(fullfile(outdir,'residual_saturation_results.mat'), ...
        'T','betas','valid','lim','resAdaptive','resFrozen', ...
        'qAdaptive','qFrozen','-v7.3');
    plot_summary(T,methods,betas,outdir,phantom_case);
    for mi=1:2
        plot_maps(resFrozen(:,mi),qFrozen(:,mi),valid,betas, ...
            x,z,lim,outdir,phantom_case,methods{mi});
    end

    fprintf('\nResidual saturation summary (same fixed pixels):\n');
    disp(T(:,{'Estimator','SpeedMode','Beta_rad_per_mm', ...
        'MedianResidual','P90Residual','MeanAbsDeltaResidual', ...
        'MeanQres','FracQresAt1'}));
    fprintf('Saved to %s\n',outdir);
    fprintf(['Interpretation: A changing PRE-clip residual alongside an ' ...
        'unchanged clipped Q with high FracQresAt1 supports clipping. ' ...
        'If raw residual is also unchanged, clipping alone does NOT ' ...
        'explain the insensitivity.\n']);
end

function [rA,rF,amp]=same_wave_residual(u,da,db,cos_t,sin_t,x,prf,cBeta,c0)
% IDENTICAL finite-difference function to analyze_rmi_speed_dependence.m.
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

function q=clip01(x)
    q=max(0,min(1,x));
end

function plot_summary(T,methods,betas,outdir,p)
    f=figure('Color','w','Name',['RMI saturation trends: ' p], ...
        'Position',[30 40 1450 860]);
    tl=tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
    for mi=1:2
        g=T(T.Estimator==string(methods{mi}),:);
        a=g(g.SpeedMode=="Adaptive",:);
        fr=g(g.SpeedMode=="Frozen",:);
        [~,ia]=sort(a.Beta_rad_per_mm);a=a(ia,:);
        [~,iff]=sort(fr.Beta_rad_per_mm);fr=fr(iff,:);
        ax=nexttile(tl,(mi-1)*3+1);hold(ax,'on');
        plot(ax,a.Beta_rad_per_mm,a.MedianResidual,'-o','LineWidth',1.6, ...
            'DisplayName','Adaptive c_beta');
        plot(ax,fr.Beta_rad_per_mm,fr.MedianResidual,'-s','LineWidth',1.6, ...
            'DisplayName','Frozen c_0');
        yline(ax,0,':');grid(ax,'on');
        ylabel(ax,'Median residual (Eq. 6; before clipping)');
        title(ax,[methods{mi} ' / residual']);
        legend(ax,'Location','best');
        ax=nexttile(tl,(mi-1)*3+2);hold(ax,'on');
        plot(ax,a.Beta_rad_per_mm,100*a.FracQresAt1,'-o', ...
            'LineWidth',1.6,'DisplayName','Adaptive');
        plot(ax,fr.Beta_rad_per_mm,100*fr.FracQresAt1,'-s', ...
            'LineWidth',1.6,'DisplayName','Frozen');
        ylabel(ax,'Fraction Q_{res}=1 [%]');ylim(ax,[0 100]);
        title(ax,[methods{mi} ' / saturation']);grid(ax,'on');
        ax=nexttile(tl,(mi-1)*3+3);hold(ax,'on');
        plot(ax,a.Beta_rad_per_mm,a.MeanQres,'-o', ...
            'LineWidth',1.6,'DisplayName','Adaptive');
        plot(ax,fr.Beta_rad_per_mm,fr.MeanQres,'-s', ...
            'LineWidth',1.6,'DisplayName','Frozen');
        ylabel(ax,'Mean clipped Q_{res}');ylim(ax,[0 1]);
        title(ax,[methods{mi} ' / clipped score']);grid(ax,'on');
        for i=1:3
            ax=nexttile(tl,(mi-1)*3+i);
            xlabel(ax,'Phase ramp beta [rad/mm]');
            xlim(ax,[min(betas)-0.025 max(betas)+0.025]);
        end
    end
    title(tl,['Residual vs clipping: ' strrep(p,'_',' ')]);
    exportgraphics(f,fullfile(outdir,'residual_saturation_trends.png'), ...
        'Resolution',160);
end

function plot_maps(resCell,qCell,valid,betas,x,z,lim,outdir,p,method)
    nB=numel(betas);
    f=figure('Color','w','Name',[method ' RMI saturation maps'], ...
        'Position',[25 30 1450 1090]);
    tl=tiledlayout(f,3,nB,'TileSpacing','compact','Padding','compact');
    % Shared color limits to avoid disguising changes between ramp cases.
    vals=[];
    for k=1:nB
        r=double(resCell{k}(valid));
        vals=[vals; log10(max(r,rmin()))]; %#ok<AGROW>
    end
    bounds=prctile(vals,[2 98]);
    if bounds(1)>=bounds(2),bounds(2)=bounds(1)+1;end
    for k=1:nB
        raw=double(resCell{k});q=double(qCell{k});
        raw(~valid)=NaN;q(~valid)=NaN;
        clips=nan(size(raw));
        clips(valid)=0;
        clips(valid & raw<=lim.res_min)=1;
        clips(valid & raw>=lim.res_max)=-1;
        ax=nexttile(tl,k);
        imagesc(ax,x*1e3,z*1e3,log10(max(raw,rmin())));
        set(ax,'YDir','normal');axis(ax,'image');clim(ax,bounds);colorbar(ax);
        title(ax,sprintf('Eq. 6 residual log10: beta=%+.2f',betas(k)));
        ax=nexttile(tl,k+nB);
        imagesc(ax,x*1e3,z*1e3,q);
        set(ax,'YDir','normal');axis(ax,'image');clim(ax,[0 1]);colorbar(ax);
        title(ax,sprintf('Clipped Q_{res}: beta=%+.2f',betas(k)));
        ax=nexttile(tl,k+2*nB);
        imagesc(ax,x*1e3,z*1e3,clips);
        set(ax,'YDir','normal');axis(ax,'image');clim(ax,[-1 1]);colorbar(ax);
        title(ax,sprintf('Clip state: beta=%+.2f',betas(k)));
        for ri=0:2
            ax=nexttile(tl,k+ri*nB);
            xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
        end
    end
    title(tl,[method ' / Frozen c0: ' strrep(p,'_',' ') ...
        ' | Clip: +1=Qres at 1, -1=Qres at 0, 0=not clipped']);
    exportgraphics(f,fullfile(outdir, ...
        ['residual_saturation_maps_' lower(method) '.png']), ...
        'Resolution',150);
end

function m=rmin()
% Positive floor for log10 visualization only; no effect on metrics.
    m=realmin('double');
end
