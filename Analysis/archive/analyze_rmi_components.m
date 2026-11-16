function T = analyze_rmi_components(phantom_case)
%ANALYZE_RMI_COMPONENTS Explain RMI's response to imposed phase ramps.
%
% THIS IS A VALIDATION/DIAGNOSTIC OF METRIC 3 (RMI-1D), NOT A NEW METRIC.
% Uses the SAME wave-equation residual, empirical 5/95 calibration bounds,
% phase-ramped displacement, PG/AC maps, and RMI implementation as
% validate_rmi_phase_ramp.m and validate_sws_rmi.m.
%
% RMI = 0.8*Q_res + 0.2*Q_amp
% Q_res = clip((res_max-res)/(res_max-res_min),0,1)
% Q_amp = clip((u_rep-u_min)/(u_max-u_min),0,1)
% Reference: Choi et al. 2015, DOI 10.1109/ULTSYM.2015.0382,
% Equations (4), (5), (6), (8); Eq. (6) is a 1D ADAPTATION.
%
% Q_res, Q_amp and RMI are ground-truth-independent calculations.
% Ground truth enters only error and regional validation statistics.
% Both components are CALIBRATION DEPENDENT. No parameter refitting here.
%
% Run (after existing RMI phase-ramp validation):
%   analyze_rmi_components('large_inclusions')
%   analyze_rmi_components('high_contrast')
%
% Place in the project's Analysis/ folder. The original scripts/files
% are NOT modified. Source phase-ramp reconstruction is reused; the RMI
% residual must be recomputed from displacement (not SWS alone).
%
% Outputs to Results/sws_measurement_quality/<case>/rmi_component_diagnostic/:
%   rmi_component_summary.csv
%   rmi_component_trends.png
%   rmi_component_maps_pg.png
%   rmi_component_maps_ac.png
%   rmi_component_results.mat
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

    fprintf('\n=== RMI COMPONENT ANALYSIS | PHASE RAMPS: %s ===\n',phantom_case);
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
    nBeta=numel(betas);
    rmaps=cell(nBeta,2);
    residualScoreMaps=cell(nBeta,2); % normalized 0..1, unweighted
    amplitudeScoreMaps=cell(nBeta,1); % same for PG and AC
    % Reconstruct identical phase-ramped DISPLACEMENT and RMI components.
    for k=1:nBeta
        beta=betas(k);
        theta=beta*(Xmm-mean(x)*1e3);
        da=a.*(cos(theta)-1)+b.*sin(theta);
        db=b.*(cos(theta)-1)-a.*sin(theta);
        fprintf('Decomposing beta=%+.2f rad/mm...\n',beta);
        cPG=double(S.maps(k).PG);
        cAC=double(S.maps(k).AC);
        cPG(~isfinite(cPG)|cPG<=0|cPG>=10)=NaN;
        cAC(~isfinite(cAC)|cAC<=0|cAC>=10)=NaN;
        [rp,ra,amp]=residual_ramped_1d(u,da,db,cos_t,sin_t,x,prf,cPG,cAC);
        amplitudeScoreMaps{k}=clip01((amp-lim.u_min)/(lim.u_max-lim.u_min));
        residualScoreMaps{k,1}=clip01((lim.res_max-rp)/(lim.res_max-lim.res_min));
        residualScoreMaps{k,2}=clip01((lim.res_max-ra)/(lim.res_max-lim.res_min));
        for mi=1:2
            rr=residualScoreMaps{k,mi};
            aa=amplitudeScoreMaps{k};
            q=0.8*rr+0.2*aa;
            q(~isfinite(rr)|~isfinite(aa))=NaN;
            rmaps{k,mi}=q;
        end
        clear da db rp ra amp;
    end
    clear u a b;

    k0=find(betas==0,1);
    baseline_mask=logical(R.basic) & ...
        isfinite(rmaps{k0,1}) & isfinite(rmaps{k0,2});
    assert(nnz(baseline_mask)>100,'Too few baseline RMI pixels.');
    refPG=double(R.RmiPG); refAC=double(R.RmiAC);
    bdiffPG=max(abs(rmaps{k0,1}(baseline_mask)-refPG(baseline_mask)));
    bdiffAC=max(abs(rmaps{k0,2}(baseline_mask)-refAC(baseline_mask)));
    fprintf('Baseline RMI agreement: PG max %.5g; AC max %.5g\n',bdiffPG,bdiffAC);
    assert(bdiffPG<2e-3 && bdiffAC<2e-3, ...
        ['Cannot reproduce original RMI baseline. ' ...
         'Verify calibration, displacement conventions and MATLAB versions.']);

    valid=logical(S.common) & logical(R.basic) & isfinite(double(S.GT));
    for k=1:nBeta
        M=S.maps(k);
        valid=valid & isfinite(M.Q) & ...
            isfinite(M.PG) & M.PG>0 & M.PG<10 & ...
            isfinite(M.AC) & M.AC>0 & M.AC<10 & ...
            isfinite(rmaps{k,1}) & isfinite(rmaps{k,2}) & ...
            isfinite(amplitudeScoreMaps{k});
    end
    fprintf('Fixed common pixels: %d\n',nnz(valid));
    assert(nnz(valid)>100, 'Too few pixels across all phase ramps.');

    % Optional check against saved RMI phase-ramp maps, for exact matching.
    priorFile=fullfile(pdir,'rmi_phase_ramp_validation', ...
        'rmi_phase_ramp_results.mat');
    if isfile(priorFile)
        Prior=load(priorFile,'rmaps');
        if isfield(Prior,'rmaps') && isequal(size(Prior.rmaps),size(rmaps))
            for k=1:nBeta
                for mi=1:2
                    oldMap=double(Prior.rmaps{k,mi});
                    discrepancy=max(abs(oldMap(valid)-rmaps{k,mi}(valid)));
                    assert(discrepancy<2e-3, ...
                        'RMI discrepancy with prior phase-ramp test (beta index %d, method %d).',k,mi);
                end
            end
            fprintf('Verified component sum against previous phase-ramp maps.\n');
        end
    end

    methods={'PG','AC'};
    regions=R.regions;
    T=table();
    for ri=1:numel(regions)
        regionMask=logical(R.masks{ri}) & valid;
        if nnz(regionMask)<30
            warning('Skipping region %s: <30 common pixels.',regions{ri});
            continue;
        end
        for mi=1:2
            for k=1:nBeta
                r=double(residualScoreMaps{k,mi}(regionMask));
                aamp=double(amplitudeScoreMaps{k}(regionMask));
                rmi=double(rmaps{k,mi}(regionMask));
                speed=double(S.maps(k).(methods{mi}));
                err=speed(regionMask)-double(S.GT(regionMask));
                qph=double(S.maps(k).Q(regionMask));
                meanR=mean(r); meanA=mean(aamp);
                meanQ=mean(rmi);
                assert(abs(meanQ-(0.8*meanR+0.2*meanA))<1e-8, ...
                    'RMI weighted component sum failed.');
                row=table(string(phantom_case),string(regions{ri}), ...
                    string(methods{mi}),betas(k),nnz(regionMask), ...
                    meanR,meanA,0.8*meanR,0.2*meanA,meanQ, ...
                    mean(qph),sqrt(mean(err.^2)),mean(abs(err)), ...
                    mean(r<=1e-12),mean(r>=1-1e-12), ...
                    mean(aamp<=1e-12),mean(aamp>=1-1e-12), ...
                    'VariableNames',{'Phantom','Region','Estimator', ...
                    'Beta_rad_per_mm','N_Pixels', ...
                    'Mean_Qres','Mean_Qamp', ...
                    'WeightedRes_0p8','WeightedAmp_0p2','Mean_RMI', ...
                    'Mean_PhaseQ','RMSE_mps','MAE_mps', ...
                    'FracResClipped0','FracResClipped1', ...
                    'FracAmpClipped0','FracAmpClipped1'});
                T=[T;row]; %#ok<AGROW>
            end
        end
    end
    T.Delta_WeightedRes=nan(height(T),1);
    T.Delta_WeightedAmp=nan(height(T),1);
    T.Delta_MeanRMI=nan(height(T),1);
    T.Delta_RMSE_mps=nan(height(T),1);
    for i=1:height(T)
        same=T.Region==T.Region(i) & T.Estimator==T.Estimator(i) & ...
            T.Beta_rad_per_mm==0;
        assert(nnz(same)==1,'Expected exactly one beta=0 per region and estimator.');
        T.Delta_WeightedRes(i)=T.WeightedRes_0p8(i)-T.WeightedRes_0p8(same);
        T.Delta_WeightedAmp(i)=T.WeightedAmp_0p2(i)-T.WeightedAmp_0p2(same);
        T.Delta_MeanRMI(i)=T.Mean_RMI(i)-T.Mean_RMI(same);
        T.Delta_RMSE_mps(i)=T.RMSE_mps(i)-T.RMSE_mps(same);
    end
    assert(all(abs(T.Delta_WeightedRes+T.Delta_WeightedAmp- ...
        T.Delta_MeanRMI)<1e-8),'Delta-RMI decomposition failed.');

    outdir=fullfile(pdir,'rmi_component_diagnostic');
    if ~isfolder(outdir),mkdir(outdir);end
    writetable(T,fullfile(outdir,'rmi_component_summary.csv'));
    save(fullfile(outdir,'rmi_component_results.mat'), ...
        'T','betas','valid','lim','rmaps','residualScoreMaps', ...
        'amplitudeScoreMaps','bdiffPG','bdiffAC','-v7.3');
    plot_component_trends(T,betas,methods,outdir,phantom_case);
    plot_component_maps(residualScoreMaps,amplitudeScoreMaps,betas, ...
        methods,x,z,outdir,phantom_case);

    fprintf('\nGLOBAL COMPONENT DECOMPOSITION (weighted scores):\n');
    disp(T(T.Region=="Global", ...
        {'Beta_rad_per_mm','Estimator','RMSE_mps','Mean_Qres', ...
         'Mean_Qamp','Mean_RMI','Delta_WeightedRes', ...
         'Delta_WeightedAmp','Delta_MeanRMI','FracResClipped0'}));
    fprintf('Saved to: %s\n',outdir);
    fprintf(['RMI-1D Eq. (6) adaptation, not full Eq. (7). ' ...
        'GT is used only for evaluating error/regions.\n']);
end

function a=clip01(a)
    a=max(0,min(1,a));
end

function plot_component_trends(T,betas,methods,outdir,p)
    f=figure('Color','w','Name',['RMI components: ' p], ...
        'Position',[40 50 1300 820]);
    tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
    for mi=1:2
        g=T(T.Region=="Global" & T.Estimator==string(methods{mi}),:);
        [~,idx]=sort(g.Beta_rad_per_mm);
        g=g(idx,:);
        ax=nexttile(tl,mi);hold(ax,'on');
        plot(ax,g.Beta_rad_per_mm,g.Mean_Qres,'-o','LineWidth',1.4, ...
            'DisplayName','Q_{res} (unweighted)');
        plot(ax,g.Beta_rad_per_mm,g.Mean_Qamp,'-s','LineWidth',1.4, ...
            'DisplayName','Q_{amp} (unweighted)');
        plot(ax,g.Beta_rad_per_mm,g.Mean_RMI,'-d','LineWidth',1.6, ...
            'DisplayName','RMI = 0.8 Q_{res} + 0.2 Q_{amp}');
        xlabel(ax,'Phase ramp [rad/mm]');ylabel(ax,'Mean normalized score');
        title(ax,[methods{mi} ' / component scores']);
        ylim(ax,[0 1]);grid(ax,'on');legend(ax,'Location','best');

        ax=nexttile(tl,mi+2);hold(ax,'on');
        bar(ax,g.Beta_rad_per_mm, ...
            [g.Delta_WeightedRes g.Delta_WeightedAmp], ...
            'grouped');
        plot(ax,g.Beta_rad_per_mm,g.Delta_MeanRMI,'k-o', ...
            'LineWidth',1.4,'MarkerFaceColor','k');
        yline(ax,0,'k:');
        xlabel(ax,'Phase ramp [rad/mm]');ylabel(ax,'Change relative to beta=0');
        title(ax,[methods{mi} ' / contributions to Delta RMI']);
        grid(ax,'on');legend(ax,{'0.8 Delta Q_{res}', ...
            '0.2 Delta Q_{amp}','Delta RMI'},'Location','best');
    end
    title(tl,['Why does RMI change? ' strrep(p,'_',' ')]);
    exportgraphics(f,fullfile(outdir,'rmi_component_trends.png'), ...
        'Resolution',170);
end

function plot_component_maps(resMaps,ampMaps,betas,methods,x,z,outdir,p)
    for mi=1:2
        f=figure('Color','w','Name',[methods{mi} ' RMI component maps'], ...
            'Position',[40 60 1450 750]);
        tl=tiledlayout(f,2,numel(betas),'TileSpacing','compact', ...
            'Padding','compact');
        for k=1:numel(betas)
            ax=nexttile(tl,k);
            imagesc(ax,x*1e3,z*1e3,resMaps{k,mi});
            set(ax,'YDir','normal');axis(ax,'image');clim(ax,[0 1]);colorbar(ax);
            title(ax,sprintf('Qres: beta=%+.2f rad/mm',betas(k)));
            xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
            ax=nexttile(tl,k+numel(betas));
            imagesc(ax,x*1e3,z*1e3,ampMaps{k});
            set(ax,'YDir','normal');axis(ax,'image');clim(ax,[0 1]);colorbar(ax);
            title(ax,sprintf('Qamp: beta=%+.2f rad/mm',betas(k)));
            xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
        end
        title(tl,[methods{mi} ' / unweighted RMI component maps: ' ...
            strrep(p,'_',' ')]);
        exportgraphics(f,fullfile(outdir, ...
            ['rmi_component_maps_' lower(methods{mi}) '.png']), ...
            'Resolution',170);
    end
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

