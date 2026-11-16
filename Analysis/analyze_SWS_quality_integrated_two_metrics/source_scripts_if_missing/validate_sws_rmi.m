function T = validate_sws_rmi(phantom_case, mode)
%VALIDATE_SWS_RMI Choi et al. (2015) RMI, 1D-equation adaptation (Eq. 6).
%
% THIS INTRODUCES A THIRD METRIC (published RMI), NOT AN EXTENSION OF
% PHASE PLANARITY Q.  We implement the 1D wave-equation alternative from
% Choi et al., IEEE IUS 2015, DOI: 10.1109/ULTSYM.2015.0382.
%
% Important: Choi et al. prefer their modified radial model (Eq. 7) for
% their point-push geometry. This code uses Eq. (6) instead because the
% corresponding source geometry and gamma are not established for our
% harmonic k-Wave simulations. Call it "RMI-1D adaptation" in the thesis,
% NOT an exact reproduction of the paper's commercial implementation.
%
% Eq. (6): res = sum_t |u_tt - c_est^2*u_xx|^2 / sum_t |u_xx|^2.
% Eq. (4,5,8): RMI = 0.8*clip((res_max-res)/(res_max-res_min),0,1)
%                  + 0.2*clip((u_rep-u_min)/(u_max-u_min),0,1).
% u_rep = max_t abs(displacement), for each pixel.
%
% c_est = previously reconstructed PG or AC SWS, NEVER ground truth.
% GT is loaded for VALIDATION ONLY, after both RMI maps are computed.
% The paper does not report numerical limits for res/u. We fix 5th/95th
% percentile limits using ONE calibration phantom, without GT, and reuse
% them unchanged for all subsequent phantoms. This is exploratory and
% NOT equivalent to the original empirical calibration.
%
% INSTALL: place this file in <repo>/Analysis/, with prior pilot scripts.
% FIRST RUN (creates calibration):
%    validate_sws_rmi('large_inclusions','calibrate')
% THEN (uses the SAME calibration):
%    validate_sws_rmi('high_contrast')
%    validate_sws_rmi('low_contrast')
%    validate_sws_rmi('chaoran_baseline')
% The initial phantom is in-sample for calibration (not independent test).
%
% INPUTS (read-only):
%   Data/<phantom>_displacement_with_quality.mat (or repo root)
%   Results/sws_measurement_quality/<phantom>/phase_consistency_diagnostic/
%       phase_quality_results.mat
%   Results/method_comparison/batch/<phantom>/method_comparison.mat
% OUTPUTS: Results/sws_measurement_quality/rmi_1d_calibration.mat
%   Results/sws_measurement_quality/<phantom>/rmi_1d_diagnostic/
%       rmi_1d_maps.png, rmi_1d_spearman.png,
%       rmi_1d_summary.csv, rmi_1d_results.mat
%
    if nargin < 1 || isempty(phantom_case)
        phantom_case = 'large_inclusions';
    end
    if nargin < 2 || isempty(mode)
        mode = 'use';
    end
    phantom_case = char(string(phantom_case));
    mode = char(lower(string(mode)));
    assert(ismember(mode, {'calibrate','use'}), ...
        'Mode must be ''calibrate'' or ''use''.');

    here = fileparts(mfilename('fullpath'));
    if isfolder(fullfile(here,'Results')), repo = here;
    else, repo = fileparts(here); end
    result_root = fullfile(repo,'Results','sws_measurement_quality');
    if ~isfolder(result_root), mkdir(result_root); end
    calibration_file = fullfile(result_root,'rmi_1d_calibration.mat');

    % USTB class must be on path BEFORE loading USTB MAT objects.
    if isempty(which('uff.beamformed_data'))
        ustb_root = fullfile(fileparts(repo), 'elastography', ...
            'ElastographyToolbox','USTB');
        if isfolder(ustb_root), addpath(genpath(ustb_root)); rehash; end
    end
    assert(~isempty(which('uff.beamformed_data')), ...
        'USTB class uff.beamformed_data not found on MATLAB path.');

    phase_file = fullfile(result_root, phantom_case, ...
        'phase_consistency_diagnostic','phase_quality_results.mat');
    compare_file = fullfile(repo,'Results','method_comparison', ...
        'batch',phantom_case,'method_comparison.mat');
    assert(isfile(phase_file), 'Missing phase results: %s', phase_file);
    assert(isfile(compare_file), 'Missing method comparison: %s', compare_file);

    A = load(phase_file,'Qold','Qphase','common','PG','AC','x','z');
    B = load(compare_file,'roi');
    needed = {'Qold','Qphase','common','PG','AC','x','z'};
    assert(all(isfield(A,needed)), 'Phase MAT is missing required maps.');
    assert(isfield(B,'roi'), 'Method comparison lacks region masks.');
    x = double(A.x(:)); z = double(A.z(:));
    [Nz,Nx] = size(A.PG);
    assert(Nx==numel(x) && Nz==numel(z), 'Scan axes do not match SWS.');
    assert(isequal(size(A.PG),size(A.AC),size(A.Qold), ...
        size(A.Qphase),size(A.common)), 'Saved map shapes do not match.');
    assert(Nx > 4 && all(diff(x)>0), 'Lateral axis must be increasing.');

    data_name = [phantom_case '_displacement_with_quality.mat'];
    candidates = {fullfile(repo,'Data',data_name), ...
        fullfile(repo,'Data',phantom_case,data_name), ...
        fullfile(repo,data_name), fullfile(repo,'Analysis',data_name)};
    data_file = '';
    for j=1:numel(candidates)
        if isfile(candidates{j}), data_file=candidates{j}; break; end
    end
    assert(~isempty(data_file),'Cannot find %s.',data_name);
    fprintf('\n=== RMI-1D VALIDATION: %s [%s] ===\n',phantom_case,mode);
    fprintf('Reading: %s\n',data_file);
    D = load(data_file,'channel_data','displacement_estimation_autocorr');
    assert(isfield(D,'channel_data') && ...
        isfield(D,'displacement_estimation_autocorr'), ...
        'Missing channel_data/displacement estimation.');
    prf = double(D.channel_data.PRF);
    assert(isscalar(prf) && isfinite(prf) && prf>0,'Invalid PRF.');
    raw = D.displacement_estimation_autocorr.data;
    assert(isnumeric(raw) && isreal(raw) && ...
        mod(numel(raw),Nz*Nx)==0, 'Displacement dimensions invalid.');
    Nt = numel(raw)/(Nz*Nx);
    assert(Nt>=5, 'Insufficient displacement frames.');
    u = reshape(single(raw),Nz,Nx,Nt);
    clear D raw;
    fprintf('PRF %.1f Hz; displacement %d x %d x %d\n',prf,Nz,Nx,Nt);

    % Self-test the numerical wave residual on a known synthetic wave.
    test_wave_operator();

    % --- METRIC CALCULATION: no GT data used in this section. ---
    cPG = double(A.PG);
    cAC = double(A.AC);
    cPG(~isfinite(cPG) | cPG<=0 | cPG>=10) = NaN;
    cAC(~isfinite(cAC) | cAC<=0 | cAC>=10) = NaN;
    [resPG, resAC, urep, denominator] = ...
        wave_residual_1d(u, x, prf, cPG, cAC);
    clear u;

    metric_support = isfinite(resPG) & isfinite(resAC) & ...
        isfinite(urep) & urep>=0;
    assert(nnz(metric_support)>100,'Too few finite RMI residual pixels.');

    if strcmp(mode,'calibrate')
        assert(~isfile(calibration_file), ...
            ['Calibration already exists: %s\n' ...
             'Run with mode ''use'' to avoid overwriting fixed limits.'], ...
            calibration_file);
        rv = [resPG(metric_support);resAC(metric_support)];
        uv = urep(metric_support);
        limits = struct();
        limits.res_min = percentile_simple(rv,5);
        limits.res_max = percentile_simple(rv,95);
        limits.u_min = percentile_simple(uv,5);
        limits.u_max = percentile_simple(uv,95);
        limits.reference_phantom = phantom_case;
        limits.percentiles = [5 95];
        limits.algorithm = 'Choi2015_Eq6_RMI_1D_adaptation';
        limits.wave_frequency_source = 'channel_data.PRF';
        assert(limits.res_max>limits.res_min && ...
            limits.u_max>limits.u_min, ...
            'Calibration ranges collapsed; inspect residual and amplitude.');
        save(calibration_file,'limits');
        fprintf('Saved FIXED calibration: %s\n',calibration_file);
    else
        assert(isfile(calibration_file), ...
            ['No calibration file found. First run:\n' ...
             '  validate_sws_rmi(''large_inclusions'',''calibrate'')']);
        K=load(calibration_file,'limits');
        limits=K.limits;
        assert(strcmp(limits.algorithm,'Choi2015_Eq6_RMI_1D_adaptation'), ...
            'Calibration file was produced by a different algorithm.');
        fprintf('Using fixed calibration from %s\n',limits.reference_phantom);
    end
    fprintf('Residual bounds [%.4g %.4g], amplitude bounds [%.4g %.4g]\n', ...
        limits.res_min,limits.res_max,limits.u_min,limits.u_max);
    RmiPG = make_rmi(resPG,urep,limits);
    RmiAC = make_rmi(resAC,urep,limits);

    % --- VALIDATION ONLY: ground truth enters only below this line. ---
    G = load(compare_file,'GT');
    GT = double(G.GT);
    assert(isequal(size(GT),[Nz Nx]),'GT has incorrect dimensions.');
    assert(all(isfield(B.roi,{'inclusion_core', ...
        'background_core','image_mask'})), 'Required ROI masks missing.');
    [X,Z] = meshgrid(x,z);
    basic = logical(A.common) & logical(B.roi.image_mask) & ...
        isfinite(GT) & isfinite(cPG) & isfinite(cAC) & ...
        isfinite(A.Qold) & isfinite(A.Qphase) & ...
        isfinite(RmiPG) & isfinite(RmiAC);
    assert(nnz(basic)>100, 'Too few common pixels across three Q metrics.');
    edge_mm=3.5;
    boundary = X <= min(x)+edge_mm*1e-3 | ...
               X >= max(x)-edge_mm*1e-3;
    gt_background = abs(GT - min(GT(isfinite(GT))))<1e-5;
    masks = {basic, ...
             basic & logical(B.roi.inclusion_core), ...
             basic & logical(B.roi.background_core) & ~boundary, ...
             basic & logical(B.roi.image_mask) & boundary & gt_background};
    regions = {'Global','InclusionCore', ...
        'BackgroundInterior','LateralBoundary'};

    % RMI depends on estimated c, so each estimator has its own RMI map.
    methods={'PG','AC'};
    c_maps={cPG,cAC};
    rmi_maps={RmiPG,RmiAC};
    qual_labels={'Correlation Q (adapted)', ...
        'Phase Planarity Q (candidate)','RMI-1D (Choi Eq. 6 adaptation)'};
    T=table();
    quality_threshold=0.8;
    block_mm=1;
    bx=floor((X-min(x))/(block_mm*1e-3));
    bz=floor((Z-min(z))/(block_mm*1e-3));
    block_ids = bz*(max(bx(:))+1) + bx + 1;
    for ri=1:numel(regions)
        mask = masks{ri};
        for mi=1:2
            error_abs = abs(c_maps{mi}-GT);
            quality_maps={double(A.Qold),double(A.Qphase),rmi_maps{mi}};
            for qi=1:3
                q=quality_maps{qi}(mask);
                e=error_abs(mask);
                ids=block_ids(mask);
                hi=q>=quality_threshold;
                [qb,eb]=block_means(ids,q,e);
                qmean=mean(q);
                if isempty(q), qmean=NaN; end
                fcr=NaN;
                if any(hi), fcr=mean(e(hi)>1); end % GT-based validation ONLY
                row=table(string(phantom_case),string(regions{ri}), ...
                    string(methods{mi}),string(qual_labels{qi}), ...
                    numel(q),qmean,safe_median(q), ...
                    safe_mean(q>=1-1e-6),safe_spearman(q,e), ...
                    safe_spearman(qb,eb),nnz(hi),nnz(~hi), ...
                    safe_mean(e(hi)),safe_mean(e(~hi)),fcr, ...
                    'VariableNames',{'Phantom','Region','Estimator', ...
                    'QualityMetric','ValidPixels','MeanQ','MedianQ', ...
                    'FractionSaturated','PixelSpearman','BlockSpearman', ...
                    'N_Qge08','N_Qlt08','MAE_Qge08_mps','MAE_Qlt08_mps', ...
                    'FCR_Qge08_ErrorGT1mps'});
                T=[T;row]; %#ok<AGROW>
            end
        end
    end
    disp(T(strcmp(T.Region,'Global'),:));
    fprintf('Common validation pixels: %d (out of %d)\n', ...
        nnz(basic),nnz(A.common));

    outdir = fullfile(result_root,phantom_case,'rmi_1d_diagnostic');
    if ~isfolder(outdir), mkdir(outdir); end
    writetable(T,fullfile(outdir,'rmi_1d_summary.csv'));
    save(fullfile(outdir,'rmi_1d_results.mat'), 'T','limits', ...
        'RmiPG','RmiAC','resPG','resAC','urep','basic','masks', ...
        'regions','x','z','GT','cPG','cAC','prf','-v7.3');

    f=figure('Color','w','Name',['RMI comparison: ' phantom_case]);
    tl=tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
    plot_map(nexttile(tl),x,z,double(A.Qold),[0 1], ...
        'Correlation Q (adapted)');
    plot_map(nexttile(tl),x,z,double(A.Qphase),[0 1], ...
        'Phase Planarity Q (candidate)');
    plot_map(nexttile(tl),x,z,RmiPG,[0 1], 'RMI-1D / PG');
    plot_map(nexttile(tl),x,z,RmiAC,[0 1], 'RMI-1D / AC');
    plot_map(nexttile(tl),x,z,abs(cPG-GT),[0 4], ...
        '|PG - GT| [m/s]');
    plot_map(nexttile(tl),x,z,abs(cAC-GT),[0 4], ...
        '|AC - GT| [m/s]');
    title(tl,['3 quality metrics: ' strrep(phantom_case,'_',' ')]);
    exportgraphics(f,fullfile(outdir,'rmi_1d_maps.png'),'Resolution',160);

    f2=figure('Color','w','Name',['Global Spearman: ' phantom_case]);
    ix=strcmp(T.Region,'Global');
    V=T.PixelSpearman(ix); W=T.BlockSpearman(ix);
    bar([V,W]);
    set(gca,'XTick',1:6, 'XTickLabel', ...
        {'Corr/PG','Phase/PG','RMI/PG','Corr/AC','Phase/AC','RMI/AC'});
    legend('Pixel Spearman','1-mm block Spearman','Location','best');
    ylabel('Spearman(Q, |SWS - GT|)');
    title(['Global association: ' strrep(phantom_case,'_',' ')]);
    grid on; yline(0,'k:');
    exportgraphics(f2,fullfile(outdir,'rmi_1d_spearman.png'), ...
        'Resolution',160);
    fprintf('Results saved: %s\n',outdir);
    fprintf(['Caution: Eq. (6) is a 1D approximation, bounds are fixed ' ...
        'empirical percentiles, and c_est enters RMI (possible shared bias).\n']);
    fprintf(['FCR is a GT-BASED VALIDATION diagnostic, NOT a quality ' ...
        'metric deployable without GT.\n']);
end

function [rp,ra,amp,den] = wave_residual_1d(u,x,prf,cPG,cAC)
% Finite differences of Choi Eq. (6), streaming over temporal frames.
    [Nz,Nx,Nt]=size(u);
    numPG=zeros(Nz,Nx); numAC=zeros(Nz,Nx);
    den=zeros(Nz,Nx);
    amp=max(abs(double(u(:,:,1))),abs(double(u(:,:,end))));
    dt2=prf^2;
    dxr = x(3:end)-x(2:end-1);
    dxl = x(2:end-1)-x(1:end-2);
    dxspan = x(3:end)-x(1:end-2);
    cp2=cPG.^2; ca2=cAC.^2;
    for it=2:Nt-1
        u0=double(u(:,:,it));
        up=double(u(:,:,it+1));
        um=double(u(:,:,it-1));
        utt=(up-2*u0+um)*dt2;
        ux2=nan(Nz,Nx);
        fwd=bsxfun(@rdivide,u0(:,3:end)-u0(:,2:end-1),dxr');
        bwd=bsxfun(@rdivide,u0(:,2:end-1)-u0(:,1:end-2),dxl');
        ux2(:,2:end-1)=bsxfun(@rdivide,2*(fwd-bwd),dxspan');
        residPG=utt(:,2:end-1)-cp2(:,2:end-1).*ux2(:,2:end-1);
        residAC=utt(:,2:end-1)-ca2(:,2:end-1).*ux2(:,2:end-1);
        numPG(:,2:end-1)=numPG(:,2:end-1)+residPG.^2;
        numAC(:,2:end-1)=numAC(:,2:end-1)+residAC.^2;
        den(:,2:end-1)=den(:,2:end-1)+ux2(:,2:end-1).^2;
        amp=max(amp,abs(u0));
    end
    good=den(isfinite(den)&den>0);
    assert(~isempty(good),'Zero spatial curvature throughout displacement field.');
    floor_den=median(good)*1e-12;
    rp=numPG./den;
    ra=numAC./den;
    rp(:,[1,end])=NaN; ra(:,[1,end])=NaN;
    bad=~isfinite(den) | den<=floor_den;
    rp(bad)=NaN;ra(bad)=NaN;
end

function q = make_rmi(res,amp,limits)
    r=(limits.res_max-res)/(limits.res_max-limits.res_min);
    a=(amp-limits.u_min)/(limits.u_max-limits.u_min);
    r=min(1,max(0,r));
    a=min(1,max(0,a));
    q=0.8*r+0.2*a;
    q(~isfinite(res)|~isfinite(amp))=NaN;
end

function p = percentile_simple(v, pct)
% Toolbox-independent sample percentile via deterministic sorted quantile.
    v=double(v(isfinite(v)));
    assert(numel(v)>=100,'Not enough values for percentile calibration.');
    v=sort(v(:)); n=numel(v);
    k=1+(n-1)*(pct/100);
    i=floor(k); j=ceil(k);
    p=v(i)+(k-i)*(v(j)-v(i));
end

function [qb,eb] = block_means(ids,q,e)
    if isempty(ids),qb=[];eb=[];return;end
    [~,~,g]=unique(double(ids(:)));
    qb=accumarray(g,double(q(:)),[],@mean);
    eb=accumarray(g,double(e(:)),[],@mean);
end

function r = safe_spearman(a,b)
    ok=isfinite(a)&isfinite(b);
    a=double(a(ok)); b=double(b(ok));
    if numel(a)<5 || numel(unique(a))<3 || numel(unique(b))<3
        r=NaN;return;
    end
    r=corr(a(:),b(:),'Type','Spearman');
end

function m = safe_mean(v)
    v=v(isfinite(v));
    if isempty(v),m=NaN;else,m=mean(v);end
end

function m = safe_median(v)
    v=v(isfinite(v));
    if isempty(v),m=NaN;else,m=median(v);end
end

function plot_map(ax,x,z,im,lims,ttl)
    imagesc(ax,x*1e3,z*1e3,im);
    set(ax,'YDir','normal');axis(ax,'image');
    clim(ax,lims);colorbar(ax);
    title(ax,ttl,'Interpreter','none');
    xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
end

function test_wave_operator()
% A known sinusoidal wave should satisfy u_tt = c^2 * u_xx better
% for the correct c than an intentionally incorrect c.
    f0=100; prf=3000; c=3; dx=0.0004;
    omega=2*pi*f0; k=omega/c;
    x=(-40:40)*dx; t=(0:80)/prf;
    u=cos(omega*t - k*x');
    tt=(u(40,43)-2*u(40,42)+u(40,41))*prf^2;
    xx=(u(41,42)-2*u(40,42)+u(39,42))/dx^2;
    err_true=abs(tt-c^2*xx);
    err_wrong=abs(tt-(c*1.5)^2*xx);
    assert(err_true < err_wrong*0.1, ...
        'Synthetic wave-equation sanity check failed.');
end
