function T = validate_phase_ramp_sws(phantom_case, mode)
%VALIDATE_PHASE_RAMP_SWS Controlled phase-ramp SWS reconstruction test.
%
% Research question:
%   Can a phase-planarity score remain high when a known artificial spatial
%   phase gradient changes reconstructed shear-wave speed (SWS)?
%
% Input (read-only):
%   Results/method_comparison/batch/<case>/method_comparison.mat
%   Results/sws_measurement_quality/<case>/phase_consistency_diagnostic/
%       phase_quality_results.mat
%   Data/<case>_displacement_with_quality.mat (or other usual locations)
%
% Output (new folder only):
%   Results/sws_measurement_quality/<case>/phase_ramp_sws_validation/
%       phase_ramp_metrics.csv, phase_ramp_metrics.mat
%       phase_ramp_maps.png, phase_ramp_trends.png
%
% Usage:
%   validate_phase_ramp_sws('large_inclusions')          % -0.2, 0, +0.2 rad/mm
%   validate_phase_ramp_sws('large_inclusions','full')   % +/-0.2, +/-0.5, 0
%   validate_phase_ramp_sws('high_contrast')
%
% Method:
%   Fit a DC + cos(2*pi*f0*t) + sin(2*pi*f0*t) model to each displacement
%   time series. Rotate its harmonic coefficients by exp(1i*beta*x_mm),
%   then replace ONLY that harmonic component in the ORIGINAL time series.
%   The residual time-domain signal remains unchanged. Run the EXACT PG and
%   AC preprocessing/reconstruction recipe from compare_SWS_methods_batch.
%   Q is computed WITHOUT ground truth; GT only evaluates SWS error.
%
% Interpretation:
%   This is a controlled, potentially nonphysical perturbation. Adding a
%   phase ramp does NOT guarantee RMSE will increase in BOTH directions;
%   one direction may partially cancel pre-existing systematic bias.
%   A low model residual / high Q does not establish a correct SWS map.
%
% Requires the project's USTB elastography toolbox and MATLAB.

    if nargin < 1 || isempty(phantom_case)
        phantom_case = 'large_inclusions';
    end
    if nargin < 2 || isempty(mode)
        mode = 'quick';
    end
    phantom_case = char(string(phantom_case));
    mode = char(lower(string(mode)));
    assert(ismember(mode, {'quick','full'}), ...
        'mode must be ''quick'' or ''full''.');
    if strcmp(mode,'quick')
        betas = [-0.2, 0, 0.2];   % radians per millimetre
    else
        betas = [-0.5, -0.2, 0, 0.2, 0.5];
    end

    here = fileparts(mfilename('fullpath'));
    if isfolder(fullfile(here,'Results'))
        repo = here;
    else
        repo = fileparts(here);
    end
    toolbox = fullfile(fileparts(repo),'elastography','ElastographyToolbox');
    assert(isfolder(toolbox), 'ElastographyToolbox missing: %s',toolbox);
    if isempty(which('uff.elastography')) || ...
       isempty(which('postprocess.elastography.autocorrelation_SWS_method'))
        addpath(genpath(toolbox));
        rehash;
    end
    assert(~isempty(which('uff.beamformed_data')) && ...
        ~isempty(which('uff.elastography')) && ...
        ~isempty(which('postprocess.elastography.phase_gradient_method')) && ...
        ~isempty(which('postprocess.elastography.autocorrelation_SWS_method')), ...
        'Required USTB classes missing from MATLAB path.');

    basefolder = fullfile(repo,'Results','method_comparison','batch',phantom_case);
    compare_file = fullfile(basefolder,'method_comparison.mat');
    phase_file = fullfile(repo,'Results','sws_measurement_quality', ...
        phantom_case,'phase_consistency_diagnostic','phase_quality_results.mat');
    assert(isfile(compare_file), 'Missing comparison file: %s',compare_file);
    assert(isfile(phase_file), 'Missing saved phase Q: %s',phase_file);

    C = load(compare_file,'PG','AC','GT','roi','x_axis','z_axis', ...
        'PRF','target_frequency','filter_speed');
    A = load(phase_file,'Qphase','window_size','min_support', ...
        'amplitude_floor');
    needC={'PG','AC','GT','roi','x_axis','z_axis','PRF', ...
        'target_frequency','filter_speed'};
    assert(all(isfield(C,needC)), 'Comparison file lacks required fields.');
    assert(all(isfield(A,{'Qphase','window_size', ...
        'min_support','amplitude_floor'})), 'Phase file lacks settings.');
    x = double(C.x_axis(:)); z = double(C.z_axis(:));
    GT = double(C.GT);
    [Nz,Nx] = size(GT);
    assert(isequal(size(C.PG),[Nz Nx]) && ...
        isequal(size(C.AC),[Nz Nx]) && ...
        isequal(size(A.Qphase),[Nz Nx]),'Saved maps differ in size.');
    assert(Nx==numel(x) && Nz==numel(z), 'Scan axes differ from map.');
    assert(all(isfield(C.roi,{'image_mask','inclusion_core','background_core'})), ...
        'Missing ROI masks.');

    filename = [phantom_case '_displacement_with_quality.mat'];
    paths = {fullfile(repo,'Data',filename), ...
        fullfile(repo,'Data',phantom_case,filename), ...
        fullfile(repo,filename),fullfile(repo,'Analysis',filename)};
    datafile = '';
    for j=1:numel(paths)
        if isfile(paths{j}), datafile=paths{j}; break; end
    end
    assert(~isempty(datafile),'Could not find displacement input: %s',filename);

    fprintf('\n=== PHASE-RAMP SWS VALIDATION: %s [%s] ===\n',phantom_case,mode);
    fprintf('Loading: %s\n',datafile);
    D = load(datafile,'channel_data','b_data', ...
        'displacement_estimation_autocorr','shear_params');
    assert(all(isfield(D,{'channel_data','b_data', ...
        'displacement_estimation_autocorr','shear_params'})), ...
        'Dataset lacks ultrasound objects or simulation parameters.');
    assert(isobject(D.displacement_estimation_autocorr), ...
        'Displacement did not load as USTB object; check toolbox path.');
    prf = double(D.channel_data.PRF);
    f0 = double(D.shear_params.source_freq);
    filter_speed = double(C.filter_speed);
    assert(isscalar(prf) && isfinite(prf) && prf>2*f0 && ...
        f0>0 && filter_speed>0,'Invalid PRF, f0, or filter speed.');
    assert(abs(prf-double(C.PRF))<1e-5 && ...
        abs(f0-double(C.target_frequency))<1e-5, ...
        'Dataset PRF/f0 differs from stored comparison.');
    fprintf('Confirmed PRF = %.1f Hz, f0 = %.1f Hz, filter speed = %.3f m/s\n', ...
        prf,f0,filter_speed);

    scan = D.b_data.scan;
    assert(max(abs(double(scan.x_axis(:))-x))<1e-9 && ...
        max(abs(double(scan.z_axis(:))-z))<1e-9, ...
        'Ultrasound scan axes do not match stored SWS maps.');
    raw=D.displacement_estimation_autocorr.data;
    assert(isnumeric(raw) && isreal(raw) && ...
        mod(numel(raw),Nz*Nx)==0,'Displacement array size is invalid.');
    Nt=numel(raw)/(Nz*Nx);
    assert(Nt>=8,'Not enough temporal displacement frames.');
    u=-reshape(double(raw),Nz,Nx,Nt); % original sign convention
    clear raw D;
    assert(all(isfinite(u(:))), ...
        'Nonfinite displacement samples: review how original estimator handles them.');
    fprintf('Displacement: %d x %d x %d\n',Nz,Nx,Nt);

    % Least squares: u(t)=a*cos(w*t)+b*sin(w*t)+dc+residual.
    tt=(0:Nt-1)'/prf;
    omega=2*pi*f0;
    cos_t=cos(omega*tt); sin_t=sin(omega*tt);
    design=[cos_t sin_t ones(Nt,1)];
    assert(rank(design)==3,'Temporal harmonic fit is rank deficient.');
    p=pinv(design);
    a=zeros(Nz,Nx); b=zeros(Nz,Nx);
    for it=1:Nt
        frame=u(:,:,it);
        a=a+p(1,it)*frame;
        b=b+p(2,it)*frame;
    end

    [Xmm,~]=meshgrid(x*1e3,z*1e3);
    x0=mean(x)*1e3;
    settings=[double(A.window_size) double(A.min_support) ...
              double(A.amplitude_floor)];
    Hbaseline=harmonic_phasor_ramp(u,prf,f0);
    [Qbaseline,~,~]=phase_planarity_quality_ramp(Hbaseline, ...
        settings(1),settings(2),settings(3));
    qmask=isfinite(Qbaseline) & isfinite(A.Qphase);
    qdiff=max(abs(Qbaseline(qmask)-double(A.Qphase(qmask))));
    fprintf('Reproduced saved Q: max absolute difference %.4g\n',qdiff);
    assert(nnz(qmask)>100 && qdiff<1e-4, ...
        'Phase Q differs from saved pilot; check implementation/parameters.');
    clear Hbaseline;

    maps=repmat(struct('beta',NaN,'PG',[],'AC',[],'Q',[]),numel(betas),1);
    for k=1:numel(betas)
        beta=betas(k);
        if beta==0
            up=u;
            Q=Qbaseline;
        else
            theta=beta*(Xmm-x0); % beta: rad/mm, Xmm: mm
            ct=cos(theta); st=sin(theta);
            % (a-i*b)*exp(i*theta) => a_new=a*ct+b*st,
            % b_new=b*ct-a*st.
            da=a.*(ct-1)+b.*st;
            db=b.*(ct-1)-a.*st;
            up=zeros(size(u));
            for it=1:Nt
                up(:,:,it)=u(:,:,it)+da*cos_t(it)+db*sin_t(it);
            end
            H=harmonic_phasor_ramp(up,prf,f0);
            [Q,~,~]=phase_planarity_quality_ramp(H, ...
                settings(1),settings(2),settings(3));
            clear H da db theta ct st;
        end
        fprintf('Reconstructing beta = %+.2f rad/mm (PG + AC)...\n',beta);
        [PG,AC]=reconstruct_both_ramp(up,scan,prf,f0,filter_speed);
        if beta~=0, clear up; end
        assert(isequal(size(PG),[Nz Nx]) && isequal(size(AC),[Nz Nx]), ...
            'Unexpected SWS map dimensions from USTB.');
        maps(k).beta=beta;
        maps(k).PG=PG;
        maps(k).AC=AC;
        maps(k).Q=Q;
    end
    clear u a b;

    k0=find(betas==0,1);
    % Diagnostic consistency: old and newly regenerated PG/AC baseline.
    pg_old=double(C.PG); ac_old=double(C.AC);
    cmp=isfinite(pg_old)&isfinite(maps(k0).PG)&logical(C.roi.image_mask);
    pg_baseline_diff=sqrt(mean((maps(k0).PG(cmp)-pg_old(cmp)).^2));
    cmp=isfinite(ac_old)&isfinite(maps(k0).AC)&logical(C.roi.image_mask);
    ac_baseline_diff=sqrt(mean((maps(k0).AC(cmp)-ac_old(cmp)).^2));
    fprintf('Regenerated baseline differences: PG %.5g, AC %.5g m/s (RMSE)\n', ...
        pg_baseline_diff,ac_baseline_diff);
    if pg_baseline_diff>1e-3 || ac_baseline_diff>1e-3
        warning(['Baseline differs from stored maps. The perturbed comparison ' ...
            'is internally consistent but inspect toolbox versions/settings.']);
    end

    % Common validity for EVERY reconstruction, BOTH estimators, and Q.
    base=logical(C.roi.image_mask) & isfinite(GT) & isfinite(Qbaseline);
    common=base;
    for k=1:numel(maps)
        pg=maps(k).PG; ac=maps(k).AC;
        common=common & isfinite(pg) & pg>0 & pg<10 & ...
            isfinite(ac) & ac>0 & ac<10 & isfinite(maps(k).Q);
    end
    fprintf('Common pixels across all ramps/estimators: %d of %d (%.2f%%)\n', ...
        nnz(common),nnz(base),100*nnz(common)/max(nnz(base),1));
    assert(nnz(common)>100,'Too few jointly valid SWS pixels.');
    if nnz(common)/nnz(base)<0.9
        warning('Joint validity below 90%%; consider reporting failure fraction explicitly.');
    end

    near_boundary=(Xmm<=min(x)*1e3+3.5) | ...
                  (Xmm>=max(x)*1e3-3.5);
    region_masks={common,common & logical(C.roi.inclusion_core), ...
        common & logical(C.roi.background_core) & ~near_boundary, ...
        common & near_boundary};
    region_names={'Global','InclusionCore','BackgroundInterior','LateralBoundary'};
    methods={'PG','AC'};
    T=table();
    for k=1:numel(maps)
        for j=1:2
            method=methods{j};
            actual=double(maps(k).(method));
            for r=1:numel(region_masks)
                mask=region_masks{r};
                if nnz(mask)<20, continue; end
                e=actual(mask)-GT(mask);
                q=double(maps(k).Q(mask));
                q0=double(Qbaseline(mask));
                valid_method=base & isfinite(actual) & actual>0 & actual<10;
                valid_fraction=nnz(valid_method)/nnz(base);
                row=table(string(phantom_case),maps(k).beta,string(method), ...
                    string(region_names{r}),nnz(mask), ...
                    sqrt(mean(e.^2)),mean(abs(e)),mean(e), ...
                    mean(q),median(q),mean(q>=1-1e-6), ...
                    mean(q-q0),max(abs(q-q0)),valid_fraction, ...
                    'VariableNames',{'Phantom','Beta_rad_per_mm','Estimator', ...
                    'Region','CommonPixels','RMSE_mps','MAE_mps','Bias_mps', ...
                    'MeanQ','MedianQ','FractionQSaturated','DeltaMeanQ', ...
                    'MaxAbsDeltaQ','ValidFraction'});
                T=[T;row]; %#ok<AGROW>
            end
        end
    end
    % Delta error on exactly the same mask, compared to reconstructed beta=0.
    T.DeltaRMSE_mps=nan(height(T),1);
    T.DeltaMAE_mps=nan(height(T),1);
    for k=1:height(T)
        ref=T.Beta_rad_per_mm==0 & T.Estimator==T.Estimator(k) & ...
            T.Region==T.Region(k);
        assert(nnz(ref)==1,'Could not identify one baseline row.');
        T.DeltaRMSE_mps(k)=T.RMSE_mps(k)-T.RMSE_mps(ref);
        T.DeltaMAE_mps(k)=T.MAE_mps(k)-T.MAE_mps(ref);
    end
    disp(T(T.Region=="Global",:));

    outdir=fullfile(repo,'Results','sws_measurement_quality',phantom_case, ...
        'phase_ramp_sws_validation');
    if ~isfolder(outdir), mkdir(outdir); end
    writetable(T,fullfile(outdir,'phase_ramp_metrics.csv'));
    save(fullfile(outdir,'phase_ramp_metrics.mat'), ...
        'T','maps','GT','common','base','region_names','betas', ...
        'pg_baseline_diff','ac_baseline_diff','prf','f0', ...
        'filter_speed','settings','x','z','-v7.3');

    chosen=[1,k0,numel(betas)];
    fig=figure('Color','w','Name','Phase-ramp Q and SWS errors', ...
        'NumberTitle','off','Position',[50 80 1500 1080]);
    tl=tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
    for j=1:3
        ix=chosen(j);
        label=sprintf('beta = %+.1f rad/mm',betas(ix));
        show_ramp_map(nexttile(tl,j),x,z,maps(ix).Q,[0 1], ...
            ['Phase Q: ' label]);
        show_ramp_map(nexttile(tl,3+j),x,z,abs(maps(ix).PG-GT),[0 4], ...
            ['PG absolute error: ' label]);
        show_ramp_map(nexttile(tl,6+j),x,z,abs(maps(ix).AC-GT),[0 4], ...
            ['AC absolute error: ' label]);
    end
    title(tl,['Phase-ramp SWS validation: ' strrep(phantom_case,'_',' ')]);
    exportgraphics(fig,fullfile(outdir,'phase_ramp_maps.png'),'Resolution',160);

    fig2=figure('Color','w','Name','Phase-ramp error and Q trends', ...
        'NumberTitle','off','Position',[80 80 1150 750]);
    tl2=tiledlayout(fig2,2,2,'TileSpacing','compact','Padding','compact');
    g=T(T.Region=="Global",:);
    ax=nexttile(tl2); hold(ax,'on');
    for j=1:2
        v=g(g.Estimator==string(methods{j}),:);
        plot(ax,v.Beta_rad_per_mm,v.RMSE_mps,'-o','DisplayName',methods{j});
    end
    hold(ax,'off');grid(ax,'on');legend(ax,'Location','best');
    xlabel(ax,'Imposed phase ramp [rad/mm]');ylabel(ax,'RMSE [m/s]');
    title(ax,'Fixed common pixels: SWS RMSE');

    ax=nexttile(tl2); hold(ax,'on');
    for j=1:2
        v=g(g.Estimator==string(methods{j}),:);
        plot(ax,v.Beta_rad_per_mm,v.DeltaRMSE_mps,'-o', ...
            'DisplayName',methods{j});
    end
    yline(ax,0,':k');hold(ax,'off');grid(ax,'on');
    xlabel(ax,'Imposed phase ramp [rad/mm]');ylabel(ax,'Delta RMSE [m/s]');
    title(ax,'Change relative to regenerated baseline');

    v=g(g.Estimator=="PG",:); % Q is shared between PG and AC.
    ax=nexttile(tl2);plot(ax,v.Beta_rad_per_mm,v.MeanQ,'-o');
    grid(ax,'on');ylim(ax,[0 1.02]);
    xlabel(ax,'Imposed phase ramp [rad/mm]');ylabel(ax,'Mean Q');
    title(ax,'Phase planarity (no GT)');

    ax=nexttile(tl2);plot(ax,v.Beta_rad_per_mm,v.MaxAbsDeltaQ,'-o');
    grid(ax,'on');
    xlabel(ax,'Imposed phase ramp [rad/mm]');ylabel(ax,'Max absolute change in Q');
    title(ax,'Is Q invariant to linear phase ramps?');
    title(tl2,['Phase-ramp metrics: ' strrep(phantom_case,'_',' ')]);
    exportgraphics(fig2,fullfile(outdir,'phase_ramp_trends.png'), ...
        'Resolution',170);

    fprintf('\nSaved outputs to:\n%s\n',outdir);
    fprintf(['CAUTION: Synthetic phase ramps are a controlled perturbation. ' ...
        'They are not physically simulated wave propagation.\n' ...
        'Compare DeltaRMSE and DeltaMeanQ on a common pixel mask. ' ...
        'A change in RMSE is NOT necessarily positive for both ramp signs.\n']);
end

function [PG,AC]=reconstruct_both_ramp(u,scan,PRF,target_frequency,filter_speed)
% Same recipe as compare_SWS_methods_batch.m estimate_both_sws.
    elasto=uff.elastography();
    elasto.PRF=PRF;
    elasto.target_frequency=target_frequency;
    elasto.scan=scan;
    elasto.autocorrelation_field=uff.data_field();
    elasto.autocorrelation_field.set_component('u', ...
        reshape(u,size(u,1),size(u,2),1,size(u,3)));
    num_angles=12;
    elasto.directional_filtering_angles= ...
        2*pi*linspace(0,1-1/num_angles,num_angles);
    extractor=postprocess.elastography.extract_freq_component(elasto);
    extractor.target_frequency=target_frequency;
    extractor.component='u';
    extractor.run();
    omega=2*pi*target_frequency;
    elasto.radial_filter=uff.spatial_filter('gaussian_charite','radial');
    elasto.radial_filter.sigma=filter_speed/omega;
    elasto.angular_filter=uff.spatial_filter('gaussian_charite','angular');
    elasto.angular_filter.deltaTheta= ...
        median(diff(elasto.directional_filtering_angles));
    filter_operator=postprocess.elastography.apply_filter(elasto);
    filter_operator.run();
    pgm=postprocess.elastography.phase_gradient_method(elasto);
    pgm.run();
    PG=squeeze(double(elasto.sws_map.u));
    acm=postprocess.elastography.autocorrelation_SWS_method(elasto);
    acm.run();
    AC=squeeze(double(elasto.sws_map.u));
end

function show_ramp_map(ax,x,z,M,limits,ttl)
    imagesc(ax,x*1e3,z*1e3,M);
    axis(ax,'image');set(ax,'YDir','normal');
    clim(ax,limits);colorbar(ax);
    xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
    title(ax,ttl,'Interpreter','none');
    try,ax.Toolbar.Visible='off';catch,end
end

function H = harmonic_phasor_ramp(u,prf,f0)
% Exact temporal harmonic lock-in used by test_sws_phase_consistency.m.
    [Nz,Nx,Nt]=size(u);
    w=exp(-1i*2*pi*f0*(0:Nt-1)/prf);
    H=complex(zeros(Nz,Nx));
    temporalSum=zeros(Nz,Nx);
    for it=1:Nt
        frame=double(u(:,:,it));
        frame(~isfinite(frame))=0;
        temporalSum=temporalSum+frame;
        H=H+frame*w(it);
    end
    H=(H-temporalSum*(sum(w)/Nt))/Nt;
end

function [Q,phaseRms,amp] = phase_planarity_quality_ramp(H,win,minSup,ampFrac)
% Exact phase-planarity method from original pilot, without using GT.
    [Nz,Nx]=size(H);
    amp=abs(H);
    goodAmp=amp(isfinite(amp)&amp>0);
    assert(~isempty(goodAmp),'Harmonic displacement is all zero.');
    threshold=ampFrac*median(goodAmp);
    good=isfinite(amp)&amp>=threshold;
    U=complex(zeros(Nz,Nx));
    U(good)=H(good)./amp(good);
    K=ones(win,win);
    pairX=complex(zeros(Nz,Nx));
    pairZ=complex(zeros(Nz,Nx));
    okX=false(Nz,Nx); okZ=false(Nz,Nx);
    okX(:,1:end-1)=good(:,1:end-1)&good(:,2:end);
    okZ(1:end-1,:)=good(1:end-1,:)&good(2:end,:);
    pairX(:,1:end-1)=U(:,2:end).*conj(U(:,1:end-1));
    pairZ(1:end-1,:)=U(2:end,:).*conj(U(1:end-1,:));
    pairX(~okX)=0;
    pairZ(~okZ)=0;
    kxStep=angle(conv2(pairX,K,'same'));
    kzStep=angle(conv2(pairZ,K,'same'));
    nx=conv2(double(okX),K,'same');
    nz=conv2(double(okZ),K,'same');
    haveSlopes=good & nx>=max(4,round(0.2*win*win)) & ...
        nz>=max(4,round(0.2*win*win));
    re=zeros(Nz,Nx); im=zeros(Nz,Nx); counts=zeros(Nz,Nx);
    radius=(win-1)/2;
    for zz=-radius:radius
        for xx=-radius:radius
            r1=max(1,1-zz);r2=min(Nz,Nz-zz);
            c1=max(1,1-xx);c2=min(Nx,Nx-xx);
            rows=r1:r2;cols=c1:c2;
            nr=rows+zz;nc=cols+xx;
            usable=haveSlopes(rows,cols)&good(nr,nc);
            phaseDifference=angle(U(nr,nc).*conj(U(rows,cols)));
            predicted=xx*kxStep(rows,cols)+zz*kzStep(rows,cols);
            residual=phaseDifference-predicted;
            re(rows,cols)=re(rows,cols)+cos(residual).*usable;
            im(rows,cols)=im(rows,cols)+sin(residual).*usable;
            counts(rows,cols)=counts(rows,cols)+double(usable);
        end
    end
    Q=hypot(re,im)./max(counts,1);
    Q(~haveSlopes|counts<minSup*win*win)=NaN;
    Q=max(0,min(1,Q));
    phaseRms=sqrt(max(0,2*(1-Q)));
end

