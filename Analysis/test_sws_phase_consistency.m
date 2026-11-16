function summary = test_sws_phase_consistency(phantom_case)
%TEST_SWS_PHASE_CONSISTENCY Ground-truth-independent phase-planarity pilot.
%
% WHY:
% The earlier displacement-correlation Q saturated for harmonic motion.
% Here we test a DIFFERENT, exploratory phase-based candidate. It is not a
% published/validated metric, and it is NOT a reproduction of Schoen et al.
%
% INPUTS: existing, read-only
% Results/sws_measurement_quality/<case>/quality_validation.mat
% <case>_displacement_with_quality.mat (in Data, Data/<case>, or repo root)
%
% OUTPUTS: newly created directory ONLY
% Results/sws_measurement_quality/<case>/phase_consistency_diagnostic/
%   phase_quality_vs_error.png
%   phase_quality_scatter.png
%   phase_quality_summary.csv
%   phase_quality_results.mat
%
% USAGE: save in PROJECT/Analysis/ and run
%   test_sws_phase_consistency('large_inclusions')
%   test_sws_phase_consistency('high_contrast')
%
% Q_phase uses displacement and excitation frequency ONLY. Ground truth and
% saved SWS maps are loaded to evaluate its relationship with SWS error.
%
% Interpretation caveats: high phase planarity does not prove SWS accuracy.
% Systematic wavefield/model errors can produce Q_phase ~ 1. Do not tune
% parameters on these phantoms and call the same test independent validation.

    if nargin < 1 || isempty(phantom_case)
        phantom_case = 'large_inclusions';
    end
    phantom_case = char(string(phantom_case));

    here = fileparts(mfilename('fullpath'));
    if isfolder(fullfile(here,'Results'))
        repo = here;
    else
        repo = fileparts(here);
    end

    % Parameters set before inspecting ground-truth validation outcomes.
    window_size = 9;         % odd-numbered spatial window, in pixels
    min_support = 0.65;      % fraction of the neighborhood with valid phase
    amplitude_floor = 0.01;  % minimum amplitude = 1% of median harmonic amp
    fprintf('\n=== PHASE CONSISTENCY PILOT: %s ===\n',phantom_case);

    % Sanity check the phase-planarity implementation with artificial data.
    self_test_phase_consistency(window_size, min_support);

    pilot_file = fullfile(repo,'Results','sws_measurement_quality', ...
        phantom_case,'quality_validation.mat');
    assert(isfile(pilot_file), ...
        'Missing previous Q pilot: %s\nRun test_sws_measurement_quality first.', ...
        pilot_file);
    % Load previous quality-metric results
    S = load(pilot_file);
    Qold = double(S.Q);
    GT   = double(S.GT);
    PG   = double(S.PG);
    AC   = double(S.AC);
    x = double(S.x_axis(:));
    z = double(S.z_axis(:));
    [Nz,Nx] = size(GT);
    assert(isequal(size(Qold),size(GT),size(PG),size(AC),size(S.valid)), ...
        'The saved maps and masks have inconsistent dimensions.');
    assert(numel(x)==Nx && numel(z)==Nz, 'Axes do not match SWS image.');
    assert(Nx >= 3*window_size && Nz >= 3*window_size, ...
        'The grid is too small for the chosen window.');

    % Previously confirmed toolbox location on this PC. The USTB class
    % must exist before MATLAB can load uff.beamformed_data objects.
    if isempty(which('uff.beamformed_data'))
        ustb_root = fullfile(fileparts(repo),'elastography', ...
            'ElastographyToolbox','USTB');
        if isfolder(ustb_root)
            addpath(genpath(ustb_root));
            rehash;
        end
    end
    assert(~isempty(which('uff.beamformed_data')), ...
        ['Could not find uff.beamformed_data. Add ' ...
         'MasterOppgave/elastography/ElastographyToolbox/USTB to MATLAB path.']);

    file_name = [phantom_case '_displacement_with_quality.mat'];
    candidates = {fullfile(repo,'Data',file_name), ...
        fullfile(repo,'Data',phantom_case,file_name), ...
        fullfile(repo,file_name), ...
        fullfile(repo,'Analysis',file_name)};
    data_file = '';
    for j=1:numel(candidates)
        if isfile(candidates{j})
            data_file = candidates{j};
            break;
        end
    end
    assert(~isempty(data_file), 'Cannot find %s in project.',file_name);
    fprintf('Loading displacement from: %s\n',data_file);
    D = load(data_file, ...
    'displacement_estimation_autocorr', ...
    'shear_params', ...
    'channel_data');
    assert(isfield(D,'displacement_estimation_autocorr') && ...
        isfield(D,'shear_params'), 'Displacement/shear_params missing.');
    assert(isobject(D.displacement_estimation_autocorr) || ...
        isstruct(D.displacement_estimation_autocorr), ...
        'Displacement object did not load. Check USTB compatibility.');

    raw = D.displacement_estimation_autocorr.data;
    assert(isnumeric(raw) && isreal(raw) && mod(numel(raw),Nz*Nx)==0, ...
        'Unexpected displacement data type/dimensions.');
    Nt = numel(raw)/(Nz*Nx);
    assert(Nt >= 8, 'Too few temporal frames.');
    u = reshape(single(raw),Nz,Nx,Nt);
    clear raw;

   %% Confirm sampling frequency

f0 = double(D.shear_params.source_freq);

prf = double(D.channel_data.PRF);

assert(isfinite(prf) && prf > 0, ...
    'Invalid PRF in channel_data.');

fprintf('Confirmed PRF: %.1f Hz\n', prf);
fprintf('Excitation frequency: %.1f Hz\n', f0);
    clear D;
    assert(isscalar(f0) && isfinite(f0) && f0>0 && f0<prf/2, ...
        'Target frequency must lie below Nyquist.');
    fprintf('Displacement array: %d x %d x %d, PRF=%.4g Hz, f0=%.4g Hz\n', ...
        Nz,Nx,Nt,prf,f0);

    % Extract complex harmonic displacement at f0, removing temporal DC.
    H = harmonic_phasor(u,prf,f0);
    clear u;
    [Qphase,phaseResidual,phaseAmp,phaseGradX,phaseGradZ] = ...
        phase_planarity_quality(H,window_size,min_support,amplitude_floor);
    clear H;

    % Ground truth enters ONLY here for validation; NOT Qphase calculation.
    common = logical(S.valid) & isfinite(GT) & isfinite(Qold) & ...
        isfinite(Qphase) & isfinite(PG) & isfinite(AC);
    assert(nnz(common)>100, 'Too few common valid phase/Q/SWS pixels.');
    fprintf('Common validation pixels: %d\n',nnz(common));

    errPG = abs(PG-GT);
    errAC = abs(AC-GT);
    methods = {'PG','AC'};
    errors = {errPG,errAC};
    qmaps = {Qold,Qphase};
    names = {'Published Q adaptation','Phase planarity (candidate)'};
    summary=table();
    for i=1:numel(qmaps)
        for j=1:2
            q=double(qmaps{i}(common));
            e=double(errors{j}(common));
            row=table(string(phantom_case),string(methods{j}), ...
                string(names{i}), nnz(common), mean(q), median(q), ...
                std(q),mean(q>=1-1e-6), ...
                correlation_spearman(q,e), ...
                mean(e(q>=0.8),'omitnan'),mean(e(q<0.8),'omitnan'), ...
                'VariableNames',{'Phantom','Method','QualityMetric', ...
                'ValidPixels','MeanQ','MedianQ','StdQ','FractionQSaturated', ...
                'SpearmanQvsAbsError','MAE_Qge08_mps','MAE_Qlt08_mps'});
            summary=[summary;row]; %#ok<AGROW>
        end
    end
    disp(summary);
    if mean(Qphase(common)>=1-1e-6)>0.8
        warning(['Phase quality also saturates strongly. ' ...
            'High phase planarity is not evidence of accurate SWS.']);
    end

    outdir=fullfile(repo,'Results','sws_measurement_quality', ...
        phantom_case,'phase_consistency_diagnostic');
    if ~isfolder(outdir), mkdir(outdir); end
    writetable(summary,fullfile(outdir,'phase_quality_summary.csv'));
    save(fullfile(outdir,'phase_quality_results.mat'), ...
        'summary','Qold','Qphase','phaseResidual','phaseAmp', ...
        'phaseGradX','phaseGradZ','common','GT','PG','AC','x','z', ...
        'prf','f0','window_size','amplitude_floor','min_support','-v7.3');

    % Comparisons on identical spatial axes and a common validation region.
    f=figure('Color','w','NumberTitle','off', ...
        'Name',['Phase quality versus SWS errors: ' phantom_case]);
    tl=tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
    x_mm=x*1e3; z_mm=z*1e3;
    plot_map(nexttile(tl),x_mm,z_mm,Qold,[0 1], ...
        'Previous displacement Q');
    plot_map(nexttile(tl),x_mm,z_mm,Qphase,[0 1], ...
        'New local phase planarity');
    plot_map(nexttile(tl),x_mm,z_mm,Qphase-Qold,[-1 1], ...
        'Phase Q minus previous Q');
    plot_map(nexttile(tl),x_mm,z_mm,errPG,[0 4], ...
        '|PG - ground truth| [m/s]');
    plot_map(nexttile(tl),x_mm,z_mm,errAC,[0 4], ...
        '|AC - ground truth| [m/s]');
    amp_db=20*log10(max(phaseAmp/median(phaseAmp(isfinite(phaseAmp)&phaseAmp>0)),1e-6));
    plot_map(nexttile(tl),x_mm,z_mm,amp_db,[-30 20], ...
        'Harmonic displacement magnitude [dB rel. median]');
    title(tl,['Phase-planarity pilot: ' strrep(phantom_case,'_',' ')]);
    exportgraphics(f,fullfile(outdir,'phase_quality_vs_error.png'), ...
        'Resolution',160);

    f2=figure('Color','w','NumberTitle','off', ...
        'Name',['Quality-error scatter: ' phantom_case]);
    t2=tiledlayout(f2,1,2,'TileSpacing','compact','Padding','compact');
    evPG=errPG(common); evAC=errAC(common);
    old=Qold(common); new=Qphase(common);
    ids=unique(round(linspace(1,nnz(common),min(nnz(common),4000))));
    scatter_plot(nexttile(t2),old(ids),evPG(ids),evAC(ids), ...
        'Previous displacement Q');
    scatter_plot(nexttile(t2),new(ids),evPG(ids),evAC(ids), ...
        'New phase planarity');
    title(t2,['Q versus absolute SWS error: ' strrep(phantom_case,'_',' ')]);
    exportgraphics(f2,fullfile(outdir,'phase_quality_scatter.png'), ...
        'Resolution',160);

    fprintf('\nSaved new outputs to:\n  %s\n',outdir);
    fprintf(['A good confidence score should generally give NEGATIVE ' ...
        'Spearman(Q, |SWS-GT|). Inspect Q saturation as well.\n']);
end

function H = harmonic_phasor(u,prf,f0)
% Temporal lock-in extraction. NO GT, no SWS map used here.
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

function [Q,phaseRms,amp,kxStep,kzStep] = ...
        phase_planarity_quality(H,win,minSup,ampFrac)
% Novel *candidate* (not validated): local complex-residual phase coherence.
% Estimate a plane-wave phase slope from nearby spatial phase increments.
% Remove that fitted phase ramp and compute coherence of residual phases:
% Q = abs(mean(exp(1i*wrapped_phase_residual))) for each local window.
% Q=1 for a perfectly planar phase; random phase tends toward zero.
    [Nz,Nx]=size(H);
    amp=abs(H);
    goodAmp=amp(isfinite(amp)&amp>0);
    assert(~isempty(goodAmp),'Harmonic displacement has no nonzero data.');
    threshold=ampFrac*median(goodAmp);
    good=isfinite(amp)&amp>=threshold;
    U=complex(zeros(Nz,Nx));
    U(good)=H(good)./amp(good);
    K=ones(win,win);

    pairX=complex(zeros(Nz,Nx));
    pairZ=complex(zeros(Nz,Nx));
    okX=false(Nz,Nx);
    okZ=false(Nz,Nx);
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

    % Compare local signal with a fitted constant-slope phase plane.
    re=zeros(Nz,Nx); im=zeros(Nz,Nx); counts=zeros(Nz,Nx);
    radius=(win-1)/2;
    for dz=-radius:radius
        for dx=-radius:radius
            r1=max(1,1-dz); r2=min(Nz,Nz-dz);
            c1=max(1,1-dx); c2=min(Nx,Nx-dx);
            rows=r1:r2; cols=c1:c2;
            nr=rows+dz; nc=cols+dx;
            usable=haveSlopes(rows,cols)&good(nr,nc);
            phaseDifference=angle(U(nr,nc).*conj(U(rows,cols)));
            predicted=dx*kxStep(rows,cols)+dz*kzStep(rows,cols);
            residual=phaseDifference-predicted;
            re(rows,cols)=re(rows,cols)+cos(residual).*usable;
            im(rows,cols)=im(rows,cols)+sin(residual).*usable;
            counts(rows,cols)=counts(rows,cols)+double(usable);
        end
    end
    Q=hypot(re,im)./max(counts,1);
    Q(~haveSlopes|counts<minSup*win*win)=NaN;
    Q=max(0,min(1,Q));
    % Circular RMS proxy, valid for small residuals. Diagnostic only.
    phaseRms=sqrt(max(0,2*(1-Q)));
end

function self_test_phase_consistency(win,minSup)
% Check expected behavior on artificial phase fields before real data.
    [xx,zz]=meshgrid(0:79,0:79);
    Hplane=exp(1i*(0.14*xx+0.09*zz));
    [qplane,~,~,~,~]=phase_planarity_quality(Hplane,win,minSup,0.01);
    core=qplane(12:end-11,12:end-11);
    assert(median(core(:),'omitnan')>0.995, ...
        'Phase implementation failed pure-planewave sanity check.');
    rng_state=rng;
    rng(7123,'twister');
    Hrnd=exp(1i*2*pi*rand(size(xx)));
    rng(rng_state);
    [qnoise,~,~,~,~]=phase_planarity_quality(Hrnd,win,minSup,0.01);
    core=qnoise(12:end-11,12:end-11);
    assert(median(core(:),'omitnan')<0.8, ...
        'Phase implementation failed random-phase sanity check.');
    fprintf('Synthetic sanity tests passed (plane wave vs random phase).\n');
end

function r=correlation_spearman(a,b)
    a=double(a(:)); b=double(b(:));
    ok=isfinite(a)&isfinite(b);
    a=a(ok); b=b(ok);
    if numel(a)<30 || std(a)<1e-10 || std(b)<1e-10
        r=NaN;
    else
        r=corr(a,b,'Type','Spearman');
    end
end

function plot_map(ax,x,z,map,lims,ttl)
    imagesc(ax,x,z,map);
    axis(ax,'image'); set(ax,'YDir','normal');
    clim(ax,lims); colorbar(ax);
    xlabel(ax,'Lateral [mm]'); ylabel(ax,'Depth [mm]');
    title(ax,ttl,'Interpreter','none');
    try, ax.Toolbar.Visible='off'; catch, end
end

function scatter_plot(ax,Q,errPG,errAC,ttl)
    plot(ax,Q,errPG,'.','MarkerSize',3); hold(ax,'on');
    plot(ax,Q,errAC,'.','MarkerSize',3); hold(ax,'off');
    grid(ax,'on'); xlim(ax,[0 1]);
    xlabel(ax,'Quality score (no GT)');
    ylabel(ax,'Absolute SWS error [m/s] (GT for validation)');
    legend(ax,{'PG','AC'},'Location','best');
    title(ax,ttl,'Interpreter','none');
end
