function T = stress_test_sws_phase_quality(phantom_case)
%STRESS_TEST_SWS_PHASE_QUALITY Controlled perturbations of harmonic displacement.
%
% PURPOSE
%   Test whether the ground-truth-independent phase-planarity Q responds to
%   known changes in its INPUT signal. This is a sensitivity test, NOT a
%   validation of reconstructed SWS accuracy under the perturbations.
%
% INPUTS (read-only)
%   Results/sws_measurement_quality/<case>/phase_consistency_diagnostic/
%       phase_quality_results.mat
%   Results/method_comparison/batch/<case>/method_comparison.mat
%   <case>_displacement_with_quality.mat
%
% OUTPUTS (only in new directory)
%   Results/sws_measurement_quality/<case>/phase_stress_test/
%       phase_stress_summary.csv
%       phase_stress_summary.mat
%       phase_stress_maps.png
%       phase_stress_trends.png
%
% EXAMPLES
%   stress_test_sws_phase_quality('large_inclusions')
%   stress_test_sws_phase_quality('high_contrast')
%
% DESIGNS
%   - AWGN on COMPLEX HARMONIC DISPLACEMENT at 10-40 dB, 3 realizations.
%     This is NOT the SNR of beamformed ultrasound/RF measurements.
%   - Spatial Gaussian smoothing of harmonic displacement (mm).
%   - Added linear lateral phase ramp (negative control: locally planar).
%   - Local Gaussian phase defect (tests localized phase inconsistency).
%
% The quality score uses only perturbed displacement data, not GT.
% Existing GT/ROI-informed masks are used only to stratify the evaluation.
% Changing the phasor without re-running SWS estimators DOES NOT provide
% new SWS errors. In particular, no conclusions about SWS accuracy can be
% drawn directly from this sensitivity-only test.

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

    phase_file = fullfile(repo,'Results','sws_measurement_quality', ...
        phantom_case,'phase_consistency_diagnostic','phase_quality_results.mat');
    comparison_file = fullfile(repo,'Results','method_comparison','batch', ...
        phantom_case,'method_comparison.mat');
    assert(isfile(phase_file), ...
        'Missing phase-quality pilot: %s\nRun test_sws_phase_consistency first.',phase_file);
    assert(isfile(comparison_file), ...
        'Missing method comparison: %s',comparison_file);

    A = load(phase_file,'Qphase','common','x','z','prf','f0', ...
        'window_size','min_support','amplitude_floor');
    assert(all(isfield(A,{'Qphase','common','x','z','prf','f0', ...
        'window_size','min_support','amplitude_floor'})), ...
        'Phase-quality pilot lacks required data.');
    B = load(comparison_file,'roi');
    assert(isfield(B,'roi') && isfield(B.roi,'image_mask') && ...
        isfield(B.roi,'inclusion_core') && isfield(B.roi,'background_core'), ...
        'Missing ROI masks in method_comparison.mat.');

    x = double(A.x(:)); z = double(A.z(:));
    [Nz,Nx] = size(A.Qphase);
    assert(numel(x)==Nx && numel(z)==Nz, 'Inconsistent map dimensions.');
    assert(isequal(size(A.common),[Nz Nx]), 'Invalid common mask size.');

    % Find USTB as in the earlier pilot.
    if isempty(which('uff.beamformed_data'))
        ustb_dir = fullfile(fileparts(repo),'elastography', ...
            'ElastographyToolbox','USTB');
        if isfolder(ustb_dir)
            addpath(genpath(ustb_dir));
            rehash;
        end
    end
    assert(~isempty(which('uff.beamformed_data')), ...
        'USTB missing. Add .../elastography/ElastographyToolbox/USTB to path.');

    name = [phantom_case '_displacement_with_quality.mat'];
    candidates = {fullfile(repo,'Data',name), ...
        fullfile(repo,'Data',phantom_case,name),fullfile(repo,name), ...
        fullfile(repo,'Analysis',name)};
    data_file = '';
    for i=1:numel(candidates)
        if isfile(candidates{i}), data_file = candidates{i}; break; end
    end
    assert(~isempty(data_file), 'Cannot find input displacement file: %s',name);

    fprintf('\n=== CONTROLLED PHASE-Q STRESS TEST: %s ===\n',phantom_case);
    fprintf('Loading displacement from: %s\n',data_file);
    D = load(data_file,'displacement_estimation_autocorr', ...
        'shear_params','channel_data');
    assert(isfield(D,'displacement_estimation_autocorr') && ...
        isfield(D,'channel_data') && isfield(D,'shear_params'), ...
        'Input dataset lacks required variables.');
    assert(isobject(D.displacement_estimation_autocorr), ...
        'Displacement did not load as a USTB object.');

    prf = double(D.channel_data.PRF);
    f0 = double(D.shear_params.source_freq);
    assert(isscalar(prf) && isfinite(prf) && prf>0 && ...
        isscalar(f0) && isfinite(f0) && f0>0 && f0<prf/2, ...
        'Invalid source frequency or PRF.');
    assert(abs(prf-double(A.prf))<1e-5 && abs(f0-double(A.f0))<1e-5, ...
        'Stored pilot PRF/f0 differs from input dataset.');
    fprintf('Confirmed PRF = %.1f Hz, f0 = %.1f Hz\n',prf,f0);

    raw = D.displacement_estimation_autocorr.data;
    assert(isnumeric(raw) && isreal(raw) && mod(numel(raw),Nz*Nx)==0, ...
        'Unexpected displacement dimensions.');
    Nt = numel(raw)/(Nz*Nx);
    u = reshape(single(raw),Nz,Nx,Nt);
    clear raw D;
    fprintf('Displacement size: %d x %d x %d\n',Nz,Nx,Nt);

    win = double(A.window_size);
    minSup = double(A.min_support);
    ampFrac = double(A.amplitude_floor);
    H0 = harmonic_phasor_stress(u,prf,f0);
    clear u;

    % Exact copy of the original phase-Q algorithm, for a FAIR test.
    [Q0,~,~] = phase_planarity_quality_stress(H0,win,minSup,ampFrac);
    same = isfinite(Q0) & isfinite(A.Qphase) & logical(A.common);
    assert(nnz(same)>100,'No common pixels to check the original Q.');
    max_difference = max(abs(Q0(same)-double(A.Qphase(same))));
    fprintf('Baseline-Q reproduction: max |Qnew-Qsaved| = %.3g\n',max_difference);
    assert(max_difference<1e-4, ...
        ['Baseline Q does not reproduce saved phase_quality_results.mat. ' ...
         'Check that this script and the phase pilot use the same Q algorithm.']);

    [Xmm,Zmm] = meshgrid(x*1e3,z*1e3);
    dx = median(diff(x))*1e3;
    dz = median(diff(z))*1e3;
    assert(dx>0 && dz>0,'Expected increasing spatial axes.');
    base = logical(A.common) & logical(B.roi.image_mask) & isfinite(Q0);
    near_boundary = (Xmm <= min(Xmm(:))+3.5) | ...
                    (Xmm >= max(Xmm(:))-3.5);
    masks = {base, base & logical(B.roi.inclusion_core), ...
        base & logical(B.roi.background_core) & ~near_boundary, ...
        base & near_boundary};
    region_names = {'Global','InclusionCore','BackgroundInterior','LateralBoundary'};
    assert(nnz(masks{1})>100,'No global Q pixels in mask.');

    % These values are fixed in advance; no calibration against SWS error.
    snr_db = [40 30 20 10];
    noise_repetitions = 3;
    blur_mm = [0.5 1.5];
    phase_ramp_rad_mm = [0.2 0.5];
    defect_amplitudes_rad = [1.5 3.0];
    defect_center_mm = [-3 15];
    defect_sigma_mm = 2.0;

    outdir = fullfile(repo,'Results','sws_measurement_quality', ...
        phantom_case,'phase_stress_test');
    if ~isfolder(outdir), mkdir(outdir); end

    % Reproducible AWGN without changing the caller's RNG state.
    old_rng = rng;
    cleanup_rng = onCleanup(@()rng(old_rng)); %#ok<NASGU>
    rng(20261110,'twister');

    T = table();
    example_maps = struct();
    example_maps.Baseline = Q0;
    ref_rms = sqrt(mean(abs(H0(base)).^2));
    assert(ref_rms>0,'Zero baseline harmonic signal.');

    % Baseline reference.
    record_case('Baseline',0,1,H0,Q0);

    % Noise is injected into the extracted HARMONIC displacement, not RF.
    for snr = snr_db
        noise_rms = ref_rms*10^(-snr/20);
        for r=1:noise_repetitions
            eta = (noise_rms/sqrt(2)) * ...
                (randn(Nz,Nx)+1i*randn(Nz,Nx));
            H = H0+eta;
            [Q,~,~] = phase_planarity_quality_stress(H,win,minSup,ampFrac);
            record_case('HarmonicAWGN_dB',snr,r,H,Q);
            if snr==10 && r==1, example_maps.Noise10dB=Q; end
            clear eta H Q;
        end
    end

    % Blur applied to harmonic displacement (not a new SWS reconstruction).
    for sigma_mm = blur_mm
        if exist('imgaussfilt','file')==2
            H = imgaussfilt(real(H0),[sigma_mm/dz sigma_mm/dx], ...
                'Padding','symmetric') + 1i* ...
                imgaussfilt(imag(H0),[sigma_mm/dz sigma_mm/dx], ...
                'Padding','symmetric');
        else
            H = gaussian_blur_fallback(H0,sigma_mm/dz,sigma_mm/dx);
        end
        [Q,~,~] = phase_planarity_quality_stress(H,win,minSup,ampFrac);
        record_case('SpatialBlur_mm',sigma_mm,1,H,Q);
        if sigma_mm==1.5, example_maps.Blur15mm=Q; end
        clear H Q;
    end

    % Negative control: a perfectly linear phase ramp is locally planar.
    for beta = phase_ramp_rad_mm
        H = H0.*exp(1i*beta*(Xmm-mean(x*1e3)));
        [Q,~,~] = phase_planarity_quality_stress(H,win,minSup,ampFrac);
        record_case('LinearPhaseRamp_radPerMm',beta,1,H,Q);
        if beta==0.5, example_maps.Ramp05=Q; end
        clear H Q;
    end

    % Local phase perturbation with a known spatial location.
    local_shape = exp(-0.5*((Xmm-defect_center_mm(1)).^2+ ...
                            (Zmm-defect_center_mm(2)).^2) / ...
                          (defect_sigma_mm^2));
    for peak = defect_amplitudes_rad
        H = H0.*exp(1i*peak*local_shape);
        [Q,~,~] = phase_planarity_quality_stress(H,win,minSup,ampFrac);
        record_case('LocalPhaseDefect_rad',peak,1,H,Q);
        if peak==3.0, example_maps.Defect3rad=Q; end
        clear H Q;
    end

    disp(T);
    writetable(T,fullfile(outdir,'phase_stress_summary.csv'));
    save(fullfile(outdir,'phase_stress_summary.mat'), ...
        'T','snr_db','noise_repetitions','blur_mm', ...
        'phase_ramp_rad_mm','defect_amplitudes_rad', ...
        'Q0','example_maps','prf','f0','win','minSup','ampFrac', ...
        'region_names','-v7.3');

    % Diagnostic comparison: same color bar, only the Q maps differ.
    f = figure('Color','w','Name',['Phase-Q stress maps: ' phantom_case], ...
        'NumberTitle','off','Position',[50 80 1600 780]);
    tl=tiledlayout(f,2,4,'TileSpacing','compact','Padding','compact');
    map_plot(nexttile(tl),x,z,Q0,[0 1],'Baseline Q');
    map_plot(nexttile(tl),x,z,example_maps.Noise10dB,[0 1], ...
        'Harmonic noise: 10 dB');
    map_plot(nexttile(tl),x,z,example_maps.Blur15mm,[0 1], ...
        'Blur: 1.5 mm');
    map_plot(nexttile(tl),x,z,example_maps.Ramp05,[0 1], ...
        'Phase ramp: 0.5 rad/mm');
    map_plot(nexttile(tl),x,z,example_maps.Defect3rad,[0 1], ...
        'Local phase defect: 3 rad');
    map_plot(nexttile(tl),x,z,example_maps.Noise10dB-Q0,[-1 1], ...
        'Noise Q minus baseline');
    map_plot(nexttile(tl),x,z,example_maps.Blur15mm-Q0,[-1 1], ...
        'Blur Q minus baseline');
    map_plot(nexttile(tl),x,z,example_maps.Defect3rad-Q0,[-1 1], ...
        'Defect Q minus baseline');
    title(tl,['Controlled perturbations of harmonic displacement: ' ...
        strrep(phantom_case,'_',' ')]);
    exportgraphics(f,fullfile(outdir,'phase_stress_maps.png'),'Resolution',150);

    % Compact, comparable per-condition summary (global ROI only).
    G = T(T.Region=="Global",:);
    [keys,~,ic] = unique(string(G.Condition)+" "+string(G.Level), 'stable');
    means = accumarray(ic,G.MedianQ,[],@mean);
    frac_low = accumarray(ic,G.FractionBelow095,[],@mean);
    labels=cellstr(keys);
    f2=figure('Color','w','Name',['Phase-Q stress trends: ' phantom_case], ...
        'NumberTitle','off','Position',[100 100 1350 650]);
    tl2=tiledlayout(f2,2,1,'TileSpacing','compact');
    ax=nexttile(tl2);
    bar(ax,means);
    ylim(ax,[0 1]); ylabel(ax,'Median Q'); grid(ax,'on');
    xticks(ax,1:numel(keys)); xticklabels(ax,labels);
    xtickangle(ax,35); title(ax,'Global median phase-planarity Q');
    ax=nexttile(tl2);
    bar(ax,100*frac_low);
    ylim(ax,[0 100]); ylabel(ax,'Q below 0.95 [%]'); grid(ax,'on');
    xticks(ax,1:numel(keys)); xticklabels(ax,labels);
    xtickangle(ax,35);
    title(ax,'Fraction of Q below 0.95 (not SWS success rate)');
    exportgraphics(f2,fullfile(outdir,'phase_stress_trends.png'), ...
        'Resolution',150);

    fprintf('\nSaved files to: %s\n',outdir);
    fprintf(['Interpretation: this tests INPUT-signal sensitivity only. ' ...
       'It does not prove Q predicts perturbed SWS reconstruction error.\n']);

    % ---------------------------------------------------------------
    % Nested helper can read local masks, baseline, variables and table.
    % ---------------------------------------------------------------
    function record_case(condition,level,rep,H,Q)
        assert(isequal(size(Q),[Nz Nx]),'Unexpected Q map dimensions.');
        for region=1:numel(masks)
            mask=masks{region} & isfinite(Q);
            n=nnz(mask);
            if n<20
                fprintf('Skipping %s %s due to too few pixels.\n', ...
                    condition,region_names{region});
                continue;
            end
            q=double(Q(mask));
            d=q-double(Q0(mask));
            h=H(mask)-H0(mask);
            h0=H0(mask);
            nrms = sqrt(mean(abs(h).^2))/max(sqrt(mean(abs(h0).^2)),eps);
            row=table(string(phantom_case),string(condition),level,rep, ...
                string(region_names{region}),n,mean(q),median(q), ...
                std(q),mean(q<0.95),mean(q>=1-1e-6), ...
                mean(d),median(d),nrms, ...
                'VariableNames',{'Phantom','Condition','Level','Replicate', ...
                'Region','ValidPixels','MeanQ','MedianQ','StdQ', ...
                'FractionBelow095','FractionSaturated','MeanDeltaQ', ...
                'MedianDeltaQ','RelativeHarmonicPerturbation'});
            T=[T;row]; %#ok<AGROW>
        end
        fprintf('%-27s level %6.2f, rep %d: Q computed\n',condition,level,rep);
    end
end

function H = harmonic_phasor_stress(u,prf,f0)
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

function [Q,phaseRms,amp] = phase_planarity_quality_stress(H,win,minSup,ampFrac)
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

function H = gaussian_blur_fallback(H,sz,sx)
% Separable convolution fallback. Note: zero padding influences image edges.
    rz=max(1,ceil(3*sz));rx=max(1,ceil(3*sx));
    kz=exp(-0.5*((-rz:rz)/sz).^2);kz=kz(:)/sum(kz);
    kx=exp(-0.5*((-rx:rx)/sx).^2);kx=kx(:).'/sum(kx);
    H = conv2(kz,kx,real(H),'same')+ ...
        1i*conv2(kz,kx,imag(H),'same');
end

function map_plot(ax,x,z,Q,limits,ttl)
    imagesc(ax,x*1e3,z*1e3,Q);
    axis(ax,'image');set(ax,'YDir','normal');
    clim(ax,limits);colorbar(ax);
    xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
    title(ax,ttl,'Interpreter','none');
    try, ax.Toolbar.Visible='off'; catch, end
end
