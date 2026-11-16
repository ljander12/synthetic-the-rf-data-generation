function T = diagnose_kwave_phase_geometry(source_name)
%DIAGNOSE_KWAVE_PHASE_GEOMETRY Check whether k-Wave SWS bias is due to wave direction.
% This is a diagnostic calculation, NOT a new quality metric or a replacement
% for the PG/AC estimators. No new k-Wave simulation or USTB computation.
%
% INSTALL: synthetic-the-rf-data-generation/Analysis/
% RUN: T = diagnose_kwave_phase_geometry();
%      T = diagnose_kwave_phase_geometry('homogeneous_ground_truth_theta0.mat');
%
% Requires the saved outputs from analyze_existing_homogeneous_kwave.m.
% Computes independent harmonic phase slopes on the native 1-mm 2-D grid,
% using the SAME 4-frame moving average and SAME interior ROI.
%
% LIMITATION: phase slopes in a 2-D slice cannot determine out-of-plane
% wave propagation. Strong interference or non-monochromatic motion can
% also bias local phase slopes. Interpret this as a diagnostic only.

    if nargin < 1 || isempty(source_name)
        source_name = 'homogeneous_ground_truth.mat';
    end
    [~,stem,ext] = fileparts(char(string(source_name)));
    if isempty(ext), ext = '.mat'; end
    filename = [stem ext];
    root = fileparts(fileparts(mfilename('fullpath')));
    assert(isfolder(fullfile(root,'Analysis')) && isfolder(fullfile(root,'Data')), ...
        'Put the diagnostic in the repository Analysis folder.');

    paths = {fullfile(root,'Data',filename),fullfile(root,filename)};
    datafile = '';
    for k = 1:numel(paths)
        if isfile(paths{k}), datafile = paths{k}; break; end
    end
    if isempty(datafile)
        F = dir(fullfile(root,'**',filename));
        assert(~isempty(F), 'Cannot find %s anywhere below repository.',filename);
        datafile = fullfile(F(1).folder,F(1).name);
    end

    priorDir = fullfile(root,'Results','sws_estimator_diagnostic', ...
        'existing_homogeneous_kwave',stem);
    priorFile = fullfile(priorDir,'existing_homogeneous_results.mat');
    assert(isfile(priorFile), ...
        'Run T = analyze_existing_homogeneous_kwave(''%s'') first.',filename);
    S = load(priorFile,'x','z','PG','AC','GT','roi','meta');
    D = load(datafile,'disp_gt','vec_T','shear_params','prf');

    assert(all(isfield(D,{'disp_gt','vec_T','shear_params','prf'})), ...
        'Missing displacement, time, material or PRF information.');
    assert(all(isfield(S,{'x','z','PG','AC','GT','roi','meta'})), ...
        'Missing prior reconstructed maps / coordinates.');

    raw = double(D.disp_gt);
    [Nz,Nx,Nraw] = size(raw);
    x = double(S.x(:)); z = double(S.z(:));
    assert(isequal(size(S.PG),[Nz,Nx]) && numel(x)==Nx && numel(z)==Nz,...
        'Dimension mismatch between saved GT and previously reconstructed SWS.');
    assert(numel(D.vec_T)==Nraw && Nraw>=12, 'Time vector mismatch.');
    assert(max(abs(diff(x)-median(diff(x))))<1e-9 && ...
           max(abs(diff(z)-median(diff(z))))<1e-9, ...
           'This diagnostic assumes regularly spaced spatial axes.');
    dx = median(diff(x)); dz = median(diff(z));
    f0 = double(D.shear_params.source_freq);
    ctrue = double(D.shear_params.c_shear_bkg);
    PRF = double(D.prf);
    w = 2*pi*f0;
    tt0 = double(D.vec_T(:));
    assert(abs(median(diff(tt0))*PRF-1)<0.02, ...
        'Time increments do not match PRF. Inspect units of vec_T and PRF.');

    % Reproduce the input used in analyze_existing_homogeneous_kwave.m.
    packet = double(S.meta.packet_size);
    assert(packet==fix(packet) && packet>=1 && packet<Nraw, ...
        'Unexpected packet size.');
    Nt = Nraw-packet+1;
    u = zeros(Nz,Nx,Nt);
    t = zeros(Nt,1);
    for j=1:Nt
        u(:,:,j) = mean(raw(:,:,j:j+packet-1),3);
        t(j) = mean(tt0(j:j+packet-1));
    end
    clear raw;

    % Least-squares fit at the KNOWN source frequency. The coefficients
    % yield a complex spatial phasor U=C_cos - i*C_sin, with phase phi.
    % Adding an intercept reduces DC leakage for a short time record.
    B = [ones(Nt,1),cos(w*t),sin(w*t)];
    Y = reshape(u,[],Nt).'; % time x pixel
    fitCoeff = B\Y;
    complexU = reshape(fitCoeff(2,:) - 1i*fitCoeff(3,:),Nz,Nx);
    amplitude = abs(complexU);
    Yfit = B*fitCoeff;
    SSerr = sum(abs(Y-Yfit).^2,1);
    SSvar = sum(abs(Y-mean(Y,1)).^2,1);
    R2 = reshape(1-SSerr./max(SSvar,realmin),Nz,Nx);
    clear Yfit fitCoeff Y u;

    % Wrapped central differences of complex harmonic phase. We take
    % the phasor product first, preventing a global 2*pi unwrap error.
    % Sampling requirement: phase change over 2dx/2dz must stay < pi.
    kx = nan(Nz,Nx);
    kz = nan(Nz,Nx);
    kx(:,2:end-1) = angle(conj(complexU(:,1:end-2)) .* ...
        complexU(:,3:end)) ./ (x(3:end)-x(1:end-2)).';
    kz(2:end-1,:) = angle(conj(complexU(1:end-2,:)) .* ...
        complexU(3:end,:)) ./ (z(3:end)-z(1:end-2));
    kabs = hypot(kx,kz);
    cLateral = w ./ abs(kx);
    cPlane = w ./ kabs;
    directionDeg = atan2d(abs(kz),abs(kx));

    roi = logical(S.roi);
    ampRef = median(amplitude(roi & isfinite(amplitude)));
    ampMin = 0.10*ampRef; % Fixed, non-optimised exclusion of near-zero phase signal.
    common = roi & amplitude>ampMin & isfinite(R2) & ...
        isfinite(cLateral) & cLateral>0 & cLateral<30 & ...
        isfinite(cPlane) & cPlane>0 & cPlane<30 & ...
        isfinite(S.PG) & S.PG>0 & S.PG<30 & ...
        isfinite(S.AC) & S.AC>0 & S.AC<30;
    assert(nnz(common)>=20, 'Too few common reliable phase pixels (%d).',nnz(common));

    names = ["PG / previous"; "AC / previous"; ...
        "Direct lateral phase"; "Direct 2D phase"];
    maps = {double(S.PG),double(S.AC),cLateral,cPlane};
    med = nan(4,1); bias = nan(4,1); rmse = nan(4,1);
    for j=1:4
        v = maps{j}(common);
        med(j) = median(v);
        bias(j) = mean(v-ctrue);
        rmse(j) = sqrt(mean((v-ctrue).^2));
    end
    T = table(names,repmat(ctrue,4,1),med,bias,rmse, ...
        100*rmse/ctrue,'VariableNames', ...
        {'Method','TrueSWS_mps','MedianSWS_mps','Bias_mps', ...
         'RMSE_mps','RelativeRMSE_percent'});

    medianAbsKx = median(abs(kx(common)));
    medianAbsKz = median(abs(kz(common)));
    medianK = median(kabs(common));
    medianDirection = median(directionDeg(common));
    medianR2 = median(R2(common));
    finiteFrac = nnz(common)/nnz(roi);

    fprintf('\n=== HOMOGENEOUS k-WAVE PHASE GEOMETRY DIAGNOSTIC ===\n');
    fprintf('Input: %s\n',datafile);
    fprintf('Ground-truth c: %.4f m/s | f0: %.1f Hz | PRF: %.1f Hz\n', ...
        ctrue,f0,PRF);
    fprintf('dx=%.3f mm dz=%.3f mm, packet=%d, %d time frames\n', ...
        1e3*dx,1e3*dz,packet,Nt);
    fprintf('Common interior pixels: %d / %d (%.1f%%)\n', ...
        nnz(common),nnz(roi),100*finiteFrac);
    fprintf('Expected wavenumber: %.2f rad/m\n',w/ctrue);
    fprintf('Median |kx| = %.2f rad/m; median |kz| = %.2f rad/m\n', ...
        medianAbsKx,medianAbsKz);
    fprintf('Median 2D |k|=%.2f rad/m; median in-plane phase angle=%.1f deg\n', ...
        medianK,medianDirection);
    fprintf('Median 200-Hz harmonic fit R^2 = %.3f\n',medianR2);
    fprintf('Caution: a single 2D slice cannot resolve out-of-plane propagation.\n');
    disp(T);

    out = fullfile(priorDir,'phase_geometry_diagnostic');
    if ~isfolder(out),mkdir(out);end
    writetable(T,fullfile(out,'phase_geometry_summary.csv'));
    save(fullfile(out,'phase_geometry_results.mat'), ...
        'T','common','kx','kz','cLateral','cPlane', ...
        'directionDeg','R2','amplitude','x','z', ...
        'ctrue','f0','PRF','medianDirection','medianR2','-v7');

    fig=figure('Color','w','Position',[60 50 1500 850]);
    tl=tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
    arrays={S.GT,S.PG,S.AC,cPlane,cLateral,directionDeg,R2,amplitude};
    labels={'Material GT','USTB PG','USTB AC','Direct 2D phase speed', ...
        'Direct lateral phase speed','Phase tilt in x-z [deg]', ...
        '200-Hz fit R^2','200-Hz phasor amplitude'};
    for j=1:8
        ax=nexttile(tl,j);
        imagesc(ax,x*1e3,z*1e3,arrays{j});
        set(ax,'YDir','normal'); axis(ax,'image');colorbar(ax);
        title(ax,labels{j}); xlabel(ax,'x [mm]'); ylabel(ax,'z [mm]');
        if ismember(j,[1 2 3 4 5]), clim(ax,[0 8]); end
        if j==6,clim(ax,[0 90]); end
        if j==7,clim(ax,[0 1]);end
    end
    title(tl,sprintf('Homogeneous k-Wave | phase direction check | %s',stem), ...
        'Interpreter','none');
    exportgraphics(fig,fullfile(out,'phase_geometry_maps.png'), ...
        'Resolution',160);
    close(fig);
    fprintf('Saved to: %s\n',out);
end
