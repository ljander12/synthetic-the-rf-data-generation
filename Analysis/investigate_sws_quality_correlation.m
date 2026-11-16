function summary = investigate_sws_quality_correlation(phantom_case)
%INVESTIGATE_SWS_QUALITY_CORRELATION  Diagnose saturation of displacement-based Q.
%
% PURPOSE
%   Test whether our exploratory displacement-correlation quality map contains
%   useful information BEFORE applying the piecewise Q transform.
%   Compare lateral separations of 1, 4, 8, 16, and 32 pixels.
%
% REFERENCE
%   Schoen et al. (2025), Ultrasound in Medicine & Biology,
%   DOI: 10.1016/j.ultrasmedbio.2025.03.003
%
% IMPORTANT
%   This is a CONTROLLED EXTENSION OF OUR PROTOTYPE, NOT an exact
%   reproduction of the published method. It uses zero-lag, mean-centered,
%   normalized temporal correlation of lateral displacement waveforms.
%   Q and abs(rho) require NO ground truth. Ground truth is used only for
%   retrospective evaluation, never in construction of a quality map.
%   Do NOT select the "best" lag from these phantoms and report its performance
%   on the same phantoms as independent validation.
%
% INPUTS (already produced by the previous pilot)
%   Results/sws_measurement_quality/<phantom>/quality_validation.mat
%   Data/<phantom>_displacement_with_quality.mat (or Data/<phantom>/...)
%
% OUTPUT (existing files are NOT overwritten)
%   Results/sws_measurement_quality/<phantom>/correlation_diagnostic/
%       correlation_lag_summary.csv
%       q_lag_comparison.png
%       lag_metric_trends.png
%       correlation_lag_summary.mat
%
% USAGE: Put in Analysis/ folder, then:
%   investigate_sws_quality_correlation('large_inclusions')
%   investigate_sws_quality_correlation('high_contrast')

    if nargin < 1 || isempty(phantom_case)
        phantom_case = 'large_inclusions';
    end
    phantom_case = char(string(phantom_case));

    script_dir = fileparts(mfilename('fullpath'));
    if isfolder(fullfile(script_dir,'Results'))
        repo_dir = script_dir;
    else
        repo_dir = fileparts(script_dir);
    end

    prior_file = fullfile(repo_dir,'Results','sws_measurement_quality', ...
        phantom_case,'quality_validation.mat');
    assert(isfile(prior_file), ...
        'Missing previous Q pilot: %s. Run test_sws_measurement_quality first.', prior_file);
    S = load(prior_file,'Q','C12dB','GT','PG','AC','valid','x_axis','z_axis');
    assert(all(isfield(S,{'Q','C12dB','GT','PG','AC','valid','x_axis','z_axis'})), ...
        'The saved Q pilot does not contain all required variables.');
    GT = double(S.GT);
    PG = double(S.PG);
    AC = double(S.AC);
    [Nz,Nx] = size(GT);
    assert(isequal(size(S.Q),size(GT),size(S.C12dB), ...
        size(PG),size(AC),size(S.valid)), 'Saved map sizes differ.');

    % The OLD Q mapping uses the lower C12 from adjacent pairs for each
    % interior pixel. The saved variable "rho" instead averaged signed
    % pair correlations, so it is NOT the right inverse of the Q mapping.
    % Recover the relevant effective magnitude from the actual saved C12.
    abs_rho_prior = 10.^(double(S.C12dB)/20);
    fprintf('\n=== QUALITY CORRELATION DIAGNOSTIC: %s ===\n',phantom_case);
    prior_valid = logical(S.valid) & isfinite(S.C12dB) & isfinite(S.Q);
    fprintf('Saved pilot: median |rho| = %.5f; Q = 1 on %.2f%% of valid pixels\n', ...
        median(abs_rho_prior(prior_valid)),100*mean(S.Q(prior_valid)>=1-1e-10));

    % Class uff.beamformed_data must exist BEFORE the MAT file is loaded.
    if isempty(which('uff.beamformed_data'))
        ustb = fullfile(fileparts(repo_dir),'elastography', ...
            'ElastographyToolbox','USTB');
        if isfolder(ustb)
            addpath(genpath(ustb));
            rehash;
        end
    end
    assert(~isempty(which('uff.beamformed_data')), ...
        ['Cannot find uff.beamformed_data. Add the USTB folder to ' ...
         'MATLAB path, e.g. elastography/ElastographyToolbox/USTB.']);

    base_name = [phantom_case '_displacement_with_quality.mat'];
    paths = {fullfile(repo_dir,'Data',base_name), ...
             fullfile(repo_dir,'Data',phantom_case,base_name), ...
             fullfile(repo_dir,base_name), ...
             fullfile(repo_dir,'Analysis',base_name)};
    displacement_file = '';
    for p=1:numel(paths)
        if isfile(paths{p})
            displacement_file=paths{p};
            break;
        end
    end
    assert(~isempty(displacement_file),'Cannot find %s',base_name);
    fprintf('Loading displacement: %s\n',displacement_file);
    D = load(displacement_file,'displacement_estimation_autocorr');
    assert(isfield(D,'displacement_estimation_autocorr') && ...
        (isobject(D.displacement_estimation_autocorr) || ...
         isstruct(D.displacement_estimation_autocorr)), ...
        'USTB displacement object could not be restored.');
    raw = D.displacement_estimation_autocorr.data;
    assert(isnumeric(raw) && isreal(raw) && mod(numel(raw),Nz*Nx)==0, ...
        'Unexpected displacement size/type.');
    Nt = numel(raw)/(Nz*Nx);
    u = reshape(single(raw),Nz,Nx,Nt);
    clear raw D
    fprintf('Displacement array: %d x %d x %d\n',Nz,Nx,Nt);

    lags = [1 4 8 16 32];
    lags = lags(lags < floor(Nx/3));
    dx_mm = median(abs(diff(double(S.x_axis(:)))))*1e3;
    Qmaps = cell(numel(lags),1);
    absRhoMaps = cell(numel(lags),1);
    for j=1:numel(lags)
        L = lags(j);
        fprintf('Calculating lateral correlation, lag %d pixels (%.3f mm)...\n', ...
            L,L*dx_mm);
        [Qmaps{j},absRhoMaps{j}] = correlate_at_lag(u,L);
    end
    clear u

    % Use EXACTLY the same spatial region and valid samples for every lag,
    % to prevent a changing ROI from explaining changes in correlation.
    maxlag = max(lags);
    common = logical(S.valid);
    common(:,1:maxlag) = false;
    common(:,Nx-maxlag+1:end) = false;
    for j=1:numel(lags)
        common = common & isfinite(Qmaps{j}) & isfinite(absRhoMaps{j});
    end
    assert(nnz(common)>=100,'Too few common valid pixels across all lags.');
    fprintf('Common evaluation pixels (same for all lags): %d\n',nnz(common));

    errPG = abs(PG-GT);
    errAC = abs(AC-GT);
    summary = table();
    for j=1:numel(lags)
        Q = Qmaps{j};
        absRho = absRhoMaps{j};
        qv = double(Q(common));
        rv = double(absRho(common));
        for method = 1:2
            if method==1
                label = "PG";
                err = double(errPG(common));
            else
                label = "AC";
                err = double(errAC(common));
            end
            row = table(string(phantom_case),label,lags(j),lags(j)*dx_mm, ...
                nnz(common),median(rv),mean(qv),mean(qv>=1-1e-10), ...
                spearman_safe(qv,err),spearman_safe(rv,err), ...
                'VariableNames',{'Phantom','Method','LagPixels','LagMM', ...
                'ValidPixels','MedianAbsRho','MeanQ','FractionQSaturated', ...
                'SpearmanQvsAbsError','SpearmanAbsRhovsAbsError'});
            summary = [summary;row]; %#ok<AGROW>
        end
    end
    disp(summary);

    out_dir = fullfile(repo_dir,'Results','sws_measurement_quality', ...
        phantom_case,'correlation_diagnostic');
    if ~isfolder(out_dir), mkdir(out_dir); end
    writetable(summary,fullfile(out_dir,'correlation_lag_summary.csv'));

    % Fig 1: directly inspect the raw correlation and the saturation.
    idx1 = 1;
    idx2 = numel(lags);
    figure('Name',['Raw correlation diagnostic: ' phantom_case], ...
        'Color','w','NumberTitle','off');
    tl=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    x_mm=double(S.x_axis(:))*1e3;
    z_mm=double(S.z_axis(:))*1e3;
    plotmap(nexttile(tl),x_mm,z_mm,absRhoMaps{idx1},[0 1], ...
        sprintf('|rho|: lag %d (%.2f mm)',lags(idx1),lags(idx1)*dx_mm));
    plotmap(nexttile(tl),x_mm,z_mm,Qmaps{idx1},[0 1], ...
        sprintf('Q: lag %d',lags(idx1)));
    plotmap(nexttile(tl),x_mm,z_mm,errPG,[0 4], ...
        '|PG - GT| [m/s]');
    plotmap(nexttile(tl),x_mm,z_mm,absRhoMaps{idx2},[0 1], ...
        sprintf('|rho|: lag %d (%.2f mm)',lags(idx2),lags(idx2)*dx_mm));
    plotmap(nexttile(tl),x_mm,z_mm,Qmaps{idx2},[0 1], ...
        sprintf('Q: lag %d',lags(idx2)));
    plotmap(nexttile(tl),x_mm,z_mm,errAC,[0 4], ...
        '|AC - GT| [m/s]');
    title(tl,['Quality correlations: ' strrep(phantom_case,'_',' ')]);
    exportgraphics(gcf,fullfile(out_dir,'q_lag_comparison.png'),'Resolution',160);

    % Fig 2: quantify if changing separation makes Q more informative.
    fig=figure('Name',['Lag sensitivity: ' phantom_case], ...
        'Color','w','NumberTitle','off');
    tiledlayout(fig,1,2,'Padding','compact');
    ax1=nexttile;
    pg=summary(summary.Method=="PG",:);
    ac=summary(summary.Method=="AC",:);
    plot(ax1,pg.LagMM,100*pg.FractionQSaturated,'-o','LineWidth',1.4);
    grid(ax1,'on'); xlabel(ax1,'Lateral separation [mm]');
    ylabel(ax1,'Pixels with Q=1 [%]'); title(ax1,'Q saturation');
    ylim(ax1,[0 100]);
    ax2=nexttile;
    plot(ax2,pg.LagMM,pg.SpearmanQvsAbsError,'-o','LineWidth',1.4);
    hold(ax2,'on');
    plot(ax2,ac.LagMM,ac.SpearmanQvsAbsError,'-s','LineWidth',1.4);
    yline(ax2,0,':'); hold(ax2,'off'); grid(ax2,'on');
    xlabel(ax2,'Lateral separation [mm]');
    ylabel(ax2,'Spearman(Q, absolute SWS error)');
    title(ax2,'Does Q rank SWS errors?'); legend(ax2,'PG','AC','Location','best');
    exportgraphics(fig,fullfile(out_dir,'lag_metric_trends.png'),'Resolution',160);

    save(fullfile(out_dir,'correlation_lag_summary.mat'), ...
        'summary','lags','dx_mm','common','Qmaps','absRhoMaps','-v7.3');
    fprintf('\nDONE. Saved new files to:\n%s\n',out_dir);
    fprintf(['Interpretation: Q=1 and Spearman near 0 imply weak predictive ' ...
        'ability of THIS prototype, not proof the published method fails.\n']);
end

function [Q,absRhoMap] = correlate_at_lag(u,L)
% Zero-lag correlation of mean-centered displacement time series at
% laterally separated positions. No GT or SWS estimator output is used.
    [Nz,Nx,Nt]=size(u);
    W=Nx-L;
    sa=zeros(Nz,W); sb=sa; saa=sa; sbb=sa; sab=sa;
    for t=1:Nt
        a=double(u(:,1:W,t));
        b=double(u(:,1+L:Nx,t));
        sa=sa+a; sb=sb+b;
        saa=saa+a.^2; sbb=sbb+b.^2; sab=sab+a.*b;
    end
    va=max(0,saa-sa.^2/Nt);
    vb=max(0,sbb-sb.^2/Nt);
    den=sqrt(va.*vb);
    pair=sab-sa.*sb/Nt;
    rho=pair./den;
    rho(den<=eps(max(den(:))))=NaN;
    rho=max(-1,min(1,rho));

    % Mask signals that are extremely weak relative to the dataset.
    pairRms=sqrt((saa+sbb)/(2*Nt));
    good=pairRms(isfinite(pairRms)&pairRms>0);
    if ~isempty(good)
        rho(pairRms<median(good)*1e-4)=NaN;
    end
    rhoAbs=abs(rho);
    cdb=20*log10(max(rhoAbs,realmin));
    q=nan(Nz,W);
    q(cdb>=-10)=1;
    m=cdb>-50 & cdb<-10;
    q(m)=cos((pi/2)*(cdb(m)+10)/40);
    q(cdb<=-50)=0;

    % Assign each pair to the pixel nearest its midpoint.
    middle=1+floor(L/2);
    cols=middle:(middle+W-1);
    Q=nan(Nz,Nx);
    absRhoMap=nan(Nz,Nx);
    Q(:,cols)=q;
    absRhoMap(:,cols)=rhoAbs;
end

function r=spearman_safe(a,b)
    a=double(a(:)); b=double(b(:));
    ok=isfinite(a)&isfinite(b);
    a=a(ok); b=b(ok);
    if numel(a)<20 || std(a)<1e-10 || std(b)<1e-10
        r=NaN;
    else
        r=corr(a,b,'Type','Spearman');
    end
end

function plotmap(ax,x,z,img,limits,ttl)
    imagesc(ax,x,z,img);
    axis(ax,'image'); set(ax,'YDir','normal');
    clim(ax,limits); colorbar(ax);
    xlabel(ax,'Lateral [mm]'); ylabel(ax,'Depth [mm]');
    title(ax,ttl,'Interpreter','none');
end
