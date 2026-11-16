function [T, top50, detail] = diagnose_high_contrast_quality_blindspots()
%DIAGNOSE_HIGH_CONTRAST_QUALITY_BLINDSPOTS
% Are high SWS errors in the HIGH_CONTRAST inclusion assigned HIGH quality?
%
% DIAGNOSTIC OF THREE EXISTING METRICS -- NOT A FOURTH METRIC.
% Correlation Q (adapted), phase-planarity Q (experimental), RMI-1D.
%
% No simulations / USTB reconstruction / Q recalculation. All source
% SWS and Q maps are READ ONLY. Ground truth is ONLY used for validation.
%
% Region score inspection:
%  (a) Mean/median Q inside inclusion vs background (descriptive).
%  (b) At Q >= 0.8 or 0.9: fraction of BAD pixels classified HIGH Q,
%      plus error rate among the HIGH-Q pixels. These are illustrative
%      UNCALIBRATED thresholds, not a fair common reliability guarantee.
%  (c) Selection of TOP 50% globally by Q, then calculate what fraction
%      of inclusion pixels (and bad inclusion pixels) get rejected.
%      Uses fractional selection at ties -- no GT-based tie-breaking.
%
% Sources (already created by prior diagnostics):
% Results/sws_measurement_quality/high_contrast/
%   phase_consistency_diagnostic/phase_quality_results.mat
%   rmi_1d_diagnostic/rmi_1d_results.mat
%
% OUTPUTS (new folder, originals untouched):
% Results/sws_measurement_quality/high_contrast/quality_blindspot_diagnostic/
%   regional_Q_summary.csv
%   fixed_Q_threshold_diagnostic.csv
%   global_top50_inclusion_detection.csv
%   quality_blindspot_maps.png
%   high_Q_wrong_SWS_maps.png
%   blindspot_results.mat
%
% Usage: [T, top50, detail] = diagnose_high_contrast_quality_blindspots();

    here = fileparts(mfilename('fullpath'));
    if isfolder(fullfile(here, 'Results'))
        repo = here;
    else
        repo = fileparts(here);
    end
    root = fullfile(repo, 'Results', 'sws_measurement_quality', 'high_contrast');
    phase_file = fullfile(root, 'phase_consistency_diagnostic', 'phase_quality_results.mat');
    rmi_file = fullfile(root, 'rmi_1d_diagnostic', 'rmi_1d_results.mat');
    assert(isfile(phase_file), 'Missing Q results: %s', phase_file);
    assert(isfile(rmi_file), 'Missing RMI results: %s', rmi_file);

    % Minimal loads; do not load large raw displacement / USTB objects.
    P = load(phase_file, 'Qold', 'Qphase', 'common', 'GT');
    R = load(rmi_file, 'RmiPG', 'RmiAC', 'GT', 'cPG', 'cAC', ...
        'basic', 'masks', 'regions', 'x', 'z');
    assert(all(isfield(P, {'Qold','Qphase','common','GT'})), ...
        'Phase file is missing one or more required maps.');
    assert(all(isfield(R, {'RmiPG','RmiAC','GT','cPG','cAC', ...
        'basic','masks','regions','x','z'})), ...
        'RMI file is missing one or more required maps.');

    gt = double(R.GT);
    qold = double(P.Qold);
    qphase = double(P.Qphase);
    rmi = {double(R.RmiPG), double(R.RmiAC)};
    speed = {double(R.cPG), double(R.cAC)};
    common = logical(R.basic) & logical(P.common);
    sz = size(gt);
    assert(numel(sz) == 2, 'Expected 2D ground-truth SWS map.');
    assert(isequal(size(P.GT), sz, size(qold), sz, size(qphase), sz, ...
        size(rmi{1}), sz, size(rmi{2}), sz, ...
        size(speed{1}), sz, size(speed{2}), sz, size(common), sz), ...
        'Maps have different sizes. Cannot compare their pixels.');
    mismatch = abs(gt(isfinite(gt) & isfinite(P.GT)) - ...
        double(P.GT(isfinite(gt) & isfinite(P.GT))));
    assert(~isempty(mismatch) && max(mismatch) < 1e-6, ...
        'Ground-truth SWS differs between saved diagnostic files.');
    assert(numel(R.regions) == numel(R.masks), ...
        'RMI region names / masks disagree.');
    region_labels = string(R.regions(:));
    ii = find(region_labels == "InclusionCore", 1);
    bi = find(region_labels == "BackgroundInterior", 1);
    assert(~isempty(ii) && ~isempty(bi), ...
        'Could not find InclusionCore / BackgroundInterior masks.');
    inc = logical(R.masks{ii}) & common;
    bg = logical(R.masks{bi}) & common;
    common = common & isfinite(gt) & isfinite(qold) & isfinite(qphase) & ...
        isfinite(rmi{1}) & isfinite(rmi{2}) & ...
        isfinite(speed{1}) & isfinite(speed{2});
    inc = inc & common;
    bg = bg & common;
    assert(nnz(inc) > 100 && nnz(bg) > 100 && nnz(common) > 1000, ...
        'Insufficient common pixels for inclusion / background analysis.');
    assert(~any(inc(:) & bg(:)), 'Inclusion and background masks overlap.');
    assert(numel(R.x) == sz(2) && numel(R.z) == sz(1), ...
        'Coordinates and map dimensions disagree.');
    assert(all(diff(R.x(:)) > 0) && all(diff(R.z(:)) > 0), ...
        'Expected increasing x/z axes.');

    % Set in advance. Values 0.8/0.9 are illustrative (not calibrated).
    qlimits = [0.8 0.9];
    error_limits = [0.5 1.0]; % Absolute SWS error in m/s.
    metric_names = ["Correlation Q", "Phase Planarity Q", "RMI-1D"];
    methods = ["PG", "AC"];
    region_names = ["Global", "InclusionCore", "BackgroundInterior"];
    region_masks = {common, inc, bg};
    T = table();
    detail = table();
    top50 = table();

    % Sanity-check tie handling for scores saturated at 1.
    [testw, testcut] = select_top_fraction([1;1;1;0], 0.5);
    assert(abs(sum(testw)-2) < 1e-12 && testcut==1 && ...
        all(abs(testw(1:3)-2/3)<1e-12), ...
        'Top-fraction plateau handling failed its self-test.');

    for mi = 1:2
        err_map = abs(speed{mi} - gt);
        qmaps = {qold, qphase, rmi{mi}};
        for qi = 1:3
            qmap = qmaps{qi};
            qin = qmap(inc);
            qbg = qmap(bg);
            for ri = 1:3
                mask = region_masks{ri};
                q = qmap(mask);
                e = err_map(mask);
                assert(~isempty(q) && all(isfinite(q)) && ...
                    all(isfinite(e)), 'Invalid common pixels detected.');
                assert(all(q>=-1e-6 & q<=1+1e-6), ...
                    'Quality score outside [0,1] for %s.', metric_names(qi));
                for ei = 1:numel(error_limits)
                    eps_mps = error_limits(ei);
                    bad = e > eps_mps;
                    for ti = 1:numel(qlimits)
                        qmin = qlimits(ti);
                        high = q >= qmin;
                        nbad = nnz(bad);
                        nhigh = nnz(high);
                        if nbad==0
                            high_among_bad = NaN;
                        else
                            high_among_bad = nnz(high & bad) / nbad;
                        end
                        if nhigh==0
                            bad_among_high = NaN;
                        else
                            bad_among_high = nnz(high & bad) / nhigh;
                        end
                        row = table(methods(mi), metric_names(qi), ...
                            region_names(ri), eps_mps, qmin, numel(q), ...
                            mean(bad), mean(high), high_among_bad, ...
                            bad_among_high, ...
                            'VariableNames', {'Estimator','Metric','Region', ...
                            'ErrorLimit_mps','Qthreshold','N_Pixels', ...
                            'ErrorPrevalence','FractionHighQ', ...
                            'FractionWrongAssignedHighQ', ...
                            'ErrorRateAmongHighQ'});
                        detail = [detail; row]; %#ok<AGROW>
                    end
                end
                row = table(methods(mi),metric_names(qi),region_names(ri), ...
                    numel(q),mean(q),median(q),local_percentile(q,10),local_percentile(q,90), ...
                    mean(q>=1-1e-6),mean(e),mean(e>0.5), ...
                    mean(q>=0.8),mean(q>=0.9), ...
                    'VariableNames',{'Estimator','Metric','Region','N_Pixels', ...
                    'MeanQ','MedianQ','P10Q','P90Q','FractionQat1', ...
                    'MAE_mps','FractionErrorGT05', ...
                    'FractionQge08','FractionQge09'});
                T = [T;row]; %#ok<AGROW>
            end

            % Threshold-free, fair cross-metric comparator:
            % Select best 50% of ALL common pixels, then ask whether the
            % inclusion pixels are REJECTED. GT never enters selection.
            w = zeros(sz);
            [selected, cutoff] = select_top_fraction(qmap(common),0.5);
            w(common) = selected;
            badinc = inc & (err_map > 0.5);
            badglobal = common & (err_map > 0.5);
            inclusion_selected = sum(w(inc)) / nnz(inc);
            background_selected = sum(w(bg)) / nnz(bg);
            wrong_inclusion_rejected = 1 - ...
                sum(w(badinc)) / max(1,nnz(badinc));
            wrong_global_rejected = 1 - ...
                sum(w(badglobal)) / max(1,nnz(badglobal));
            row = table(methods(mi),metric_names(qi),0.5,cutoff, ...
                mean(qin),mean(qbg),nnz(inc),nnz(badinc), ...
                inclusion_selected,background_selected, ...
                1-inclusion_selected,wrong_inclusion_rejected, ...
                wrong_global_rejected, ...
                'VariableNames',{'Estimator','Metric','GlobalCoverage', ...
                'GlobalQCutoff','MeanQ_Inclusion','MeanQ_Background', ...
                'N_Inclusion','N_WrongInclusion', ...
                'FractionInclusionAccepted','FractionBackgroundAccepted', ...
                'FractionInclusionRejected', ...
                'FractionWrongInclusionRejected', ...
                'FractionWrongGlobalRejected'});
            top50 = [top50; row]; %#ok<AGROW>
        end
    end

    out = fullfile(root,'quality_blindspot_diagnostic');
    if ~isfolder(out),mkdir(out);end
    writetable(T,fullfile(out,'regional_Q_summary.csv'));
    writetable(detail,fullfile(out,'fixed_Q_threshold_diagnostic.csv'));
    writetable(top50,fullfile(out,'global_top50_inclusion_detection.csv'));
    save(fullfile(out,'blindspot_results.mat'),'T','detail','top50', ...
        'qlimits','error_limits','methods','metric_names','-v7');

    make_quality_maps(R.x, R.z,gt,speed,qold,qphase,rmi,common, ...
        methods,metric_names,fullfile(out,'quality_blindspot_maps.png'));
    make_error_high_q_maps(R.x,R.z,gt,speed,qold,qphase,rmi,common, ...
        methods,metric_names,0.5,0.8, ...
        fullfile(out,'high_Q_wrong_SWS_maps.png'));

    fprintf('\n=== HIGH CONTRAST: QUALITY BLIND-SPOT CHECK ===\n');
    fprintf('Same valid pixels: global %d; inclusion %d; background %d.\n', ...
        nnz(common),nnz(inc),nnz(bg));
    fprintf('Error > 0.5 m/s in inclusion: PG %.2f%%, AC %.2f%%.\n', ...
        100*mean(abs(speed{1}(inc)-gt(inc))>0.5), ...
        100*mean(abs(speed{2}(inc)-gt(inc))>0.5));
    fprintf('\n--- Mean and median Q INSIDE the faulty inclusion ---\n');
    disp(T(T.Region=="InclusionCore", ...
        {'Estimator','Metric','MeanQ','MedianQ','FractionQge08', ...
         'FractionQge09','FractionErrorGT05'}));
    fprintf('\n--- Bad pixels assigned Q >= 0.8, epsilon=0.5 m/s ---\n');
    subset = detail(detail.Region=="InclusionCore" & ...
        detail.ErrorLimit_mps==0.5 & detail.Qthreshold==0.8,:);
    disp(subset(:,{'Estimator','Metric','FractionHighQ', ...
        'FractionWrongAssignedHighQ','ErrorRateAmongHighQ'}));
    fprintf('\n--- GLOBAL top-50%% selection: does it reject inclusion? ---\n');
    disp(top50(:,{'Estimator','Metric','FractionInclusionRejected', ...
        'FractionWrongInclusionRejected','FractionBackgroundAccepted'}));
    fprintf('Saved outputs: %s\n',out);
    fprintf(['IMPORTANT: Q>=0.8 is only a descriptive, uncalibrated cutoff. ' ...
        'Global top-50%% is the fairer relative comparison.\n']);
    fprintf(['If ALL inclusion pixels exceed the error threshold, ' ...
        'regional ranking cannot reduce their error prevalence. ' ...
        'Evaluate rejection/absolute scores instead.\n']);
    fprintf(['GT is used ONLY to define validation regions, errors and ' ...
        'figure contours; NOT to compute Q or choose accepted pixels.\n']);
    fprintf(['RMI-1D used previous calibration from large_inclusions. ' ...
        'Results are exploratory; pixels are spatially correlated.\n']);
end

function [w,cutoff] = select_top_fraction(q,fraction)
% Fractionally accept cutoff Q-ties, independently of ground truth.
    q = double(q(:));
    n = numel(q);
    assert(n>0 && all(isfinite(q)) && fraction>0 && fraction<=1, ...
        'Invalid scores/selection fraction.');
    qsorted = sort(q,'descend');
    cutoff = qsorted(max(1,ceil(fraction*n)));
    better = q>cutoff;
    tied = q==cutoff;
    frac_ties = (fraction*n-nnz(better))/nnz(tied);
    assert(frac_ties>-1e-10 && frac_ties<1+1e-10, ...
        'Internal tie handling failed.');
    w = double(better)+max(0,min(1,frac_ties))*double(tied);
    assert(abs(sum(w)-fraction*n)<max(1e-7,1e-10*n), ...
        'Cannot achieve target accepted fraction.');
end

function make_quality_maps(x,z,gt,speed,qold,qphase,rmi,mask, ...
        methods,metric_names,target)
    xv=double(x(:)')*1e3; zv=double(z(:)')*1e3;
    gtbin = double(gt > min(gt(isfinite(gt))) + 1e-5);
    fig=figure('Color','w','Visible','off','Position',[100 100 1500 740], ...
        'Name','High contrast: existing Q maps versus SWS error');
    tl=tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
    for mi=1:2
        Q = {qold,qphase,rmi{mi}};
        err = abs(speed{mi}-gt);
        ax=nexttile(tl);
        image_panel(ax,xv,zv,err,mask,gtbin,[0 4], ...
            sprintf('%s |SWS - GT| [m/s]',methods(mi)));
        for qi=1:3
            ax=nexttile(tl);
            image_panel(ax,xv,zv,Q{qi},mask,gtbin,[0 1], ...
                sprintf('%s: %s',methods(mi),metric_names(qi)));
        end
    end
    title(tl,'HIGH CONTRAST: SWS errors and saved quality maps | white = GT inclusion boundary');
    exportgraphics(fig,target,'Resolution',145);
    close(fig);
end

function make_error_high_q_maps(x,z,gt,speed,qold,qphase,rmi,mask, ...
        methods,metric_names,eps_mps,qmin,target)
    xv=double(x(:)')*1e3; zv=double(z(:)')*1e3;
    gtbin=double(gt > min(gt(isfinite(gt))) + 1e-5);
    fig=figure('Color','w','Visible','off','Position',[100 100 1200 750], ...
        'Name','High-Q but incorrect SWS');
    tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    for mi=1:2
        Q={qold,qphase,rmi{mi}};
        bad=abs(speed{mi}-gt)>eps_mps;
        for qi=1:3
            ax=nexttile(tl);
            wrong_hi = double(bad & Q{qi}>=qmin);
            image_panel(ax,xv,zv,wrong_hi,mask,gtbin,[0 1], ...
                sprintf('%s %s | high Q + wrong SWS',methods(mi),metric_names(qi)));
        end
    end
    title(tl,sprintf(['High-Q wrong-SWS pixels: Q >= %.1f AND ' ...
        '|SWS-GT| > %.1f m/s | 1 = false reassurance'],qmin,eps_mps));
    exportgraphics(fig,target,'Resolution',145);
    close(fig);
end

function image_panel(ax,x,z,vals,mask,gtbin,limits,label)
    vals = double(vals);
    vals(~mask)=NaN;
    imagesc(ax,x,z,vals);
    set(ax,'YDir','normal');
    axis(ax,'image');
    caxis(ax,limits);
    colorbar(ax);
    hold(ax,'on');
    contour(ax,x,z,gtbin,[0.5 0.5], 'w','LineWidth',1.35);
    hold(ax,'off');
    title(ax,label,'Interpreter','none','FontSize',10);
    xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
end

function v = local_percentile(x,pct)
% Linear percentile interpolation without Statistics Toolbox.
    x=sort(double(x(:)));
    assert(~isempty(x) && pct>=0 && pct<=100);
    pos=1+(numel(x)-1)*pct/100;
    ilo=floor(pos); ihi=ceil(pos);
    v=x(ilo)+(pos-ilo)*(x(ihi)-x(ilo));
end
