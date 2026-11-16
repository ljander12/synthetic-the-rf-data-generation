function [T, ranking] = benchmark_sws_quality_false_confidence(phantom_cases)
%BENCHMARK_SWS_QUALITY_FALSE_CONFIDENCE
% Compare three EXISTING ground-truth-independent SWS quality metrics.
% This is a VALIDATION EXPERIMENT, not a fourth quality metric.
%
% Metrics: correlation Q (adapted), experimental phase planarity Q,
% RMI-1D (Choi et al. 2015, Eq. 6 adaptation).
% Reuses precomputed PG/AC SWS and Q maps: NO reconstruction, re-simulation,
% filtering, or changes to saved prior results.
%
% Primary endpoint (chosen before inspecting this benchmark):
%   P(abs(c_est-c_GT) > 0.5 m/s | selected TOP 50% by Q)
% Comparators: random acceptor at same coverage, other coverages 20%,80%,
% error threshold 1.0 m/s, and fixed inclusion/background regions.
% IMPORTANT: GT is used ONLY for offline evaluation of the selected pixels.
%
% Critical fair-comparison feature: ties are weighted fractionally.
% If a Q plateau crosses the coverage threshold, expected event rate under
% unbiased random tie-breaking is computed (essential for correlation Q).
% Selectors require Q values only; NEVER use GT to select pixels.
%
% Source data under Results/sws_measurement_quality/<phantom>/:
%   phase_consistency_diagnostic/phase_quality_results.mat
%   rmi_1d_diagnostic/rmi_1d_results.mat
% Optionally reads: risk_coverage_diagnostic/risk_coverage_summary.csv
%
% Run from a file placed in Analysis/:
%   [T, ranking] = benchmark_sws_quality_false_confidence();
%
% Writes to Results/sws_measurement_quality/final_quality_benchmark/:
%   false_confidence_all.csv
%   false_confidence_primary.csv
%   false_confidence_inclusion_primary.csv
%   primary_metric_ranking.csv
%   false_confidence_global.png
%   false_confidence_inclusion.png
%   benchmark_results.mat
%
% MATLAB R2020b+ recommended.

    if nargin < 1 || isempty(phantom_cases)
        phantom_cases = ["large_inclusions"; "high_contrast"; ...
            "low_contrast"; "chaoran_baseline"];
    end
    phantom_cases = string(phantom_cases(:));
    here = fileparts(mfilename('fullpath'));
    if isfolder(fullfile(here,'Results')), repo = here;
    else, repo = fileparts(here); end
    root = fullfile(repo,'Results','sws_measurement_quality');
    out = fullfile(root,'final_quality_benchmark');
    assert(isfolder(root),'Missing results folder: %s',root);
    if ~isfolder(out), mkdir(out); end

    covers = [0.20, 0.50, 0.80];
    errlimits = [0.5, 1.0]; % Absolute SWS error [m/s]. FIXED for all cases.
    methods = ["PG", "AC"];
    metrics = ["Correlation Q", "Phase Planarity Q", "RMI-1D"];
    units = ["Pixel", "Block1mm"];
    min_units = 20;
    T = table();
    self_test_ties();
    fprintf('\n=== FINAL QUALITY BENCHMARK: false confidence ===\n');
    fprintf('Primary endpoint: 50%% coverage, |SWS-GT| > 0.5 m/s.\n');
    fprintf('Fixed filter/reference: use EXISTING saved SWS/Q maps.\n');
    fprintf('Fractional weighting for equal Q values; spatial units NOT independent.\n');

    for pi=1:numel(phantom_cases)
        p = char(phantom_cases(pi));
        base = fullfile(root,p);
        rmi_file = fullfile(base,'rmi_1d_diagnostic','rmi_1d_results.mat');
        phase_file = fullfile(base,'phase_consistency_diagnostic', ...
            'phase_quality_results.mat');
        assert(isfile(rmi_file),'Missing RMI maps: %s',rmi_file);
        assert(isfile(phase_file),'Missing phase/correlation maps: %s',phase_file);

        R = load(rmi_file,'RmiPG','RmiAC','GT','cPG','cAC', ...
            'basic','masks','regions','x','z');
        P = load(phase_file,'Qold','Qphase');
        reqR = {'RmiPG','RmiAC','GT','cPG','cAC', ...
            'basic','masks','regions','x','z'};
        assert(all(isfield(R,reqR)) && all(isfield(P,{'Qold','Qphase'})), ...
            'Input MAT has missing maps for %s.',p);
        sz = size(R.GT);
        qall = {double(P.Qold),double(P.Qphase)};
        rmi = {double(R.RmiPG),double(R.RmiAC)};
        speeds = {double(R.cPG),double(R.cAC)};
        common = logical(R.basic);
        assert(isequal(size(qall{1}),sz,size(qall{2}),sz, ...
            size(rmi{1}),sz,size(rmi{2}),sz, ...
            size(speeds{1}),sz,size(speeds{2}),sz,size(common),sz), ...
            'Size mismatch in %s.',p);
        assert(numel(R.regions)==numel(R.masks), ...
            'Mismatch region/mask counts in %s.',p);
        assert(numel(R.x)==sz(2) && numel(R.z)==sz(1), ...
            'Coordinate dimension mismatch in %s.',p);
        gt = double(R.GT);
        [X,Z] = meshgrid(double(R.x(:)),double(R.z(:)));
        xbin = floor((X-min(X(:)))/0.001); % 1-mm blocks, SI axes
        zbin = floor((Z-min(Z(:)))/0.001);
        blockID = zbin * (max(xbin(:))+1) + xbin + 1;
        clear X Z xbin zbin;

        fprintf('\n%s: %d common valid pixels\n',p,nnz(common));
        for ri=1:numel(R.regions)
            region = string(R.regions{ri});
            mask = logical(R.masks{ri}) & common;
            if nnz(mask)<min_units
                fprintf('Skipping %s / %s: < %d pixels.\n',p,region,min_units);
                continue;
            end
            for mi=1:2
                e = abs(speeds{mi}(mask)-gt(mask));
                qs = {qall{1}(mask),qall{2}(mask),rmi{mi}(mask)};
                ids = blockID(mask);
                assert(all(isfinite(e)) && all(e>=0) && ...
                    all(cellfun(@(v)all(isfinite(v)),qs)), ...
                    'Nonfinite/invalid data in common %s/%s/%s.', ...
                    p,region,char(methods(mi)));
                for qi=1:3
                    for ui=1:2
                        if ui==1
                            q = double(qs{qi}(:));
                            er = double(e(:));
                        else
                            [q,er] = average_blocks(ids,qs{qi},e);
                        end
                        if numel(q)<min_units, continue; end
                        frac_sat = mean(q>=1-1e-6);
                        for ci=1:numel(covers)
                            cov = covers(ci);
                            [w, cutoff, tie_fraction] = top_fraction(q,cov);
                            accepted_mae = sum(w.*er)/sum(w);
                            for ei=1:numel(errlimits)
                                eps_mps = errlimits(ei);
                                error_event = double(er>eps_mps);
                                prevalence = mean(error_event);
                                false_confidence = sum(w.*error_event)/sum(w);
                                if prevalence>0
                                    rate_ratio = false_confidence/prevalence;
                                    gain = 100*(1-rate_ratio);
                                else
                                    rate_ratio = NaN;
                                    gain = NaN; % no large-error events exist
                                end
                                row = table(string(p),region,methods(mi), ...
                                    metrics(qi),units(ui),100*cov,eps_mps, ...
                                    numel(er),nnz(mask),mean(er), ...
                                    accepted_mae,prevalence,false_confidence, ...
                                    rate_ratio,gain,frac_sat,cutoff,tie_fraction, ...
                                    'VariableNames',{'Phantom','Region', ...
                                    'Estimator','Metric','Unit','Coverage_pct', ...
                                    'ErrorThreshold_mps','N_Units','N_Pixels', ...
                                    'FullMAE_mps','AcceptedMAE_mps', ...
                                    'RandomErrorRate','FalseConfidenceRate', ...
                                    'FC_over_Random','FC_Reduction_pct', ...
                                    'Q_SaturatedFraction','ScoreCutoff', ...
                                    'BoundaryTieFraction'});
                                T=[T;row]; %#ok<AGROW>
                            end
                        end
                    end
                end
            end
        end
    end

    % Exact same primary endpoint for every metric, estimator and phantom.
    primary = T(T.Region=="Global" & T.Unit=="Pixel" & ...
        T.Coverage_pct==50 & T.ErrorThreshold_mps==0.5,:);
    assert(~isempty(primary),'No primary benchmark rows were produced.');
    writetable(T,fullfile(out,'false_confidence_all.csv'));
    writetable(primary,fullfile(out,'false_confidence_primary.csv'));
    inclusion_primary = T(T.Region=="InclusionCore" & ...
        T.Unit=="Pixel" & T.Coverage_pct==50 & ...
        T.ErrorThreshold_mps==0.5,:);
    writetable(inclusion_primary, ...
        fullfile(out,'false_confidence_inclusion_primary.csv'));

    % Exploratory across-phantom consistency, NOT a statistical validation.
    % large_inclusions was used to calibrate RMI and is reported separately;
    % no tuning or refitting is permitted here.
    ranking = table();
    for qi=1:numel(metrics)
        metric_rows = primary(primary.Metric==metrics(qi) & ...
            primary.Phantom~="large_inclusions",:);
        good = isfinite(metric_rows.FC_over_Random);
        tested = metric_rows(good,:);
        if isempty(tested)
            mr = NaN; worst = NaN; mgain = NaN;
        else
            mr = mean(tested.FC_over_Random);
            worst = max(tested.FC_over_Random);
            mgain = mean(tested.FC_Reduction_pct);
        end
        incl = inclusion_primary(inclusion_primary.Metric==metrics(qi) & ...
            inclusion_primary.Phantom~="large_inclusions",:);
        incl = incl(isfinite(incl.FC_over_Random),:);
        if isempty(incl), incl_mean=NaN; incl_worst=NaN;
        else
            incl_mean=mean(incl.FC_over_Random);
            incl_worst=max(incl.FC_over_Random);
        end
        r = table(metrics(qi),height(tested),mr,worst,mgain, ...
            height(incl),incl_mean,incl_worst, ...
            'VariableNames',{'Metric','N_NonCalibrationComparisons', ...
            'MeanFC_over_Random','WorstFC_over_Random', ...
            'MeanFC_Reduction_pct', ...
            'N_InclusionComparisons','InclusionMeanFC_over_Random', ...
            'InclusionWorstFC_over_Random'});
        ranking=[ranking;r]; %#ok<AGROW>
    end
    ranking=sortrows(ranking,'MeanFC_over_Random','ascend');
    writetable(ranking,fullfile(out,'primary_metric_ranking.csv'));

    save(fullfile(out,'benchmark_results.mat'), ...
        'T','primary','inclusion_primary','ranking','covers','errlimits','metrics', ...
        'methods','units','phantom_cases','-v7.3');
    plot_overview(primary,phantom_cases,metrics,methods, ...
        'Global',fullfile(out,'false_confidence_global.png'));
    plot_overview(T(T.Unit=="Pixel" & T.Coverage_pct==50 & ...
        T.ErrorThreshold_mps==0.5,:), ...
        phantom_cases,metrics,methods,'InclusionCore', ...
        fullfile(out,'false_confidence_inclusion.png'));

    fprintf('\n=== PRIMARY: top 50%%, error > 0.5 m/s, GLOBAL PIXELS ===\n');
    disp(primary(:,{'Phantom','Estimator','Metric', ...
        'RandomErrorRate','FalseConfidenceRate','FC_over_Random', ...
        'FC_Reduction_pct'}));
    fprintf('\n=== EXPLORATORY METRIC RANKING (check inclusion columns) ===\n');
    disp(ranking);
    fprintf('Saved to: %s\n',out);
    fprintf(['Caution: RMI calibration used large_inclusions; ' ...
        'its performance there is IN-SAMPLE.\n']);
    fprintf(['Cases have been inspected during development: ' ...
        'ranking is EXPLORATORY, not independent validation.\n']);
    fprintf(['Threshold comparisons are per-case equal coverage. ' ...
        'Unchanged Q>=0.8 is NOT a fair cross-metric operating point.\n']);
    fprintf(['Tied score selections are fractional expectations; ' ...
        'pixels and spatial blocks are NOT independent observations.\n']);
    fprintf(['DO NOT select a winner from Global ranking alone: ' ...
        'inspect InclusionCore and previously saved AURC as well.\n']);
end

function [w,cutoff,tie_fraction] = top_fraction(q,cov)
% No GT input. Unbiased fractional inclusion of cutoff Q plateau.
    q=double(q(:));
    n=numel(q);
    assert(n>0 && cov>0 && cov<=1 && all(isfinite(q)), ...
        'Invalid Q or coverage.');
    k=max(1,min(n,ceil(cov*n)));
    sq=sort(q,'descend');
    cutoff=sq(k);
    better=(q>cutoff);
    tied=(q==cutoff);
    target=cov*n;
    nt=nnz(tied);
    assert(nt>=1,'No cutoff ties, internal error.');
    tie_fraction=max(0,min(1,(target-nnz(better))/nt));
    w=double(better)+tie_fraction*double(tied);
    assert(abs(sum(w)-target)<max(1e-7,1e-10*n), ...
        'Could not achieve exact target coverage.');
end

function [qb,eb] = average_blocks(blockid,q,e)
% Equal-weight 1-mm spatial blocks. Block event is mean MAE > epsilon.
    [~,~,g]=unique(double(blockid(:)));
    n=accumarray(g,1);
    qb=accumarray(g,double(q(:)))./n;
    eb=accumarray(g,double(e(:)))./n;
end

function self_test_ties()
% All tied scores must behave exactly like random acceptance.
    [w,~,tf]=top_fraction(ones(10,1),0.2);
    assert(all(abs(w-0.2)<1e-12) && abs(tf-0.2)<1e-12, ...
        'Tie weighting self-test failed.');
    [w,~,~]=top_fraction((10:-1:1)',0.3);
    assert(isequal(w,[ones(3,1);zeros(7,1)]), ...
        'Non-tied ranking self-test failed.');
end

function plot_overview(V, cases, metrics, methods, region, filename)
    f=figure('Color','w','Name',['False confidence: ' region], ...
        'Position',[80 80 1450 630]);
    t=tiledlayout(f,2,numel(cases),'TileSpacing','compact', ...
        'Padding','compact');
    for mi=1:2
        for pi=1:numel(cases)
            ax=nexttile(t);
            r=V(V.Phantom==cases(pi) & V.Region==string(region) & ...
                V.Estimator==methods(mi),:);
            rates=nan(1,numel(metrics));
            for qi=1:numel(metrics)
                rr=r(r.Metric==metrics(qi),:);
                if ~isempty(rr), rates(qi)=100*rr.FalseConfidenceRate(1); end
            end
            bar(ax,1:numel(metrics),rates);
            hold(ax,'on');
            if ~isempty(r)
                yline(ax,100*r.RandomErrorRate(1),'k--', ...
                    'Random selection','LabelHorizontalAlignment','left');
            end
            set(ax,'XTick',1:3,'XTickLabel',{'Corr','Phase','RMI'});
            ylabel(ax,'False confidence [%]');
            title(ax,sprintf('%s / %s',strrep(char(cases(pi)),'_',' '), ...
                char(methods(mi))),'Interpreter','none');
            ylim(ax,[0 100]); grid(ax,'on');
        end
    end
    title(t,sprintf(['False confidence: %s | top 50%% by Q | ' ...
        '|SWS-GT| > 0.5 m/s | lower is better'],region), ...
        'Interpreter','none');
    exportgraphics(f,filename,'Resolution',155);
    close(f);
end
