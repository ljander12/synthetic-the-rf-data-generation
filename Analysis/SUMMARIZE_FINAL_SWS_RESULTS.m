function T = SUMMARIZE_FINAL_SWS_RESULTS()
% SUMMARIZE_FINAL_SWS_RESULTS
% One-run READ-ONLY summary of existing PG/AC elastography results.
%
% NO k-Wave, Field II, USTB, displacement estimation, or SWS reconstruction.
% NO ground-truth-independent metrics; this summarizes existing GT metrics.
%
% Preferred input from earlier compare_SWS_methods_batch.m:
%   Results/method_comparison/batch/batch_summary.csv
% Alternative input: batch_summary.mat or case-by-case method_metrics.csv
%
% Run:
%   T = SUMMARIZE_FINAL_SWS_RESULTS();
%
% Results:
%   Results/final_sws_evaluation/final_sws_comparison.csv
%   Results/final_sws_evaluation/global_regional_rmse.png
%   Results/final_sws_evaluation/contrast_metrics.png
%   Results/final_sws_evaluation/final_sws_table.tex
%   Results/final_sws_evaluation/key_findings.txt
%
% The selected four phantoms are the previously analysed cases.

    cases = ["chaoran_baseline", "low_contrast", ...
             "high_contrast", "large_inclusions"];
    methods = ["PG", "AC"];
    analysis_dir = fileparts(mfilename('fullpath'));
    repo_dir = fileparts(analysis_dir);
    results_root = fullfile(repo_dir, 'Results');
    input_root = fullfile(results_root, 'method_comparison', 'batch');
    outdir = fullfile(results_root, 'final_sws_evaluation');
    if ~isfolder(outdir), mkdir(outdir); end

    %% Load existing results ONLY (never recompute or overwrite original data)
    csv_file = fullfile(input_root, 'batch_summary.csv');
    mat_file = fullfile(input_root, 'batch_summary.mat');
    if isfile(csv_file)
        Traw = readtable(csv_file, 'TextType', 'string');
        source_file = csv_file;
    elseif isfile(mat_file)
        S = load(mat_file, 'batch_summary');
        assert(isfield(S,'batch_summary'), 'Missing batch_summary in MAT file.');
        Traw = S.batch_summary;
        source_file = mat_file;
    else
        tables = cell(numel(cases),1);
        for i = 1:numel(cases)
            fn = fullfile(input_root, char(cases(i)), 'method_metrics.csv');
            if isfile(fn)
                tables{i} = readtable(fn,'TextType','string');
            else
                fm = fullfile(input_root, char(cases(i)), 'method_comparison.mat');
                if isfile(fm)
                    D = load(fm, 'this_table');
                    if isfield(D,'this_table'), tables{i} = D.this_table; end
                end
            end
        end
        if any(cellfun(@isempty,tables))
            missing = cases(cellfun(@isempty,tables));
            fprintf('\nCannot locate existing PG/AC comparison results for:\n');
            fprintf('  %s\n',missing);
            fprintf('\nExpected: %s\n',csv_file);
            fprintf('You can generate the missing method-comparison results from\n');
            fprintf('SAVED displacement data (no k-Wave/Field II) using:\n');
            fprintf('  compare_SWS_methods_batch(["chaoran_baseline", "low_contrast", ...\n');
            fprintf('      "high_contrast", "large_inclusions"]);\n');
            error('Missing existing method-comparison data. No results were invented.');
        end
        Traw = vertcat(tables{:});
        source_file = input_root;
    end

    required = { 'Phantom', 'Method', 'GlobalRMSE_mps', ...
        'InclusionRMSE_mps', 'BackgroundRMSE_mps', ...
        'InclusionBias_mps', 'BackgroundBias_mps', ...
        'CNR', 'DirectKDE_gCNR', 'ContrastRecovery', ...
        'EdgeAccuracy_mm', 'EdgeSharpness_mps_per_mm' };
    absent = setdiff(required, Traw.Properties.VariableNames);
    assert(isempty(absent), 'Missing result columns: %s', strjoin(absent, ', '));

    phantom_names = strtrim(string(Traw.Phantom));
    method_names = strtrim(lower(string(Traw.Method)));
    method_key = strings(height(Traw),1);
    method_key(method_names == "phase gradient" | method_names == "pg") = "PG";
    method_key(method_names == "spatial autocorrelation" | ...
        method_names == "autocorrelation" | method_names == "ac") = "AC";
    assert(~any(method_key == "" & ismember(phantom_names,cases)), ...
        'Unrecognized PG/AC method label.');

    % Fixed ordering prevents mismatched PG/AC comparisons.
    order = zeros(8,1);
    p = strings(8,1); m = strings(8,1);
    k = 0;
    for i = 1:numel(cases)
        for j = 1:numel(methods)
            k = k+1;
            hits = find(phantom_names == cases(i) & method_key == methods(j));
            assert(numel(hits) == 1, ...
                'Expected exactly one row for %s / %s, found %d.', ...
                cases(i),methods(j),numel(hits));
            order(k) = hits;
            p(k) = cases(i); m(k) = methods(j);
        end
    end
    T = Traw(order,:);
    T.Phantom = p;
    T.Method = m;

    if ismember('Frequency_Hz',T.Properties.VariableNames)
        f = double(T.Frequency_Hz);
        if any(isfinite(f) & abs(f-200)>1e-6)
            warning('Some rows have frequencies different from 200 Hz. Verify comparability.');
        end
    end
    if ismember('GlobalPixels',T.Properties.VariableNames)
        for i=1:2:height(T)
            if T.GlobalPixels(i) ~= T.GlobalPixels(i+1)
                warning('PG/AC evaluated on different pixel counts in %s.', T.Phantom(i));
            end
        end
    end

    % In the existing comparison scripts ContrastRecovery is a fraction.
    contrast = double(T.ContrastRecovery);
    if all(abs(contrast(isfinite(contrast))) <= 2)
        T.ContrastRecovery_percent = 100*contrast;
    else
        error(['ContrastRecovery appears not to be stored as a fraction. ' ...
               'Check units before summarizing; no automatic mixed-unit conversion.']);
    end
    % Edge metrics with NaN mean not evaluated (e.g. boundary undetected).
    if any(isnan(double(T.EdgeAccuracy_mm)))
        fprintf('Edge accuracy contains NaN: these values must NOT be interpreted as zero.\n');
    end

    csv_out = fullfile(outdir,'final_sws_comparison.csv');
    writetable(T, csv_out);
    fprintf('\n=== FINAL SWS COMPARISON (reported existing data) ===\n');
    disp(T(:,{'Phantom','Method','GlobalRMSE_mps','InclusionRMSE_mps', ...
        'BackgroundRMSE_mps','InclusionBias_mps','BackgroundBias_mps', ...
        'CNR','DirectKDE_gCNR','ContrastRecovery_percent', ...
        'EdgeAccuracy_mm','EdgeSharpness_mps_per_mm'}));

    %% Plots using common y-axes across both estimators
    x = 1:4;
    legends = {'Phase gradient','Autocorrelation'};
    ticknames = {'Baseline','Low contrast','High contrast','Large inclusions'};
    f1=figure('Color','w','Name','Final SWS RMSE comparison', ...
        'Position',[100 100 1280 430]);
    metrics = {'GlobalRMSE_mps','InclusionRMSE_mps','BackgroundRMSE_mps'};
    titles = {'Global RMSE','Inclusion RMSE','Background RMSE'};
    for col = 1:3
        subplot(1,3,col);
        vals = reshape(double(T.(metrics{col})),[2,4]).';
        bar(x,vals,'grouped');
        grid on;
        set(gca,'XTick',x,'XTickLabel',ticknames,'XTickLabelRotation',25);
        ylabel('RMSE [m/s]'); title(titles{col});
        if col==1, legend(legends,'Location','northwest'); end
    end
    sgtitle('PG versus AC | Four phantom cases | Existing results');
    exportgraphics(f1,fullfile(outdir,'global_regional_rmse.png'),'Resolution',180);

    f2=figure('Color','w','Name','Final SWS contrast comparison', ...
        'Position',[100 100 1270 430]);
    names = {'CNR','DirectKDE_gCNR','ContrastRecovery_percent'};
    labels = {'CNR','Direct-KDE gCNR','Contrast recovery [%]'};
    for col=1:3
        subplot(1,3,col);
        vals = reshape(double(T.(names{col})),[2,4]).';
        bar(x,vals,'grouped'); grid on;
        set(gca,'XTick',x,'XTickLabel',ticknames,'XTickLabelRotation',25);
        ylabel(labels{col}); title(labels{col});
        if col==1, legend(legends,'Location','northwest'); end
    end
    sgtitle('Contrast-related SWS metrics | Four phantom cases');
    exportgraphics(f2,fullfile(outdir,'contrast_metrics.png'),'Resolution',180);

    %% Export LaTeX summary table, without external toolbox dependencies
    fid = fopen(fullfile(outdir,'final_sws_table.tex'),'w');
    assert(fid > 0,'Cannot write LaTeX output.');
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%% Automatically generated from %s\n',source_file);
    fprintf(fid,'%% Requires booktabs and graphicx packages.\n');
    fprintf(fid,'\\begin{table}[htbp]\\centering\\scriptsize\n');
    fprintf(fid,'\\caption{Comparison of PG and AC using shared phantom-specific evaluation masks. RMSE in m/s; contrast recovery in percent.}\n');
    fprintf(fid,'\\label{tab:final_sws_comparison}\n');
    fprintf(fid,'\\resizebox{\\textwidth}{!}{%%\n');
    fprintf(fid,'\\begin{tabular}{llrrrrrrrrrr}\\toprule\n');
    fprintf(fid,'Phantom & Est. & Global RMSE & Incl. RMSE & Bg. RMSE & Incl. bias & Bg. bias & CNR & gCNR & CR [\\%%] & Edge acc. [mm] & Edge sharp. [m/s/mm] \\\\ \\midrule\n');
    for i=1:height(T)
        label = strrep(char(T.Phantom(i)),'_',' ');
        fprintf(fid,'%s & %s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.2f & %s & %s \\\\ \n', ...
            label, char(T.Method(i)), ...
            T.GlobalRMSE_mps(i), T.InclusionRMSE_mps(i), T.BackgroundRMSE_mps(i), ...
            T.InclusionBias_mps(i), T.BackgroundBias_mps(i), T.CNR(i), ...
            T.DirectKDE_gCNR(i), T.ContrastRecovery_percent(i), ...
            fmt_na(T.EdgeAccuracy_mm(i), '%.3f'), ...
            fmt_na(T.EdgeSharpness_mps_per_mm(i), '%.3f'));
    end
    fprintf(fid,'\\bottomrule\\end{tabular}}\n');
    fprintf(fid,'\\end{table}\n');

    %% Data-grounded qualitative comparisons; no fabricated physical cause
    ix_pg = 1:2:8; ix_ac = 2:2:8;
    d_inc = double(T.InclusionRMSE_mps(ix_pg))-double(T.InclusionRMSE_mps(ix_ac));
    d_bg = double(T.BackgroundRMSE_mps(ix_pg))-double(T.BackgroundRMSE_mps(ix_ac));
    fnotes = fopen(fullfile(outdir,'key_findings.txt'),'w');
    assert(fnotes>0,'Cannot write findings output.');
    fprintf(fnotes,'SWS COMPARISON - INTERPRETATION OF EXISTING RESULTS\n');
    fprintf(fnotes,'Input: %s\n\n',source_file);
    fprintf(fnotes,'AC lower inclusion RMSE in %d / 4 cases.\n',sum(d_inc>0));
    fprintf(fnotes,'PG lower background RMSE in %d / 4 cases.\n',sum(d_bg<0));
    fprintf(fnotes,'Inclusion RMSE max: %.4f m/s in %s (%s).\n', ...
        max(double(T.InclusionRMSE_mps)), ...
        char(T.Phantom(find(T.InclusionRMSE_mps==max(T.InclusionRMSE_mps),1))), ...
        char(T.Method(find(T.InclusionRMSE_mps==max(T.InclusionRMSE_mps),1))));
    fprintf(fnotes,['Interpretation: No universal estimator winner follows from these cases.\n' ...
        'Distinguish reconstruction accuracy (RMSE/bias) from contrast separability (CNR/gCNR).\n' ...
        'NaN edge values denote unmeasurable edge metrics, not zero error.\n' ...
        'Causes of systematic bias require separate physical verification.\n']);
    fclose(fnotes);
    fprintf('\nSaved outputs to: %s\n',outdir);
end

function txt = fmt_na(x, form)
    if ~isfinite(x)
        txt = '--';
    else
        txt = sprintf(form,x);
    end
end
