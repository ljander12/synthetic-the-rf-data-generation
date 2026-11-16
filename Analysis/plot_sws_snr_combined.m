function [fig, combined] = plot_sws_snr_combined()
%PLOT_SWS_SNR_COMBINED  Combined SWS SNR results for two digital phantoms.
%
% Run from anywhere after adding this .m file to the MATLAB path:
%   plot_sws_snr_combined
%
% Reads FULL-STUDY results saved by run_sws_snr_comparison.m:
%   Results/snr_method_comparison/large_inclusions/full/measurements.csv
%   Results/snr_method_comparison/high_contrast/full/measurements.csv
%
% Creates ONE figure: 4 metrics (rows) x 2 phantoms (columns), PG vs AC.
% Error bars: sample SD across independent noise realizations (n = 10).
% Clean reference (SNR = Inf) is a separate marker, NOT a numeric 50 dB point.
% Does not rerun reconstruction or modify any source files.

script_dir = fileparts(mfilename('fullpath'));
repo_dir = script_dir;
for attempt = 1:8
    if isfolder(fullfile(repo_dir,'Results','snr_method_comparison'))
        break;
    end
    parent = fileparts(repo_dir);
    if strcmp(parent,repo_dir)
        error('Repository not found above %s.',script_dir);
    end
    repo_dir = parent;
end
assert(isfolder(fullfile(repo_dir,'Results','snr_method_comparison')), ...
    'Cannot locate Results/snr_method_comparison above script.');

cases = {'large_inclusions','high_contrast'};
case_titles = {'Large inclusions','High contrast'};
methods = ["Phase gradient", "Spatial autocorrelation"];
levels = [10 20 30 40];
fields = {'GlobalRMSE_mps','InclusionRMSE_mps','CNR','ContrastRecovery'};
labels = {'Global RMSE [m/s]', 'Inclusion RMSE [m/s]', ...
    'Contrast-to-noise ratio (CNR)', 'Contrast recovery [%]'};

T = cell(1,2);
combined = table();
for c = 1:2
    csv_path = fullfile(repo_dir,'Results','snr_method_comparison', ...
        cases{c},'full','measurements.csv');
    assert(isfile(csv_path),'Missing full-study CSV: %s',csv_path);
    T{c} = readtable(csv_path,'TextType','string');
    needed = [{'SNR_dB','Method'},fields];
    assert(all(ismember(needed,T{c}.Properties.VariableNames)), ...
        'Required columns missing from %s',csv_path);
    assert(nnz(isinf(T{c}.SNR_dB)) == 2, ...
        'Expected exactly two clean-reference rows: %s',csv_path);
    if isempty(combined)
        combined = T{c};
    else
        combined = [combined;T{c}]; %#ok<AGROW>
    end
end

fig = figure('Name','Combined SWS noise validation', ...
    'NumberTitle','off','Color','w', ...
    'Position',[90 70 1350 1600]);
layout = tiledlayout(fig,4,2, ...
    'Padding','compact','TileSpacing','compact');

title(layout,'SWS reconstruction sensitivity to beamformed-data noise', ...
    'FontSize',17,'FontWeight','bold');

for row = 1:4
    field = fields{row};
    factor = 1;
    if strcmp(field,'ContrastRecovery'), factor = 100; end

    for c = 1:2
        ax = nexttile(layout,(row-1)*2+c);
        hold(ax,'on');
        colors = get(ax,'ColorOrder');
        for j = 1:2
            selected_method = methods(j);
            avg = nan(1,numel(levels));
            sd = nan(1,numel(levels));
            for ii = 1:numel(levels)
                mask = T{c}.SNR_dB == levels(ii) & ...
                    string(T{c}.Method) == selected_method;
                values = double(T{c}.(field)(mask));
                assert(numel(values)==10, ...
                    ['Expected 10 realizations for %s, %s, %d dB. ' ...
                     'Found %d.'],cases{c},selected_method,levels(ii),numel(values));
                assert(all(isfinite(values)), ...
                    'Nonfinite metric values for %s, %s, %g dB.', ...
                    cases{c},selected_method,levels(ii));
                avg(ii) = mean(values)*factor;
                sd(ii) = std(values,0)*factor;  % Sample SD (N-1).
            end

            mask_clean = isinf(T{c}.SNR_dB) & ...
                string(T{c}.Method) == selected_method;
            clean = double(T{c}.(field)(mask_clean));
            assert(numel(clean)==1, ...
                'Expected one clean value for %s, %s.', cases{c},selected_method);
            clean = clean*factor;

            linestyle = '-'; marker = 'o';
            if j==2, linestyle = '--'; marker = 's'; end
            offset = -0.065 + (j-1)*0.13;
            h = errorbar(ax,(1:4)+offset,avg,sd, ...
                'LineStyle',linestyle,'Marker',marker, ...
                'LineWidth',1.5,'MarkerSize',5,'CapSize',4, ...
                'Color',colors(j,:), ...
                'DisplayName',char(selected_method));

            % Clean is shown separately, NOT treated as an SNR of 50 dB.
            plot(ax,5+offset,clean,'d','LineStyle','none', ...
                'MarkerSize',7,'MarkerFaceColor',h.Color, ...
                'MarkerEdgeColor',h.Color,'HandleVisibility','off');
        end

        xline(ax,4.5,':','HandleVisibility','off');
        xlim(ax,[0.55 5.45]);
        xticks(ax,1:5);
        xticklabels(ax,{'10','20','30','40','Clean'});
        ylabel(ax,labels{row});
        grid(ax,'on');
        ax.GridAlpha = 0.15;
        ax.FontSize = 10;
        if row==1
            title(ax,case_titles{c},'FontSize',13);
            legend(ax,'Location','best','FontSize',9);
        end
        if row==4
            xlabel(ax,'Beamformed-data SNR [dB] / clean reference');
        end
    end
end

output_dir = fullfile(repo_dir,'Results','snr_method_comparison','combined');
if ~isfolder(output_dir), mkdir(output_dir); end
png_file = fullfile(output_dir,'SNR_comparison_two_phantoms.png');
fig_file = fullfile(output_dir,'SNR_comparison_two_phantoms.fig');
combined_csv = fullfile(output_dir,'measurements_two_phantoms.csv');
exportgraphics(fig,png_file,'Resolution',300);
savefig(fig,fig_file);
writetable(combined,combined_csv);
fprintf('\nCombined SNR figure saved:\n%s\n',png_file);
fprintf('Editable MATLAB figure:\n%s\n',fig_file);
fprintf('Combined raw measurements:\n%s\n',combined_csv);
fprintf('Error bars: sample SD of 10 realizations. Clean: one reference.\n');
end
