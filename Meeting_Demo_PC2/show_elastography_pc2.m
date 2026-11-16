function show_elastography_pc2(mode)
%SHOW_ELASTOGRAPHY_PC2  Standalone presentation viewer for computer 2.
% No external toolboxes or datasets required. Uses PNG and CSV export only.
% Run without arguments for all sections, or choose:
%   show_elastography_pc2('overview')
%   show_elastography_pc2('high_contrast')
%   show_elastography_pc2('large_inclusions')
%   show_elastography_pc2('metrics')
%   show_elastography_pc2('stress')
%   show_elastography_pc2('diagnostic')

if nargin<1 || isempty(mode), mode='all'; end
mode=lower(char(string(mode)));
root=fileparts(mfilename('fullpath'));

fprintf('\nTIME-HARMONIC ELASTOGRAPHY - PORTABLE RESULTS DEMO\n');
fprintf('Folder: %s\n', root);
fprintf('Data: four simulated cases; small_inclusions not yet simulated.\n\n');

cases={'chaoran_baseline','low_contrast','high_contrast','large_inclusions'};
if any(strcmp(mode,cases))
    show_case(root,mode);
    return;
end

switch mode
    case 'all'
        show_overview(root);
        show_metrics(root);
        show_stress(root);
        show_diagnostic(root);
    case 'overview'
        show_overview(root);
    case 'metrics'
        show_metrics(root);
    case 'stress'
        show_stress(root);
    case 'diagnostic'
        show_diagnostic(root);
    otherwise
        error('Unsupported mode: %s',mode);
end
end

function show_overview(root)
show_image(fullfile(root,'all_phantoms_overview.png'), ...
    'All four digital phantoms | Ground truth / Phase gradient / Spatial autocorrelation');
end

function show_case(root,phantom)
dir_case=fullfile(root,'phantoms',phantom);
show_image(fullfile(dir_case,'method_maps.png'), ...
    [strrep(phantom,'_',' ') ' | Ground truth and both estimators']);
show_image(fullfile(dir_case,'method_errors.png'), ...
    [strrep(phantom,'_',' ') ' | SWS error maps']);
metric_path=fullfile(dir_case,'method_metrics.csv');
if isfile(metric_path)
    fprintf('\n%s: Method metrics\n', phantom);
    disp(readtable(metric_path));
end
end

function show_metrics(root)
filename=fullfile(root,'batch_summary.csv');
if ~isfile(filename)
    warning('Metrics CSV missing: %s',filename); return;
end
T=readtable(filename);
fields={'Phantom','Method','GlobalRMSE_mps','InclusionRMSE_mps', ...
    'BackgroundBias_mps','CNR','DirectKDE_gCNR','ContrastRecovery'};
if ~all(ismember(fields,T.Properties.VariableNames))
    warning('Unexpected columns in batch_summary.csv');
    disp(T); return;
end
fprintf('\nFOUR-PHANTOM QUANTITATIVE METRICS\n');
disp(T(:,fields));

phantoms=["chaoran_baseline","low_contrast","high_contrast","large_inclusions"];
methods=["Phase gradient","Spatial autocorrelation"];
metricnames={'InclusionRMSE_mps','CNR','ContrastRecovery','BackgroundBias_mps'};
metriclabels={'Inclusion RMSE [m/s] (lower better)', 'CNR (higher better)', ...
    'Contrast recovery [%] (higher better)', 'Background bias [m/s] (near zero better)'};
fig=figure('Name','THE | Metrics across four phantoms', ...
    'NumberTitle','off','Color','w','Position',[100 100 1200 850]);
L=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for m=1:4
    vals=nan(4,2);
    for i=1:4
        for j=1:2
            ix=find(string(T.Phantom)==phantoms(i) & ...
                    string(T.Method)==methods(j),1);
            if ~isempty(ix), vals(i,j)=double(T.(metricnames{m})(ix));end
        end
    end
    if m==3, vals=100*vals;end
    ax=nexttile(L);
    bar(ax,vals,'grouped');
    xticks(ax,1:4);
    xticklabels(ax,{'Baseline','Low contrast','High contrast','Large incl.'});
    title(ax,metriclabels{m});
    grid(ax,'on');
    if m==1
        legend(ax,{'Phase gradient','Spatial autocorrelation'},'Location','best');
    end
end
title(L,'Same ground truth and ROI masks | Both SWS methods');
end

function show_stress(root)
stress_dir=fullfile(root,'stress_tests');
files={'global_bias_absolute_error.png','map_noise_global_rmse.png', ...
    'blur_sharpness_metric_comparison.png','sharpness_noise_sensitivity.png'};
labels={'Global SWS bias','Map-domain Gaussian noise', ...
    'Gaussian blur','Sharpness response to noise'};
fig=figure('Name','THE | Controlled metric stress tests', ...
    'NumberTitle','off','Color','w','Position',[100 100 1400 900]);
L=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for k=1:4
    ax=nexttile(L);
    f=fullfile(stress_dir,files{k});
    if isfile(f)
        image(ax,imread(f)); axis(ax,'image'); axis(ax,'off');
    else
        axis(ax,'off');
        text(ax,.5,.5,'Figure unavailable','Units','normalized', ...
            'HorizontalAlignment','center');
    end
    title(ax,labels{k});
end
title(L,'Stress tests applied to reconstructed SWS maps (not input-SNR sweep)');
end

function show_diagnostic(root)
show_image(fullfile(root,'diagnostic','gt_vs_estimated_maps.png'), ...
    'Exploratory diagnostic | GT vs estimated displacement (alignment limitations)');
f=fullfile(root,'diagnostic','gt_vs_estimated_metrics.csv');
if isfile(f)
    fprintf('\nGT versus estimated displacement diagnostic (exploratory)\n');
    disp(readtable(f));
end
end

function show_image(filename,windowtitle)
if ~isfile(filename)
    warning('Figure not available: %s',filename); return;
end
im=imread(filename);
fig=figure('Name',windowtitle,'NumberTitle','off','Color','w', ...
    'Units','normalized','Position',[.06 .06 .88 .85]);
ax=axes('Parent',fig,'Position',[.02 .02 .96 .91]);
image(ax,im); axis(ax,'image'); axis(ax,'off');
title(ax,windowtitle,'Interpreter','none');
end
