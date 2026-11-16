function output_dir = elastography_meeting_demo(mode)
%ELASTOGRAPHY_MEETING_DEMO One-button portable MSc elastography presentation.
%
% PC 1 (computer with the repository and saved Results):
%   elastography_meeting_demo           % exports + opens ALL figures
%   elastography_meeting_demo('export') % exports only
%
% Copy the created Meeting_Demo_PC2 directory in its entirety to PC 2.
% PC 2 (with no simulation code or USTB installed):
%   Open Meeting_Demo_PC2/elastography_meeting_demo.m and press Run.
%   elastography_meeting_demo('overview')
%   elastography_meeting_demo('large_inclusions')
%   elastography_meeting_demo('metrics')
%   elastography_meeting_demo('stress')
%   elastography_meeting_demo('diagnostic')
%
% Only reads previous results; it never starts k-Wave or Field II.
% Compatible with MATLAB R2023b. Script has not been executed in MATLAB here.

if nargin < 1 || isempty(mode), mode = 'all'; end
mode = lower(char(string(mode)));
valid_modes = {'all','export','overview','metrics','stress','diagnostic', ...
    'chaoran_baseline','low_contrast','high_contrast','large_inclusions'};
assert(ismember(mode,valid_modes), 'Unknown mode: %s',mode);

script_dir = fileparts(mfilename('fullpath'));
marker = fullfile(script_dir,'EXPORT_COMPLETE.txt');

if isfile(marker)
    % On PC 2: the MATLAB file lives inside the previously exported folder.
    output_dir = script_dir;
    fprintf('\nPORTABLE MODE: reading exported results from:\n  %s\n', output_dir);
    if strcmp(mode,'export')
        fprintf('This directory is already exported. Copy the whole folder to PC 2.\n');
        return;
    end
else
    % On PC 1: locate saved Results from the location of THIS .m file.
    repo_dir = locate_project_root(script_dir);
    output_dir = fullfile(repo_dir,'Meeting_Demo_PC2');
    fprintf('\nPC 1: Found project:\n  %s\n',repo_dir);
    fprintf('Exporting previously computed data for PC 2...\n');
    export_saved_results(repo_dir,output_dir,mfilename('fullpath'));
    fprintf('\nEXPORT COMPLETE:\n  %s\n',output_dir);
    fprintf('Copy the entire Meeting_Demo_PC2 folder to PC 2.\n');
    if strcmp(mode,'export'), return; end
end

% This viewer works on PC 2 with standard MATLAB only: PNG + CSV.
switch mode
    case 'all'
        show_overview(output_dir);
        show_metrics(output_dir);
        show_case(output_dir,'large_inclusions',false);
        show_stress(output_dir);
        show_diagnostic(output_dir);
    case 'overview'
        show_overview(output_dir);
    case 'metrics'
        show_metrics(output_dir);
    case 'stress'
        show_stress(output_dir);
    case 'diagnostic'
        show_diagnostic(output_dir);
    otherwise
        show_case(output_dir,mode,true);
end
fprintf('\nDemo complete. No simulations were run.\n');
end

function repo_dir = locate_project_root(start_dir)
repo_dir = start_dir;
for k = 1:12
    check = fullfile(repo_dir,'Results','method_comparison','batch');
    if isfolder(check), return; end
    parent = fileparts(repo_dir);
    if strcmp(parent,repo_dir), break; end
    repo_dir = parent;
end
error(['Cannot locate Results/method_comparison/batch above this script.\n' ...
    'Put elastography_meeting_demo.m in the project Analysis folder or one of its subfolders.\n' ...
    'Script location: %s'],start_dir);
end

function export_saved_results(repo_dir,out_dir,source_script)
results_dir = fullfile(repo_dir,'Results');
batch_dir = fullfile(results_dir,'method_comparison','batch');
assert(isfolder(batch_dir),'Missing batch folder: %s',batch_dir);
if ~isfolder(out_dir), mkdir(out_dir); end

cases = {'chaoran_baseline','low_contrast','high_contrast','large_inclusions'};
labels = {'Baseline','Low contrast','High contrast','Large inclusions'};

% Export individual phantom figures, generating a PNG from saved maps if
% a previous PNG is missing. No inversion or simulation is performed.
num_cases = 0;
for k = 1:numel(cases)
    in_dir = fullfile(batch_dir,cases{k});
    case_out = fullfile(out_dir,'phantoms',cases{k});
    if ~isfolder(case_out), mkdir(case_out); end
    map_in = fullfile(in_dir,'method_maps.png');
    err_in = fullfile(in_dir,'method_errors.png');
    mat_in = fullfile(in_dir,'method_comparison.mat');
    if ~isfile(map_in) && isfile(mat_in)
        fprintf('Rendering missing method_maps.png from %s MAT file...\n',cases{k});
        try
            render_case_maps(mat_in,fullfile(case_out,'method_maps.png'));
        catch ME
            warning('Could not render %s: %s',cases{k},ME.message);
        end
    else
        copy_if_available(map_in,fullfile(case_out,'method_maps.png'));
    end
    copy_if_available(err_in,fullfile(case_out,'method_errors.png'));
    copy_if_available(fullfile(in_dir,'method_metrics.csv'), ...
        fullfile(case_out,'method_metrics.csv'));
    if isfile(fullfile(case_out,'method_maps.png')), num_cases=num_cases+1; end
end
assert(num_cases>0, ['No phantom figures or reconstruction MAT files were found.\n' ...
    'First run compare_SWS_methods_batch on PC 1.']);

copy_if_available(fullfile(batch_dir,'batch_summary.csv'), ...
    fullfile(out_dir,'batch_summary.csv'));

% Generate a genuine 4-row x 3-column overview with common colour scale.
% Individual pre-rendered PNG images are included regardless, as fallback.
overview_file = fullfile(out_dir,'all_phantoms_overview.png');
try
    fig=figure('Visible','off','Color','w','Position',[80 70 1700 1250]);
    L=tiledlayout(fig,4,3,'TileSpacing','compact','Padding','compact');
    captions={'Ground truth','Phase gradient','Spatial autocorrelation'};
    n_valid=0;
    for k=1:4
        mat_file=fullfile(batch_dir,cases{k},'method_comparison.mat');
        S=struct();
        if isfile(mat_file)
            S=load(mat_file,'GT','PG','AC','x_axis','z_axis');
        end
        complete=all(isfield(S,{'GT','PG','AC','x_axis','z_axis'}));
        if complete, n_valid=n_valid+1; end
        for j=1:3
            ax=nexttile(L,(k-1)*3+j);
            if complete
                A={S.GT,S.PG,S.AC};
                imagesc(ax, double(S.x_axis(:))*1e3, ...
                    double(S.z_axis(:))*1e3, double(A{j}));
                axis(ax,'image'); set(ax,'YDir','normal');
                clim(ax,[0,8]);
                if j==3, colorbar(ax); end
                if k==4, xlabel(ax,'Lateral [mm]'); end
                if j==1, ylabel(ax,'Depth [mm]'); end
            else
                axis(ax,'off');
                text(ax,.5,.5,'MAT map unavailable','Units','normalized', ...
                    'HorizontalAlignment','center');
            end
            title(ax,sprintf('%s | %s',labels{k},captions{j}), ...
                'Interpreter','none','FontSize',10);
        end
    end
    title(L,'Digital phantoms | SWS [m/s] | 200 Hz | Common colour scale');
    exportgraphics(fig,overview_file,'Resolution',170);
    close(fig);
    fprintf('Overview: %d/4 complete phantom rows.\n',n_valid);
catch ME
    if exist('fig','var') && isgraphics(fig), close(fig); end
    warning('Overview creation failed: %s. Individual phantom PNGs are still available.',ME.message);
end

% Export stress-test plots if genuinely available, never fabricate figures.
stress_in = fullfile(results_dir,'metric_stress_test');
stress_out = fullfile(out_dir,'stress_tests');
if ~isfolder(stress_out), mkdir(stress_out); end
stress_count = copy_png_csv(stress_in,stress_out);
if stress_count == 0
    warning(['No stress-test files at %s. This is optional: the demo will ' ...
        'skip its stress-test figure, not show empty panels.'],stress_in);
end

% Export exploratory k-Wave GT vs displacement-estimate diagnostic.
diag_in = fullfile(results_dir,'gt_displacement_diagnostic','large_inclusions');
diag_out = fullfile(out_dir,'diagnostic');
if ~isfolder(diag_out), mkdir(diag_out); end
copy_png_csv(diag_in,diag_out);

% Copy this same file to the portable output directory. Do this before the
% completion marker so 'Run' on PC 2 always means viewing existing data.
copyfile(source_script,fullfile(out_dir,'elastography_meeting_demo.m'));
marker=fullfile(out_dir,'EXPORT_COMPLETE.txt');
fd=fopen(marker,'w');
assert(fd>=0,'Could not create export completion marker: %s',marker);
fprintf(fd,'This directory contains exported results. Run elastography_meeting_demo in MATLAB.\n');
fclose(fd);

fprintf('Exported %d of 4 phantom comparison image sets.\n',num_cases);
fprintf('Exported %d stress-test PNG/CSV files.\n',stress_count);
end

function count=copy_png_csv(src,dst)
count=0;
if ~isfolder(src), return; end
patterns={'*.png','*.csv'};
for j=1:numel(patterns)
    listing=dir(fullfile(src,patterns{j}));
    for k=1:numel(listing)
        copyfile(fullfile(src,listing(k).name),fullfile(dst,listing(k).name));
        count=count+1;
    end
end
end

function ok=copy_if_available(src,dst)
ok=false;
if isfile(src)
    copyfile(src,dst); ok=true;
else
    fprintf('Optional file not found: %s\n',src);
end
end

function render_case_maps(mat_file,png_file)
S=load(mat_file,'GT','PG','AC','x_axis','z_axis');
assert(all(isfield(S,{'GT','PG','AC','x_axis','z_axis'})), ...
    'MAT missing the saved GT/PG/AC arrays.');
fig=figure('Visible','off','Color','w','Position',[100 100 1450 450]);
L=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
labels={'Ground truth','Phase gradient','Spatial autocorrelation'};
A={S.GT,S.PG,S.AC};
for j=1:3
    ax=nexttile(L);
    imagesc(ax,double(S.x_axis(:))*1e3,double(S.z_axis(:))*1e3,double(A{j}));
    axis(ax,'image');set(ax,'YDir','normal');clim(ax,[0 8]);colorbar(ax);
    title(ax,labels{j});xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
end
exportgraphics(fig,png_file,'Resolution',180);close(fig);
end

function show_overview(root)
png=fullfile(root,'all_phantoms_overview.png');
if isfile(png)
    show_image(png,'ALL PHANTOMS | Ground truth / Phase gradient / Spatial autocorrelation');
else
    fprintf('\nNo combined overview figure; showing saved phantom figures individually.\n');
    cases={'chaoran_baseline','low_contrast','high_contrast','large_inclusions'};
    for k=1:numel(cases), show_case(root,cases{k},false); end
end
end

function show_case(root,phantom,with_errors)
folder=fullfile(root,'phantoms',phantom);
show_image(fullfile(folder,'method_maps.png'), ...
    ['SWS | ' strrep(phantom,'_',' ')]);
if with_errors
    show_image(fullfile(folder,'method_errors.png'), ...
        ['SWS errors | ' strrep(phantom,'_',' ')]);
end
csv=fullfile(folder,'method_metrics.csv');
if with_errors && isfile(csv)
    fprintf('\n%s - method metrics:\n',phantom);
    disp(readtable(csv));
end
end

function show_metrics(root)
file=fullfile(root,'batch_summary.csv');
if ~isfile(file)
    fprintf('\n[SKIPPED] Metrics: batch_summary.csv was not exported.\n');
    return;
end
T=readtable(file);
fprintf('\n=== METRICS: FOUR PHANTOMS AND TWO ESTIMATORS ===\n');
disp(T);
columns={'Phantom','Method','InclusionRMSE_mps','CNR','ContrastRecovery','BackgroundBias_mps'};
if ~all(ismember(columns,T.Properties.VariableNames))
    fprintf('Metrics table exists, but its columns do not match this plot layout.\n');
    return;
end
cases={'chaoran_baseline','low_contrast','high_contrast','large_inclusions'};
methods={'Phase gradient','Spatial autocorrelation'};
metrics={'InclusionRMSE_mps','CNR','ContrastRecovery','BackgroundBias_mps'};
names={'Inclusion RMSE [m/s] (lower is better)','CNR (higher is better)', ...
    'Contrast recovery [%] (higher is better)', ...
    'Background bias [m/s] (closer to zero is better)'};
fig=figure('Name','ELASTOGRAPHY | Quantitative method comparison', ...
    'NumberTitle','off','Color','w','Position',[100 100 1200 850]);
L=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for m=1:4
    v=nan(4,2);
    for i=1:4
        for j=1:2
            ix=find(strcmp(string(T.Phantom),cases{i}) & ...
                strcmp(string(T.Method),methods{j}),1);
            if ~isempty(ix), v(i,j)=double(T.(metrics{m})(ix)); end
        end
    end
    if m==3, v=100*v; end
    ax=nexttile(L);bar(ax,v,'grouped');grid(ax,'on');
    xticks(ax,1:4);xticklabels(ax,{'Baseline','Low contrast','High contrast','Large incl.'});
    title(ax,names{m});
    if m==1
        legend(ax,methods,'Location','best');
    end
end
title(L,'Same SWS ground truth and common evaluation masks');
end

function show_stress(root)
folder=fullfile(root,'stress_tests');
files={'global_bias_absolute_error.png','map_noise_global_rmse.png', ...
    'blur_sharpness_metric_comparison.png','sharpness_noise_sensitivity.png'};
titles={'Global SWS bias','Gaussian noise added to SWS map', ...
    'Gaussian blur','Sharpness sensitivity to noise'};
keep=find(cellfun(@(f)isfile(fullfile(folder,f)),files));
if isempty(keep)
    fprintf(['\n[SKIPPED] Stress-test figures were not exported.\n' ...
        'Check Results/metric_stress_test on PC 1.\n']);
    return;
end
fig=figure('Name','ELASTOGRAPHY | Controlled metric stress tests', ...
    'NumberTitle','off','Color','w','Position',[100 100 1300 800]);
n=numel(keep); ncols=min(2,n); nrows=ceil(n/ncols);
L=tiledlayout(fig,nrows,ncols,'TileSpacing','compact','Padding','compact');
for k=1:n
    idx=keep(k);
    ax=nexttile(L); image(ax,imread(fullfile(folder,files{idx})));
    axis(ax,'image');axis(ax,'off');title(ax,titles{idx});
end
title(L,'Controlled perturbations of reconstructed SWS maps (not input SNR)');
fprintf('Stress plots: %d of 4 available.\n',n);
end

function show_diagnostic(root)
file=fullfile(root,'diagnostic','gt_vs_estimated_maps.png');
if ~isfile(file)
    fprintf('\n[SKIPPED] GT-displacement diagnostic PNG was not exported.\n');
    return;
end
show_image(file,'Exploratory: GT displacement vs estimated displacement SWS');
csv=fullfile(root,'diagnostic','gt_vs_estimated_metrics.csv');
if isfile(csv)
    fprintf('\nExploratory GT displacement diagnostic:\n');disp(readtable(csv));
end
end

function show_image(path,name)
if ~isfile(path)
    fprintf('[SKIPPED] Figure not found: %s\n',path);
    return;
end
fig=figure('Name',name,'NumberTitle','off','Color','w', ...
    'Units','normalized','Position',[.06 .06 .88 .84]);
ax=axes('Parent',fig,'Position',[.02 .03 .96 .90]);
image(ax,imread(path)); axis(ax,'image');axis(ax,'off');
title(ax,name,'Interpreter','none');
end
