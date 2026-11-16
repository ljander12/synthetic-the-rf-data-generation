function report = diagnose_gt_geometry(phantom_case)
%DIAGNOSE_GT_GEOMETRY Check whether stored SWS ground truth uses the
% physical imaging plane specified by SIM_PHANTOM_MAIN/get_GT_kwave.
%
% Usage:
%   diagnose_gt_geometry('large_inclusions')
%   diagnose_gt_geometry('high_contrast')
%
% This diagnostic NEVER overwrites original data or SWS results.
% It requires k-Wave, create_phantom.m, and coor_transformation.m.

if nargin < 1 || isempty(phantom_case)
    phantom_case = 'chaoran_baseline';
end
phantom_case = char(string(phantom_case));

% Find project root independent of location of this .m file.
repo_dir = fileparts(mfilename('fullpath'));
found = false;
for k = 1:8
    if isfolder(fullfile(repo_dir,'Data')) && ...
            isfolder(fullfile(repo_dir,'Analysis'))
        found = true;
        break;
    end
    parent_dir = fileparts(repo_dir);
    if strcmp(parent_dir, repo_dir), break; end
    repo_dir = parent_dir;
end
if ~found
    error(['Could not find project root containing Data and Analysis. ' ...
        'Place diagnose_gt_geometry.m in the project Analysis folder.']);
end

% Add only project helper folders; external k-Wave must already be installed.
for d = {'kwave_fun','utils','FIELD_fun'}
    p = fullfile(repo_dir,d{1});
    if isfolder(p), addpath(p); end
end

file_name = [phantom_case '_displacement_with_quality.mat'];
candidates = { ...
    fullfile(repo_dir,'Data',phantom_case,file_name), ...
    fullfile(repo_dir,'Data',file_name), ...
    fullfile(repo_dir,'Analysis',file_name), ...
    fullfile(repo_dir,file_name)};
idx = find(cellfun(@isfile,candidates),1,'first');
if isempty(idx)
    fprintf('Dataset locations checked:\n');
    fprintf('  %s\n',candidates{:});
    error('Dataset not found for %s.',phantom_case);
end
data_file = candidates{idx};
fprintf('\n=== GROUND TRUTH GEOMETRY DIAGNOSTIC ===\n');
fprintf('Phantom: %s\nData: %s\n',phantom_case,data_file);

% Load ONLY numerical ground-truth data and metadata. Loading b_data
% requires USTB's uff.beamformed_data class; it is NOT needed for this test.
% Query names first so old datasets with missing geometry do not emit noise.
var_names = who('-file',data_file);
wanted = {'sws_gt','sws_gt_info','shear_params','prb_center','prb_theta'};
selected = wanted(ismember(wanted,var_names));
S = load(data_file,selected{:});
if ~isfield(S,'sws_gt') || ~isfield(S,'shear_params')
    error('Dataset lacks sws_gt or shear_params. No data was changed.');
end

% Retrieve the actual ultrasound coordinate axes from saved metadata or
% SWS-analysis results, without requiring any USTB objects.
axes_source = '';
if isfield(S,'sws_gt_info') && isstruct(S.sws_gt_info) && ...
        isfield(S.sws_gt_info,'x_axis') && ...
        isfield(S.sws_gt_info,'z_axis')
    x_axis = double(S.sws_gt_info.x_axis(:));
    z_axis = double(S.sws_gt_info.z_axis(:));
    axes_source = 'sws_gt_info';
else
    res_name = [phantom_case '_SWS_quality_results.mat'];
    result_candidates = { ...
        fullfile(repo_dir,'Results',phantom_case,res_name), ...
        fullfile(repo_dir,'Results',res_name)};
    res_idx = find(cellfun(@isfile,result_candidates),1,'first');
    if ~isempty(res_idx)
        R = load(result_candidates{res_idx},'results');
        if isfield(R,'results') && isfield(R.results,'x_axis') && ...
                isfield(R.results,'z_axis')
            x_axis = double(R.results.x_axis(:));
            z_axis = double(R.results.z_axis(:));
            axes_source = 'saved SWS quality results';
        end
    end
    if isempty(axes_source)
        batch_file = fullfile(repo_dir,'Results','method_comparison', ...
            'batch',phantom_case,'method_comparison.mat');
        if isfile(batch_file)
            B = load(batch_file,'x_axis','z_axis');
            if isfield(B,'x_axis') && isfield(B,'z_axis')
                x_axis = double(B.x_axis(:));
                z_axis = double(B.z_axis(:));
                axes_source = 'saved method-comparison results';
            end
        end
    end
end
if isempty(axes_source)
    error(['Could not find ultrasound x_axis/z_axis in sws_gt_info, ' ...
        'SWS quality results or batch method results. ' ...
        'Do not load b_data without USTB.']);
end
fprintf('Scan axes loaded from %s.\n',axes_source);
[X,Z] = meshgrid(x_axis,z_axis);
N = numel(X);
gt_stored = double(S.sws_gt);
assert(isequal(size(gt_stored),size(X)), ...
    'Stored GT size %s differs from scan grid size %s.', ...
    mat2str(size(gt_stored)),mat2str(size(X)));

[kgrid, medium] = create_phantom(S.shear_params);
if isfield(S,'prb_center') && isfield(S,'prb_theta')
    prb_center = double(S.prb_center);
    prb_theta = double(S.prb_theta);
elseif strcmpi(phantom_case,'large_inclusions')
    % Same documented legacy fallback as analyze_SWS_quality.m.
    warning(['Probe geometry missing: using documented legacy geometry ' ...
        'for large_inclusions. Verify against original simulation.']);
    prb_center = [0 -20e-3 double(kgrid.z_vec(1))];
    prb_theta = [deg2rad(90) deg2rad(0)];
else
    error(['Probe geometry not stored in data. ' ...
        'Cannot independently validate this phantom.']);
end

% The GENERATING SIMULATION uses TF (not TF_rev), and maps [x;0;z;1].
[TF,~] = coor_transformation(prb_center,prb_theta);
F = griddedInterpolant( ...
    {double(kgrid.x_vec(:)),double(kgrid.y_vec(:)),double(kgrid.z_vec(:))}, ...
    double(medium.sound_speed_shear),'nearest','none');

% Reproduce all three hypothetical plane constructions from the analysis.
plane_maps = cell(1,3);
local_points = { ...
    [X(:).';zeros(1,N);Z(:).';ones(1,N)], ...
    [X(:).';Z(:).';zeros(1,N);ones(1,N)], ...
    [zeros(1,N);X(:).';Z(:).';ones(1,N)]};

for mode = 1:3
    P = TF * local_points{mode};
    v = F(P(1,:).',P(2,:).',P(3,:).');
    plane_maps{mode} = reshape(v,size(X));
end
% Mode 1 is defined by get_GT_kwave, regardless of inclusion counts.
gt_physical = plane_maps{1};

valid = isfinite(gt_stored) & isfinite(gt_physical);
diff_map = gt_stored - gt_physical;
abs_diff = abs(diff_map(valid));
if isempty(abs_diff)
    error('No spatial overlap between stored and physically defined GT.');
end
n_changed = nnz(abs_diff > 1e-4);
n_valid = nnz(valid);
c_bkg = double(S.shear_params.c_shear_bkg);
c_incl = double(S.shear_params.c_shear_incl);
threshold = (c_bkg + c_incl)/2;
if c_incl >= c_bkg
    stored_inc = gt_stored > threshold;
    physical_inc = gt_physical > threshold;
else
    stored_inc = gt_stored < threshold;
    physical_inc = gt_physical < threshold;
end
n_class_changed = nnz(valid & xor(stored_inc,physical_inc));

fprintf('Probe orientation: [%.1f %.1f] degrees\n',rad2deg(prb_theta));
fprintf('Imaging plane: [x, 0, z] transformed with TF (mode 1)\n');
if isfield(S,'sws_gt_info') && isstruct(S.sws_gt_info) && ...
        isfield(S.sws_gt_info,'plane_mode')
    fprintf('Stored GT metadata plane_mode: %d\n', ...
        S.sws_gt_info.plane_mode);
end
for mode = 1:3
    m = plane_maps{mode};
    incl = isfinite(m) & ((c_incl >= c_bkg & m > threshold) | ...
                            (c_incl < c_bkg & m < threshold));
    fprintf('Mode %d: finite %5.1f%%, inclusion pixels %d\n', ...
        mode,100*nnz(isfinite(m))/numel(m),nnz(incl));
end
fprintf('\nStored GT vs physically defined GT:\n');
fprintf('Compared pixels: %d\n',n_valid);
fprintf('Differing SWS pixels: %d (%.4f%%)\n', ...
    n_changed,100*n_changed/n_valid);
fprintf('Differing inclusion labels: %d (%.4f%%)\n', ...
    n_class_changed,100*n_class_changed/n_valid);
fprintf('Max absolute SWS difference: %.6f m/s\n',max(abs_diff));
fprintf('RMSE of GT map difference: %.6f m/s\n', ...
    sqrt(mean(abs_diff.^2)));
if n_changed == 0 && nnz(xor(isfinite(gt_stored),isfinite(gt_physical))) == 0
    fprintf('PASS: Stored GT agrees with the simulation-defined image plane.\n');
else
    fprintf(['ATTENTION: Stored GT and simulation-defined plane differ.\n' ...
        'Do not change original files yet: inspect the figure.\n']);
end

report = struct('phantom_case',phantom_case,'data_file',data_file, ...
    'n_compared',n_valid,'n_changed',n_changed, ...
    'n_inclusion_labels_changed',n_class_changed, ...
    'gt_map_difference_rmse_mps',sqrt(mean(abs_diff.^2)), ...
    'max_abs_gt_difference_mps',max(abs_diff), ...
    'prb_center',prb_center,'prb_theta',prb_theta);

fig = figure('Color','w','Name',['GT geometry diagnostic: ' phantom_case], ...
    'Position',[100 90 1600 850]);
t = tiledlayout(fig,2,3,'Padding','compact','TileSpacing','compact');
full_range = [min(c_bkg,c_incl),max(c_bkg,c_incl)];
if full_range(1) == full_range(2)
    full_range = full_range + [-1 1]*0.1;
end
for k=1:3
    ax=nexttile(t,k);
    imagesc(ax,x_axis*1e3,z_axis*1e3,plane_maps{k});
    axis(ax,'image');set(ax,'YDir','normal');
    clim(ax,full_range);colorbar(ax);
    title(ax,sprintf('Candidate mode %d%s',k, ...
        ternary(k==1,' (simulation)','')));
    xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
end
ax=nexttile(t,4);
imagesc(ax,x_axis*1e3,z_axis*1e3,gt_stored);
axis(ax,'image');set(ax,'YDir','normal');clim(ax,full_range);colorbar(ax);
title(ax,'Stored GT');xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
ax=nexttile(t,5);
imagesc(ax,x_axis*1e3,z_axis*1e3,diff_map);
axis(ax,'image');set(ax,'YDir','normal');colorbar(ax);
err_lim=max(abs(diff_map(valid)));
if err_lim==0,err_lim=1;end
clim(ax,[-err_lim err_lim]);
title(ax,'Stored - simulation-defined GT [m/s]');
xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
ax=nexttile(t,6);
imagesc(ax,x_axis*1e3,z_axis*1e3,double(xor(stored_inc,physical_inc) & valid));
axis(ax,'image');set(ax,'YDir','normal');clim(ax,[0 1]);colorbar(ax);
title(ax,'Different inclusion classifications');
xlabel(ax,'Lateral [mm]');ylabel(ax,'Depth [mm]');
title(t,['GT geometry validation | ' strrep(phantom_case,'_',' ')]);

output_dir = fullfile(repo_dir,'Results','gt_geometry_diagnostic',phantom_case);
if ~isfolder(output_dir),mkdir(output_dir);end
exportgraphics(fig,fullfile(output_dir,'gt_geometry_diagnostic.png'),'Resolution',160);
save(fullfile(output_dir,'gt_geometry_diagnostic.mat'), ...
    'report','gt_stored','gt_physical','plane_maps','x_axis','z_axis');
fprintf('\nDiagnostic saved to: %s\n',output_dir);
fprintf('Original dataset and analysis outputs were NOT modified.\n\n');
end

function v = ternary(cond,a,b)
if cond,v=a;else,v=b;end
end
