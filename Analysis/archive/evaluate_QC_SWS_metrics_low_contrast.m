clear;
clc;
close all;

%% ========================================================================
% GT-BASED SWS EVALUATION - LOW CONTRAST
%
% Compare:
%
%       ORIGINAL SWS
%           vs.
%       Q_C-SCREENED SWS
%
% Ground truth is reconstructed directly from the SAME digital phantom
% definition used by k-Wave:
%
%       shear_params
%           ->
%       create_phantom()
%           ->
%       medium.sound_speed_shear
%           ->
%       sample into ultrasound imaging plane
%
% This avoids assumptions about cx/cy/cz coordinate interpretation.
%
% NO k-Wave propagation simulation is performed.
% NO Field II simulation is performed.
%
% Metrics:
%
%   Inclusion:
%       mean
%       median
%       bias
%       MAE
%       RMSE
%       SD
%
%   Background:
%       mean
%       median
%       bias
%       MAE
%       RMSE
%       SD
%
%   Image:
%       CNR
%       global bias
%       global MAE
%       global RMSE
%
% INPUT:
%
%   low_contrast_displacement_with_quality.mat
%   low_contrast_QC_screened_SWS.mat
%   low_contrast_ground_truth_theta90.mat
%
% OUTPUT:
%
%   low_contrast_QC_SWS_metrics.mat
%   low_contrast_QC_SWS_metrics.csv
%
% ========================================================================


%% ========================================================================
% 1. Paths
% ========================================================================

repo_dir = ...
    'C:\Users\louis\Documents\GitHub\MasterOppgave\synthetic-the-rf-data-generation';

toolbox_dir = ...
    'C:\Users\louis\Documents\GitHub\MasterOppgave\elastography\ElastographyToolbox';


addpath(genpath(repo_dir));
addpath(genpath(toolbox_dir));


%% ------------------------------------------------------------------------
% Add k-Wave if stored in HOME
% -------------------------------------------------------------------------

kwave_home = ...
    fullfile( ...
        getenv('HOME'), ...
        'k-Wave');


if isfolder(kwave_home)

    addpath(genpath(kwave_home));

end


rehash toolboxcache;


fprintf('\n');
fprintf('Checking required functions/classes...\n');

fprintf('uff.beamformed_data: %s\n', ...
    which('uff.beamformed_data'));

fprintf('create_phantom:      %s\n', ...
    which('create_phantom'));

fprintf('coor_transformation: %s\n', ...
    which('coor_transformation'));

fprintf('kWaveGrid:           %s\n', ...
    which('kWaveGrid'));


if isempty(which('uff.beamformed_data'))

    error('UFF toolbox was not found.');

end


if isempty(which('create_phantom'))

    error('create_phantom.m was not found.');

end


if isempty(which('coor_transformation'))

    error('coor_transformation.m was not found.');

end


if isempty(which('kWaveGrid'))

    error( ...
        ['k-Wave was not found on the MATLAB path. ' ...
         'Add the k-Wave folder before running this script.']);

end


%% ========================================================================
% 2. Input files
% ========================================================================

DATA_FILE = ...
    'low_contrast_displacement_with_quality.mat';

SWS_FILE = ...
    'low_contrast_QC_screened_SWS.mat';

GT_FILE = ...
    'low_contrast_ground_truth_theta90.mat';


if ~isfile(DATA_FILE)

    error('Missing file: %s',DATA_FILE);

end


if ~isfile(SWS_FILE)

    error('Missing file: %s',SWS_FILE);

end


if ~isfile(GT_FILE)

    error('Missing file: %s',GT_FILE);

end


%% ========================================================================
% 3. Load only required variables
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('GT-BASED SWS EVALUATION - LOW CONTRAST\n');
fprintf('============================================================\n');


data = ...
    load( ...
        DATA_FILE, ...
        'b_data', ...
        'shear_params', ...
        'prb_center', ...
        'prb_theta');


sws_data = ...
    load( ...
        SWS_FILE, ...
        'sws_original', ...
        'sws_screened');


gt_data = ...
    load( ...
        GT_FILE, ...
        'prb_center', ...
        'prb_theta', ...
        'TF');


%% ========================================================================
% 4. Extract variables
% ========================================================================

b_data = ...
    data.b_data;


shear_params = ...
    data.shear_params;


sws_original = ...
    double( ...
        sws_data.sws_original);


sws_screened = ...
    double( ...
        sws_data.sws_screened);


fprintf('\nLoaded b_data class: %s\n', ...
    class(b_data));


if ~isa(b_data,'uff.beamformed_data')

    error( ...
        'b_data did not load as uff.beamformed_data.');

end


%% ========================================================================
% 5. Ground-truth material values
% ========================================================================

c_bkg = ...
    double( ...
        shear_params.c_shear_bkg);


c_incl = ...
    double( ...
        shear_params.c_shear_incl);


fprintf('\n');

fprintf('Ground-truth background SWS: %.3f m/s\n', ...
    c_bkg);

fprintf('Ground-truth inclusion SWS:  %.3f m/s\n', ...
    c_incl);


%% ========================================================================
% 6. Ultrasound scan coordinates
% ========================================================================

x_scan = ...
    double( ...
        b_data.scan.x_axis(:).');


z_scan = ...
    double( ...
        b_data.scan.z_axis(:));


[X_scan, Z_scan] = ...
    meshgrid( ...
        x_scan, ...
        z_scan);


Nz = ...
    size(X_scan,1);

Nx = ...
    size(X_scan,2);


fprintf('\n');

fprintf('Scan dimensions: %d x %d\n', ...
    Nz, ...
    Nx);

fprintf('Lateral range: %.2f to %.2f mm\n', ...
    min(x_scan)*1e3, ...
    max(x_scan)*1e3);

fprintf('Depth range:   %.2f to %.2f mm\n', ...
    min(z_scan)*1e3, ...
    max(z_scan)*1e3);


%% ========================================================================
% 7. Verify SWS dimensions
% ========================================================================

if ~isequal( ...
        size(sws_original), ...
        [Nz Nx])

    error( ...
        'Original SWS dimensions do not match scan.');

end


if ~isequal( ...
        size(sws_screened), ...
        [Nz Nx])

    error( ...
        'Screened SWS dimensions do not match scan.');

end


%% ========================================================================
% 8. Recreate digital phantom geometry
% ========================================================================
%
% IMPORTANT:
%
% This only builds:
%
%       kgrid
%       medium
%
% It does NOT call pstdElastic3D.
%
% Therefore this is NOT a new shear-wave simulation.
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('RECREATING DIGITAL PHANTOM GEOMETRY\n');
fprintf('============================================================\n');


[kgrid, medium] = ...
    create_phantom( ...
        shear_params);


fprintf('Phantom created.\n');

fprintf('Medium dimensions: %s\n', ...
    mat2str( ...
        size( ...
            medium.sound_speed_shear)));


fprintf('Medium SWS range: %.3f to %.3f m/s\n', ...
    min( ...
        medium.sound_speed_shear(:)), ...
    max( ...
        medium.sound_speed_shear(:)));


%% ========================================================================
% 9. Probe geometry
% ========================================================================

if isfield(gt_data,'prb_center')

    prb_center = ...
        double( ...
            gt_data.prb_center);

elseif isfield(data,'prb_center')

    prb_center = ...
        double( ...
            data.prb_center);

else

    error('prb_center not found.');

end


if isfield(gt_data,'prb_theta')

    prb_theta = ...
        double( ...
            gt_data.prb_theta);

elseif isfield(data,'prb_theta')

    prb_theta = ...
        double( ...
            data.prb_theta);

else

    error('prb_theta not found.');

end


fprintf('\nProbe center [mm]:\n');

disp( ...
    prb_center * 1e3);


fprintf('Probe theta [deg]:\n');

disp( ...
    rad2deg(prb_theta));


%% ========================================================================
% 10. Obtain exact transformation matrices
% ========================================================================

[TF, TF_rev] = ...
    coor_transformation( ...
        prb_center, ...
        prb_theta);


%% ========================================================================
% 11. Determine transformation convention automatically
% ========================================================================
%
% Different helper functions sometimes store homogeneous transformations
% using row-vector or column-vector convention.
%
% Instead of assuming one, test:
%
%   TF * p
%   p * TF
%   TF_rev * p
%   p * TF_rev
%
% The correct local -> world transformation must map:
%
%       local origin [0 0 0]
%
% to:
%
%       prb_center
%
% ========================================================================

origin_col = ...
    [0;0;0;1];


origin_row = ...
    [0 0 0 1];


candidate_origin = ...
    nan(4,3);


%% Candidate 1
tmp = ...
    TF * origin_col;

candidate_origin(1,:) = ...
    tmp(1:3).';


%% Candidate 2
tmp = ...
    origin_row * TF;

candidate_origin(2,:) = ...
    tmp(1:3);


%% Candidate 3
tmp = ...
    TF_rev * origin_col;

candidate_origin(3,:) = ...
    tmp(1:3).';


%% Candidate 4
tmp = ...
    origin_row * TF_rev;

candidate_origin(4,:) = ...
    tmp(1:3);


origin_error = ...
    vecnorm( ...
        candidate_origin - ...
        prb_center(:).', ...
        2, ...
        2);


[best_origin_error, transform_mode] = ...
    min(origin_error);


fprintf('\n');
fprintf('Transformation convention test:\n');

for q = 1:4

    fprintf( ...
        'Mode %d origin error: %.6f mm\n', ...
        q, ...
        origin_error(q)*1e3);

end


fprintf( ...
    'Selected transform mode: %d\n', ...
    transform_mode);


fprintf( ...
    'Origin mapping error: %.6f mm\n', ...
    best_origin_error*1e3);


if best_origin_error > 1e-3

    warning( ...
        ['No transformation convention mapped the local origin ' ...
         'closely to prb_center.']);

end


%% ========================================================================
% 12. k-Wave coordinate vectors
% ========================================================================

kx = ...
    double( ...
        kgrid.x_vec(:));


ky = ...
    double( ...
        kgrid.y_vec(:));


kz = ...
    double( ...
        kgrid.z_vec(:));


medium_sws = ...
    double( ...
        medium.sound_speed_shear);


%% ========================================================================
% 13. Create interpolant of exact digital phantom
% ========================================================================

F_sws = ...
    griddedInterpolant( ...
        {kx,ky,kz}, ...
        medium_sws, ...
        'nearest', ...
        'none');


%% ========================================================================
% 14. Determine probe-plane coordinate convention
% ========================================================================
%
% Expected Field-II convention is:
%
%       local x = lateral
%       local y = elevation
%       local z = depth
%
% i.e.
%
%       [X_scan, 0, Z_scan]
%
% However, we test the three possible placements of the zero/elevation
% coordinate and select the one that:
%
%   1. lies predominantly inside the k-Wave domain
%   2. intersects the known inclusion material
%
% ========================================================================

Npts = ...
    numel(X_scan);


local_candidates = ...
    cell(3,1);


%% Expected convention
local_candidates{1} = ...
    [ ...
    X_scan(:).'; ...
    zeros(1,Npts); ...
    Z_scan(:).'; ...
    ones(1,Npts)];


%% Alternative convention
local_candidates{2} = ...
    [ ...
    X_scan(:).'; ...
    Z_scan(:).'; ...
    zeros(1,Npts); ...
    ones(1,Npts)];


%% Alternative convention
local_candidates{3} = ...
    [ ...
    zeros(1,Npts); ...
    X_scan(:).'; ...
    Z_scan(:).'; ...
    ones(1,Npts)];


candidate_GT = ...
    cell(3,1);


candidate_finite = ...
    zeros(3,1);


candidate_inclusion = ...
    zeros(3,1);


fprintf('\n');
fprintf('============================================================\n');
fprintf('TESTING IMAGING-PLANE MAPPING\n');
fprintf('============================================================\n');


for plane_mode = 1:3

    local_h = ...
        local_candidates{plane_mode};


    world_h = ...
        apply_transform( ...
            local_h, ...
            TF, ...
            TF_rev, ...
            transform_mode);


    wx = ...
        world_h(1,:);


    wy = ...
        world_h(2,:);


    wz = ...
        world_h(3,:);


    gt_values = ...
        F_sws( ...
            wx(:), ...
            wy(:), ...
            wz(:));


    gt_map = ...
        reshape( ...
            gt_values, ...
            Nz, ...
            Nx);


    candidate_GT{plane_mode} = ...
        gt_map;


    valid_gt = ...
        isfinite(gt_map);


    inclusion_gt = ...
        valid_gt & ...
        abs(gt_map-c_incl) < 1e-6;


    candidate_finite(plane_mode) = ...
        mean(valid_gt(:));


    candidate_inclusion(plane_mode) = ...
        nnz(inclusion_gt);


    fprintf( ...
        ['Plane mode %d: finite = %.2f %% | ' ...
         'inclusion pixels = %d\n'], ...
        plane_mode, ...
        100*candidate_finite(plane_mode), ...
        candidate_inclusion(plane_mode));

end


%% ========================================================================
% 15. Select best imaging-plane mapping
% ========================================================================
%
% Prefer a mapping with:
%
%       high domain coverage
%       AND
%       actual inclusion pixels
%
% ========================================================================

valid_plane_candidate = ...
    candidate_finite > 0.80;


score = ...
    candidate_inclusion;


score(~valid_plane_candidate) = ...
    -1;


[best_score, plane_mode] = ...
    max(score);


if best_score <= 0

    fprintf('\n');
    fprintf('shear_params values:\n');
    disp(shear_params);

    error( ...
        ['Unable to find the inclusion in any plausible probe-plane ' ...
         'mapping. Send this output before continuing.']);

end


sws_gt = ...
    candidate_GT{plane_mode};


fprintf('\nSelected plane mode: %d\n', ...
    plane_mode);


fprintf( ...
    'Finite GT coverage: %.2f %%\n', ...
    100*candidate_finite(plane_mode));


fprintf( ...
    'GT inclusion pixels: %d\n', ...
    candidate_inclusion(plane_mode));


%% ========================================================================
% 16. Ground-truth masks
% ========================================================================

TOL = ...
    1e-6;


mask_inclusion_full = ...
    isfinite(sws_gt) & ...
    abs(sws_gt-c_incl) < TOL;


mask_background_full = ...
    isfinite(sws_gt) & ...
    abs(sws_gt-c_bkg) < TOL;


fprintf('\n');
fprintf('============================================================\n');
fprintf('GROUND-TRUTH MATERIAL MAP\n');
fprintf('============================================================\n');


fprintf('Inclusion pixels: %d\n', ...
    nnz(mask_inclusion_full));


fprintf('Background pixels: %d\n', ...
    nnz(mask_background_full));


fprintf('GT map minimum: %.3f m/s\n', ...
    min(sws_gt(:),[],'omitnan'));


fprintf('GT map maximum: %.3f m/s\n', ...
    max(sws_gt(:),[],'omitnan'));


if nnz(mask_inclusion_full) == 0

    error('Ground-truth map contains no inclusion.');

end


%% ========================================================================
% 17. Build conservative inclusion/background ROIs
% ========================================================================
%
% Remove approximately 2 mm around the material transition.
%
% This is done using convolution rather than Image Processing Toolbox.
%
% Inclusion core:
%   local neighbourhood must contain ONLY inclusion pixels.
%
% Background:
%   local neighbourhood must contain NO inclusion pixels.
%
% ========================================================================

ROI_MARGIN = ...
    2e-3;


dx = ...
    median( ...
        diff(x_scan));


dz = ...
    median( ...
        diff(z_scan));


half_x = ...
    max( ...
        1, ...
        ceil(ROI_MARGIN/dx));


half_z = ...
    max( ...
        1, ...
        ceil(ROI_MARGIN/dz));


kernel = ...
    ones( ...
        2*half_z+1, ...
        2*half_x+1);


kernel_N = ...
    numel(kernel);


%% ------------------------------------------------------------------------
% Inclusion core
% -------------------------------------------------------------------------

local_inclusion_count = ...
    conv2( ...
        double(mask_inclusion_full), ...
        kernel, ...
        'same');


mask_inclusion_core = ...
    local_inclusion_count == kernel_N;


%% ------------------------------------------------------------------------
% Background outside transition zone
% -------------------------------------------------------------------------

mask_background_core = ...
    mask_background_full & ...
    local_inclusion_count == 0;


%% ========================================================================
% 18. Remove image boundary
% ========================================================================

IMAGE_MARGIN = ...
    3e-3;


mask_image_interior = ...
    X_scan >= min(x_scan)+IMAGE_MARGIN & ...
    X_scan <= max(x_scan)-IMAGE_MARGIN & ...
    Z_scan >= min(z_scan)+IMAGE_MARGIN & ...
    Z_scan <= max(z_scan)-IMAGE_MARGIN;


mask_inclusion_core = ...
    mask_inclusion_core & ...
    mask_image_interior;


mask_background_core = ...
    mask_background_core & ...
    mask_image_interior;


%% ========================================================================
% 19. Common valid SWS pixels
% ========================================================================

valid_original = ...
    isfinite(sws_original) & ...
    sws_original > 0 & ...
    sws_original < 10;


valid_screened = ...
    isfinite(sws_screened) & ...
    sws_screened > 0 & ...
    sws_screened < 10;


valid_common = ...
    valid_original & ...
    valid_screened;


mask_inclusion_eval = ...
    mask_inclusion_core & ...
    valid_common;


mask_background_eval = ...
    mask_background_core & ...
    valid_common;


fprintf('\n');
fprintf('============================================================\n');
fprintf('EVALUATION ROIS\n');
fprintf('============================================================\n');


fprintf('ROI margin: %.2f mm\n', ...
    ROI_MARGIN*1e3);


fprintf('Inclusion core pixels: %d\n', ...
    nnz(mask_inclusion_eval));


fprintf('Background pixels:     %d\n', ...
    nnz(mask_background_eval));


if nnz(mask_inclusion_eval) < 20

    error( ...
        'Too few pixels in inclusion core ROI.');

end


if nnz(mask_background_eval) < 20

    error( ...
        'Too few pixels in background ROI.');

end


%% ========================================================================
% 20. Original SWS metrics
% ========================================================================

metrics_original = ...
    calculate_SWS_metrics( ...
        sws_original, ...
        mask_inclusion_eval, ...
        mask_background_eval, ...
        c_incl, ...
        c_bkg);


%% ========================================================================
% 21. Q_C-screened SWS metrics
% ========================================================================

metrics_screened = ...
    calculate_SWS_metrics( ...
        sws_screened, ...
        mask_inclusion_eval, ...
        mask_background_eval, ...
        c_incl, ...
        c_bkg);


%% ========================================================================
% 22. Global GT comparison
% ========================================================================

global_mask = ...
    mask_image_interior & ...
    valid_common & ...
    isfinite(sws_gt);


error_original_global = ...
    sws_original(global_mask) - ...
    sws_gt(global_mask);


error_screened_global = ...
    sws_screened(global_mask) - ...
    sws_gt(global_mask);


global_bias_original = ...
    mean( ...
        error_original_global);


global_bias_screened = ...
    mean( ...
        error_screened_global);


global_MAE_original = ...
    mean( ...
        abs(error_original_global));


global_MAE_screened = ...
    mean( ...
        abs(error_screened_global));


global_RMSE_original = ...
    sqrt( ...
        mean( ...
            error_original_global.^2));


global_RMSE_screened = ...
    sqrt( ...
        mean( ...
            error_screened_global.^2));


%% ========================================================================
% 23. Print quantitative results
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('GROUND-TRUTH SWS METRICS\n');
fprintf('============================================================\n');


fprintf('\n');
fprintf('                         ORIGINAL      Q_C SCREENED\n');
fprintf('------------------------------------------------------------\n');


fprintf('INCLUSION\n');


fprintf('Mean SWS [m/s]           %8.4f       %8.4f\n', ...
    metrics_original.inclusion_mean, ...
    metrics_screened.inclusion_mean);


fprintf('Median SWS [m/s]         %8.4f       %8.4f\n', ...
    metrics_original.inclusion_median, ...
    metrics_screened.inclusion_median);


fprintf('Bias [m/s]               %8.4f       %8.4f\n', ...
    metrics_original.inclusion_bias, ...
    metrics_screened.inclusion_bias);


fprintf('MAE [m/s]                %8.4f       %8.4f\n', ...
    metrics_original.inclusion_MAE, ...
    metrics_screened.inclusion_MAE);


fprintf('RMSE [m/s]               %8.4f       %8.4f\n', ...
    metrics_original.inclusion_RMSE, ...
    metrics_screened.inclusion_RMSE);


fprintf('SD [m/s]                 %8.4f       %8.4f\n', ...
    metrics_original.inclusion_std, ...
    metrics_screened.inclusion_std);


fprintf('\nBACKGROUND\n');


fprintf('Mean SWS [m/s]           %8.4f       %8.4f\n', ...
    metrics_original.background_mean, ...
    metrics_screened.background_mean);


fprintf('Median SWS [m/s]         %8.4f       %8.4f\n', ...
    metrics_original.background_median, ...
    metrics_screened.background_median);


fprintf('Bias [m/s]               %8.4f       %8.4f\n', ...
    metrics_original.background_bias, ...
    metrics_screened.background_bias);


fprintf('MAE [m/s]                %8.4f       %8.4f\n', ...
    metrics_original.background_MAE, ...
    metrics_screened.background_MAE);


fprintf('RMSE [m/s]               %8.4f       %8.4f\n', ...
    metrics_original.background_RMSE, ...
    metrics_screened.background_RMSE);


fprintf('SD [m/s]                 %8.4f       %8.4f\n', ...
    metrics_original.background_std, ...
    metrics_screened.background_std);


fprintf('\nIMAGE QUALITY\n');


fprintf('CNR                       %8.4f       %8.4f\n', ...
    metrics_original.CNR, ...
    metrics_screened.CNR);


fprintf('\nGLOBAL INTERIOR\n');


fprintf('Bias [m/s]               %8.4f       %8.4f\n', ...
    global_bias_original, ...
    global_bias_screened);


fprintf('MAE [m/s]                %8.4f       %8.4f\n', ...
    global_MAE_original, ...
    global_MAE_screened);


fprintf('RMSE [m/s]               %8.4f       %8.4f\n', ...
    global_RMSE_original, ...
    global_RMSE_screened);


fprintf('============================================================\n');


%% ========================================================================
% 24. Relative effect of Q_C screening
% ========================================================================

change_inclusion_RMSE = ...
    100 * ...
    ( ...
    metrics_screened.inclusion_RMSE - ...
    metrics_original.inclusion_RMSE) ...
    / metrics_original.inclusion_RMSE;


change_background_RMSE = ...
    100 * ...
    ( ...
    metrics_screened.background_RMSE - ...
    metrics_original.background_RMSE) ...
    / metrics_original.background_RMSE;


change_global_RMSE = ...
    100 * ...
    ( ...
    global_RMSE_screened - ...
    global_RMSE_original) ...
    / global_RMSE_original;


change_CNR = ...
    100 * ...
    ( ...
    metrics_screened.CNR - ...
    metrics_original.CNR) ...
    / metrics_original.CNR;


fprintf('\n');
fprintf('============================================================\n');
fprintf('EFFECT OF Q_C SCREENING\n');
fprintf('============================================================\n');


fprintf('Inclusion RMSE change:   %+7.2f %%\n', ...
    change_inclusion_RMSE);


fprintf('Background RMSE change:  %+7.2f %%\n', ...
    change_background_RMSE);


fprintf('Global RMSE change:      %+7.2f %%\n', ...
    change_global_RMSE);


fprintf('CNR change:              %+7.2f %%\n', ...
    change_CNR);


fprintf('\nNegative RMSE change = improvement.\n');
fprintf('Positive CNR change  = improvement.\n');


fprintf('============================================================\n');


%% ========================================================================
% 25. Results table
% ========================================================================

Metric = ...
    [ ...
    "Inclusion mean SWS"; ...
    "Inclusion median SWS"; ...
    "Inclusion bias"; ...
    "Inclusion MAE"; ...
    "Inclusion RMSE"; ...
    "Inclusion SD"; ...
    "Background mean SWS"; ...
    "Background median SWS"; ...
    "Background bias"; ...
    "Background MAE"; ...
    "Background RMSE"; ...
    "Background SD"; ...
    "CNR"; ...
    "Global bias"; ...
    "Global MAE"; ...
    "Global RMSE"];


Original = ...
    [ ...
    metrics_original.inclusion_mean; ...
    metrics_original.inclusion_median; ...
    metrics_original.inclusion_bias; ...
    metrics_original.inclusion_MAE; ...
    metrics_original.inclusion_RMSE; ...
    metrics_original.inclusion_std; ...
    metrics_original.background_mean; ...
    metrics_original.background_median; ...
    metrics_original.background_bias; ...
    metrics_original.background_MAE; ...
    metrics_original.background_RMSE; ...
    metrics_original.background_std; ...
    metrics_original.CNR; ...
    global_bias_original; ...
    global_MAE_original; ...
    global_RMSE_original];


QC_screened = ...
    [ ...
    metrics_screened.inclusion_mean; ...
    metrics_screened.inclusion_median; ...
    metrics_screened.inclusion_bias; ...
    metrics_screened.inclusion_MAE; ...
    metrics_screened.inclusion_RMSE; ...
    metrics_screened.inclusion_std; ...
    metrics_screened.background_mean; ...
    metrics_screened.background_median; ...
    metrics_screened.background_bias; ...
    metrics_screened.background_MAE; ...
    metrics_screened.background_RMSE; ...
    metrics_screened.background_std; ...
    metrics_screened.CNR; ...
    global_bias_screened; ...
    global_MAE_screened; ...
    global_RMSE_screened];


results_table = ...
    table( ...
        Metric, ...
        Original, ...
        QC_screened);


fprintf('\n');

disp(results_table);


%% ========================================================================
% 26. Coordinates for plotting
% ========================================================================

x_mm = ...
    x_scan * 1e3;


z_mm = ...
    z_scan * 1e3;


%% ========================================================================
% 27. Ground-truth map
% ========================================================================

figure;

imagesc( ...
    x_mm, ...
    z_mm, ...
    sws_gt);


set(gca,'YDir','normal');

axis image;

clim([1.5 4.0]);

colorbar;


xlabel('Lateral position [mm]');
ylabel('Depth [mm]');

title('Ground-truth SWS in ultrasound imaging plane');


%% ========================================================================
% 28. Original SWS with GT boundary
% ========================================================================

figure;

imagesc( ...
    x_mm, ...
    z_mm, ...
    sws_original);


set(gca,'YDir','normal');

axis image;

clim([1.5 4.0]);

colorbar;

hold on;


contour( ...
    X_scan*1e3, ...
    Z_scan*1e3, ...
    mask_inclusion_full, ...
    [1 1], ...
    'k--', ...
    'LineWidth', ...
    1.5);


contour( ...
    X_scan*1e3, ...
    Z_scan*1e3, ...
    mask_inclusion_core, ...
    [1 1], ...
    'k-', ...
    'LineWidth', ...
    1.25);


xlabel('Lateral position [mm]');
ylabel('Depth [mm]');

title('Original SWS with ground-truth inclusion');

hold off;


%% ========================================================================
% 29. Q_C-screened SWS with GT boundary
% ========================================================================

figure;

imagesc( ...
    x_mm, ...
    z_mm, ...
    sws_screened);


set(gca,'YDir','normal');

axis image;

clim([1.5 4.0]);

colorbar;

hold on;


contour( ...
    X_scan*1e3, ...
    Z_scan*1e3, ...
    mask_inclusion_full, ...
    [1 1], ...
    'k--', ...
    'LineWidth', ...
    1.5);


contour( ...
    X_scan*1e3, ...
    Z_scan*1e3, ...
    mask_inclusion_core, ...
    [1 1], ...
    'k-', ...
    'LineWidth', ...
    1.25);


xlabel('Lateral position [mm]');
ylabel('Depth [mm]');

title('Q_C-screened SWS with ground-truth inclusion');

hold off;


%% ========================================================================
% 30. Evaluation ROI map
% ========================================================================

ROI_map = ...
    zeros(size(X_scan));


ROI_map(mask_background_eval) = ...
    1;


ROI_map(mask_inclusion_eval) = ...
    2;


figure;

imagesc( ...
    x_mm, ...
    z_mm, ...
    ROI_map);


set(gca,'YDir','normal');

axis image;

colorbar;


xlabel('Lateral position [mm]');
ylabel('Depth [mm]');

title('Evaluation ROIs: background = 1, inclusion = 2');


%% ========================================================================
% 31. Error maps
% ========================================================================

error_original = ...
    sws_original - ...
    sws_gt;


error_screened = ...
    sws_screened - ...
    sws_gt;


figure;

imagesc( ...
    x_mm, ...
    z_mm, ...
    error_original);


set(gca,'YDir','normal');

axis image;

colorbar;


xlabel('Lateral position [mm]');
ylabel('Depth [mm]');

title('Original SWS error');


figure;

imagesc( ...
    x_mm, ...
    z_mm, ...
    error_screened);


set(gca,'YDir','normal');

axis image;

colorbar;


xlabel('Lateral position [mm]');
ylabel('Depth [mm]');

title('Q_C-screened SWS error');


%% ========================================================================
% 32. Save
% ========================================================================

results = ...
    struct();


results.transform_mode = ...
    transform_mode;


results.plane_mode = ...
    plane_mode;


results.metrics_original = ...
    metrics_original;


results.metrics_screened = ...
    metrics_screened;


results.global_bias_original = ...
    global_bias_original;


results.global_bias_screened = ...
    global_bias_screened;


results.global_MAE_original = ...
    global_MAE_original;


results.global_MAE_screened = ...
    global_MAE_screened;


results.global_RMSE_original = ...
    global_RMSE_original;


results.global_RMSE_screened = ...
    global_RMSE_screened;


results.change_inclusion_RMSE_percent = ...
    change_inclusion_RMSE;


results.change_background_RMSE_percent = ...
    change_background_RMSE;


results.change_global_RMSE_percent = ...
    change_global_RMSE;


results.change_CNR_percent = ...
    change_CNR;


results.c_bkg = ...
    c_bkg;


results.c_incl = ...
    c_incl;


results.ROI_MARGIN = ...
    ROI_MARGIN;


results.IMAGE_MARGIN = ...
    IMAGE_MARGIN;


save( ...
    'low_contrast_QC_SWS_metrics.mat', ...
    'results', ...
    'results_table', ...
    'sws_gt', ...
    'mask_inclusion_full', ...
    'mask_inclusion_core', ...
    'mask_background_core', ...
    'mask_inclusion_eval', ...
    'mask_background_eval', ...
    'error_original', ...
    'error_screened', ...
    'X_scan', ...
    'Z_scan', ...
    '-v7.3');


writetable( ...
    results_table, ...
    'low_contrast_QC_SWS_metrics.csv');


fprintf('\n');
fprintf('============================================================\n');
fprintf('GT SWS EVALUATION COMPLETE\n');
fprintf('============================================================\n');


fprintf('Saved:\n');

fprintf('  low_contrast_QC_SWS_metrics.mat\n');
fprintf('  low_contrast_QC_SWS_metrics.csv\n');


fprintf('============================================================\n');


%% ========================================================================
% LOCAL FUNCTION 1
% ========================================================================

function world_h = ...
    apply_transform( ...
        local_h, ...
        TF, ...
        TF_rev, ...
        mode)

    switch mode

        case 1

            world_h = ...
                TF * ...
                local_h;


        case 2

            world_h = ...
                ( ...
                local_h.' * ...
                TF ...
                ).';


        case 3

            world_h = ...
                TF_rev * ...
                local_h;


        case 4

            world_h = ...
                ( ...
                local_h.' * ...
                TF_rev ...
                ).';


        otherwise

            error('Unknown transformation mode.');

    end

end


%% ========================================================================
% LOCAL FUNCTION 2
% ========================================================================

function metrics = ...
    calculate_SWS_metrics( ...
        sws, ...
        mask_inclusion, ...
        mask_background, ...
        c_incl, ...
        c_bkg)


    %% --------------------------------------------------------------------
    % Inclusion
    % ---------------------------------------------------------------------

    incl = ...
        double( ...
            sws(mask_inclusion));


    incl = ...
        incl( ...
            isfinite(incl));


    incl_error = ...
        incl - ...
        c_incl;


    metrics.inclusion_N = ...
        numel(incl);


    metrics.inclusion_mean = ...
        mean(incl);


    metrics.inclusion_median = ...
        median(incl);


    metrics.inclusion_std = ...
        std(incl);


    metrics.inclusion_bias = ...
        mean(incl_error);


    metrics.inclusion_MAE = ...
        mean( ...
            abs(incl_error));


    metrics.inclusion_RMSE = ...
        sqrt( ...
            mean( ...
                incl_error.^2));


    %% --------------------------------------------------------------------
    % Background
    % ---------------------------------------------------------------------

    bkg = ...
        double( ...
            sws(mask_background));


    bkg = ...
        bkg( ...
            isfinite(bkg));


    bkg_error = ...
        bkg - ...
        c_bkg;


    metrics.background_N = ...
        numel(bkg);


    metrics.background_mean = ...
        mean(bkg);


    metrics.background_median = ...
        median(bkg);


    metrics.background_std = ...
        std(bkg);


    metrics.background_bias = ...
        mean(bkg_error);


    metrics.background_MAE = ...
        mean( ...
            abs(bkg_error));


    metrics.background_RMSE = ...
        sqrt( ...
            mean( ...
                bkg_error.^2));


    %% --------------------------------------------------------------------
    % CNR
    % ---------------------------------------------------------------------

    denominator = ...
        sqrt( ...
            metrics.inclusion_std^2 + ...
            metrics.background_std^2);


    if denominator > 0

        metrics.CNR = ...
            abs( ...
                metrics.inclusion_mean - ...
                metrics.background_mean) ...
            / denominator;

    else

        metrics.CNR = ...
            NaN;

    end

end