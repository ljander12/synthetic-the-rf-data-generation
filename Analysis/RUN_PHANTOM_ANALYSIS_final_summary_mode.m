%% ========================================================================
% RUN_PHANTOM_ANALYSIS
%
% Main controller for quantitative digital-phantom analysis.
% Original displacement / SWS analyses are retained.
%
% Optional literature-based SWS reliability analysis:
%   - Correlation Q (Schoen Jr. et al., 2025; adapted)
%   - RMI-1D (Choi et al., 2015; adapted)
%   - Original and quality-masked stiffness maps for PG and AC
%
% Phase Planarity Q and Log-RMI are NOT included in the outputs.
%
% Requires the UPDATED analyze_SWS_quality.m (integrated version).
% ========================================================================

clear;
clc;
close all;

%% 1. Select phantom case
phantom_case = 'large_inclusions';

% Other examples:
% phantom_case = 'small_inclusions';
% phantom_case = 'high_contrast';
% phantom_case = 'low_contrast';
% phantom_case = 'chaoran_baseline';

%% FINAL EVALUATION MODE (recommended for thesis wrap-up)
% true: read precomputed PG/AC batch results only; NO reconstructions.
% Requires SUMMARIZE_FINAL_SWS_RESULTS.m in this Analysis folder.
RUN_FINAL_SWS_SUMMARY_ONLY = true;
if RUN_FINAL_SWS_SUMMARY_ONLY
    analysis_dir = fileparts(mfilename('fullpath'));
    addpath(analysis_dir);
    final_sws_results = SUMMARIZE_FINAL_SWS_RESULTS(); %#ok<NASGU>
    fprintf('\nFinal SWS evaluation complete. Original analyses were not rerun.\n');
    return;
end

%% 2. Select analyses
RUN_DISPLACEMENT_ANALYSIS = true;
RUN_NOISE_SWEEP           = false;

RUN_SWS_ANALYSIS          = true;
RUN_SWS_NOISE_TEST        = false;

% New: quality and stiffness comparison is part of analyze_SWS_quality.
RUN_PUBLISHED_QUALITY_METRICS = false; % Paused following supervisor feedback

% 'saved' = use existing Q/RMI-1D maps, do not rerun quality algorithms.
% 'auto'  = calculate missing Q/RMI-1D data with original source scripts.
%           Note: old RMI script can require historical phase data in auto.
QUALITY_MODE = 'saved';

% Illustrative quality thresholds, NOT validated accuracy guarantees.
QUALITY_CORR_THRESHOLD = 0.63;
QUALITY_RMI_THRESHOLD  = 0.85;

% Stiffness: E ~= 3*rho*c_s^2 (Young's modulus, near-incompressible)
QUALITY_RHO_KG_M3 = 1000;
QUALITY_MODULUS = 'young';    % 'young' or 'shear'
QUALITY_SHOW_FIGURES = true;

%% 3. Resolve repository paths
analysis_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(analysis_dir);
masteroppgave_dir = fileparts(repo_dir);

addpath(genpath(repo_dir));

toolbox_dir = fullfile( ...
    masteroppgave_dir, ...
    'elastography', ...
    'ElastographyToolbox');

if ~isfolder(toolbox_dir)
    error('ElastographyToolbox was not found at:\n%s', toolbox_dir);
end

addpath(genpath(toolbox_dir));

%% 4. Data and result directories
data_root = fullfile(repo_dir, 'Data');
results_root = fullfile(repo_dir, 'Results');
results_dir = fullfile(results_root, phantom_case);

if ~isfolder(data_root), mkdir(data_root); end
if ~isfolder(results_root), mkdir(results_root); end
if ~isfolder(results_dir), mkdir(results_dir); end

%% 5. Locate dataset
switch lower(phantom_case)
    case 'low_contrast'
        data_name = 'low_contrast_displacement_with_quality.mat';
    case 'large_inclusions'
        data_name = 'large_inclusions_displacement_with_quality.mat';
    case 'small_inclusions'
        data_name = 'small_inclusions_displacement_with_quality.mat';
    case 'high_contrast'
        data_name = 'high_contrast_displacement_with_quality.mat';
    case 'chaoran_baseline'
        data_name = 'chaoran_baseline_displacement_with_quality.mat';
    otherwise
        error('Unknown phantom case: %s', phantom_case);
end

candidate_files = { ...
    fullfile(data_root, phantom_case, data_name), ...
    fullfile(data_root, data_name), ...
    fullfile(analysis_dir, data_name), ...
    fullfile(repo_dir, data_name)};

data_file = '';

for k = 1:numel(candidate_files)
    if isfile(candidate_files{k})
        data_file = candidate_files{k};
        break;
    end
end

if isempty(data_file)
    fprintf('\nDataset not found. Locations checked:\n');
    for k = 1:numel(candidate_files)
        fprintf('  %s\n', candidate_files{k});
    end
    error('Could not locate dataset for %s.', phantom_case);
end

%% 6. Configuration
config = struct();

config.phantom_case = phantom_case;
config.repo_dir = repo_dir;
config.analysis_dir = analysis_dir;
config.data_root = data_root;
config.data_file = data_file;
config.results_dir = results_dir;

% Displacement analysis (unchanged)
config.run_noise_sweep = RUN_NOISE_SWEEP;
config.edge_margin_samples = 3;
config.noise_levels_dB = [40 30 20 10];
config.n_noise_realizations = 10;

% SWS analysis (unchanged)
config.run_sws_noise_test = RUN_SWS_NOISE_TEST;
config.sws_noise_dB = 20;
config.n_sws_noise_realizations = 10;
config.sws_filter_speed = 3.2;
config.sws_roi_margin = 2e-3;
config.sws_image_margin = 3e-3;

% New: two published quality-metric adaptations + stiffness masking.
% analyze_SWS_quality(config) reads these fields directly.
config.run_published_quality_metrics = RUN_PUBLISHED_QUALITY_METRICS;
config.quality_mode = QUALITY_MODE;
config.quality_options = struct( ...
    'rho', QUALITY_RHO_KG_M3, ...
    'modulus', QUALITY_MODULUS, ...
    'corrThreshold', QUALITY_CORR_THRESHOLD, ...
    'rmiThreshold', QUALITY_RMI_THRESHOLD, ...
    'showFigures', QUALITY_SHOW_FIGURES, ...
    'dpi', 220);

% Check unsupported switches before running the expensive SWS analysis.
if RUN_PUBLISHED_QUALITY_METRICS && ~RUN_SWS_ANALYSIS
    error(['RUN_PUBLISHED_QUALITY_METRICS requires RUN_SWS_ANALYSIS = true, ' ...
           'because it is invoked inside analyze_SWS_quality(config).']);
end

supported_quality_cases = { ...
    'large_inclusions', 'low_contrast', ...
    'chaoran_baseline', 'high_contrast'};
if RUN_PUBLISHED_QUALITY_METRICS && ...
        ~any(strcmpi(phantom_case, supported_quality_cases))
    error(['No combined quality-metric workflow is configured for %s. ' ...
           'Set RUN_PUBLISHED_QUALITY_METRICS = false or create the ' ...
           'required Q/RMI-1D results first.'], phantom_case);
end

if RUN_PUBLISHED_QUALITY_METRICS && ...
        ~ismember(lower(QUALITY_MODE), {'saved', 'auto'})
    error('QUALITY_MODE must be ''saved'' or ''auto''.');
end

% In saved mode, fail early if Q/RMI maps have not been computed yet.
% This avoids re-running the PG baseline before discovering missing inputs.
if RUN_PUBLISHED_QUALITY_METRICS && strcmpi(QUALITY_MODE, 'saved')
    q_root = fullfile(results_root, 'sws_measurement_quality', phantom_case);
    corr_file = fullfile(q_root, 'quality_validation.mat');
    rmi_file = fullfile(q_root, 'rmi_1d_diagnostic', 'rmi_1d_results.mat');
    if ~isfile(corr_file) || ~isfile(rmi_file)
        fprintf('\nSaved quality inputs missing (required by QUALITY_MODE = ''saved''):\n');
        if ~isfile(corr_file), fprintf('  Correlation Q: %s\n', corr_file); end
        if ~isfile(rmi_file), fprintf('  RMI-1D:      %s\n', rmi_file); end
        error(['Missing saved Q/RMI-1D results. Run the corresponding ' ...
            'original analysis scripts, or use QUALITY_MODE = ''auto'' ' ...
            'to compute missing stages (may invoke legacy phase prerequisite).']);
    end
end

% Ensure the updated analysis function is actually available on the path.
if RUN_SWS_ANALYSIS && isempty(which('analyze_SWS_quality'))
    error('analyze_SWS_quality.m is not on the MATLAB path.');
end

%% 7. Summary
fprintf('\n');
fprintf('========================================\n');
fprintf('PHANTOM ANALYSIS\n');
fprintf('========================================\n');
fprintf('Phantom case:          %s\n', phantom_case);
fprintf('Dataset:               %s\n', data_file);
fprintf('Displacement analysis: %s\n', on_off(RUN_DISPLACEMENT_ANALYSIS));
fprintf('Displacement noise:    %s\n', on_off(RUN_NOISE_SWEEP));
fprintf('SWS analysis:          %s\n', on_off(RUN_SWS_ANALYSIS));
fprintf('SWS noise test:        %s\n', on_off(RUN_SWS_NOISE_TEST));
fprintf('Correlation Q + RMI:  %s\n', on_off(RUN_PUBLISHED_QUALITY_METRICS));
if RUN_PUBLISHED_QUALITY_METRICS
    fprintf('Quality source:        %s\n', QUALITY_MODE);
    fprintf('Quality thresholds:    Corr %.2f | RMI %.2f\n', ...
        QUALITY_CORR_THRESHOLD, QUALITY_RMI_THRESHOLD);
    fprintf('Stiffness modulus:     %s | rho %.0f kg/m^3\n', ...
        QUALITY_MODULUS, QUALITY_RHO_KG_M3);
end
fprintf('SWS filter speed:      %.2f m/s\n', config.sws_filter_speed);
fprintf('========================================\n\n');

%% 8. Run displacement analysis
if RUN_DISPLACEMENT_ANALYSIS
    displacement_results = analyze_displacement_quality(config); %#ok<NASGU>
end

%% 9. Run SWS analysis
% The UPDATED analyze_SWS_quality.m automatically performs the optional
% Correlation Q + RMI-1D + stiffness-mapping stage when its config flag is ON.
if RUN_SWS_ANALYSIS
    sws_results = analyze_SWS_quality(config); %#ok<NASGU>
end

%% 10. Final summary
fprintf('\n');
fprintf('========================================\n');
fprintf('PHANTOM ANALYSIS COMPLETE\n');
fprintf('========================================\n');
fprintf('Baseline results directory:\n%s\n', results_dir);
if RUN_PUBLISHED_QUALITY_METRICS
    fprintf('Quality/stiffness figures directory:\n%s\n', ...
        fullfile(results_root, 'sws_measurement_quality', ...
        'combined_metrics_published'));
end
fprintf('========================================\n\n');

%% Local function
function txt = on_off(value)
if value
    txt = 'ON';
else
    txt = 'OFF';
end
end
