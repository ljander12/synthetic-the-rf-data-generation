%% ========================================================================
% RUN_PHANTOM_ANALYSIS
%
% Main controller for quantitative phantom analysis.
% Uses established/reference-based metrics only.
% ========================================================================

clear;
clc;
close all;

%% 1. Select phantom case
phantom_case = 'large_inclusions';

% Other examples:
% phantom_case = 'large_inclusions';
% phantom_case = 'small_inclusions';
% phantom_case = 'high_contrast';
% phantom_case = 'low_contrast';
% phantom_case = 'chaoran_baseline'

%% 2. Select analyses
RUN_DISPLACEMENT_ANALYSIS = true;
RUN_NOISE_SWEEP           = false;

RUN_SWS_ANALYSIS          = false;
RUN_SWS_NOISE_TEST        = false;

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
    error(['ElastographyToolbox was not found at:\n%s'], toolbox_dir);
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

% Displacement analysis
config.run_noise_sweep = RUN_NOISE_SWEEP;
config.edge_margin_samples = 3;
config.noise_levels_dB = [40 30 20 10];
config.n_noise_realizations = 10;

% SWS analysis
config.run_sws_noise_test = RUN_SWS_NOISE_TEST;
config.sws_noise_dB = 20;
config.n_sws_noise_realizations = 10;
config.sws_filter_speed = 3.2;
config.sws_roi_margin = 2e-3;
config.sws_image_margin = 3e-3;

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
fprintf('========================================\n\n');

%% 8. Run displacement analysis
if RUN_DISPLACEMENT_ANALYSIS
    displacement_results = analyze_displacement_quality(config); %#ok<NASGU>
end

%% 9. Run SWS analysis
if RUN_SWS_ANALYSIS
    sws_results = analyze_SWS_quality(config); %#ok<NASGU>
end

%% 10. Final summary
fprintf('\n');
fprintf('========================================\n');
fprintf('PHANTOM ANALYSIS COMPLETE\n');
fprintf('========================================\n');
fprintf('Results directory:\n%s\n', results_dir);
fprintf('========================================\n\n');

%% Local function
function txt = on_off(value)
if value
    txt = 'ON';
else
    txt = 'OFF';
end
end
