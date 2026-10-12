clear;
clc;
close all;

%% ========================================================================
% RUN_PHANTOM_ANALYSIS
%
% Main controller for synthetic THE phantom analysis.
%
% Change ONLY phantom_case and the RUN_* switches below.
%
% Supported cases:
%
%   'low_contrast'
%   'large_inclusions'
%   'small_inclusions'
%
% ========================================================================


%% ========================================================================
% 1. Select phantom
% ========================================================================

phantom_case = ...
    'low_contrast';

% Examples:
%
% phantom_case = 'large_inclusions';
% phantom_case = 'small_inclusions';


%% ========================================================================
% 2. Select analyses
% ========================================================================

RUN_DISPLACEMENT_ANALYSIS = ...
    false;

RUN_NOISE_SWEEP = ...
    false;

RUN_SWS_ANALYSIS = ...
    true;

RUN_SWS_NOISE_TEST = ...
    false;


%% ========================================================================
% 3. Locate repository
% ========================================================================
%
% Expected location:
%
% synthetic-the-rf-data-generation/
%     Analysis/
%         run_phantom_analysis.m
%
% ========================================================================

analysis_dir = ...
    fileparts( ...
        mfilename('fullpath'));


repo_dir = ...
    fileparts( ...
        analysis_dir);


masteroppgave_dir = ...
    fileparts( ...
        repo_dir);


fprintf('\n');
fprintf('============================================================\n');
fprintf('PHANTOM ANALYSIS CONTROLLER\n');
fprintf('============================================================\n');

fprintf('Phantom:        %s\n', ...
    phantom_case);

fprintf('Repository:     %s\n', ...
    repo_dir);

fprintf('Analysis dir:   %s\n', ...
    analysis_dir);


%% ========================================================================
% 4. Add repository code
% ========================================================================

addpath( ...
    genpath( ...
        repo_dir));


%% ========================================================================
% 5. Add elastography toolbox / USTB
% ========================================================================

toolbox_dir = ...
    fullfile( ...
        masteroppgave_dir, ...
        'elastography', ...
        'ElastographyToolbox');


if ~isfolder(toolbox_dir)

    error( ...
        ['Elastography toolbox not found:\n%s'], ...
        toolbox_dir);

end


addpath( ...
    genpath( ...
        toolbox_dir));


rehash toolboxcache;


fprintf('Elastography:   %s\n', ...
    toolbox_dir);


%% ========================================================================
% 6. Data and result directories
% ========================================================================

data_root = ...
    fullfile( ...
        repo_dir, ...
        'Data');


results_root = ...
    fullfile( ...
        repo_dir, ...
        'Results');


if ~isfolder(data_root)

    mkdir(data_root);

end


if ~isfolder(results_root)

    mkdir(results_root);

end


case_results_dir = ...
    fullfile( ...
        results_root, ...
        phantom_case);


if ~isfolder(case_results_dir)

    mkdir(case_results_dir);

end


%% ========================================================================
% 7. Build configuration structure
% ========================================================================

config = ...
    struct();


config.phantom_case = ...
    phantom_case;


config.repo_dir = ...
    repo_dir;


config.analysis_dir = ...
    analysis_dir;


config.data_root = ...
    data_root;


config.results_dir = ...
    case_results_dir;


%% ========================================================================
% 8. Common displacement-analysis settings
% ========================================================================

config.run_noise_sweep = ...
    RUN_NOISE_SWEEP;


config.Q_threshold = ...
    0.05;


config.edge_margin_samples = ...
    3;


config.noise_levels_dB = ...
    [40 30 20 10];


config.n_noise_realizations = ...
    10;


%% ========================================================================
% 9. SWS-analysis settings
% ========================================================================

config.run_sws_noise_test = ...
    RUN_SWS_NOISE_TEST;


config.sws_noise_dB = ...
    20;


config.n_sws_noise_realizations = ...
    10;


config.sws_filter_speed = ...
    3.2;


config.sws_roi_margin = ...
    2e-3;


config.sws_image_margin = ...
    3e-3;


%% ========================================================================
% 10. Resolve input dataset
% ========================================================================
%
% Search order:
%
%   Data/<phantom_case>/
%   Data/
%   Analysis/
%   repository root
%
% ========================================================================

switch lower(phantom_case)

    case 'low_contrast'

        dataset_name = ...
            'low_contrast_displacement_with_quality.mat';


    case 'large_inclusions'

        dataset_name = ...
            'large_inclusions_displacement_with_quality.mat';


    case 'small_inclusions'

        dataset_name = ...
            'small_inclusions_displacement_with_quality.mat';


    otherwise

        error( ...
            'Unknown phantom case: %s', ...
            phantom_case);

end


candidate_files = ...
    { ...
    fullfile( ...
        data_root, ...
        phantom_case, ...
        dataset_name), ...
    fullfile( ...
        data_root, ...
        dataset_name), ...
    fullfile( ...
        analysis_dir, ...
        dataset_name), ...
    fullfile( ...
        repo_dir, ...
        dataset_name) ...
    };


config.data_file = ...
    '';


for k = 1:numel(candidate_files)

    if isfile(candidate_files{k})

        config.data_file = ...
            candidate_files{k};

        break;

    end

end


if isempty(config.data_file)

    fprintf('\n');
    fprintf('Dataset was not found.\n');
    fprintf('Searched locations:\n');

    for k = 1:numel(candidate_files)

        fprintf('  %s\n', ...
            candidate_files{k});

    end


    error( ...
        'Could not find dataset: %s', ...
        dataset_name);

end


fprintf('\n');
fprintf('Input dataset:\n');
fprintf('%s\n', ...
    config.data_file);


fprintf('\n');
fprintf('Results directory:\n');
fprintf('%s\n', ...
    config.results_dir);


%% ========================================================================
% 11. Print selected analysis configuration
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('ANALYSIS SETTINGS\n');
fprintf('============================================================\n');


fprintf('Displacement analysis:   %s\n', ...
    on_off( ...
        RUN_DISPLACEMENT_ANALYSIS));


fprintf('Displacement noise:      %s\n', ...
    on_off( ...
        RUN_NOISE_SWEEP));


fprintf('SWS analysis:            %s\n', ...
    on_off( ...
        RUN_SWS_ANALYSIS));


fprintf('SWS noise test:          %s\n', ...
    on_off( ...
        RUN_SWS_NOISE_TEST));


fprintf('Q_C threshold:           %.3f\n', ...
    config.Q_threshold);


if RUN_SWS_NOISE_TEST

    fprintf('SWS noise level:         %d dB\n', ...
        config.sws_noise_dB);

    fprintf('SWS realizations:        %d\n', ...
        config.n_sws_noise_realizations);

end


fprintf('============================================================\n');


%% ========================================================================
% 12. Run displacement quality analysis
% ========================================================================

if RUN_DISPLACEMENT_ANALYSIS

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('STARTING DISPLACEMENT QUALITY ANALYSIS\n');
    fprintf('============================================================\n');


    displacement_results = ...
        analyze_displacement_quality( ...
            config);


    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('DISPLACEMENT QUALITY ANALYSIS COMPLETE\n');
    fprintf('============================================================\n');

end


%% ========================================================================
% 13. Run SWS quality analysis
% ========================================================================

if RUN_SWS_ANALYSIS

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('STARTING SWS QUALITY ANALYSIS\n');
    fprintf('============================================================\n');


    sws_results = ...
        analyze_SWS_quality( ...
            config);


    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('SWS QUALITY ANALYSIS COMPLETE\n');
    fprintf('============================================================\n');

end


%% ========================================================================
% 14. Final summary
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('PHANTOM ANALYSIS FINISHED\n');
fprintf('============================================================\n');


fprintf('Phantom:      %s\n', ...
    phantom_case);


fprintf('Results:      %s\n', ...
    config.results_dir);


if RUN_DISPLACEMENT_ANALYSIS

    fprintf('Displacement: complete\n');

end


if RUN_SWS_ANALYSIS

    fprintf('SWS:          complete\n');

end


fprintf('============================================================\n');


%% ========================================================================
% LOCAL FUNCTION
% ========================================================================

function txt = ...
    on_off(value)

if value

    txt = ...
        'ON';

else

    txt = ...
        'OFF';

end

end