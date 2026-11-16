%% RUN_SWS_QUALITY.m -- Open this file in MATLAB and click RUN.
% Runs ONLY Correlation Q and RMI-1D. Also generates stiffness maps with
% the same SWS estimate, before and after each quality mask.
% No Phase Planarity Q or experimental Log-RMI in the reported results.
clearvars; clc;

%% Select phantom(s) and settings
phantoms = 'high_contrast';    % 'all' or a case name
mode = 'auto';                % 'auto' fills missing maps; 'saved' reads only

opts = struct();
opts.rho = 1000;              % kg/m^3
opts.modulus = 'young';       % 'young' -> E~3*rho*c^2; 'shear' -> mu=rho*c^2
opts.corrThreshold = 0.63;   % illustrative PIXELWISE Correlation Q cutoff
opts.rmiThreshold = 0.85;    % exploratory RMI cutoff; not universal
opts.stiffnessColorMaxKPa = [];  % auto: shared 99th percentile for PG + AC
opts.showFigures = true;     % false for HPC / unattended batch mode
opts.dpi = 220;

[T, globalSummary, maskSummary] = ...
    run_all_sws_quality_metrics(phantoms, mode, opts);
