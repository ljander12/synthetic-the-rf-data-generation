%% ========================================================================
%
%   CONTINUE_LOW_CONTRAST_FIELDII
%
%   Continue the low-contrast synthetic-data pipeline from the saved
%   k-Wave ground-truth checkpoint.
%
%   IMPORTANT:
%       pstdElastic3D is NOT run in this script.
%
%   Input:
%       low_contrast_ground_truth_theta90.mat
%
%   Output:
%       low_contrast_RF_checkpoint.mat
%
% ========================================================================

clear;
clc;
close all;


%% ========================================================================
% 1. Add dependencies
% ========================================================================

addpath(fullfile(getenv('HOME'),'k-Wave'));

% Field II installation
fieldii_root = ...
    'C:\Users\louis\Documents\MATLAB\Field_II_ver_3_30_windows';

addpath(fieldii_root);

% Start parallel pool if necessary
p = gcp('nocreate');

if isempty(p)
    parpool;
end

% Make Field II available on all workers
pctRunOnAll addpath('C:\Users\louis\Documents\MATLAB\Field_II_ver_3_30_windows')

% Project functions
addpath("FIELD_fun");
addpath("kwave_fun");
addpath("utils");

%% ========================================================================
% 2. Load saved k-Wave ground truth
% ========================================================================

GT_FILE = 'low_contrast_ground_truth_theta90.mat';

fprintf('\n========================================\n');
fprintf('LOADING K-WAVE CHECKPOINT\n');
fprintf('%s\n', GT_FILE);
fprintf('========================================\n\n');

load(GT_FILE);

fprintf('Ground-truth size: %s\n', mat2str(size(disp_gt)));
fprintf('Motion intervals:  %d\n', size(disp_gt,3));

fprintf('Background SWS: %.3f m/s\n', ...
    shear_params.c_shear_bkg);

fprintf('Inclusion SWS:  %.3f m/s\n', ...
    shear_params.c_shear_incl);

fprintf('Probe angle:    %.1f deg\n', ...
    rad2deg(prb_theta(1)));


%% ========================================================================
% 3. Ultrasound parameters
% ========================================================================

c0 = 1540;

f0 = 5.1333e6;

lambda = c0/f0;

pulse_duration = 2.5;


%% ------------------------------------------------------------------------
% Probe
% -------------------------------------------------------------------------

prb.element_height = 5e-3;

prb.pitch = 0.300e-3;

prb.kerf = 0.03e-3;

prb.element_width = ...
    prb.pitch - prb.kerf;

prb.N = 128;


%% ------------------------------------------------------------------------
% Sampling
% -------------------------------------------------------------------------

fs = 100e6;

dt = 1/fs;

attenfreq = 0.55;


%% ========================================================================
% 4. Recreate k-Wave geometry
% ========================================================================
%
% We only recreate the grid and sensor geometry required to map the
% ALREADY SAVED motion onto the scatterers.
%
% pstdElastic3D is NOT called.
% ========================================================================

fprintf('\nRecreating k-Wave grid and sensor geometry...\n');

[kgrid, medium] = ...
    create_phantom(shear_params);


%% ------------------------------------------------------------------------
% Recreate coordinate transforms
% -------------------------------------------------------------------------

[TF, TF_rev] = ...
    coor_transformation( ...
        prb_center, ...
        prb_theta);


%% ------------------------------------------------------------------------
% Recreate sensor
% -------------------------------------------------------------------------

sensor = ...
    create_sensor( ...
        kgrid, ...
        array_length, ...
        prb.element_height, ...
        scene_depth, ...
        TF);


fprintf('Geometry recreated.\n');


%% ========================================================================
% 5. Create initial scatterer map
% ========================================================================

fprintf('\nCreating scatterer distribution...\n');

sca_per_cell = 50;

scatterers = ...
    create_scatterers( ...
        prb, ...
        scene_depth, ...
        pulse_duration, ...
        lambda, ...
        sca_per_cell);


%% ========================================================================
% 6. Apply saved k-Wave motion to scatterers
% ========================================================================

fprintf('Applying saved motion to scatterers...\n');

sca_mesh = ...
    time_vary_sca_kwave( ...
        kgrid, ...
        sensor, ...
        motion, ...
        scatterers, ...
        TF);

fprintf('Time-varying scatterer mesh created.\n');


%% ========================================================================
% 7. Define scatterer amplitudes
% ========================================================================

fprintf('Defining scatterer amplitudes...\n');

amp = ...
    define_phantom_scatterer_amplitudes( ...
        sca_mesh, ...
        shear_params);


%% ========================================================================
% 8. Plane-wave transmit sequence
% ========================================================================

alpha_max = deg2rad(15);

Na = 1;

if Na == 1

    alpha = 0;

else

    alpha = ...
        linspace( ...
            -alpha_max, ...
             alpha_max, ...
             Na);

end


%% ========================================================================
% 9. Ultrasound pulse
% ========================================================================

[impulse_response, excitation, lag] = ...
    pulse_setup( ...
        f0, ...
        pulse_duration, ...
        0.65, ...
        fs);


%% ========================================================================
% 10. Run Field II
% ========================================================================

fprintf('\n========================================\n');
fprintf('STARTING FIELD II RF SIMULATION\n');
fprintf('========================================\n\n');

tic;

RF = ...
    FIELD_calc_RF( ...
        prb, ...
        sca_mesh, ...
        amp, ...
        TF_rev, ...
        alpha, ...
        scene_depth, ...
        fs, ...
        attenfreq, ...
        f0, ...
        c0, ...
        impulse_response, ...
        excitation);

fieldII_time = toc;


fprintf('\n========================================\n');
fprintf('FIELD II COMPLETE\n');
fprintf('========================================\n');

fprintf('RF size: %s\n', ...
    mat2str(size(RF)));

fprintf('Field II runtime: %.2f min\n', ...
    fieldII_time/60);


%% ========================================================================
% 11. Save RF checkpoint
% ========================================================================
%
% After this file is saved, Field II does not need to be rerun for
% beamforming/autocorrelation experiments.
% ========================================================================

RF_SAVE_FILE = ...
    'low_contrast_RF_checkpoint.mat';


fprintf('\nSaving RF checkpoint...\n');

save( ...
    RF_SAVE_FILE, ...
    'RF', ...
    'disp_gt', ...
    'motion', ...
    'vec_T', ...
    'shear_params', ...
    'prb', ...
    'prf', ...
    'array_length', ...
    'scene_depth', ...
    'prb_center', ...
    'prb_theta', ...
    'TF', ...
    'TF_rev', ...
    'f0', ...
    'fs', ...
    'c0', ...
    'dt', ...
    'lag', ...
    'alpha', ...
    'Na', ...
    'pulse_duration', ...
    'attenfreq', ...
    'fieldII_time', ...
    '-v7.3');


fprintf('\n========================================\n');
fprintf('RF CHECKPOINT SAVED\n');
fprintf('%s\n', RF_SAVE_FILE);
fprintf('========================================\n\n');


%% ========================================================================
% 12. Stop before USTB processing
% ========================================================================

fprintf('Field II stage finished successfully.\n');
fprintf('Stopping before beamforming/autocorrelation.\n');

return