clear;
clc;
close all;

%% ================================================================
% Load completed Field II simulation
% ================================================================

load('low_contrast_RF_checkpoint.mat');

phantom_case = 'low_contrast';

fprintf('RF size:      %s\n', mat2str(size(RF)));
fprintf('GT disp size: %s\n', mat2str(size(disp_gt)));


%% ================================================================
% Add USTB / ElastographyToolbox
% ================================================================

toolbox_dir = ...
    'C:\Users\louis\Documents\GitHub\MasterOppgave\elastography\ElastographyToolbox';

addpath(genpath(toolbox_dir));


%% ================================================================
% Create USTB probe
% ================================================================

probe = uff.linear_array();

probe.element_height = prb.element_height;
probe.pitch          = prb.pitch;
probe.element_width  = prb.element_width;
probe.N              = prb.N;


%% ================================================================
% Create transmit sequence
% ================================================================

seq(Na) = uff.wave();

for n = 1:Na

    seq(n) = uff.wave();

    seq(n).probe = probe;

    seq(n).source.azimuth = alpha(n);
    seq(n).source.distance = Inf;

    seq(n).sound_speed = c0;

    seq(n).delay = -lag*dt;

end


%% ================================================================
% Create pulse
% ================================================================

pulse = uff.pulse();

pulse.fractional_bandwidth = 0.65;
pulse.center_frequency = f0;


%% ================================================================
% Create channel data
% ================================================================

channel_data = uff.channel_data();

channel_data.sampling_frequency = fs;
channel_data.sound_speed = c0;
channel_data.initial_time = 0;

channel_data.pulse = pulse;
channel_data.probe = probe;
channel_data.sequence = seq;

channel_data.PRF = prf;

channel_data.data = RF;


%% ================================================================
% Define scan
% ================================================================

sca = uff.linear_scan( ...
    'x_axis', ...
    linspace( ...
        -array_length/2, ...
         array_length/2, ...
         512).', ...
    'z_axis', ...
    linspace( ...
        5e-3, ...
        scene_depth, ...
        1024).');


%% ================================================================
% DAS beamforming
% ================================================================

fprintf('\n');
fprintf('========================================\n');
fprintf('STARTING USTB BEAMFORMING\n');
fprintf('========================================\n\n');

F_number = 1.2;

das = midprocess.das();

das.channel_data = channel_data;
das.scan = sca;

das.transmit_apodization.window = uff.window.tukey25;
das.transmit_apodization.f_number = F_number;

das.receive_apodization.window = uff.window.tukey25;
das.receive_apodization.f_number = F_number;

b_data = das.go();

fprintf('\nBeamforming complete.\n');
fprintf('Beamformed data size: %s\n', ...
    mat2str(size(b_data.data)));


%% ================================================================
% Autocorrelation displacement estimation
% ================================================================

fprintf('\n');
fprintf('========================================\n');
fprintf('STARTING AUTOCORRELATION ESTIMATION\n');
fprintf('========================================\n\n');

autocorr_params.z_gate = 4;
autocorr_params.x_gate = 2;
autocorr_params.packet_size = 4;

displacement_estimator = ...
    postprocess.autocorrelation_displacement_estimation();

displacement_estimator.channel_data = channel_data;
displacement_estimator.input = b_data;

displacement_estimator.z_gate = autocorr_params.z_gate;
displacement_estimator.x_gate = autocorr_params.x_gate;
displacement_estimator.packet_size = autocorr_params.packet_size;

displacement_estimation_autocorr = ...
    displacement_estimator.go();


%% ================================================================
% Extract autocorrelation quality indicator C
% ================================================================

quality_indicator_autocorr = ...
    displacement_estimator.quality_indicator;


%% ================================================================
% Check outputs
% ================================================================

fprintf('\n');
fprintf('========================================\n');
fprintf('AUTOCORRELATION OUTPUT\n');
fprintf('========================================\n');

fprintf('Ground-truth displacement size: ');
disp(size(disp_gt));

fprintf('Estimated displacement size:    ');
disp(size(displacement_estimation_autocorr.data));

fprintf('Quality indicator size:         ');
disp(size(quality_indicator_autocorr));

expected_estimator_frames = ...
    size(disp_gt,3) - ...
    (autocorr_params.packet_size - 1) + 1;

actual_estimator_frames = ...
    size(displacement_estimation_autocorr.data,4);

fprintf('Expected estimator frames: %d\n', ...
    expected_estimator_frames);

fprintf('Actual estimator frames:   %d\n', ...
    actual_estimator_frames);

if actual_estimator_frames ~= expected_estimator_frames
    warning('Unexpected number of estimator frames.');
end

if isequal( ...
        size(displacement_estimation_autocorr.data), ...
        size(quality_indicator_autocorr))

    fprintf('Quality-indicator dimensions match displacement output.\n');

else
    warning('Quality-indicator dimensions do not match.');
end

fprintf('========================================\n\n');


%% ================================================================
% Save complete low-contrast dataset
% ================================================================

SAVE_FILE = ...
    'low_contrast_displacement_with_quality.mat';

fprintf('Saving complete low-contrast dataset...\n');

save( ...
    SAVE_FILE, ...
    'b_data', ...
    'channel_data', ...
    'disp_gt', ...
    'displacement_estimation_autocorr', ...
    'quality_indicator_autocorr', ...
    'shear_params', ...
    'autocorr_params', ...
    'prb_center', ...
    'prb_theta', ...
    'scene_depth', ...
    'prf', ...
    'array_length', ...
    'phantom_case', ...
    '-v7.3');

fprintf('\n');
fprintf('========================================\n');
fprintf('LOW-CONTRAST DATASET SAVED\n');
fprintf('%s\n', SAVE_FILE);
fprintf('========================================\n');