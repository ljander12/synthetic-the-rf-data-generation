%% ========================================================================
%
%   SIM_PHANTOM_MAIN
%
%   Generate synthetic ultrasound RF channel data for the low-contrast
%   phantom used in the time-harmonic elastography case study.
%
%   Pipeline:
%
%       Digital phantom
%           ->
%       k-Wave shear-wave simulation
%           ->
%       PRF-sampled ground-truth displacement
%           ->
%       Time-varying scatterer distribution
%           ->
%       Field II RF simulation
%           ->
%       USTB beamforming
%           ->
%       Autocorrelation displacement estimation
%           ->
%       Autocorrelation quality indicator C
%
%   Two files are saved:
%
%       1) low_contrast_ground_truth_theta90.mat
%
%          Contains the k-Wave ground-truth simulation. This allows the
%          expensive k-Wave part to be recovered later if required.
%
%       2) low_contrast_displacement_with_quality.mat
%
%          Contains beamformed data, displacement estimates, ground truth
%          and the autocorrelation quality indicator required for the
%          quality-metric experiments.
%
%   Dependencies:
%       - k-Wave
%       - Field II
%       - USTB
%
% ========================================================================


clear;
clc;
close all;


%% ========================================================================
% 1. Add dependencies
% ========================================================================

addpath(fullfile(getenv('HOME'),'k-Wave'));
addpath(fullfile(getenv('HOME'),'Field_II'));

addpath("FIELD_fun");
addpath("kwave_fun");
addpath("utils");


%% ========================================================================
% 2. Parameters for shear-wave simulation
% ========================================================================

phantom_case = 'low_contrast';

shear_params = get_phantom_params(phantom_case);


%% ------------------------------------------------------------------------
% Verify selected phantom
% -------------------------------------------------------------------------

fprintf('\n');
fprintf('========================================\n');
fprintf('PHANTOM CASE: %s\n', upper(phantom_case));
fprintf('Background SWS: %.3f m/s\n', shear_params.c_shear_bkg);
fprintf('Inclusion SWS:  %.3f m/s\n', shear_params.c_shear_incl);
fprintf('Source freq:    %.1f Hz\n', shear_params.source_freq);
fprintf('========================================\n\n');


%% ========================================================================
% 3. Parameters for ultrasound pulse-echo imaging simulation
% ========================================================================

c0 = 1540;                         % Compression-wave speed [m/s]

f0 = 5.1333e6;                     % Transducer center frequency [Hz]

lambda = c0/f0;                    % Ultrasound wavelength [m]


%% ------------------------------------------------------------------------
% Probe
% -------------------------------------------------------------------------

prb.element_height = 5e-3;         % Element height [m]

prb.pitch = 0.300e-3;              % Element pitch [m]

prb.kerf = 0.03e-3;                % Gap between elements [m]

prb.element_width = ...
    prb.pitch - prb.kerf;          % Element width [m]

prb.N = 128;                       % Number of probe elements


%% ------------------------------------------------------------------------
% Acquisition
% -------------------------------------------------------------------------

pulse_duration = 2.5;              % Pulse duration [cycles]

prf = 3000;                        % Pulse repetition frequency [Hz]

fs = 100e6;                        % Sampling frequency [Hz]

dt = 1/fs;                         % Sampling interval [s]

attenfreq = 0.55;                  % Attenuation [dB/MHz/cm]


%% ------------------------------------------------------------------------
% Imaging geometry
% -------------------------------------------------------------------------

array_length = ...
    (prb.N-1)*prb.pitch;

scene_depth = ...
    50e-3;


%% ========================================================================
% 4. Create k-Wave phantom
% ========================================================================

fprintf('\nCreating k-Wave phantom...\n');

[kgrid, medium] = ...
    create_phantom(shear_params);


%% ------------------------------------------------------------------------
% Verify actual shear-wave-speed distribution
% -------------------------------------------------------------------------

fprintf('\n');
fprintf('----------------------------------------\n');
fprintf('K-WAVE PHANTOM CHECK\n');
fprintf('Phantom case:   %s\n', phantom_case);

fprintf( ...
    'Background cs: %.3f m/s\n', ...
    shear_params.c_shear_bkg);

fprintf( ...
    'Inclusion cs:  %.3f m/s\n', ...
    shear_params.c_shear_incl);

fprintf( ...
    'Minimum cs:    %.3f m/s\n', ...
    min(medium.sound_speed_shear(:)));

fprintf( ...
    'Maximum cs:    %.3f m/s\n', ...
    max(medium.sound_speed_shear(:)));

fprintf('----------------------------------------\n\n');


%% ========================================================================
% 5. Define k-Wave shear-wave source
% ========================================================================

source = ...
    create_source( ...
        kgrid, ...
        shear_params.source_freq);


%% ========================================================================
% 6. Define k-Wave sensor / probe geometry
% ========================================================================

% Probe center
prb_center = ...
    [0 -20e-3 kgrid.z_vec(1)];


% Probe orientation:
%
% 90 degrees was selected based on previous propagation-geometry
% validation. This aligns the dominant shear-wave propagation direction
% with the ultrasound imaging plane.

prb_theta = ...
    [deg2rad(90), deg2rad(0)];


%% ------------------------------------------------------------------------
% Coordinate transformation
% -------------------------------------------------------------------------

[TF, TF_rev] = ...
    coor_transformation( ...
        prb_center, ...
        prb_theta);


%% ------------------------------------------------------------------------
% k-Wave sensor
% -------------------------------------------------------------------------

sensor = ...
    create_sensor( ...
        kgrid, ...
        array_length, ...
        prb.element_height, ...
        scene_depth, ...
        TF);


fprintf('Probe angle: %.1f degrees\n', ...
    rad2deg(prb_theta(1)));


%% ========================================================================
% 7. Visualize shear-wave simulation setup
% ========================================================================

viewer = viewer3d( ...
    BackgroundColor="white", ...
    BackgroundGradient="off");


%% ------------------------------------------------------------------------
% Shear-wave-speed phantom
% -------------------------------------------------------------------------

volshow( ...
    permute( ...
        medium.sound_speed_shear, ...
        [2,1,3]), ...
    Parent=viewer, ...
    Colormap=cool, ...
    Alphamap=[0.005,0.1]);


%% ------------------------------------------------------------------------
% Source
% -------------------------------------------------------------------------

volshow( ...
    permute( ...
        source.s_mask, ...
        [2,1,3]), ...
    Parent=viewer, ...
    Colormap=hot, ...
    Alphamap=[0.00,0.2]);


%% ------------------------------------------------------------------------
% Sensor
% -------------------------------------------------------------------------

volshow( ...
    permute( ...
        sensor.mask, ...
        [2,1,3]), ...
    Parent=viewer, ...
    Colormap=parula, ...
    Alphamap=[0.00,0.1]);


%% ========================================================================
% 8. Run k-Wave shear-wave simulation
% ========================================================================

fprintf('\n');
fprintf('========================================\n');
fprintf('STARTING K-WAVE SIMULATION\n');
fprintf('========================================\n\n');


display_mask = ...
    source.s_mask;


DATA_CAST = ...
    'single';


input_args_shear = { ...
    'DisplayMask', display_mask, ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'PlotSim', false, ...
    'DataCast', DATA_CAST};


sensor_data_shear = ...
    pstdElastic3D( ...
        kgrid, ...
        medium, ...
        source, ...
        sensor, ...
        input_args_shear{:});


fprintf('\n');
fprintf('========================================\n');
fprintf('K-WAVE SIMULATION COMPLETE\n');
fprintf('========================================\n\n');


%% ========================================================================
% 9. Sample motion at ultrasound PRF
% ========================================================================

fprintf('Sampling k-Wave motion at ultrasound PRF...\n');


t_start = ...
    sqrt( ...
        shear_params.xrange^2 + ...
        shear_params.yrange^2 + ...
        shear_params.zrange^2) ...
    / shear_params.c_shear_bkg;


t_start = ...
    t_start + 50*kgrid.dt;


[motion, vec_T] = ...
    sample_motion_at_prf( ...
        "kwave", ...
        sensor_data_shear, ...
        kgrid, ...
        prf, ...
        t_start);


fprintf('Number of sampled motion intervals: %d\n', ...
    length(vec_T));


%% ========================================================================
% 10. Calculate ground-truth displacement
% ========================================================================

fprintf('\nCalculating ground-truth displacement...\n');


disp_gt = ...
    get_GT_kwave( ...
        kgrid, ...
        sensor, ...
        array_length, ...
        scene_depth, ...
        motion, ...
        prb_theta, ...
        TF);


fprintf( ...
    'Ground-truth displacement size: %s\n', ...
    mat2str(size(disp_gt)));


fprintf( ...
    'Ground-truth displacement range: %.4f to %.4f um\n', ...
    min(disp_gt(:))*1e6, ...
    max(disp_gt(:))*1e6);


%% ========================================================================
% 11. Save k-Wave ground-truth checkpoint
% ========================================================================
%
% Saving here means that the expensive k-Wave simulation does not need
% to be repeated if something fails later in Field II or USTB.
% ========================================================================

GT_SAVE_FILE = ...
    'low_contrast_ground_truth_theta90.mat';


save( ...
    GT_SAVE_FILE, ...
    'disp_gt', ...
    'motion', ...
    'vec_T', ...
    'shear_params', ...
    'prf', ...
    'array_length', ...
    'scene_depth', ...
    'prb_center', ...
    'prb_theta', ...
    'TF', ...
    '-v7.3');


fprintf('\n');
fprintf('========================================\n');
fprintf('GROUND TRUTH SAVED\n');
fprintf('%s\n', GT_SAVE_FILE);
fprintf('========================================\n\n');

disp('Saved low_contrast_ground_truth_theta90.mat');
return

%% ========================================================================
% 12. Create initial scatterer map
% ========================================================================

fprintf('Creating initial scatterer distribution...\n');


sca_per_cell = ...
    50;


scatterers = ...
    create_scatterers( ...
        prb, ...
        scene_depth, ...
        pulse_duration, ...
        lambda, ...
        sca_per_cell);


%% ========================================================================
% 13. Compute time-varying scatterer positions
% ========================================================================

fprintf('Applying k-Wave motion to scatterers...\n');


sca_mesh = ...
    time_vary_sca_kwave( ...
        kgrid, ...
        sensor, ...
        motion, ...
        scatterers, ...
        TF);


%% ========================================================================
% 14. Define scatterer amplitudes
% ========================================================================

amp = ...
    define_phantom_scatterer_amplitudes( ...
        sca_mesh, ...
        shear_params);


%% ========================================================================
% 15. Plane-wave transmit sequence
% ========================================================================

alpha_max = ...
    deg2rad(15);


Na = ...
    1;


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
% 16. Ultrasound pulse setup
% ========================================================================

[impulse_response, excitation, lag] = ...
    pulse_setup( ...
        f0, ...
        pulse_duration, ...
        0.65, ...
        fs);


%% ========================================================================
% 17. Field II RF simulation
% ========================================================================

fprintf('\n');
fprintf('========================================\n');
fprintf('STARTING FIELD II RF SIMULATION\n');
fprintf('========================================\n\n');


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


fprintf('\n');
fprintf('========================================\n');
fprintf('FIELD II SIMULATION COMPLETE\n');
fprintf('RF data size: %s\n', mat2str(size(RF)));
fprintf('========================================\n\n');


%% ========================================================================
% 18. Add USTB
% ========================================================================

addpath( ...
    genpath( ...
        fullfile( ...
            getenv('HOME'), ...
            'USTB')));


%% ========================================================================
% 19. Create USTB probe
% ========================================================================

probe = ...
    uff.linear_array();


probe.element_height = ...
    prb.element_height;


probe.pitch = ...
    prb.pitch;


probe.element_width = ...
    prb.element_width;


probe.N = ...
    prb.N;


%% ========================================================================
% 20. Create USTB transmit sequence
% ========================================================================

seq(Na) = ...
    uff.wave();


for n = 1:Na

    seq(n) = ...
        uff.wave();

    seq(n).probe = ...
        probe;

    seq(n).source.azimuth = ...
        alpha(n);

    seq(n).source.distance = ...
        Inf;

    seq(n).sound_speed = ...
        c0;

    seq(n).delay = ...
        -lag*dt;

end


%% ========================================================================
% 21. Create ultrasound pulse
% ========================================================================

pulse = ...
    uff.pulse();


pulse.fractional_bandwidth = ...
    0.65;


pulse.center_frequency = ...
    f0;


%% ========================================================================
% 22. Create USTB channel data
% ========================================================================

channel_data = ...
    uff.channel_data();


channel_data.sampling_frequency = ...
    fs;


channel_data.sound_speed = ...
    c0;


channel_data.initial_time = ...
    0;


channel_data.pulse = ...
    pulse;


channel_data.probe = ...
    probe;


channel_data.sequence = ...
    seq;


channel_data.PRF = ...
    prf;


channel_data.data = ...
    RF;


%% ========================================================================
% 23. Define beamforming scan
% ========================================================================

sca = ...
    uff.linear_scan( ...
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


%% ========================================================================
% 24. Delay-and-sum beamforming
% ========================================================================

fprintf('\n');
fprintf('========================================\n');
fprintf('STARTING USTB BEAMFORMING\n');
fprintf('========================================\n\n');


F_number = ...
    1.2;


das = ...
    midprocess.das();


das.channel_data = ...
    channel_data;


das.scan = ...
    sca;


%% ------------------------------------------------------------------------
% Transmit apodization
% -------------------------------------------------------------------------

das.transmit_apodization.window = ...
    uff.window.tukey25;


das.transmit_apodization.f_number = ...
    F_number;


%% ------------------------------------------------------------------------
% Receive apodization
% -------------------------------------------------------------------------

das.receive_apodization.window = ...
    uff.window.tukey25;


das.receive_apodization.f_number = ...
    F_number;


%% ------------------------------------------------------------------------
% Run beamformer
% -------------------------------------------------------------------------

b_data = ...
    das.go();


fprintf('\nBeamforming complete.\n');

fprintf( ...
    'Beamformed data size: %s\n', ...
    mat2str(size(b_data.data)));


%% ------------------------------------------------------------------------
% Display B-mode image
% -------------------------------------------------------------------------

b_data.plot( ...
    [], ...
    [], ...
    40);


%% ========================================================================
% 25. Autocorrelation displacement estimation
% ========================================================================

fprintf('\n');
fprintf('========================================\n');
fprintf('STARTING AUTOCORRELATION ESTIMATION\n');
fprintf('========================================\n\n');


%% ------------------------------------------------------------------------
% Estimator parameters
% -------------------------------------------------------------------------

autocorr_params.z_gate = ...
    4;


autocorr_params.x_gate = ...
    2;


autocorr_params.packet_size = ...
    4;


%% ------------------------------------------------------------------------
% Create estimator
% -------------------------------------------------------------------------

displacement_estimator = ...
    postprocess.autocorrelation_displacement_estimation();


displacement_estimator.channel_data = ...
    channel_data;


displacement_estimator.input = ...
    b_data;


displacement_estimator.z_gate = ...
    autocorr_params.z_gate;


displacement_estimator.x_gate = ...
    autocorr_params.x_gate;


displacement_estimator.packet_size = ...
    autocorr_params.packet_size;


%% ------------------------------------------------------------------------
% Run displacement estimator
% -------------------------------------------------------------------------

displacement_estimation_autocorr = ...
    displacement_estimator.go();


%% ========================================================================
% 26. Extract autocorrelation quality indicator
% ========================================================================
%
% The modified autocorrelation_displacement_estimation class exposes the
% normalized autocorrelation quantity C for each displacement estimate.
%
% Candidate quality metric used later:
%
%                       Q_C = |C - 1|
%
% ========================================================================

quality_indicator_autocorr = ...
    displacement_estimator.quality_indicator;


%% ------------------------------------------------------------------------
% Display estimated displacement
% -------------------------------------------------------------------------

displacement_estimation_autocorr.plot( ...
    [], ...
    "Autocorrelation displacement", ...
    [], ...
    "none");


%% ========================================================================
% 27. Validate output dimensions
% ========================================================================

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


%% ------------------------------------------------------------------------
% Check number of temporal estimator frames
% -------------------------------------------------------------------------

expected_estimator_frames = ...
    size(disp_gt,3) - ...
    (autocorr_params.packet_size - 1) + 1;


actual_estimator_frames = ...
    size(displacement_estimation_autocorr.data,4);


fprintf( ...
    'Expected estimator frames: %d\n', ...
    expected_estimator_frames);


fprintf( ...
    'Actual estimator frames:   %d\n', ...
    actual_estimator_frames);


if actual_estimator_frames ~= expected_estimator_frames

    warning( ...
        ['Unexpected number of autocorrelation frames. ' ...
         'Expected %d but obtained %d.'], ...
        expected_estimator_frames, ...
        actual_estimator_frames);

end


%% ------------------------------------------------------------------------
% Check quality-indicator dimensions
% -------------------------------------------------------------------------

if ~isequal( ...
        size(displacement_estimation_autocorr.data), ...
        size(quality_indicator_autocorr))

    warning( ...
        ['The quality-indicator dimensions do not match the ' ...
         'displacement-estimation dimensions.']);

else

    fprintf( ...
        'Quality-indicator dimensions match displacement output.\n');

end


fprintf('========================================\n\n');


%% ========================================================================
% 28. Save complete low-contrast dataset
% ========================================================================

SAVE_FILE = ...
    'low_contrast_displacement_with_quality.mat';


fprintf('\n');
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
fprintf('========================================\n\n');


%% ========================================================================
% 29. Final summary
% ========================================================================

fprintf('SIMULATION SUMMARY\n');
fprintf('----------------------------------------\n');

fprintf( ...
    'Phantom:             %s\n', ...
    phantom_case);

fprintf( ...
    'Background SWS:      %.3f m/s\n', ...
    shear_params.c_shear_bkg);

fprintf( ...
    'Inclusion SWS:       %.3f m/s\n', ...
    shear_params.c_shear_incl);

fprintf( ...
    'Source frequency:    %.1f Hz\n', ...
    shear_params.source_freq);

fprintf( ...
    'Probe angle:         %.1f deg\n', ...
    rad2deg(prb_theta(1)));

fprintf( ...
    'Packet size:         %d\n', ...
    autocorr_params.packet_size);

fprintf( ...
    'z-gate:              %d\n', ...
    autocorr_params.z_gate);

fprintf( ...
    'x-gate:              %d\n', ...
    autocorr_params.x_gate);

fprintf( ...
    'GT file:             %s\n', ...
    GT_SAVE_FILE);

fprintf( ...
    'Complete data file:  %s\n', ...
    SAVE_FILE);

fprintf('----------------------------------------\n');

fprintf('\nSIM_PHANTOM_MAIN COMPLETE.\n');