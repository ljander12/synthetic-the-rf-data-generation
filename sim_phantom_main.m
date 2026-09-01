%% ========================================================================
%
%   SIM_PHANTOM_MAIN
%
%   Generate synthetic ultrasound RF channel data for the phantom
%   time-harmonic elastography case study.
%
%   This script creates the phantom-mimicking shear-wave model, runs the
%   k-Wave elastic simulation, encodes the simulated particle motion into a
%   time-varying scatterer distribution, and uses Field II to simulate
%   pulse-echo RF channel data. Optional sections compute a ground-truth
%   motion field on the ultrasound image grid and beamform the RF data with
%   USTB for downstream processing examples.
%
%   Inputs:
%       No function inputs. Simulation, transducer, acquisition and output
%       parameters are set directly in this script.
%
%   Outputs:
%       RF channel data, PRF-sampled scatterer motion, ground-truth motion
%       and optional USTB beamformed data in the MATLAB workspace and saved
%       MAT-file.
%
%   Dependencies:
%       - k-Wave (http://www.k-wave.org/)
%       - Field II (https://field-ii.dk/)
%       - USTB (https://ultrasoundtoolbox.com/) for optional beamforming
%
% =========================================================================

clear; clc; close all;

% Add dependencies
addpath('/path/to/field_II/')

addpath("FIELD_fun")
addpath("kwave_fun")
addpath("utils")

%% ---------------- Parameters for shear-wave simulation ------------------

phantom_case = 'large_inclusions';

shear_params = get_phantom_params(phantom_case);
%% -------- Parameters for ultrasound pulse-echo imaging simulation -------

c0                      = 1540;                   % Speed of sound compression wave [m/s]
f0                      = 5.1333e+06;             % Transducer center frequency [Hz]
lambda                  = c0/f0;                  % Wavelength [m]
prb.element_height      = 5e-3;                   % Height of element [m]
prb.pitch               = 0.300e-3;               % Pitch [m]
prb.kerf                = 0.03e-03;               % Gap between elements [m]
prb.element_width       = prb.pitch-prb.kerf;     % Width of element [m]
prb.N                   = 128;                    % Number of elements
pulse_duration          = 2.5;                    % Pulse duration [cycles]
prf                     = 3000;                   % Pulse repetition frequency [Hz]

fs = 100e6;    % Sampling frequency [Hz]
dt = 1/fs;     % Sampling step [s]

attenfreq = 0.55;  % Attenuation [dB/MHz/cm]

array_length = (prb.N-1)*prb.pitch;
scene_depth = 50e-3;

%% Define k-Wave objects for phantom model

[kgrid, medium] = create_phantom(shear_params);
fprintf('Phantom case: %s\n', phantom_case);
fprintf('Background cs: %.2f m/s\n', shear_params.c_shear_bkg);
fprintf('Inclusion cs:  %.2f m/s\n', shear_params.c_shear_incl);

fprintf('Min cs: %.2f m/s\n', min(medium.sound_speed_shear(:)));
fprintf('Max cs: %.2f m/s\n', max(medium.sound_speed_shear(:)));

%% Define wave source k-Wave object

source = create_source(kgrid, shear_params.source_freq);

%% Define k-Wave sensor object

prb_center = [0 -20e-3 kgrid.z_vec(1)];
prb_theta = [deg2rad(45), deg2rad(0)];

[TF, TF_rev] = coor_transformation(prb_center, prb_theta);
sensor = create_sensor(kgrid, array_length, prb.element_height, scene_depth, TF);

%% Visualize shear-wave simulation setup

viewer = viewer3d(BackgroundColor="white", BackgroundGradient="off");
volshow(permute(medium.sound_speed_shear, [2, 1, 3]), ...
    Parent=viewer, Colormap=cool, Alphamap=[0.005, 0.1]);
volshow(permute(source.s_mask, [2, 1, 3]), ...
    Parent=viewer, Colormap=hot, Alphamap=[0.00, 0.2]);
volshow(permute(sensor.mask, [2, 1, 3]), ...
    Parent=viewer, Colormap=parula, Alphamap=[0.00, 0.1]);

%% Shear wave simulation

display_mask = source.s_mask;
DATA_CAST = 'single';
input_args_shear = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false, 'DataCast', DATA_CAST};
sensor_data_shear = pstdElastic3D(kgrid, medium, source, sensor, input_args_shear{:});

%% Sample motion at the ultrasound PRF

t_start = sqrt(shear_params.xrange^2+shear_params.yrange^2+shear_params.zrange^2)/shear_params.c_shear_bkg;
t_start = t_start + 50 * kgrid.dt;
[motion, vec_T] = sample_motion_at_prf("kwave", sensor_data_shear, kgrid, prf, t_start);

%% Create initial scatterer map

sca_per_cell = 50;
scatterers = create_scatterers(prb, scene_depth, pulse_duration, lambda, sca_per_cell);

%% Compute time-varying scatterer positions

sca_mesh = time_vary_sca_kwave(kgrid, sensor, motion, scatterers, TF);

%% Define scatterers intensity

amp = define_phantom_scatterer_amplitudes(sca_mesh, shear_params);

%% Plane wave transmit sequence

alpha_max = deg2rad(15);
Na = 1;

if Na == 1
    alpha = 0;
else
    alpha = linspace(-alpha_max, alpha_max, Na);
end

%% Ultrasound pulse setup

[impulse_response, excitation, lag] = pulse_setup(f0, pulse_duration, 0.65, fs);

%% Field II compute RF signals

RF = FIELD_calc_RF(prb, sca_mesh, amp, TF_rev, alpha, ...
    scene_depth, fs, attenfreq, f0, c0, impulse_response, excitation);

%% (Optional) Get simulation ground truth

disp_gt = get_GT_kwave(kgrid, sensor, array_length, scene_depth, motion, prb_theta, TF);

%% (Optional) Beamforming using USTB toolbox
addpath('/path/to/ustb')

probe = uff.linear_array();
probe.element_height    = prb.element_height;
probe.pitch             = prb.pitch;
probe.element_width     = prb.element_width;
probe.N                 = prb.N;

seq(Na) = uff.wave();
for n = 1:Na
    seq(n) = uff.wave();
    seq(n).probe = probe;
    seq(n).source.azimuth = alpha(n);
    seq(n).source.distance = Inf;
    seq(n).sound_speed = c0;
    seq(n).delay = -lag*dt;
end

pulse = uff.pulse();
pulse.fractional_bandwidth = 0.65;
pulse.center_frequency = f0;

channel_data = uff.channel_data();
channel_data.sampling_frequency = fs;
channel_data.sound_speed = c0;
channel_data.initial_time = 0;
channel_data.pulse = pulse;
channel_data.probe = probe;
channel_data.sequence = seq;
channel_data.PRF = prf;
channel_data.data = RF;

sca = uff.linear_scan('x_axis', linspace(-array_length/2, array_length/2, 512).', ...
    'z_axis', linspace(5e-3, scene_depth, 1024).');

F_number = 1.2;
das = midprocess.das();
das.channel_data = channel_data;
das.scan = sca;
das.transmit_apodization.window = uff.window.tukey25;
das.transmit_apodization.f_number = F_number;

das.receive_apodization.window = uff.window.tukey25;
das.receive_apodization.f_number = F_number;

b_data = das.go();

b_data.plot([], [], 40);

%% (Optional) Displacement estimation using USTB autocorrelation

autocorr_params.z_gate = 4;
autocorr_params.x_gate = 2;
autocorr_params.packet_size = 4;

displacement_estimator = postprocess.autocorrelation_displacement_estimation();
displacement_estimator.channel_data = channel_data;
displacement_estimator.input = b_data;
displacement_estimator.z_gate = autocorr_params.z_gate;
displacement_estimator.x_gate = autocorr_params.x_gate;
displacement_estimator.packet_size = autocorr_params.packet_size;

displacement_estimation_autocorr = displacement_estimator.go();
displacement_estimation_autocorr.plot([], "Autocorrelation displacement", [], "none");

%% (Optional) Save useful variables

SAVE_DIR = strcat('/path/to/data.mat');
save(SAVE_DIR, "RF", "channel_data", "disp_gt", "motion", "vec_T", ...
    "autocorr_params", "displacement_estimation_autocorr", ...
    "prb_center", "prb_theta", "scene_depth", '-v7.3')
