%% ========================================================================
%
%   SIM_LIVER_MAIN
%
%   Generate synthetic ultrasound RF channel data for the liver
%   time-harmonic elastography case study.
%
%   This script loads the tetrahedral liver FEM mesh and simulated harmonic
%   shear-wave motion from an HDF5 file, samples the motion at the
%   ultrasound PRF, encodes it into a time-varying scatterer distribution,
%   and uses Field II to simulate curvilinear-probe pulse-echo RF channel
%   data. Optional sections compute a ground-truth motion field on the
%   ultrasound image grid and beamform the RF data with USTB.
%
%   Inputs:
%       No function inputs. The liver mesh, material labels and FEM motion
%       are read from Liver-sim-data.h5; transducer and acquisition
%       parameters are set directly in this script.
%
%   Outputs:
%       RF channel data, PRF-sampled scatterer motion, ground-truth motion
%       and optional USTB beamformed data in the MATLAB workspace and saved
%       MAT-file.
%
%   Dependencies:
%       - Field II (https://field-ii.dk/)
%       - USTB (https://ultrasoundtoolbox.com/) for optional beamforming
%
% =========================================================================

clear; clc; close all;

% Add dependencies
addpath('C:\Users\louis\Documents\GitHub\MasterOppgave\synthetic-the-rf-data-generation')

addpath("FIELD_fun")
addpath("FEM_fun")
addpath("utils")

%% Read mesh and simulated shear wave displacement field

h5file = 'Liver-sim-data.h5';

[model, mesh, Uxyz, TimeVector] = load_liver_sim_data(h5file);
nodes = mesh.nodes;

%% -------- Parameters for ultrasound pulse-echo imaging simulation -------

c0                      = 1540;        % Speed of sound compression wave [m/s]
f0                      = 3570000;     % Transducer center frequency [Hz]
lambda                  = c0/f0;       % Wavelength [m]
prb.element_height      = 0.0135;      % Height of element [m]
prb.pitch               = 5.0800e-04;  % Pitch [m]
prb.kerf                = 4.8000e-05;  % Gap between elements [m]
prb.element_width       = 4.6000e-04;  % Width of element [m]
prb.radius              = 0.0496;      % Probe radius [m]
prb.N                   = 128;         % Number of elements
pulse_duration          = 2.5;         % Pulse duration [cycles]
prf                     = 3000;        % Pulse repetition frequency [Hz]

fs = 100e6;    % Sampling frequency [Hz]
dt = 1/fs;     % Sampling step [s]

attenfreq = 0;  % Old liver Field II simulation did not use attenuation.

scene_depth = 120e-3;

%% Sample motion at the ultrasound PRF

t_start = 0;
[motion, vec_T] = sample_motion_at_prf("fem", Uxyz, TimeVector, prf, t_start);

%% Create initial scatterer map

sca_per_cell = 10;
scatterers = create_scatterers(prb, scene_depth, pulse_duration, lambda, sca_per_cell);

%% Coordinate transfer matrix

prb_center = [-0.05 0.02 min(nodes(3,:))];
prb_theta = [deg2rad(-20), deg2rad(0)];

[TF, TF_rev] = coor_transformation(prb_center, prb_theta);

%% Compute time-varying scatterer positions

sca_mesh = time_vary_sca_fem(model, motion, vec_T, scatterers, TF);

%% Define scatterers intensity

[amp, liver_idx] = define_liver_scatterer_amplitudes(model, sca_mesh);

%% Visualize the scene of probe

plot_liver_probe_scene(model, sca_mesh, liver_idx);

%% Diverging wave transmit sequence

Na = 1;
if Na == 1
    focal = [0, 0, -prb.radius];
else
    focal = [linspace(-30e-3, 30e-3, Na)', zeros(Na, 1), ...
        -prb.radius*ones(Na, 1)];
end

%% Ultrasound pulse setup

[impulse_response, excitation, lag] = pulse_setup(f0, pulse_duration, 0.65, fs);

%% Field II compute RF signals

RF = FIELD_calc_RF(prb, sca_mesh, amp, TF_rev, focal, ...
    scene_depth, fs, attenfreq, f0, c0, impulse_response, excitation);

%% (Optional) Get simulation ground truth

disp_gt = get_GT_fem(model, motion, vec_T, prb, scene_depth, TF);

%% (Optional) Beamforming using USTB toolbox

addpath('/path/to/ustb')

probe = uff.curvilinear_array();
probe.element_height    = prb.element_height;
probe.pitch             = prb.pitch;
probe.element_width     = prb.element_width;
probe.radius            = prb.radius;
probe.N                 = prb.N;

seq(Na) = uff.wave();
for n = 1:Na
    seq(n) = uff.wave();
    seq(n).probe = probe;
    seq(n).source.xyz = focal(n,:);
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

scan = uff.sector_scan();
scan.azimuth_axis = linspace(-probe.maximum_angle, probe.maximum_angle, 512).';
scan.depth_axis = linspace(probe.radius, probe.radius+scene_depth, 768).';
scan.origin = uff.point('xyz', [0, 0, -probe.radius]);

mid = midprocess.das();
mid.dimension = dimension.both;
mid.channel_data = channel_data;
mid.scan = scan;

F_number = 1.7;
das.transmit_apodization.window = uff.window.tukey25;
das.transmit_apodization.f_number = F_number;
das.receive_apodization.window = uff.window.tukey25;
das.receive_apodization.f_number = F_number;

b_data = mid.go();

b_data.plot();

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
save(SAVE_DIR, "RF", "disp_gt", "channel_data", "motion", "vec_T", ...
    "autocorr_params", "displacement_estimation_autocorr", ...
    "prb_center", "prb_theta", "scene_depth", '-v7.3')
