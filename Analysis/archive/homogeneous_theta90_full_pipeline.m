clear;
clc;
close all;

%% Load completed theta90 k-Wave simulation

load('homogeneous_ground_truth_theta90.mat');

fprintf('Loaded theta90 simulation\n');
fprintf('Probe angle: %.1f deg\n', rad2deg(prb_theta(1)));
fprintf('Ground-truth SWS: %.3f m/s\n', shear_params.c_shear_bkg);
fprintf('Motion frames: %d\n', length(vec_T));

%% Recreate k-Wave grid and sensor geometry

[kgrid, ~] = create_phantom(shear_params);

sensor = create_sensor( ...
    kgrid, ...
    array_length, ...
    5e-3, ...
    scene_depth, ...
    TF);