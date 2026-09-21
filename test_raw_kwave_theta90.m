clear; clc;

%% Load 0-degree simulation

load('homogeneous_ground_truth_theta90.mat');

%% Recreate k-Wave grid and sensor geometry
% This does NOT run a new k-Wave simulation

[kgrid, ~] = create_phantom(shear_params);

sensor = create_sensor( ...
    kgrid, ...
    array_length, ...
    5e-3, ...          % probe element height
    scene_depth, ...
    TF);

%% Check correspondence

mask = sensor.mask ~= 0;

fprintf('Sensor points: %d\n', sum(mask(:)));
fprintf('Motion rows:   %d\n', size(motion.z,1));

%% Sensor coordinates in probe coordinate system

x_s = kgrid.x(mask);
y_s = kgrid.y(mask);
z_s = kgrid.z(mask);

P_probe = TF \ [ ...
    x_s(:).'; ...
    y_s(:).'; ...
    z_s(:).'; ...
    ones(1,length(x_s))];

xp = P_probe(1,:).';
yp = P_probe(2,:).';
zp = P_probe(3,:).';

%% Extract 200 Hz directly from raw motion

f = shear_params.source_freq;
t = vec_T(:).';

Uraw = sum( ...
    double(motion.z) .* ...
    exp(-1i*2*pi*f*t), ...
    2);

%% Select horizontal line around z = 25 mm

z_target = 25e-3;

[~, iz] = min(abs(zp-z_target));
z_sel = zp(iz);

tol_y = kgrid.dy/2 + 1e-9;
tol_z = kgrid.dz/2 + 1e-9;

line_mask = ...
    abs(yp) <= tol_y & ...
    abs(zp-z_sel) <= tol_z & ...
    abs(xp) < 15e-3;

x_line = xp(line_mask);
U_line = Uraw(line_mask);

[x_line, order] = sort(x_line);
U_line = U_line(order);

phase_line = unwrap(angle(U_line));

%% Phase slope

p = polyfit(x_line, phase_line, 1);

k_raw = abs(p(1));

omega = 2*pi*f;
c_raw = omega/k_raw;

k_true = omega/shear_params.c_shear_bkg;

fprintf('\nRaw theta0 phase-slope test:\n');
fprintf('Probe angle: %.2f deg\n', rad2deg(prb_theta(1)));
fprintf('Number of points: %d\n', length(x_line));
fprintf('Measured k: %.2f rad/m\n', k_raw);
fprintf('Expected k: %.2f rad/m\n', k_true);
fprintf('SWS from raw motion: %.3f m/s\n', c_raw);
fprintf('Ground truth SWS:    %.3f m/s\n', shear_params.c_shear_bkg);