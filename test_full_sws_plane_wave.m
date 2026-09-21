clear; clc; close all;

%% Setup

toolbox_dir = ...
    'C:\Users\louis\Documents\GitHub\MasterOppgave\elastography\ElastographyToolbox';

addpath(genpath(toolbox_dir));

%% Use same GT grid

load(fullfile(tempdir, 'sws_from_gt_disp.mat'), ...
    'x_gt_disp', 'z_gt_disp');

scan = uff.linear_scan( ...
    'x_axis', x_gt_disp, ...
    'z_axis', z_gt_disp);

%% Parameters

PRF = 3000;      % Hz
f = 200;         % Hz
c_true = 2.32;   % m/s
N_frames = 62;

omega = 2*pi*f;
k = omega/c_true;

[X, Z] = meshgrid(x_gt_disp, z_gt_disp);

t = (0:N_frames-1)/PRF;

%% Perfect travelling shear wave

disp_data = zeros( ...
    size(X,1), ...
    size(X,2), ...
    1, ...
    N_frames);

for n = 1:N_frames

    disp_data(:,:,1,n) = ...
        cos(k*Z - omega*t(n));

end

%% Elastography object

elasto = uff.elastography();

elasto.PRF = PRF;
elasto.target_frequency = f;
elasto.scan = scan;

elasto.autocorrelation_field = uff.data_field();

elasto.autocorrelation_field.set_component( ...
    'u', disp_data);

%% Frequency extraction

extractor = ...
    postprocess.elastography.extract_freq_component(elasto);

extractor.target_frequency = f;
extractor.component = 'u';

extractor.run();

%% Directional/spatial filtering

c0 = c_true;
wj = 2*pi*f;

elasto.radial_filter = ...
    uff.spatial_filter('gaussian_charite','radial');

elasto.radial_filter.sigma = c0/wj;

numOfAngles = 12;

elasto.directional_filtering_angles = ...
    2*pi*linspace(0,1-1/numOfAngles,numOfAngles);

elasto.angular_filter = ...
    uff.spatial_filter('gaussian_charite','angular');

elasto.angular_filter.deltaTheta = ...
    median(diff(elasto.directional_filtering_angles));

filterOp = ...
    postprocess.elastography.apply_filter(elasto);

filterOp.run();

Udir = elasto.directional_filtered_field_xz.u;

angles_deg = rad2deg(elasto.directional_filtering_angles);
energy_angle = nan(1,size(Udir,4));

for a = 1:size(Udir,4)

    slice = squeeze(Udir(:,:,1,a));
    energy_angle(a) = sum(abs(slice(:)).^4);

end

energy_relative = energy_angle / sum(energy_angle);

fprintf('\nDirection weights for ideal X-wave:\n');

for a = 1:length(energy_relative)

    fprintf('%6.1f deg: relative weight = %.3f\n', ...
        angles_deg(a), energy_relative(a));

end

%% Phase-gradient SWS

elastoPG = uff.elastography();
elastoPG.copy(elasto);

pgm = ...
    postprocess.elastography.phase_gradient_method(elastoPG);

pgm.run();

sws = squeeze(elastoPG.sws_map.u);

%% Ignore boundaries

sws_inner = sws(4:end-3,4:end-3);

fprintf('\nFull ideal plane-wave pipeline:\n');
fprintf('True SWS:   %.4f m/s\n', c_true);
fprintf('Mean SWS:   %.4f m/s\n', mean(sws_inner(:),'omitnan'));
fprintf('Median SWS: %.4f m/s\n', median(sws_inner(:),'omitnan'));
fprintf('Std SWS:    %.4f m/s\n', std(sws_inner(:),'omitnan'));