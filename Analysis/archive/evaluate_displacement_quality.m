clear; clc; close all;

%% Setup

repo_dir = fileparts(mfilename("fullpath"));
cd(repo_dir);
addpath(fullfile(getenv('HOME'),'USTB'));

%% Load data

data_file = fullfile(repo_dir, 'large_inclusions_displacement.mat');

load(data_file);

%% Prepare ground-truth displacement

% Add initial zero-displacement frame
gt_abs = cat(3, zeros(size(disp_gt(:,:,1))), disp_gt);

% Ground- truth frame-to-frame displacement
gt_delta = diff(gt_abs, 1, 3);

% Match autocorrelation packet size
packet_size = 4;
N_frames = size(gt_delta,3) - packet_size + 2;

gt_packet = zeros(size(gt_delta,1),size(gt_delta,2),N_frames);

for k = 1:N_frames
    gt_packet(:,:,k) = mean(gt_delta(:,:,k:k+packet_size-2), 3);
end    

%% Reshape estimated displacement 

disp_est = reshape(displacement_estimation_autocorr.data, 1024, 512, []);

%% Sanity checks

fprintf('gt_abs size: %s\n', mat2str(size(gt_abs)));
fprintf('gt_delta size: %s\n', mat2str(size(gt_delta)));
fprintf('disp_gt size: %s\n', mat2str(size(disp_gt)));

fprintf('Ground truth size: %s\n', mat2str(size(gt_packet)));
fprintf('Estimated displacement size: %s\n', mat2str(size(disp_est)));

%% Spatial alignment

% USTB beamforming grid
x_est = b_data.scan.x_axis(:).'; %1 x 512
z_est = b_data.scan.z_axis(:).'; %1024 x 1

% Ground truth grid spacing used in create_phantom / get__GT_kwave
d_gt = 1 / shear_params.source_freq / 5;

% Physical width of the probe
array_length = (channel_data.probe.N - 1) * channel_data.probe.pitch;

%Reconstruct the spatial grid used for disp_gt
x_gt = -array_length/2 + (0:size(gt_packet,2)-1) * d_gt;
z_gt = 5e-3 + (0:size(gt_packet,1)-1) * d_gt;

X_gt = repmat(x_gt, size(gt_packet,1), 1);
Z_gt = repmat(z_gt(:),1, size(gt_packet,2));

% Allocate aligned displacement estimate 
disp_est_gtgrid = zeros(size(gt_packet));

%Interpolate each estimated displacement frame onto the GT grid

for k = 1:size(disp_est, 3)

    disp_est_gtgrid(:,:,k) = interp2(x_est,z_est,double(disp_est(:,:,k)), X_gt,Z_gt,'linear',NaN);
end

%% Spatial alignment sanity checks
fprintf('Ground truth size: %s\n', mat2str(size(gt_packet)));
fprintf('Aligned estimate size: %s\n', mat2str(size(disp_est_gtgrid)));
fprintf('NaN values after interpolation: %s\n', sum(isnan(disp_est_gtgrid(:))));


%% Visual comparison

frame = 30;

gt_frame = gt_packet(:,:,frame);
est_frame = disp_est_gtgrid(:,:,frame);
err_frame = est_frame - gt_frame;

% Color scale for GT and estimate
max_disp = max(abs([gt_frame(:); est_frame(:)]));

figure;

subplot(1,3,1)
imagesc(x_gt*1e3, z_gt*1e3, gt_frame*1e6)
axis image
set(gca, 'YDir', 'normal')
caxis([-max_disp max_disp]*1e6)
colorbar
title('Ground truth')
xlabel('x [mm]')
ylabel('z [mm]')

subplot(1,3,2)
imagesc(x_gt*1e3, z_gt*1e3, est_frame*1e6)
axis image
set(gca, 'YDir', 'normal')
caxis([-max_disp max_disp]*1e6)
colorbar
title('Estimated displacement')
xlabel('x [mm]')
ylabel('z [mm]')

subplot(1,3,3)
imagesc(x_gt*1e3, z_gt*1e3, err_frame*1e6)
axis image
set(gca, 'YDir', 'normal')
colorbar
title('Error')
xlabel('x [mm]')
ylabel('z [mm]')


%% Quaøity metrics

error = disp_est_gtgrid - gt_packet;

RMSE = sqrt(mean(error(:).*error(:)));
MAE = mean(abs(error(:)));
Bias = mean(error(:));

disp('Displacement quality metrics:');
fprintf('RMSE: %.4f um\n', RMSE*1e6);
fprintf('MAE: %.4f um\n', MAE*1e6);
fprintf('Bias: %.4f um\n', Bias*1e6);
