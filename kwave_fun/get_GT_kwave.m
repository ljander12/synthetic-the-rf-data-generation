function disp_gt = get_GT_kwave(kgrid, sensor, array_length, scene_depth, motion, prb_theta, TF)
%GET_GT_KWAVE Interpolate phantom ground-truth motion onto the image grid.
%
% The k-Wave displacement field is known on the sensor mask. This function
% reconstructs the displacement components on the full k-Wave grid,
% interpolates them onto the ultrasound image plane, and projects the
% result along the probe beam direction. The output is used as the
% ground-truth particle motion for downstream THE evaluation.
%
% Inputs:
%   kgrid        : kWaveGrid object used for the shear-wave simulation.
%   sensor       : k-Wave sensor structure defining sampled grid nodes.
%   array_length : Lateral width of the linear-probe image grid [m].
%   scene_depth  : Axial depth of the ultrasound image grid [m].
%   motion       : PRF-sampled motion structure with fields x, y and z.
%   prb_theta    : Probe rotation [yaw, roll] in radians.
%   TF           : 4-by-4 transform from probe coordinates to mesh
%                  coordinates.
%
% Outputs:
%   disp_gt : Ground-truth displacement projected onto the probe direction,
%             arranged as [depth samples, lateral samples, time frames].

prb_theta_x = prb_theta(1);
prb_theta_z = prb_theta(2);

%% Scene normal and probe direction without Symbolic Math Toolbox

N_scene = [ ...
    sin(prb_theta_x) * cos(prb_theta_z), ...
   -cos(prb_theta_x) * cos(prb_theta_z), ...
    sin(prb_theta_z)];

PrbDir = [ ...
   -tan(prb_theta_z) * sin(prb_theta_x), ...
    tan(prb_theta_z) * cos(prb_theta_x), ...
    1];

PrbDir = PrbDir / norm(PrbDir);

Nx = kgrid.Nx;
Ny = kgrid.Ny;
Nz = kgrid.Nz;

[X, Z] = meshgrid(-array_length/2:kgrid.dx:array_length/2, 5e-3:kgrid.dz:scene_depth);
gt_grid = TF*[X(:), zeros(size(X(:))), Z(:), ones(size(X(:)))]';
num_updates = size(motion.z, 2);
disp_gt = zeros([size(X), num_updates]);

for i = 1:num_updates
    disp_whole = zeros(Nx, Ny, Nz);
    disp_whole(logical(sensor.mask)) = motion.x(:,i);
    Dx_gt = interpn(kgrid.x, kgrid.y, kgrid.z, disp_whole, ...
        gt_grid(1,:), gt_grid(2,:), gt_grid(3,:), 'spline');

    disp_whole = zeros(Nx, Ny, Nz);
    disp_whole(logical(sensor.mask)) = motion.y(:,i);
    Dy_gt = interpn(kgrid.x, kgrid.y, kgrid.z, disp_whole, ...
        gt_grid(1,:), gt_grid(2,:), gt_grid(3,:), 'spline');

    disp_whole = zeros(Nx, Ny, Nz);
    disp_whole(logical(sensor.mask)) = motion.z(:,i);
    Dz_gt = interpn(kgrid.x, kgrid.y, kgrid.z, disp_whole, ...
        gt_grid(1,:), gt_grid(2,:), gt_grid(3,:), 'spline');

    Dx_gt(isnan(Dx_gt)) = 0;
    Dy_gt(isnan(Dy_gt)) = 0;
    Dz_gt(isnan(Dz_gt)) = 0;

    disp_gt(:,:,i) = reshape(Dx_gt.*-PrbDir(1) + ...
        Dy_gt.*-PrbDir(2) + Dz_gt.*-PrbDir(3), size(X));
end
