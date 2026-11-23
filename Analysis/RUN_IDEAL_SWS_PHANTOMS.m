%% RUN_IDEAL_SWS_PHANTOMS
% Generate exact, ideal 2-D SWS ground-truth maps with known inclusions.
% No k-Wave, Field II, USTB, beamforming, or external toolbox required.
% The maps are deliberately *not* reconstructed measurements.
% Authors' synthetic reference phantoms, not a copy of existing k-Wave GT.
%
% Run from the Analysis/ directory (or anywhere within this repository).

clear; clc; close all;

%% SETTINGS --------------------------------------------------------------
c_background = 2.32;    % shear wave speed [m/s]
rho = 1000;             % density [kg/m^3] for optional Young modulus
x_mm = -20:0.20:20;     % lateral axis [mm]
z_mm =   5:0.20:50;     % depth axis [mm]
show_stiffness = false; % true: also show ideal Young-modulus overview
color_limits_sws = [1.5 6.2];  % same color limits across all cases

% Prefer Results/ in the repository when this script is in Analysis/.
script_dir = fileparts(mfilename('fullpath'));
if strcmpi(get_last_folder(script_dir),'Analysis')
    repo_dir = fileparts(script_dir);
else
    repo_dir = script_dir;
end
out_dir = fullfile(repo_dir,'Results','artificial_sws_validation','ideal_phantoms');
if ~isfolder(out_dir), mkdir(out_dir); end

[X,Z] = meshgrid(x_mm,z_mm); % coordinates in mm
N = 8;
names = {'homogeneous', 'small_inclusion', 'large_inclusion', ...
    'high_contrast', 'low_contrast', 'soft_inclusion', ...
    'two_inclusions', 'elliptical_inclusion'};
titles = {'Homogeneous', 'Small inclusion', 'Large inclusion', ...
    'High contrast', 'Low contrast', 'Soft inclusion', ...
    'Two inclusions', 'Elliptical inclusion'};

% Each speed is known exactly at each sampled pixel; no random terms.
phantoms = repmat(struct('name','','sws_gt',[], ...
    'inclusion_mask',[],'x_axis',[],'z_axis',[], ...
    'description',''),1,N);

for k = 1:N
    sws_gt = c_background*ones(size(X));
    inclusion_mask = false(size(X));

    switch names{k}
        case 'homogeneous'
            desc = 'No inclusion; background 2.32 m/s';
        case 'small_inclusion'
            mask = circle_mask(X,Z,0,26,3);       % radius 3 mm
            sws_gt(mask) = 4.67;
            inclusion_mask = mask;
            desc = 'Radius 3 mm, inclusion 4.67 m/s';
        case 'large_inclusion'
            mask = circle_mask(X,Z,0,26,8);       % radius 8 mm
            sws_gt(mask) = 4.67;
            inclusion_mask = mask;
            desc = 'Radius 8 mm, inclusion 4.67 m/s';
        case 'high_contrast'
            mask = circle_mask(X,Z,0,26,4.5);
            sws_gt(mask) = 6.00;
            inclusion_mask = mask;
            desc = 'Radius 4.5 mm, inclusion 6.00 m/s';
        case 'low_contrast'
            mask = circle_mask(X,Z,0,26,4.5);
            sws_gt(mask) = 2.80; % illustrative value; not repository GT
            inclusion_mask = mask;
            desc = 'Radius 4.5 mm, inclusion 2.80 m/s (illustrative)';
        case 'soft_inclusion'
            mask = circle_mask(X,Z,0,26,5);
            sws_gt(mask) = 1.60;
            inclusion_mask = mask;
            desc = 'Radius 5 mm, soft inclusion 1.60 m/s';
        case 'two_inclusions'
            a = circle_mask(X,Z,-8,20,4);
            b = circle_mask(X,Z, 8,34,5);
            sws_gt(a) = 4.67;
            sws_gt(b) = 1.60;
            inclusion_mask = a | b;
            desc = 'Two circles: 4.67 m/s (stiff), 1.60 m/s (soft)';
        case 'elliptical_inclusion'
            mask = ((X-0)/4).^2 + ((Z-26)/10).^2 <= 1;
            sws_gt(mask) = 4.67;
            inclusion_mask = mask;
            desc = 'Ellipse: lateral radius 4 mm, axial radius 10 mm';
    end

    phantoms(k).name = names{k};
    phantoms(k).sws_gt = sws_gt;
    phantoms(k).inclusion_mask = inclusion_mask;
    phantoms(k).x_axis = x_mm(:)*1e-3; % meters, for pipeline reuse
    phantoms(k).z_axis = z_mm(:)*1e-3;
    phantoms(k).description = desc;

    % Preserve descriptive metadata alongside ground-truth arrays.
    case_info = struct('case_name',names{k},'description',desc, ...
        'background_speed_mps',c_background,'is_exact_gt',true, ...
        'source','Analytically defined 2-D geometry (not k-Wave simulation)');
    x_axis = phantoms(k).x_axis; %#ok<NASGU>
    z_axis = phantoms(k).z_axis; %#ok<NASGU>
    E_gt_kPa = 3*rho*sws_gt.^2/1000; %#ok<NASGU>
    save(fullfile(out_dir,[names{k} '_GT.mat']), ...
        'sws_gt','inclusion_mask','x_axis','z_axis','case_info', ...
        'E_gt_kPa','rho');
end

%% OVERVIEW: eight ideal SWS ground-truth maps ---------------------------
f = figure('Color','w','Position',[70 70 1600 780], ...
    'Name','Ideal SWS phantom ground truth');
t = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
for k=1:N
    ax = nexttile(t);
    imagesc(ax,x_mm,z_mm,phantoms(k).sws_gt);
    set(ax,'YDir','normal');
    axis(ax,'image');
    clim_compat(ax,color_limits_sws);
    colormap(ax,parula(256));
    title(ax,titles{k},'Interpreter','none');
    xlabel(ax,'x [mm]'); ylabel(ax,'z [mm]');
    cb=colorbar(ax);
    cb.Label.String='SWS [m/s]';
end
title(t,'Ideal SWS maps | Exact ground truth | No wave simulation');
saveas(f,fullfile(out_dir,'ideal_sws_phantoms_overview.png'));
savefig(f,fullfile(out_dir,'ideal_sws_phantoms_overview.fig'));

%% OPTIONAL: analogous stiffness maps (same patterns, different units) ---
if show_stiffness
    f2 = figure('Color','w','Position',[70 70 1600 780], ...
        'Name','Ideal Young modulus phantom ground truth');
    t2 = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
    for k=1:N
        ax = nexttile(t2);
        E_kPa = 3*rho*phantoms(k).sws_gt.^2/1000;
        imagesc(ax,x_mm,z_mm,E_kPa);
        set(ax,'YDir','normal'); axis(ax,'image');
        clim_compat(ax,3*rho*color_limits_sws.^2/1000);
        colormap(ax,parula(256));
        title(ax,titles{k},'Interpreter','none');
        xlabel(ax,'x [mm]'); ylabel(ax,'z [mm]');
        cb=colorbar(ax);
        cb.Label.String='Young''s modulus [kPa]';
    end
    title(t2,'Ideal Young''s modulus (E = 3 rho c_s^2, approximation)');
    saveas(f2,fullfile(out_dir,'ideal_youngs_modulus_overview.png'));
    savefig(f2,fullfile(out_dir,'ideal_youngs_modulus_overview.fig'));
end

save(fullfile(out_dir,'ALL_IDEAL_SWS_PHANTOMS.mat'), ...
    'phantoms','c_background','x_mm','z_mm','rho','-v7.3');

fprintf('\nCreated %d exact synthetic SWS ground-truth maps.\n',N);
fprintf('Overview: %s\n',fullfile(out_dir,'ideal_sws_phantoms_overview.png'));
fprintf('Individual MAT files and combined data saved to:\n  %s\n',out_dir);
fprintf('NOTE: Correlation Q and RMI-1D require displacement u(x,z,t).\n');
fprintf('SWS maps alone are valid for GT-based metrics, but cannot directly test Q/RMI.\n\n');

%% LOCAL HELPERS ---------------------------------------------------------
function m = circle_mask(X,Z,x0,z0,r_mm)
    m = (X-x0).^2+(Z-z0).^2 <= r_mm.^2;
end

function folder = get_last_folder(path_string)
    [~,folder] = fileparts(path_string);
end

function clim_compat(ax,lims)
    caxis(ax,lims);
end
