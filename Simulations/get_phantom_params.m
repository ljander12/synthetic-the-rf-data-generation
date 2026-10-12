function shear_params = get_phantom_params(case_name)

%% =========================================================
% Common parameters based on Chaoran's synthetic example
% ==========================================================

shear_params.xrange = 100e-3;
shear_params.zrange = 95e-3;
shear_params.yrange = 120e-3;

shear_params.cx1 = 0;
shear_params.cy1 = -20e-3;
shear_params.cz1 = 15e-3;
shear_params.cr1 = 5e-3;

shear_params.cx2 = 0;
shear_params.cy2 = 20e-3;
shear_params.cz2 = 35e-3;
shear_params.cr2 = 10e-3;

shear_params.c_shear_bkg  = 2.32;
shear_params.c_shear_incl = 4.67;

shear_params.source_freq = 200;
shear_params.rho0 = 1079;

%% =========================================================
% Select phantom
% ==========================================================

switch lower(case_name)

    case 'chaoran_baseline'
        % Original synthetic phantom.
        % Nothing is changed.

    case 'homogeneous'
        % No mechanical contrast
        shear_params.c_shear_incl = ...
            shear_params.c_shear_bkg;

    case 'low_contrast'
        shear_params.c_shear_incl = 3.0;

    case 'high_contrast'
        shear_params.c_shear_incl = 6.0;

    case 'small_inclusions'
        shear_params.cr1 = 3e-3;
        shear_params.cr2 = 5e-3;

    case 'large_inclusions'
        shear_params.cr1 = 8e-3;
        shear_params.cr2 = 15e-3;

    otherwise
        error("Unknown phantom case: %s", case_name);

end

end