clear;
clc;
close all;

%% ========================================================================
% Q_C-SCREENED SWS RECONSTRUCTION - LOW CONTRAST
%
% Compare:
%
%   1) Original autocorrelation displacement
%          ->
%      Phase-gradient SWS
%
%   2) Q_C-screened autocorrelation displacement
%          ->
%      Temporal interpolation of rejected samples
%          ->
%      Phase-gradient SWS
%
% Candidate displacement reliability metric:
%
%                       Q_C = |C - 1|
%
% Preliminary threshold:
%
%                       Q_C <= 0.05
%
% IMPORTANT:
%   - Samples with non-finite Q_C are NOT counted as rejected.
%   - Only samples with finite Q_C > threshold are screened.
%   - Rejected displacement samples are interpolated temporally.
%   - No new k-Wave or Field II simulation is performed.
%
% Input:
%   low_contrast_displacement_with_quality.mat
%
% Output:
%   low_contrast_QC_screened_SWS.mat
%
% ========================================================================


%% ========================================================================
% 1. Load data
% ========================================================================

DATA_FILE = ...
    'low_contrast_displacement_with_quality.mat';

load(DATA_FILE);


fprintf('\n');
fprintf('====================================================\n');
fprintf('Q_C-SCREENED SWS RECONSTRUCTION - LOW CONTRAST\n');
fprintf('====================================================\n');


%% ========================================================================
% 2. Add elastography toolbox
% ========================================================================

toolbox_dir = ...
    ['C:\Users\louis\Documents\GitHub\MasterOppgave\' ...
     'elastography\ElastographyToolbox'];

addpath(genpath(toolbox_dir));


%% ========================================================================
% 3. Extract dimensions
% ========================================================================

Nz = ...
    length(b_data.scan.z_axis);

Nx = ...
    length(b_data.scan.x_axis);

Nt = ...
    size( ...
        displacement_estimation_autocorr.data, ...
        4);


fprintf('Beamformed scan:       %d x %d\n', ...
    Nz, ...
    Nx);

fprintf('Temporal frames:       %d\n', ...
    Nt);


%% ========================================================================
% 4. Extract autocorrelation displacement
% ========================================================================

u_original = ...
    reshape( ...
        double( ...
            displacement_estimation_autocorr.data(:,1,1,:)), ...
        Nz, ...
        Nx, ...
        Nt);


%% ------------------------------------------------------------------------
% Flip sign
%
% The autocorrelation estimator uses the opposite displacement sign
% convention relative to the k-Wave ground truth.
%
% The global sign does not change the SWS magnitude, but keeping the
% convention consistent simplifies later comparisons.
% -------------------------------------------------------------------------

u_original = ...
    -u_original;


%% ========================================================================
% 5. Extract normalized autocorrelation quantity C
% ========================================================================

C = ...
    reshape( ...
        double( ...
            quality_indicator_autocorr(:,1,1,:)), ...
        Nz, ...
        Nx, ...
        Nt);


%% ========================================================================
% 6. Candidate quality metric
% ========================================================================

Q_C = ...
    abs(C - 1);


fprintf('\n');

fprintf('Displacement size:     %s\n', ...
    mat2str(size(u_original)));

fprintf('Q_C size:              %s\n', ...
    mat2str(size(Q_C)));


%% ========================================================================
% 7. Define quality threshold
% ========================================================================

Q_threshold = ...
    0.05;


%% ========================================================================
% 8. Classify Q_C samples
% ========================================================================
%
% IMPORTANT:
%
% valid_Q:
%     Q_C exists and can therefore be evaluated.
%
% accepted:
%     finite Q_C <= threshold.
%
% rejected:
%     finite Q_C > threshold.
%
% Samples with NaN/Inf Q_C are therefore NOT automatically classified as
% poor-quality samples.
% ========================================================================

valid_Q = ...
    isfinite(Q_C);

accepted = ...
    valid_Q & ...
    Q_C <= Q_threshold;

rejected = ...
    valid_Q & ...
    Q_C > Q_threshold;


%% ========================================================================
% 9. Quality statistics
% ========================================================================

N_total = ...
    numel(Q_C);

N_valid_Q = ...
    nnz(valid_Q);

N_accepted = ...
    nnz(accepted);

N_rejected = ...
    nnz(rejected);


valid_Q_percent = ...
    100 * N_valid_Q / N_total;

retained_percent = ...
    100 * N_accepted / N_valid_Q;

rejected_percent = ...
    100 * N_rejected / N_valid_Q;


fprintf('\n');
fprintf('----------------------------------------------------\n');
fprintf('Q_C SCREENING STATISTICS\n');
fprintf('----------------------------------------------------\n');

fprintf('Q_C threshold:         %.3f\n', ...
    Q_threshold);

fprintf('Total samples:         %d\n', ...
    N_total);

fprintf('Finite Q_C samples:    %d (%.2f %%)\n', ...
    N_valid_Q, ...
    valid_Q_percent);

fprintf('Accepted samples:      %d (%.2f %% of finite Q_C)\n', ...
    N_accepted, ...
    retained_percent);

fprintf('Rejected samples:      %d (%.2f %% of finite Q_C)\n', ...
    N_rejected, ...
    rejected_percent);

fprintf('----------------------------------------------------\n');


%% ========================================================================
% 10. Create Q_C-screened displacement field
% ========================================================================
%
% Only samples explicitly classified as rejected are removed.
%
% Samples where Q_C is undefined are left unchanged.
% ========================================================================

u_screened = ...
    u_original;

u_screened(rejected) = ...
    NaN;


%% ========================================================================
% 11. Temporally interpolate rejected samples
% ========================================================================
%
% Zero-filling should not be used here because abrupt zeros can introduce
% artificial temporal frequency components.
%
% Linear temporal interpolation is therefore used instead.
% ========================================================================

fprintf('\nInterpolating rejected displacement samples...\n');


u_screened = ...
    fillmissing( ...
        u_screened, ...
        'linear', ...
        3, ...
        'EndValues', ...
        'nearest');


%% ------------------------------------------------------------------------
% Safety fallback
%
% If a complete temporal trace contains no finite values, fillmissing may
% not be able to reconstruct it. In that case restore the corresponding
% original displacement estimate.
% -------------------------------------------------------------------------

remaining_nan = ...
    ~isfinite(u_screened);

N_remaining_nan = ...
    nnz(remaining_nan);


if N_remaining_nan > 0

    u_screened(remaining_nan) = ...
        u_original(remaining_nan);

end


fprintf( ...
    'Remaining NaNs after interpolation: %d\n', ...
    nnz(~isfinite(u_screened)));


%% ========================================================================
% 12. Quantify how much displacement data changed
% ========================================================================

changed_mask = ...
    rejected & ...
    isfinite(u_original) & ...
    isfinite(u_screened);


if any(changed_mask(:))

    displacement_change = ...
        u_screened(changed_mask) - ...
        u_original(changed_mask);

    displacement_change_RMSE_um = ...
        sqrt( ...
            mean( ...
                displacement_change.^2)) ...
        * 1e6;

    displacement_change_MAE_um = ...
        mean( ...
            abs(displacement_change)) ...
        * 1e6;

else

    displacement_change_RMSE_um = ...
        NaN;

    displacement_change_MAE_um = ...
        NaN;

end


fprintf('\n');

fprintf( ...
    'RMSE change in screened displacement: %.4f um\n', ...
    displacement_change_RMSE_um);

fprintf( ...
    'MAE change in screened displacement:  %.4f um\n', ...
    displacement_change_MAE_um);


%% ========================================================================
% 13. Run phase-gradient SWS on ORIGINAL displacement
% ========================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('RUNNING ORIGINAL SWS RECONSTRUCTION\n');
fprintf('====================================================\n');


elasto_original = ...
    run_phase_gradient_pipeline( ...
        u_original, ...
        b_data.scan, ...
        prf, ...
        shear_params.source_freq);


sws_original = ...
    double( ...
        elasto_original.sws_map.u);


fprintf( ...
    'Original SWS reconstruction complete.\n');


%% ========================================================================
% 14. Run phase-gradient SWS on Q_C-SCREENED displacement
% ========================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('RUNNING Q_C-SCREENED SWS RECONSTRUCTION\n');
fprintf('====================================================\n');


elasto_screened = ...
    run_phase_gradient_pipeline( ...
        u_screened, ...
        b_data.scan, ...
        prf, ...
        shear_params.source_freq);


sws_screened = ...
    double( ...
        elasto_screened.sws_map.u);


fprintf( ...
    'Q_C-screened SWS reconstruction complete.\n');


%% ========================================================================
% 15. Verify SWS dimensions
% ========================================================================

fprintf('\n');

fprintf('Original SWS size:     %s\n', ...
    mat2str(size(sws_original)));

fprintf('Screened SWS size:     %s\n', ...
    mat2str(size(sws_screened)));


if ~isequal( ...
        size(sws_original), ...
        size(sws_screened))

    error( ...
        ['Original and screened SWS maps have ' ...
         'different dimensions.']);

end


%% ========================================================================
% 16. Common valid SWS mask
% ========================================================================
%
% The <10 m/s criterion is only used for descriptive comparison and avoids
% extreme numerical phase-gradient values dominating the summary.
%
% The SWS maps themselves are NOT modified or clipped.
% ========================================================================

valid_SWS = ...
    isfinite(sws_original) & ...
    isfinite(sws_screened) & ...
    sws_original > 0 & ...
    sws_screened > 0 & ...
    sws_original < 10 & ...
    sws_screened < 10;


N_valid_SWS = ...
    nnz(valid_SWS);


fprintf( ...
    'Common valid SWS pixels: %d\n', ...
    N_valid_SWS);


%% ========================================================================
% 17. Overall SWS statistics
% ========================================================================

orig_vals = ...
    sws_original(valid_SWS);

screen_vals = ...
    sws_screened(valid_SWS);


%% ------------------------------------------------------------------------
% Original
% -------------------------------------------------------------------------

orig_mean = ...
    mean(orig_vals);

orig_median = ...
    median(orig_vals);

orig_std = ...
    std(orig_vals);


%% ------------------------------------------------------------------------
% Q_C-screened
% -------------------------------------------------------------------------

screen_mean = ...
    mean(screen_vals);

screen_median = ...
    median(screen_vals);

screen_std = ...
    std(screen_vals);


%% ========================================================================
% 18. Difference between SWS maps
% ========================================================================

map_difference = ...
    sws_screened - ...
    sws_original;


difference_vals = ...
    map_difference(valid_SWS);


difference_RMSE = ...
    sqrt( ...
        mean( ...
            difference_vals.^2));


difference_MAE = ...
    mean( ...
        abs(difference_vals));


difference_bias = ...
    mean( ...
        difference_vals);


difference_max_abs = ...
    max( ...
        abs(difference_vals));


%% ========================================================================
% 19. Correlation between original and screened SWS
% ========================================================================

sws_map_correlation = ...
    corr( ...
        orig_vals, ...
        screen_vals, ...
        'Rows', ...
        'complete');


%% ========================================================================
% 20. Print SWS comparison
% ========================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('SWS MAP COMPARISON\n');
fprintf('====================================================\n');


fprintf('\nORIGINAL\n');

fprintf('Mean SWS:              %.4f m/s\n', ...
    orig_mean);

fprintf('Median SWS:            %.4f m/s\n', ...
    orig_median);

fprintf('Std SWS:               %.4f m/s\n', ...
    orig_std);


fprintf('\nQ_C SCREENED\n');

fprintf('Mean SWS:              %.4f m/s\n', ...
    screen_mean);

fprintf('Median SWS:            %.4f m/s\n', ...
    screen_median);

fprintf('Std SWS:               %.4f m/s\n', ...
    screen_std);


fprintf('\nMAP DIFFERENCE\n');

fprintf('RMSE difference:       %.4f m/s\n', ...
    difference_RMSE);

fprintf('MAE difference:        %.4f m/s\n', ...
    difference_MAE);

fprintf('Mean difference:       %.4f m/s\n', ...
    difference_bias);

fprintf('Maximum abs. diff.:    %.4f m/s\n', ...
    difference_max_abs);

fprintf('Map correlation:       %.6f\n', ...
    sws_map_correlation);


fprintf('\nQUALITY SCREENING\n');

fprintf('Finite Q_C:            %.2f %%\n', ...
    valid_Q_percent);

fprintf('Q_C retention:         %.2f %%\n', ...
    retained_percent);

fprintf('Q_C rejection:         %.2f %%\n', ...
    rejected_percent);


fprintf('====================================================\n\n');


%% ========================================================================
% 21. Scan coordinates
% ========================================================================

x_mm = ...
    double( ...
        b_data.scan.x_axis(:)) ...
    * 1e3;

z_mm = ...
    double( ...
        b_data.scan.z_axis(:)) ...
    * 1e3;


%% ========================================================================
% 22. Plot original SWS map
% ========================================================================

figure;

imagesc( ...
    x_mm, ...
    z_mm, ...
    sws_original);

axis image;

set( ...
    gca, ...
    'YDir', ...
    'normal');

xlabel( ...
    'Lateral position [mm]');

ylabel( ...
    'Depth [mm]');

title( ...
    'Original phase-gradient SWS');

colorbar;

clim( ...
    [1.5 4.0]);


%% ========================================================================
% 23. Plot Q_C-screened SWS map
% ========================================================================

figure;

imagesc( ...
    x_mm, ...
    z_mm, ...
    sws_screened);

axis image;

set( ...
    gca, ...
    'YDir', ...
    'normal');

xlabel( ...
    'Lateral position [mm]');

ylabel( ...
    'Depth [mm]');

title( ...
    'Q_C-screened phase-gradient SWS');

colorbar;

clim( ...
    [1.5 4.0]);


%% ========================================================================
% 24. Plot SWS difference map
% ========================================================================

figure;

imagesc( ...
    x_mm, ...
    z_mm, ...
    map_difference);

axis image;

set( ...
    gca, ...
    'YDir', ...
    'normal');

xlabel( ...
    'Lateral position [mm]');

ylabel( ...
    'Depth [mm]');

title( ...
    'Q_C-screened minus original SWS');

colorbar;


%% ========================================================================
% 25. Rejection fraction per spatial pixel
% ========================================================================
%
% IMPORTANT:
%
% Rejection fraction is calculated only over temporal samples where Q_C is
% finite.
%
% For each spatial location:
%
%   rejection fraction =
%
%       number of finite Q_C samples above threshold
%       ---------------------------------------------
%       total number of finite Q_C samples
%
% ========================================================================

N_valid_temporal = ...
    sum( ...
        valid_Q, ...
        3);

N_rejected_temporal = ...
    sum( ...
        rejected, ...
        3);


rejection_fraction = ...
    nan( ...
        Nz, ...
        Nx);


spatial_Q_available = ...
    N_valid_temporal > 0;


rejection_fraction(spatial_Q_available) = ...
    N_rejected_temporal(spatial_Q_available) ./ ...
    N_valid_temporal(spatial_Q_available);


%% ========================================================================
% 26. Plot rejection fraction
% ========================================================================

figure;

imagesc( ...
    x_mm, ...
    z_mm, ...
    rejection_fraction);

axis image;

set( ...
    gca, ...
    'YDir', ...
    'normal');

xlabel( ...
    'Lateral position [mm]');

ylabel( ...
    'Depth [mm]');

title( ...
    'Fraction of finite Q_C samples rejected');

colorbar;

clim( ...
    [0 1]);


%% ========================================================================
% 27. Median Q_C map
% ========================================================================

median_Q_map = ...
    median( ...
        Q_C, ...
        3, ...
        'omitnan');


figure;

imagesc( ...
    x_mm, ...
    z_mm, ...
    median_Q_map);

axis image;

set( ...
    gca, ...
    'YDir', ...
    'normal');

xlabel( ...
    'Lateral position [mm]');

ylabel( ...
    'Depth [mm]');

title( ...
    'Temporal median Q_C');

colorbar;


%% ========================================================================
% 28. Summary structure
% ========================================================================

results = ...
    struct();


%% ------------------------------------------------------------------------
% Quality screening
% -------------------------------------------------------------------------

results.Q_threshold = ...
    Q_threshold;

results.N_total_Q = ...
    N_total;

results.N_valid_Q = ...
    N_valid_Q;

results.N_accepted = ...
    N_accepted;

results.N_rejected = ...
    N_rejected;

results.valid_Q_percent = ...
    valid_Q_percent;

results.retained_percent = ...
    retained_percent;

results.rejected_percent = ...
    rejected_percent;


%% ------------------------------------------------------------------------
% Displacement modification
% -------------------------------------------------------------------------

results.displacement_change_RMSE_um = ...
    displacement_change_RMSE_um;

results.displacement_change_MAE_um = ...
    displacement_change_MAE_um;


%% ------------------------------------------------------------------------
% Original SWS
% -------------------------------------------------------------------------

results.original_mean_SWS = ...
    orig_mean;

results.original_median_SWS = ...
    orig_median;

results.original_std_SWS = ...
    orig_std;


%% ------------------------------------------------------------------------
% Screened SWS
% -------------------------------------------------------------------------

results.screened_mean_SWS = ...
    screen_mean;

results.screened_median_SWS = ...
    screen_median;

results.screened_std_SWS = ...
    screen_std;


%% ------------------------------------------------------------------------
% SWS map comparison
% -------------------------------------------------------------------------

results.SWS_difference_RMSE = ...
    difference_RMSE;

results.SWS_difference_MAE = ...
    difference_MAE;

results.SWS_difference_bias = ...
    difference_bias;

results.SWS_difference_max_abs = ...
    difference_max_abs;

results.SWS_map_correlation = ...
    sws_map_correlation;


%% ========================================================================
% 29. Save results
% ========================================================================

SAVE_FILE = ...
    'low_contrast_QC_screened_SWS.mat';


save( ...
    SAVE_FILE, ...
    'sws_original', ...
    'sws_screened', ...
    'map_difference', ...
    'u_original', ...
    'u_screened', ...
    'C', ...
    'Q_C', ...
    'valid_Q', ...
    'accepted', ...
    'rejected', ...
    'rejection_fraction', ...
    'median_Q_map', ...
    'retained_percent', ...
    'rejected_percent', ...
    'Q_threshold', ...
    'results', ...
    'elasto_original', ...
    'elasto_screened', ...
    '-v7.3');


fprintf('\n');
fprintf('====================================================\n');
fprintf('Q_C-SCREENED SWS ANALYSIS COMPLETE\n');
fprintf('====================================================\n');

fprintf('Saved:\n');
fprintf('%s\n', ...
    SAVE_FILE);

fprintf('====================================================\n');


%% ========================================================================
% LOCAL FUNCTION
% ========================================================================

function elasto = ...
    run_phase_gradient_pipeline( ...
        displacement, ...
        scan, ...
        PRF, ...
        target_frequency)


    %% ====================================================================
    % Create elastography object
    % =====================================================================

    elasto = ...
        uff.elastography();


    elasto.PRF = ...
        PRF;


    elasto.target_frequency = ...
        target_frequency;


    elasto.scan = ...
        scan;


    %% ====================================================================
    % Insert autocorrelation displacement field
    % =====================================================================

    elasto.autocorrelation_field = ...
        uff.data_field();


    elasto.autocorrelation_field.set_component( ...
        'u', ...
        reshape( ...
            displacement, ...
            size(displacement,1), ...
            size(displacement,2), ...
            1, ...
            size(displacement,3)));


    %% ====================================================================
    % Directional filtering angles
    % =====================================================================

    num_angles = ...
        12;


    elasto.directional_filtering_angles = ...
        2*pi*linspace( ...
            0, ...
            1-1/num_angles, ...
            num_angles);


    %% ====================================================================
    % Extract target-frequency component
    % =====================================================================

    extractor = ...
        postprocess.elastography.extract_freq_component( ...
            elasto);


    extractor.target_frequency = ...
        target_frequency;


    extractor.component = ...
        'u';


    extractor.run();


    %% ====================================================================
    % Spatial filtering parameters
    % =====================================================================

    filter_speed = ...
        3.2;


    omega = ...
        2*pi*target_frequency;


    %% --------------------------------------------------------------------
    % Radial filter
    % ---------------------------------------------------------------------

    elasto.radial_filter = ...
        uff.spatial_filter( ...
            'gaussian_charite', ...
            'radial');


    elasto.radial_filter.sigma = ...
        filter_speed / omega;


    %% --------------------------------------------------------------------
    % Angular filter
    % ---------------------------------------------------------------------

    elasto.angular_filter = ...
        uff.spatial_filter( ...
            'gaussian_charite', ...
            'angular');


    elasto.angular_filter.deltaTheta = ...
        median( ...
            diff( ...
                elasto.directional_filtering_angles));


    %% ====================================================================
    % Apply filtering
    % =====================================================================

    filter_operator = ...
        postprocess.elastography.apply_filter( ...
            elasto);


    filter_operator.run();


    %% ====================================================================
    % Phase-gradient SWS estimation
    % =====================================================================

    pgm = ...
        postprocess.elastography.phase_gradient_method( ...
            elasto);


    pgm.run();

end