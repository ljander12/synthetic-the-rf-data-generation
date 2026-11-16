INTEGRATED ANALYZE_SWS_QUALITY (PG baseline + two published Q metrics)
==========================================================================

What changed
------------
The ORIGINAL analyze_SWS_quality(config) from your pasted file is preserved.
Only a new optional SECTION 16 was inserted and one local helper group was
appended below its existing local functions. This gives ONE entry point.

Existing baseline analysis: PG SWS, GT, ROI, bias/RMSE/MAE, CNR/gCNR,
contrast recovery, Dong edge metrics, optional repeated noise experiment.

New optional analysis, enabled by config.run_published_quality_metrics:
  - Correlation Q (Schoen et al., 2025, adaptation)
  - RMI-1D (Choi et al., 2015 Eq. 6 adaptation)
  - PG and AC global/regional AUC, false-confidence checks
  - PG and AC original Young modulus / shear modulus, and two masks each:
        original | Correlation Q masked | RMI-1D masked

EXCLUDED: Phase Planarity Q and experimental Log-RMI.
The original mathematical implementations are not rewritten.

Installation
------------
1. Make a backup of your repository Analysis/analyze_SWS_quality.m.
2. Replace it with the supplied analyze_SWS_quality.m (same name).
3. Keep these original source functions in Analysis/ (only needed if
   the corresponding saved Q or RMI map does not yet exist):
      test_sws_measurement_quality.m
      validate_sws_rmi.m
      test_sws_phase_consistency.m [historical RMI dependency only]
4. Run the SAME main controller that already calls analyze_SWS_quality(config).
   Directly before that call, add:

      config.run_published_quality_metrics = true;
      config.quality_mode = 'saved'; % or 'auto' for missing cached Q/RMI
      config.quality_options = struct( ...
          'rho',1000,'modulus','young', ...
          'corrThreshold',0.63,'rmiThreshold',0.85, ...
          'showFigures',true,'dpi',220);

IMPORTANT: the original analyze_SWS_quality requires a config with the
existing fields used by your controller: phantom_case, data_file,
results_dir, sws_filter_speed, sws_roi_margin, sws_image_margin,
run_sws_noise_test, etc. This isn't a new standalone Run script.

Default behavior
----------------
If you do not set run_published_quality_metrics, the ORIGINAL function
behaves as before. It does not unexpectedly compute or overwrite Q/RMI.
When enabled, 'saved' is the default to avoid recalculating the metrics.
If cached Q/RMI files are missing, choose 'auto' explicitly.

Output
------
Original PG baseline:
    config.results_dir/<phantom>_SWS_quality_results.mat
    config.results_dir/<phantom>_SWS_baseline.csv

New quality analysis, PG and AC:
    Results/sws_measurement_quality/combined_metrics_published/
      quality_maps_<phantom>.png
      stiffness_maps_<phantom>.png
      stiffness_comparison_<phantom>_PG.png / .fig
      stiffness_comparison_<phantom>_AC.png / .fig
      stiffness_data_<phantom>.mat
      all_metrics_comparison.csv
      global_metrics_comparison.csv
      global_auc_comparison.png
      global_false_confidence.png
      stiffness_mask_coverage.csv

The returned 'results' also contains results.reliability when enabled.

Physical interpretation / caveats
---------------------------------
E = 3*rho*c_s^2, mu = rho*c_s^2 (Pa); /1000 for kPa.
E ~ 3mu assumes isotropic, approximately incompressible elasticity.
Masks are examples of thresholds, NOT error corrections.
Correlation Q 0.63: illustrative pixelwise use, not full region-based rule
from Schoen Jr. 2025. RMI-1D 0.85 failed to generalize between phantoms.
GT is used for evaluation only, not to build Q or masked stiffness.

The baseline function computes PG directly from the source dataset.
The Q/RMI helper uses the SAVED original PG and AC maps from research
scripts. A new consistency check compares new PG and GT against cached
PG and GT and stops if they disagree materially. This avoids mixing
results from different filter/dataset versions. AC is not reconstructed
in the old baseline function; it comes from the saved quality MAT files.

The original RMI code has a legacy dependency on phase_quality_results.mat
when the RMI result has not yet been computed. In mode 'auto', the phase
script may be run as a data prerequisite; Phase Planarity Q is NOT
reported, plotted or used for masking.

TESTING: Checked structural consistency and package integrity, but cannot
execute MATLAB+USTB against your datasets in this environment.
