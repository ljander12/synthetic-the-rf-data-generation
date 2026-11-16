PHASE-RAMP SWS VALIDATION -- MATLAB

Place validate_phase_ramp_sws.m into the Analysis/ folder of your
synthetic-the-rf-data-generation repository.

Run:
  validate_phase_ramp_sws('large_inclusions')

Default (quick): beta = -0.2, 0, +0.2 radians/mm.
Extended:
  validate_phase_ramp_sws('large_inclusions','full')
  --> beta = -0.5, -0.2, 0, +0.2, +0.5 radians/mm.

Inputs (READ ONLY):
  - Results/method_comparison/batch/<phantom>/method_comparison.mat
  - Results/sws_measurement_quality/<phantom>/phase_consistency_diagnostic/phase_quality_results.mat
  - Dataset <phantom>_displacement_with_quality.mat

Output ONLY in a newly created folder:
  Results/sws_measurement_quality/<phantom>/phase_ramp_sws_validation/
      phase_ramp_metrics.csv
      phase_ramp_metrics.mat
      phase_ramp_maps.png
      phase_ramp_trends.png

Notes:
  * Requires the existing USTB toolbox (it attempts the project's usual path).
  * Uses the ORIGINAL time-series residual; modifies only the fitted
    f0 harmonic phase and reconstructs PG and AC using the existing recipe.
  * Q uses no ground truth. Ground truth is used only for post hoc errors.
  * All scenarios use exactly the SAME valid pixels for RMSE comparisons.
  * Since baseline maps can be systemically biased, a positive OR negative
    phase ramp might improve or worsen RMSE. Do not interpret the sign alone.
  * This is an intentionally synthetic perturbation, not a new k-Wave run.
  * We cannot execute MATLAB here. If MATLAB reports an error, send the
    exact output and line number so the script can be adjusted.
