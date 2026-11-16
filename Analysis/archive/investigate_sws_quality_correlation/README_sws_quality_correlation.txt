SWS quality correlation diagnostic (MATLAB)
==========================================

Place investigate_sws_quality_correlation.m in your project's Analysis folder.
Requires that you previously ran test_sws_measurement_quality('large_inclusions')
and/or ('high_contrast'), which created quality_validation.mat.

Run in MATLAB Command Window:
    investigate_sws_quality_correlation('large_inclusions')
    investigate_sws_quality_correlation('high_contrast')

What it does:
- Recomputes mean-centered, zero-lag displacement correlation using lateral
  separations of 1, 4, 8, 16 and 32 pixels.
- Compares raw correlation magnitude with the clipped quality score Q.
- Uses stored ground truth only to evaluate correlation with SWS error.
- Uses exactly the same ROI/mask across all spatial separations.
- Does not change or rerun k-Wave / Field II / USTB reconstructions.

Outputs:
Results/sws_measurement_quality/<phantom>/correlation_diagnostic/

Research caution:
The implementation is exploratory, not an exact recreation of Schoen et al.
(2025). The quality score for harmonically excited wavefields can saturate.
Spatial separation affects physical phase offset, so the highest observed
correlation with error is not necessarily the best method. Avoid tuning a
metric on the same cases used to claim independent validation.
