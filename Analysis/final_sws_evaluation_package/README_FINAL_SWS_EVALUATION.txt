FINAL SWS EVALUATION: PG VS AC, FOUR EXISTING DIGITAL PHANTOMS

Purpose
-------
Summarize ALREADY COMPUTED results from four phantoms:
  chaoran_baseline, low_contrast, high_contrast, large_inclusions.
No new k-Wave, Field II, USTB, displacement or SWS reconstruction is executed
by SUMMARIZE_FINAL_SWS_RESULTS.m. Ground-truth-independent metrics are excluded.

Install into your MATLAB repo
-----------------------------
1. Put SUMMARIZE_FINAL_SWS_RESULTS.m in the existing Analysis/ folder.
2. (Option A) Run directly in MATLAB:
       T = SUMMARIZE_FINAL_SWS_RESULTS();
   (Option B) If you want a single entry point, replace your controller
   RUN_PHANTOM_ANALYSIS.m with the supplied
   RUN_PHANTOM_ANALYSIS_final_summary_mode.m (rename to RUN_PHANTOM_ANALYSIS.m
   and keep a backup). It has RUN_FINAL_SWS_SUMMARY_ONLY = true.
   Then just press Run in RUN_PHANTOM_ANALYSIS.m.

INPUTS (existing)
-----------------
Preferred:
  Results/method_comparison/batch/batch_summary.csv
Or:
  Results/method_comparison/batch/batch_summary.mat
Or one file per phantom:
  Results/method_comparison/batch/<phantom>/method_metrics.csv
  Results/method_comparison/batch/<phantom>/method_comparison.mat

If batch outputs are missing, use your EXISTING code (no new simulations):
  compare_SWS_methods_batch(["chaoran_baseline", "low_contrast", ...
      "high_contrast", "large_inclusions"]);
This may reconstruct PG/AC from previously saved displacement using USTB.
It requires your earlier phantom-specific baseline MAT results and dataset MAT.
The summary script itself does NOT invoke that batch pipeline.

OUTPUTS
-------
  Results/final_sws_evaluation/final_sws_comparison.csv
  Results/final_sws_evaluation/global_regional_rmse.png
  Results/final_sws_evaluation/contrast_metrics.png
  Results/final_sws_evaluation/final_sws_table.tex
  Results/final_sws_evaluation/key_findings.txt

Interpretation warnings
-----------------------
- RMSE, bias: evaluate speed accuracy against known digital-phantom GT.
- CNR, direct-KDE gCNR: evaluate separability, not correctness of speed.
- ContrastRecovery input is a FRACTION and output percent is multiplied by 100.
- NaN in edge metrics means no valid edge metric. It is NOT 0 error.
- Compare PG/AC on the SAME valid pixels (as prior batch scripts do).
- If results have different frequencies, review comparability.

PRELIMINARY REFERENCE
---------------------
preliminary_sws_comparison_from_report.csv is transcribed from the EXISTING
LaTeX draft comparison table. This is a convenient reference for values
already reported, NOT a replacement for rerunning/reloading the actual batch
result files. It includes only those six previously documented columns.

MATLAB not installed in this environment: static checks only.
