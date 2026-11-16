ELASTOGRAPHY MEETING DEMO - ONE MATLAB FILE (PC 1 AND PC 2)
==========================================================

1. PC 1 (THE COMPUTER WITH THE PROJECT DATA)

Extract elastography_meeting_demo.m into your project's Analysis folder:
  synthetic-the-rf-data-generation\Analysis\

Start MATLAB and run:
  elastography_meeting_demo

It will locate your actual project Results directory, EXPORT the
saved phantom figures, comparison metrics and any stress / diagnostic
figures to the project root, then OPEN the results immediately.

The exported folder will be:
  synthetic-the-rf-data-generation\Meeting_Demo_PC2

No new k-Wave, Field II or USTB calculation is performed.

2. PC 2

Copy ENTIRE Meeting_Demo_PC2 folder to PC 2.
Open Meeting_Demo_PC2/elastography_meeting_demo.m in MATLAB and press RUN.
The same MATLAB file automatically uses its nearby PNG/CSV exports.

3. OPTIONS

  elastography_meeting_demo('overview')          4 phantoms / 3 SWS maps each
  elastography_meeting_demo('large_inclusions')  detail maps + error maps
  elastography_meeting_demo('metrics')           result table + plots
  elastography_meeting_demo('stress')            stress test figures, if saved
  elastography_meeting_demo('diagnostic')        GT diagnostic, if saved
  elastography_meeting_demo('export')            PC 1: export only

WHY THE PRIOR VIEWER SHOWED BLANK FIGURES

The previous show_elastography_pc2.m was run from the SETUP folder
instead of from a previously exported Meeting_Demo_PC2 data folder.
The new single-file entry point makes this much harder to do by mistake.

PREREQUISITES
- PC 1 must have Results/method_comparison/batch from the four-case batch.
- For the stress-test section, saved PNGs should exist under
  Results/metric_stress_test. Missing optional outputs are skipped.
- The GT vs estimated displacement diagnostic is exploratory.

Tested by static inspection, not in the MATLAB interpreter.
