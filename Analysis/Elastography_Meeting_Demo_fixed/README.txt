ELASTOGRAPHY MEETING DEMO - FIXED VERSION

Fix: MATLAB mfilename('fullpath') does not include the .m extension.

PC 1:
1. Put elastography_meeting_demo.m in the Analysis folder (or a subfolder).
2. Run in MATLAB:
   clear elastography_meeting_demo
   elastography_meeting_demo
3. Wait for EXPORT COMPLETE. This also opens the saved figures.
4. Copy the ENTIRE folder synthetic-the-rf-data-generation/Meeting_Demo_PC2 to PC 2.

PC 2:
1. Open Meeting_Demo_PC2/elastography_meeting_demo.m.
2. Run elastography_meeting_demo, or elastography_meeting_demo('overview').

No k-Wave or Field II simulation is run. Optional stress figures are shown
only if their existing PNG files can be found during export.
