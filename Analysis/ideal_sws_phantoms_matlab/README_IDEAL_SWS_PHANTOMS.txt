RUN_IDEAL_SWS_PHANTOMS.m

Purpose: Generate exact, ideal (artifact-free) 2-D SWS ground-truth maps in MATLAB, without k-Wave, Field II, USTB or any simulation toolbox.

1. Put RUN_IDEAL_SWS_PHANTOMS.m in the project's Analysis/ folder.
2. Open it in MATLAB and click Run.
3. The figure shows 8 phantoms: homogeneous, small/large inclusion, high/low contrast, soft inclusion, two inclusions, ellipse.
4. The output folder is Results/artificial_sws_validation/ideal_phantoms/.
5. Ideal SWS maps are saved as individual *_GT.mat files, and all cases in ALL_IDEAL_SWS_PHANTOMS.mat. Each *_GT.mat contains sws_gt (m/s), inclusion_mask, x_axis/z_axis (m), case_info, E_gt_kPa, rho.
6. Set show_stiffness = true to also generate ideal Young's modulus maps (approximation for isotropic incompressible material).

IMPORTANT: These are exact user-defined ground-truth SWS maps, NOT reconstructed SWS measurements and NOT physically solved wavefields. Their inclusion properties are illustrative, NOT identical to prior project phantom configurations. Ground-truth metrics (RMSE, CNR, etc.) can compare a deliberately corrupted estimated SWS map against sws_gt. Correlation Q requires measured/synthetic displacement signal pairs; RMI-1D requires a displacement wavefield and an estimated SWS map. Neither should be interpreted as evaluated by these ground-truth maps alone.
