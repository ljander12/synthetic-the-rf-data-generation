clear;
clc;

%% Simple test

gt = [2.32 2.32 2.32 2.32];

est = [2.30 2.40 2.50 2.25];

metrics = compute_quality_metrics(est, gt);

disp(metrics);