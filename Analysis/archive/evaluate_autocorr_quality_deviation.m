clear;
clc;
close all;

%% Load previous quality analysis

load('large_inclusions_autocorr_quality_analysis.mat');

C = quality_results.C_values;
err = quality_results.error_values;

%% Deviation from ideal coherence value

C_deviation = abs(C - 1);

%% Correlation

R = corrcoef(C_deviation, err);
pearson_dev = R(1,2);

spearman_dev = corr( ...
    C_deviation, ...
    err, ...
    'Type', ...
    'Spearman');

fprintf('\n============================================================\n');
fprintf(' QUALITY INDICATOR DEVIATION ANALYSIS\n');
fprintf('============================================================\n');

fprintf('Pearson corr |C-1| vs |error|:  %.4f\n', ...
    pearson_dev);

fprintf('Spearman corr |C-1| vs |error|: %.4f\n', ...
    spearman_dev);


%% -------------------------------------------------------------
% Two-sided acceptance thresholds
% -------------------------------------------------------------

deviation_thresholds = ...
    [0.01 0.02 0.03 0.05 0.075 0.10 0.15 0.20];

fprintf('\n--- Two-sided C threshold analysis ---\n');

for d = deviation_thresholds

    keep = C_deviation <= d;

    if any(keep)

        err_keep = err(keep);

        fprintf([ ...
            '|C-1| <= %.3f | retained = %6.2f %% | ' ...
            'RMSE = %.4f um | median = %.4f um | ' ...
            'P95 = %.4f um\n'], ...
            d, ...
            100*mean(keep), ...
            sqrt(mean(err_keep.^2))*1e6, ...
            median(err_keep)*1e6, ...
            prctile(err_keep,95)*1e6);

    end

end


%% -------------------------------------------------------------
% Compare samples close to and far from C = 1
% -------------------------------------------------------------

close_mask = C_deviation <= 0.05;
far_mask   = C_deviation > 0.05;

fprintf('\n--- Close to C=1 versus far from C=1 ---\n');

fprintf('Close samples: %.2f %%\n', ...
    100*mean(close_mask));

fprintf('Close median error: %.4f um\n', ...
    median(err(close_mask))*1e6);

fprintf('Close RMSE:         %.4f um\n', ...
    sqrt(mean(err(close_mask).^2))*1e6);

fprintf('\nFar samples: %.2f %%\n', ...
    100*mean(far_mask));

fprintf('Far median error: %.4f um\n', ...
    median(err(far_mask))*1e6);

fprintf('Far RMSE:         %.4f um\n', ...
    sqrt(mean(err(far_mask).^2))*1e6);


%% -------------------------------------------------------------
% Scatter plot
% -------------------------------------------------------------

figure;

scatter( ...
    C_deviation, ...
    err*1e6, ...
    5, ...
    'filled');

xlabel('|C - 1|');
ylabel('Absolute displacement error [\mum]');

title('Deviation from ideal autocorrelation coherence vs error');

grid on;


%% -------------------------------------------------------------
% Binned median error
% -------------------------------------------------------------

edges = 0:0.01:0.20;

centers = ...
    edges(1:end-1) + diff(edges)/2;

median_error = nan(size(centers));
p95_error = nan(size(centers));

for k = 1:length(centers)

    idx = ...
        C_deviation >= edges(k) & ...
        C_deviation < edges(k+1);

    if any(idx)

        median_error(k) = ...
            median(err(idx))*1e6;

        p95_error(k) = ...
            prctile(err(idx),95)*1e6;

    end

end

figure;

plot( ...
    centers, ...
    median_error, ...
    '-o', ...
    'LineWidth', ...
    1.5);

xlabel('|C - 1|');
ylabel('Median absolute displacement error [\mum]');

title('Displacement error versus deviation from C=1');

grid on;


%% Save

quality_deviation_results = struct();

quality_deviation_results.C_deviation = ...
    C_deviation;

quality_deviation_results.pearson = ...
    pearson_dev;

quality_deviation_results.spearman = ...
    spearman_dev;

save( ...
    'large_inclusions_autocorr_quality_deviation.mat', ...
    'quality_deviation_results');

fprintf('\nDeviation analysis saved.\n');