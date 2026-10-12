function metrics = compute_quality_metrics(estimated, ground_truth, mask)
% COMPUTE_QUALITY_METRICS
% Calculates ground-truth-based quality metrics for displacement or SWS.
%
% INPUTS:
%   estimated     - Estimated displacement field or SWS map
%   ground_truth  - Corresponding ground-truth data
%   mask          - Optional logical ROI mask
%
% OUTPUT:
%   metrics       - Structure containing quality metrics
%
% Example:
%   metrics = compute_quality_metrics(sws_est, sws_gt, mask_bkg);
%
% Notes:
%   - RMSE, MAE, bias, NRMSE and error standard deviation can be used
%     for both displacement and SWS.
%
%   - Relative bias, relative error and coefficient of variation are
%     only calculated when the evaluated quantity does not change sign.
%     This avoids meaningless percentage values for oscillating
%     displacement signals centred around zero.
%
%   - Correlation is returned as NaN when either input has zero variance.


    %% -------------------------------------------------------------
    % 1. Input validation
    % --------------------------------------------------------------

    if ~isequal(size(estimated), size(ground_truth))
        error(['Estimated data and ground truth must have the same size. ' ...
               'Estimated: [%s], Ground truth: [%s].'], ...
               num2str(size(estimated)), ...
               num2str(size(ground_truth)));
    end

    if nargin < 3 || isempty(mask)
        mask = true(size(ground_truth));
    end

    if ~isequal(size(mask), size(ground_truth))
        error('Mask must have the same size as the input data.');
    end


    %% -------------------------------------------------------------
    % 2. Select valid samples
    % -------------------------------------------------------------

    valid = logical(mask) & ...
            isfinite(estimated) & ...
            isfinite(ground_truth);

    est = double(estimated(valid));
    gt  = double(ground_truth(valid));

    if isempty(est)
        error('No valid samples found inside the selected ROI.');
    end

    error_values = est - gt;


    %% -------------------------------------------------------------
    % 3. General information
    % -------------------------------------------------------------

    metrics.N = numel(est);

    metrics.mean_estimated = mean(est);
    metrics.mean_ground_truth = mean(gt);


    %% -------------------------------------------------------------
    % 4. Absolute accuracy metrics
    % -------------------------------------------------------------

    % Mean signed error
    metrics.bias = mean(error_values);

    % Mean absolute error
    metrics.MAE = mean(abs(error_values));

    % Root-mean-square error
    metrics.RMSE = sqrt(mean(error_values.^2));

    % Standard deviation of the reconstruction error
    metrics.error_std = std(error_values);


    %% -------------------------------------------------------------
    % 5. Normalized RMSE
    % -------------------------------------------------------------
    %
    % NRMSE is normalized by the RMS value of the ground truth.
    % This is suitable for oscillating displacement signals where
    % the mean ground-truth value may be close to zero.

    metrics.ground_truth_RMS = sqrt(mean(gt.^2));

    if metrics.ground_truth_RMS > eps
        metrics.NRMSE_percent = ...
            100 * metrics.RMSE / metrics.ground_truth_RMS;
    else
        metrics.NRMSE_percent = NaN;
    end


    %% -------------------------------------------------------------
    % 6. Relative metrics
    % -------------------------------------------------------------
    %
    % Relative percentage metrics are meaningful for quantities such
    % as SWS, which remain positive.
    %
    % They are not meaningful for oscillating displacement data,
    % because the signal crosses zero and has a mean close to zero.

    gt_one_sided = ...
        all(gt >= 0) || all(gt <= 0);

    if gt_one_sided && abs(metrics.mean_ground_truth) > eps

        % Relative bias
        metrics.relative_bias_percent = ...
            100 * metrics.bias / metrics.mean_ground_truth;

        % Mean absolute percentage error (MAPE)
        nonzero_gt = abs(gt) > eps;

        if all(nonzero_gt)
            metrics.relative_error_percent = ...
                100 * mean(abs(error_values ./ gt));
        else
            metrics.relative_error_percent = NaN;
        end

    else

        metrics.relative_bias_percent = NaN;
        metrics.relative_error_percent = NaN;

    end


    %% -------------------------------------------------------------
    % 7. Precision / variability
    % -------------------------------------------------------------

    % Standard deviation of estimated values
    metrics.std = std(est);

    % Coefficient of variation.
    % Only meaningful for one-sided quantities such as SWS.
    est_one_sided = ...
        all(est >= 0) || all(est <= 0);

    if est_one_sided && abs(metrics.mean_estimated) > eps
        metrics.CV_percent = ...
            100 * metrics.std / abs(metrics.mean_estimated);
    else
        metrics.CV_percent = NaN;
    end


    %% -------------------------------------------------------------
    % 8. Correlation
    % -------------------------------------------------------------
    %
    % Pearson correlation requires variation in both signals.

    if numel(est) > 1 && std(est) > eps && std(gt) > eps

        R = corrcoef(est, gt);

        metrics.correlation = R(1,2);

    else

        metrics.correlation = NaN;

    end

end