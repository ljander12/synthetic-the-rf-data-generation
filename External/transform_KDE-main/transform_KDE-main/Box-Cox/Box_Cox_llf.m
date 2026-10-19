function llf = Box_Cox_llf(x, alpha)
    % Parameters
    % ----------
    % x : array
    %     Data to calculate Box-Cox log-likelihood for.
    % alpha : scalar
    %     Parameter for Box-Cox transformation. 
    % 
    % Returns
    % -------
    % llf : float
    %     Box-Cox log-likelihood of `x` given `alpha`. 
    % 
    % Notes
    % -----
    % The Box-Cox log-likelihood function is defined here as
    % 
    %     llf = (\alpha - 1) \sum_i(\log(x_i)) -
    %           N/2 \log(\sum_i (y_i - \bar{y})^2 / N),
    % 
    % where `y` is the Box-Cox transformed input data `x`.

    
    n_samples = length(x);
    if n_samples == 0; llf = NaN; return; end

    % logdata
    logdata = log(x);

    % variance of transformed data
    trans = Box_Cox(x, alpha);
    trans_var = var(trans);
    
    % log-likelihood function
    llf = (alpha - 1) * sum(logdata) - n_samples/2 * log(trans_var);
end
