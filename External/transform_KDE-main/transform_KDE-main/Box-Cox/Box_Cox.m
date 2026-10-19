function [y, alpha] = Box_Cox(x, alpha)
    % Return a dataset transformed by a Box-Cox power transformation.
    % Adapted from SciPy (v1.11.4) implementation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.boxcox.html
    %
    %    Parameters
    %    ----------
    %    x : Input array to be transformed.
    % 
    %    alpha : scalar (optional)
    %        If `alpha` is None (default), find the value of `alpha` that maximizes
    %        the log-likelihood function and return it as the second output argument.
    % 
    %        If `alpha` is not None, do the transformation for that value.
    % 
    %    Notes
    %    -----
    %    The Box-Cox transform is given by::
    % 
    %        y = (x.^alpha - 1) / alpha,  for alpha ~= 0
    %            log(x),                  for alpha == 0
    % 
    %    `Box_Cox` requires the input data to be positive.  Sometimes a Box-Cox
    %    transformation provides a shift parameter to achieve this; `Box_Cox` does
    %    not.  Such a shift parameter is equivalent to adding a positive constant to
    %    `x` before calling `Box_Cox`.
    %
    %    This transform has the opposite order of argumnts from Matlab's
    %    "boxcox" function, which is based on nargin instead of the
    %    positional parameters.
    % 
    %    References
    %    ----------
    %    G.E.P. Box and D.R. Cox, "An Analysis of Transformations", Journal of the
    %    Royal Statistical Society B, 26, 211-252 (1964).
    
    arguments
        x;
        alpha = [];
    end
    
    x = x(:);
    if isempty(x); y=0; return; end
    if ~isreal(x); error('Box-Cox transformation is not defined for complex numbers.'); end
    if any(x<=0);  error('Box-Cox transformation is not defined values <= 0.'); end

    if isempty(alpha) % alpha is not set. 
        alpha = Box_Cox_normmax(x);
    end
            
    % This is the Box-Cox transformation
    if alpha == 0
        y = log(x);
    else 
        y = (x.^alpha - 1) ./ alpha;
    end

end
    
