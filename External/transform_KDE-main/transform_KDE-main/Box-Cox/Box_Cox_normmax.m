function estimated_alpha = Box_Cox_normmax(x)
    % Compute optimal Box_Cox transform parameter.
    % 
    % Compute optimal Box_Cox transform parameter for input data, using
    % maximum likelihood estimation.
    % 
    % Parameters
    % ----------
    % x : array
    %
    % Note
    % ----------
    % Some implementations use a bracket for searching. We have so far
    % found unbounded search with the starting value of 0 to be robust.
   
    arguments
        x
    end

    neg_llf = @(alpha) -Box_Cox_llf(x, alpha);
    
    estimated_alpha = fminsearch( neg_llf, 0 );
end