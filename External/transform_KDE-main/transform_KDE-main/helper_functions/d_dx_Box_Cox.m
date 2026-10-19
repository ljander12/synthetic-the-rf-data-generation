function output = d_dx_Box_Cox(x, alpha)
    
    % Compute a small regularization term to avoid possibility of diverging density near 0.
    % The term is set at the 5th percentile.
    x_sort = sort(x);
    reg_val = x_sort( max(ceil(0.05 * length(x_sort)), 1) );

    output = (x+reg_val).^(alpha-1);
end

