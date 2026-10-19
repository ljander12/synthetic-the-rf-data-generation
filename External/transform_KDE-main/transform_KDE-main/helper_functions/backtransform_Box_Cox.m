function pdf_backtransformed = backtransform_Box_Cox(x_BC, pdf_BC, x, alpha)
    % x_BC = grid in transformed Box-Cox domain, e.g., 500 samples +/- 5 std around mean
    % pdf_BC = density (KDE) in transformed Box-Cox domain
    % x = target grid
    % alpha = Box-Cox transform parameter

    x_BC = x_BC(:);
    pdf_BC = pdf_BC(:);
    x = x(:);
    g_x = Box_Cox(x, alpha); % target grid mapping to Box-Cox domain for interpolation

    chain_rule = abs( d_dx_Box_Cox(x, alpha) ); % a small regularization term is included here
    pdf_backtransformed = interp1(x_BC, pdf_BC, g_x, 'linear', 0) .* chain_rule;
end