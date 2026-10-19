function [pdf, x, data_struct, f] = pdf_tKDE(pixel_values, gamma, kernel_type, return_plot)
    
    %% PDF estimation using the transform KDE (tKDE) method
    % Returns the pdf of the pixel_values at evaluation points x.
    % 
    % Also returns a data_struct with auxiliary information.
    % Variables ending with "_BC" are only relevant for the Box-Cox transformed domain.
    % f is a figure handle.
    %
    % pixel_values = pixels from US image (Matrix or vector. Complex values will be converted via abs).
    % return_plot = true or false (default).
    % kernel_type = [], e.g. "normal" (default) or "epanechnikov".
    
    arguments
        pixel_values
        gamma   = [];
        kernel_type = [];
        return_plot = [];
    end
   
    if isempty(gamma);   gamma   = 1;        end 
    if isempty(kernel_type); kernel_type = "normal"; end
    if isempty(return_plot); return_plot = false;    end

    %% Initialize and check
    
    data_struct = {};
    pixel_values = abs(pixel_values(:));

    if anynan(pixel_values); error('NaN in pixel_values'); end
   
    if std(pixel_values) == 0; error("PDF is a delta distribution: std(pixel_values)==0"); end

    if any(pixel_values == 0)
        N_values = length(pixel_values);
        N_zeros = sum(pixel_values==0);
        warning( strcat( num2str(N_zeros), " of ", num2str(N_values), " values (", num2str(100*N_zeros/N_values, 3) ,"%) are true zeros! Adding eps to offset, but density may spike near zero.")) 
        pixel_values(pixel_values == 0) = eps;
    end

    %% MAIN PART: Transform data and perform KDE

    % Box-Cox
    [pixel_values_BC, alpha_BC] = Box_Cox(pixel_values);

    % 500 samples between -5std and 5std seems very robust.
    x_BC = linspace( -5*std(pixel_values_BC), +5*std(pixel_values_BC), 500);

    % Silverman's rule of thumb for bandwidth selection
    bw_BC = gamma * 0.9*min(std(pixel_values_BC), iqr(pixel_values_BC)/1.34) * length(pixel_values_BC).^(-1/5);

    % Kernel density estimation (KDE)
    [pdf_BC] = ksdensity(pixel_values_BC, x_BC, 'kernel', kernel_type, 'bandwidth', bw_BC); 

    % Get the samples in the original domain
    x = inverse_Box_Cox(x_BC, alpha_BC);
    %x = x(x>10*eps); % to better avoid (potential) numerical spikes near zero
    
    % Backtransform the PDF to the original domain
    pdf = backtransform_Box_Cox(x_BC, pdf_BC, x, alpha_BC);
    
    % Normalize PDF just incase
    pdf = pdf ./ trapz(x, pdf); % might be removed if following line is taken out above: x = x(x>10*eps);
    

    %% Store auxilary data in struct
    % Linear domain
    data_struct.x   = x; 
    data_struct.pdf = pdf;
    
    % Transformed Box-Cox domain
    data_struct.x_BC     = x_BC; 
    data_struct.pdf_BC   = pdf_BC;
    data_struct.alpha_BC = alpha_BC;
    data_struct.bw_BC    = bw_BC;

   
    %% Plot if 'return_plot' is true
    if return_plot == true
        f = figure();
    
        subplot(1,2,1)
        hhi = area(x, pdf, 'LineStyle','none', 'displayname', 'PDF estimate'); hold on; grid on;
        hhi.FaceColor = [0.1 0.6, 0.3];     
        ylabel('Probability density');
        xlabel('\itx');
        title("Linear domain")

        subplot(1,2,2)
        hhi = area(x_BC, pdf_BC, 'LineStyle','none', 'displayname', 'PDF estimate'); hold on; grid on;
        hhi.FaceColor = 'r';     
        ylabel('Probability density');
        xlabel('\ity');
        title("Transformed domain (Box-Cox)")
    else
        f=0;
    end
end




