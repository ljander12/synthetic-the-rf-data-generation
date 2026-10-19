%% PDF estimation
% Draw n_samples from a Rayleigh distribution with parameter sigma.
% The Rayleigh distribution is skewed by a power transform Y = X.^exponent.

% rng(5) % A random seed can be specified.

% GENERATE DATA
n_samples = 400;                                     % number of samples
sigma = 1;                                           % Rayleigh parameter
exponent = 1.5;                                      % exponent for EPT (skews distribution)
                                                     % note that the density diverges for exponent >= 2.
samples = raylinv(rand(n_samples, 1), 1 ).^exponent; % sample Rayleigh


% PROPOSED METHOD: pdf_tKDE(pixel_values, gamma, kernel_type, return_plot)
[estimated_pdf, x, data_struct] = pdf_tKDE(samples, [], "normal", true);
% NOTE1: since "return_plot" is true, the estimated_pdf will be plotted as a function of x.
% NOTE2: "data_struct" contains auxiliary data, such as the Box-Cox domain alpha, bandwidth, axis, etc. 
% NOTE3: "f" is a figure handle.
% NOTE4: gamma = 1 is optimal for density estimation, and set as the default.


% PLOT TRUE PDF
x_true = 0:0.001:6;
pdf_true = x_true./sigma^2 .* exp( -x_true.^2 / (2*sigma.^2) ); % true Rayleigh pdf
pdf_true = transform(pdf_true, x_true, x_true.^exponent);       % transform
subplot(1,2,1)
plot(x_true, pdf_true,'-.','displayname', 'True PDF ~ Rayl.(1)', 'linewidth', 3, 'color', 'k'); hold on; grid on;
xlim([0,6])
ylim([0, 1.1*max(estimated_pdf)]);
legend



%% gCNR estimation
% Estimates gCNR from two Rayleigh distributions and compares with ground truth.

n_samples = 400; % number of samples
sigma1 = 1;
samples1 = raylinv(rand(n_samples, 1), sigma1);
sigma2 = 2;
samples2 = raylinv(rand(n_samples, 1), sigma2);


% PROPOSED METHOD: gCNR_tKDE(pixels_i, pixels_o, return_plot, gamma, kernel_type)
[gCNR_value, data_struct, f] = gCNR_tKDE(samples1, samples2, true, 0.7, []);
% NOTE: gamma = 0.7 is generally better for estimating gCNR and is set as the default.


% plot
x = data_struct.x.';
pdf_i = x./ sigma1^2 .* exp(- x.^2 / (2*sigma1^2) );
pdf_o = x./ sigma2^2 .* exp(- x.^2 / (2*sigma2^2) );
true_gCNR = 1-trapz(x.', min([pdf_i; pdf_o]));
plot(x, pdf_i,'-.','displayname', 'True PDF ~ Rayl.(1)', 'linewidth', 3, 'color', 'k'); hold on;
plot(x, pdf_o,'--','displayname', 'True PDF ~ Rayl.(2)', 'linewidth', 3, 'color', 'k'); hold on; 
grid on;

disp("true gCNR: " + num2str(true_gCNR,3))
disp("est. gCNR: " + num2str(gCNR_value,3))



%% FUNCTION FOR TRANSFORMING TRUE PDF FOR VARIOUS EXPONENTS

function f_tr = transform(f, x, g_x)
    x = x(:);
    f = f(:);
    g_x = g_x(:);

    g_inv =  interp1(g_x, x, x,'linear', 0);

    chain_rule = abs( gradient(g_inv) ./ gradient(x) );
    f_tr = interp1(x, f, g_inv,'linear', 0) .* chain_rule;
end
