
% Font size
fSize_title = 20;
fSize_legend = 16;
plot_dx = 0.01;

xlims = 0.75 * region_width/ 2 * [-1, 1];

% Mean and std dev of positions
mu_xo = 0;
s_xo = 0.00;

% Mean and std dev of diameters
mu_d = 10;
s_d = 5;

% Region width
region_width = 200;

% Number of trials
num_trials = 1E4;

% Number of points
num_points = 5E3;

% Make the random variables
xo = mu_xo * ones(num_trials, 1) + s_xo * randn(num_trials, 1);
d  = mu_d * ones(num_trials, 1) + s_d * randn(num_trials, 1);

% Make the domain
x = linspace(-region_width/2, region_width/2, num_points);
xg = repmat(x, [num_trials, 1]);

plot_colors = zeros(size(xg));

ps2 = 0;
pu2= 0;
p2u2 = 0;

p_uniform = 1 / num_trials;

subplot(1, 2, 1);
hold off;
for k = 1 : num_trials
    plot_colors(k, :) = p_uniform * exp(-(x - xo(k)).^2 / (2 * d(k)^2));
    
    ps2 = ps2 + p_uniform * (d(k))^2;
    pu2 = pu2 + p_uniform * xo(k)^2;
    p2u2 = p2u2 + p_uniform * xo(k);
    
%     
    if mod(k, 250) == 0
        ax1 = subplot(1, 2, 1);
        plot(x, plot_colors(k, :) ./ max(plot_colors(k, :)), '-k', 'linewidth', 2);

        hold on
        
        
 
        
    end
%     

    

    
end



% p = sqrt(ps2 + pu2 - p2u2^2);
% 
g_sum = sum(plot_colors, 1);
% 
g_sum_norm = g_sum ./ max(g_sum(:));
% 
% g_analytical = exp(-(x - mu_xo).^2 / (2 * p^2));
% 
[A, sx_fit, xo_fit, B, ARRAY_gaussian] = fit_gaussian_1D(g_sum_norm, x);
% 
[A, B, XO, OFFSET, ARRAY_laplace] = fit_laplace_1D(g_sum_norm,x);
% 

yt = linspace(0, 1, 6);

plot_colors = get_plot_colors;

axes(ax1);
title_str = sprintf(...
    '$\\textrm{Individual Gaussians,} \\, \\, \\sigma_{\\sigma} = %0.0f$', s_d);
title(title_str, 'FontSize', fSize_title, 'interpreter', 'latex');
set(gca, 'xticklabel', '');
set(gca, 'yticklabel', '');
grid on
p = get(gca, 'position');
p_left = p(1);
p_bottom = p(2);
p_width = p(3);
p_height = p(4);
ylim([0, 1.1]);
set(gca, 'ytick', yt);
xlim(xlims)
box on

subplot(1, 2, 2);
plot(x, g_sum_norm, '-', 'linewidth', 10, 'color', plot_colors(3, :));
hold on
plot(x, ARRAY_gaussian, '-', 'linewidth', 6, 'color', plot_colors(2, :));
plot(x, ARRAY_laplace, '-', 'linewidth', 3, 'color', plot_colors(1, :));
% xlim(region_width * [-1, 1]);
xlim(xlims);
ylim([0, 1.1]);
hold off
set(gca, 'yticklabel', '');
set(gca, 'xticklabel', '');
grid on
title('$\textrm{Sum of Gaussians}$', 'interpreter', 'latex', 'fontsize', fSize_title);
h = legend('\textrm{Sum}', '\textrm{Gaussian fit}', '\textrm{Laplace Fit}');
set(h, 'fontsize', fSize_legend, 'interpreter', 'latex');
p = get(gca, 'position');
p(1) = p_left + p_width + plot_dx;
set(gca, 'position', p);
set(gca, 'ytick', yt);
box on

set(gcf, 'color', 'white');




% % fprintf('Measured std dev: %0.2f\nAnalytical std dev: %0.2f\n\n', sx_fit, p);
% 
% h = legend('Simulation', 'Fit to simulation', 'Analytical');
% set(h, 'fontsize', 16);
