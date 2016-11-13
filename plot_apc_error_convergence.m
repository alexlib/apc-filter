% This is the path of the error results to load

diffusion_std = 4.5;

load(sprintf('/Users/matthewgiarra/Documents/School/VT/Research/Aether/piv_test_images/poiseuille_diffusion_%0.2f/error/poiseuille_error_h128_w128_diff_std_%0.2f_000001_002000.mat', diffusion_std, diffusion_std));

dp = particle_diameter;

dp_ratio = diffusion_std / (dp / 4);

% This is the number of elements over which the
% moving mean error is calculated
kernel_length = 25;

% This is the convergence criterion
conv_criterion = 0.05;

% Calculate SCC residuals
[conv_iter_scc, diff_x_scc_mean_smoothed, diff_y_scc_mean_smoothed, image_nums] = ...
    calculate_vector_convergence(tx_err_scc_spatial_mat,...
    ty_err_scc_spatial_mat, conv_criterion, kernel_length);

% Calculate RPC residuals
[conv_iter_rpc, diff_x_rpc_mean_smoothed, diff_y_rpc_mean_smoothed] = ...
    calculate_vector_convergence(tx_err_rpc_spatial_mat,...
    ty_err_rpc_spatial_mat, conv_criterion, kernel_length);

% Calculate APC residuals
[conv_iter_apc, diff_x_apc_mean_smoothed, diff_y_apc_mean_smoothed] = ...
    calculate_vector_convergence(tx_err_apc_mat,...
    ty_err_apc_mat, conv_criterion, kernel_length);


% Plot colors
plot_colors = get_plot_colors(3);
c_blue = plot_colors(2, :);
c_red = plot_colors(1, :);
c_black = plot_colors(3, :);

% Make a plot of the residuals
semilogy(image_nums, diff_x_scc_mean_smoothed, '-', 'color', c_blue, 'linewidth', 3);
hold on
semilogy(image_nums, diff_x_rpc_mean_smoothed, '-', 'color', c_red, 'linewidth', 3);
semilogy(image_nums, diff_x_apc_mean_smoothed, '-', 'color', c_black, 'linewidth', 3);

semilogy(image_nums, diff_y_scc_mean_smoothed, '--', 'color', c_blue, 'linewidth', 3);
semilogy(image_nums, diff_y_rpc_mean_smoothed, '--', 'color', c_red, 'linewidth', 3);
semilogy(image_nums, diff_y_apc_mean_smoothed, '--', 'color', c_black, 'linewidth', 3);
plot(conv_iter_scc * [1, 1], [10^-3 10^1], '-', 'color', c_blue, 'linewidth', 1);
plot(conv_iter_rpc * [1, 1], [10^-3 10^1], '-', 'color', c_red, 'linewidth', 1);
plot(conv_iter_apc * [1, 1], [10^-3 10^1], '-', 'color', c_black, 'linewidth', 1);
hold off
axis square
grid on
xlim([0, 1000]);

legend_str{1} = '$R_{\Delta x} \, \left(\textrm{SCC}\right)$';
legend_str{2} = '$R_{\Delta x} \, \left(\textrm{RPC}\right)$';
legend_str{3} = '$R_{\Delta x} \, \left(\textrm{APC}\right)$';
legend_str{4} = '$R_{\Delta y} \, \left(\textrm{SCC}\right)$';
legend_str{5} = '$R_{\Delta y} \, \left(\textrm{RPC}\right)$';
legend_str{6} = '$R_{\Delta y} \, \left(\textrm{APC}\right)$';

h = legend(legend_str);
set(h, 'fontsize', 24);
set(h, 'interpreter', 'latex');
xlabel('$\textrm{Number of ensemble pairs}$', 'fontsize', 24, 'interpreter', 'latex');
ylabel('$\textrm{Displacement estimate residual} \, \left(\mathbf{R} = \langle \Delta\mathbf{x}_{n+1}\rangle - \langle \Delta\mathbf{x}_n\rangle \right)$', 'fontsize', 24, 'interpreter', 'latex');
title({'$\textrm{Convergence of displacement estimates}$', ...
    sprintf('$\\textrm{Poiseuille flow synthetic images,} \\, \\sigma_\\mathbf{d} / \\sigma_\\tau =  %0.1f$', dp_ratio), ...
    sprintf('$\\textrm{Convergence: } \\, R < %0.2f$', conv_criterion)}, 'interpreter', 'latex', 'fontsize', 24) 
set(gca, 'fontsize', 24);
set(gcf, 'color', 'white');
set(gcf, 'position', [-1159         409         878         796]);












