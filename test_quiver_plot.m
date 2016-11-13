
% Diffusion std dev
diffusion_std_dev = 4.5;

% Convergence criterion
conv_criterion = 0.05;

% Kernel length
kernel_length = 25;

% This is the directory in which to save the plot
plot_save_dir = fullfile('..');

% This is the name of the saved plot
plot_save_name = sprintf('poiseuille_error_h128_w128_diff_std_%0.2f_profile_plot.png', diffusion_std_dev);


% % % 

% test_quiver_plot
image_repo = get_image_repo();

% Data directory
data_dir = fullfile(image_repo, sprintf('poiseuille_diffusion_%0.2f', diffusion_std_dev), 'error');

% Data name
data_name = sprintf('poiseuille_error_h128_w128_diff_std_%0.2f_000001_002000.mat', diffusion_std_dev);

% Data path
data_path = fullfile(data_dir, data_name);

% Load the data file
load(data_path);

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

% Iteration numbers to plot
iteration_numbers = [10, conv_iter_apc, conv_iter_rpc, conv_iter_scc];

% Diffusion ratio
diffusion_ratio = diffusion_std_dev / (particle_diameter / 4);

% Make the plots
draw_apc_quiver_plots(tx_scc_spatial_mat, ty_scc_spatial_mat, ...
                    tx_rpc_spatial_mat, ty_rpc_spatial_mat, ...
                    tx_apc_mat, ty_apc_mat, ...
                    grid_x, grid_y, image_size, dx_max, ...
                    iteration_numbers, diffusion_ratio);

                
% This is the path to the saved plot
plot_save_path = fullfile(plot_save_dir, plot_save_name);
export_fig('-r300', plot_save_path);
                
% % Plot save name
% save_name = sprintf('profile_plot_diffusion_%0.1f_num_pairs_%d.eps', diffusion_std_dev, num_pairs);
% 
% % Plot save directory
% save_dir = '~/Desktop/profile_plots';
% save_path = fullfile(save_dir, save_name);

% print(1, '-depsc', save_path);































