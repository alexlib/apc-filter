function plot_apc_error_analysis_results(results_path_list, particle_diameter)

if nargin < 1
   results_path_list{1} =  '/Users/matthewgiarra/Desktop/piv_test_images/poiseuille_diffusion_1.00/error/poiseuille_error_h128_w128_diff_std_1.00_000001_002000.mat';
   results_path_list{2} =  '/Users/matthewgiarra/Desktop/piv_test_images/poiseuille_diffusion_3.00/error/poiseuille_error_h128_w128_diff_std_3.00_000001_002000.mat';
   results_path_list{3} =  '/Users/matthewgiarra/Desktop/piv_test_images/poiseuille_diffusion_5.00/error/poiseuille_error_h128_w128_diff_std_5.00_000001_002000.mat';
end

if nargin < 2
    particle_diameter = 3;
end

if iscell(results_path_list)
% Load the results
    results_paths = results_path_list;
    num_jobs = length(results_paths);
else
    num_jobs = 1;
    results_paths{1} = results_path_list;
end

% Diffusion standard deviations
d_ratio = zeros(num_jobs, 1);

% Markers list
marker_list = {'-', '--', '-', '.-.'};

% Plot font size
fSize_labels_x = 25;
fSize_labels_y = 25;
fSize_legend = 16;
fSize_axes = 20;
fSize_title = 16;

% Plot spacing (horizontal)
dx_plot = 0.01;

% Get plot colors
color_list = get_plot_colors();

% List of line widths
lw_list = [1.5, 4, 4];

% Make the first plots
for n = 1 : num_jobs
    
    % Load results
    load(results_paths{n});
    line_marker = marker_list{n};
    pairs_vect = 1 : length(mean_err_apc);
    
    % Std dev ratio
    d_ratio(n) = diffusion_std / (particle_diameter / 4);
    
    ax1 = subplot(1, 2, 1);
    hold on;
    plot(pairs_vect, (mean_err_scc_spatial), line_marker, 'color', color_list(2, :), 'linewidth', lw_list(n));
    plot(pairs_vect, (mean_err_rpc_spatial), line_marker, 'color', color_list(1, :), 'linewidth', lw_list(n));
    plot(pairs_vect, (mean_err_apc), line_marker, 'color', color_list(3, :), 'linewidth', lw_list(n));
    
    set(gca, 'yscale', 'log');
     set(gca, 'fontsize', fSize_axes);
    axis square
    grid on
    box on
     title('$\textrm{Average Translation Errors}$', ...
    'interpreter', 'latex');    
    xlabel('Number of ensemble pairs', ...
        'FontSize', fSize_labels_x, 'interpreter', 'latex')
    
    ylabel('Error magnitude (pixels)',...
        'FontSize', fSize_labels_x, 'interpreter', 'latex');
    ylim([1E-2, 1E1]);
  
    
    ax2 = subplot(1, 2, 2);
    hold on;
    load(results_paths{n});
    pairs_vect = 1 : length(mean_err_apc);
    
    plot(pairs_vect, (std_err_scc_spatial), line_marker, 'color', color_list(2, :), 'linewidth', lw_list(n));
    plot(pairs_vect, (std_err_rpc_spatial), line_marker, 'color', color_list(1, :), 'linewidth', lw_list(n));
    plot(pairs_vect, (std_err_apc), line_marker, 'color', color_list(3, :), 'linewidth', lw_list(n));
    set(gca, 'yscale', 'log');
    xlim([0, 1000]); 
    axis square
    grid on
    box on
    set(gca, 'fontsize', fSize_axes);
    title('$\textrm{Standard Deviation of Translation Errors}$', ...
    'interpreter', 'latex');    
    xlabel('Number of ensemble pairs', ...
        'FontSize', fSize_labels_x, 'interpreter', 'latex');
    ylabel('');
    ylim([1E-2, 1E1]);
 
    set(gcf, 'color', 'white');

    
end

f_position = [-1595, 160, 1468, 1124];

set(gcf, 'position', f_position);

xt = linspace(0, 1000, 6);
xtl= {};
xtl = num2cell(linspace(0, 800, 5));
xtl{end + 1} = '';

axes(ax1);
p = get(gca, 'position');
p_left = p(1);
p_bottom = p(2);
p_width = p(3);
p_height = p(4);
xt = linspace(0, 1000, 6);
set(gca, 'xtick', xt);
set(ax1, 'xticklabel', xtl);


axes(ax2);
p = get(gca, 'position');
p(1) = p_left + p_width + dx_plot;
set(gca, 'position', p);
set(gca, 'yticklabel', '');
set(gca, 'xtick', xt);
set(gca, 'xticklabel', xtl);

% Allocate the legend string
legend_str = {};

% Format the legend string
legend_str{1} = '$\textrm{SCC}$';
legend_str{2} = '$\textrm{RPC}$';
legend_str{3} = '$\textrm{APC}$';

for n = 1 : num_jobs
    legend_str{3 + n} = sprintf('$\\sigma_{\\mathbf{d}} / \\sigma_\\tau = %0.1f$', d_ratio(n));
end

%
axes(ax1);
tightfig;
% Make the legend
[legend_h, obj_h] = columnlegend(2, legend_str);
h_pos = get(legend_h, 'position');
h_left = h_pos(1);
h_bottom = h_pos(2);
h_width = h_pos(3);
h_height = h_pos(4);

h_pos(1) = p_left + p_width - h_width + 0.04;
h_pos(2) = (p_bottom + p_height) - 0.28;
set(legend_h, 'position', h_pos);

set(obj_h(7), 'linewidth', lw_list(end));
set(obj_h(9), 'linewidth', lw_list(end));
set(obj_h(11), 'linewidth', lw_list(end));
set(obj_h(13), 'linestyle', marker_list{1}, 'linewidth', lw_list(1), 'color', 'k');
set(obj_h(15), 'linestyle', marker_list{2}, 'linewidth', lw_list(2), 'color', 'k');
set(obj_h(17), 'linestyle', marker_list{3}, 'linewidth', lw_list(3), 'color', 'k');

% Release the hold.
hold off;



end