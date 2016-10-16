function plot_apc_error_analysis_results(results_path_list)


if nargin < 1
    
   results_path_list{1} =  '/Users/matthewgiarra/Desktop/piv_test_images/poiseuille_diffusion_1.00/error/poiseuille_error_h128_w128_diff_std_1.00_000001_002000.mat';
   results_path_list{2} =  '/Users/matthewgiarra/Desktop/piv_test_images/poiseuille_diffusion_3.00/error/poiseuille_error_h128_w128_diff_std_3.00_000001_002000.mat';
   results_path_list{3} =  '/Users/matthewgiarra/Desktop/piv_test_images/poiseuille_diffusion_5.00/error/poiseuille_error_h128_w128_diff_std_5.00_000001_002000.mat';
end



if iscell(results_path_list)
% Load the results
    results_paths = results_path_list;
    num_jobs = length(results_paths);
else
    num_jobs = 1;
    results_paths{1} = results_path_list;
end


d_std = zeros(num_jobs, 1);

% Markers list
marker_list = {'-', '--', ':', '-.'};

% Line width
lw = 3;

% Plot font size
fSize_labels_x = 25;
fSize_labels_y = 25;
fSize_legend = 16;
fSize_axes = 20;
fSize_title = 16;

% Color specs
c_red = 1 / 255 * [178,34,34];
c_gray = 0.5 * [1, 1, 1];
c_blue = 	1 / 255 * [30, 144, 255];



% Make a new figure
figure;


% Make the first plots
for n = 1 : num_jobs
    
   
    
    
    load(results_paths{n});
    line_marker = marker_list{n};
    pairs_vect = 1 : length(mean_err_scc);
    
    d_std(n) = diffusion_std;
    
    ax1 = subplot(1, 2, 1);
    hold on;
    plot(pairs_vect, (mean_err_scc), line_marker, 'color', c_blue, 'linewidth', lw);
    plot(pairs_vect, (mean_err_rpc), line_marker, 'color', c_red, 'linewidth', lw);
    plot(pairs_vect, (mean_err_apc), line_marker, 'color', c_gray, 'linewidth', lw);
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
    line_marker = marker_list{n};
    pairs_vect = 1 : length(mean_err_scc);
    
    plot(pairs_vect, (std_err_scc), line_marker, 'color', c_blue, 'linewidth', lw);
    plot(pairs_vect, (std_err_rpc), line_marker, 'color', c_red, 'linewidth', lw);
    plot(pairs_vect, (std_err_apc), line_marker, 'color', c_gray, 'linewidth', lw);
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

xt = linspace(0, 1000, 6);
xtl= {};
xtl = num2cell(linspace(0, 800, 5));
xtl{end + 1} = '';

axis(ax1)
xt = linspace(0, 1000, 6);
set(gca, 'xtick', xt);
set(ax1, 'xticklabel', xtl);


axis(ax2)

p = get(gca, 'position');
p(1) = .98 * p(1);
set(gca, 'position', p);
set(gca, 'yticklabel', '');
set(gca, 'xtick', xt);
set(gca, 'xticklabel', xtl);

% Allocate the legend string
legend_str = {};

legend_str{1} = '$\textrm{SCC}$';
legend_str{2} = '$\textrm{RPC}$';
legend_str{3} = '$\textrm{APC}$';
legend_str{4} = '$\sigma_{\mathbf{d}} = 1.0$';
legend_str{5} = '$\sigma_{\mathbf{d}} = 3.0$';
legend_str{6} = '$\sigma_{\mathbf{d}} = 5.0$';


[legend_h, obj_h, plot_h,] = columnlegend(2, legend_str);
% legend_h = legend(legend_str);
% h = legend(legend_str);
%     set(h, 'FontSize', fSize_legend)
%     set(h, 'interpreter', 'latex');

set(obj_h(13), 'color', 'k');
set(obj_h(13), 'linestyle', marker_list{1});
set(obj_h(15), 'color', 'k');
set(obj_h(15), 'linestyle', marker_list{2});
set(obj_h(17), 'color', 'k');
set(obj_h(17), 'linestyle', marker_list{3});

% Release the hold.
hold off

% legend_pos = [0.2905    0.6089    2 * 0.1594    0.5 * 0.1427];
% 
% set(legend_h, 'position', legend_pos);


tightfig;
% set(gcf, 'position', 1.0e+03 * [ -1.2590    0.3080    0.8896    0.8676]);

% export_fig -r300 ~/Desktop/test.png;





end