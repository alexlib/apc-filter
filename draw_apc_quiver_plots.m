function draw_apc_quiver_plots(tx_scc_spatial_mat, ty_scc_spatial_mat, ...
                    tx_rpc_spatial_mat, ty_rpc_spatial_mat, ...
                    tx_apc_spatial_mat, ty_apc_spatial_mat, ...
                    grid_x, grid_y, image_size, u_max, ...
                    iteration_numbers, diffusion_ratio)
% Some plotting preferences
Scale = 36;
fSize = 16;
fSize_titles = 20;
fSize_textbox = 24;
fSize_labels = 24;
profile_x_scale = 1;
Skip_x = 8;
Skip_y = 1;
profile_y_min = -5;
profile_y_max = 14;

% Line width for profile true solution
lw = 1;

% Quiver arrow width
arrow_width = 1.1;

% Textbox spacing
dx_textbox = -85;
dy_textbox = 128;

% Title spacing
titles_dx = -80;
titles_dy = 28;

% Axes height
axes_height = 175;
axes_width = axes_height/2;

% Subplot pair spacing
dx_pair = axes_width * 2.2;
dx_sub = axes_width + 7;
dy_pair = axes_height * 1.05;

% Where the top left plot goes
x_left_anchor = 100;
y_bottom_anchor = 400;

% Number of ticks for the velocity profile
n_ticks_profile = 6;

% Labels for quiver
n_ticks_y_quiv = 11;

% Image size
image_height = image_size(1);
image_width = image_size(2);

% Number of grid points
nx = length(unique(grid_x(:)));
ny = length(unique(grid_y(:)));
gy = unique(grid_y);

% Indices to plot
x_inds = 1 : Skip_x : nx;
y_inds = 1 : Skip_y : ny;

% Number of plot columns and rows
num_plot_cols = 6;
num_plot_rows = length(iteration_numbers);

% Ticks for profiles
yt_prof = linspace(profile_y_min, profile_y_max, n_ticks_profile);

% Y ticks for the quiver plot
yt_quiv = linspace(1, image_height, n_ticks_y_quiv);

% Y tick labelsfor the quiver plot
ytl_quiv = num2cell(linspace(0, 1, n_ticks_y_quiv));
for n = 1 : n_ticks_y_quiv
   if mod(ytl_quiv{n}, 0.2)
        ytl_quiv{n} = '';
   end  
end

% X ticks for the velocity profile plot
xt_prof = linspace(1 , image_height, 10);
% X tick labels for the profile plot
xtl_prof = num2cell(0.1 : 0.1 : 1);

% X ticks for the quiver plot
xt_quiv = linspace(1, image_width, 5);
% Make x tick labels for the quiver plot
xtl_quiv = {'0', [], '0.5', [], '1'};

% X tick labels for the profile plot
for p = 1 : 2 : 9
    xtl_prof{p} = '';
end

% High resolution y coordinate
y_highres = linspace(min(gy), max(gy), 1000);

% Center of the image in the height direction
yc = image_height / 2;
    
% Radial coordinate in the height direction
r = abs((y_highres - yc) / (image_height/2));

% Velocity profile exact solution.
u_exact = - u_max * (r.^2 - 1);

% Make the grid coordinates into arrays
gx_mat = reshape(grid_x, ny, nx);
gy_mat = reshape(grid_y, ny, nx);


% These are the left edges of the quiver and profile plots
x_left_quiver = x_left_anchor + (0:2) * dx_pair;
x_left_profile = x_left_quiver + dx_sub;

% Y positions of subplots
y_bottom = y_bottom_anchor - (0 : num_plot_rows - 1) * dy_pair;

% Colors
plot_colors = get_plot_colors(3);

c_red = plot_colors(1, :);
c_blue = plot_colors(2, :);
c_black = plot_colors(3, :);

% Associate colors with correlations
scc_color = c_blue;
rpc_color = c_red;
apc_color = c_black;

% Line specs
linespec_scc_profile = {'-', 'color', scc_color};
linespec_rpc_profile = {'-', 'color', rpc_color};
linespec_apc_profile = {'-', 'color', apc_color};

% Textbox strings
textbox_str{1} = sprintf('$%d \\, \\textrm{pairs}$', iteration_numbers(1));
textbox_str{2} = {'$\textrm{APC}$', '$\textrm{Convergence}$'  ...
    sprintf('$\\left( %d \\, \\textrm{pairs} \\right)$', iteration_numbers(2))};
textbox_str{3} = {'$\textrm{RPC}$', '$\textrm{Convergence}$'  ...
    sprintf('$\\left( %d \\, \\textrm{pairs} \\right)$', iteration_numbers(3))};
textbox_str{4} = {'$\textrm{SCC}$', '$\textrm{Convergence}$'  ...
    sprintf('$\\left( %d \\, \\textrm{pairs} \\right)$', iteration_numbers(4))};

% Title strings
title_str{1} = sprintf(...
    '$\\textrm{SCC}, \\, \\sigma_\\mathbf{d} / \\sigma_\\tau = %0.1f$', ...
    diffusion_ratio);
title_str{2} = sprintf(...
    '$\\textrm{RPC}, \\, \\sigma_\\mathbf{d} / \\sigma_\\tau = %0.1f$', ...
    diffusion_ratio);
title_str{3} = sprintf(...
    '$\\textrm{APC}, \\, \\sigma_\\mathbf{d} / \\sigma_\\tau = %0.1f$', ...
    diffusion_ratio);

% Delete any existing textboxes
delete(findall(gcf,'Tag','text_box'))

% Loop over the plots
for n = 1 : num_plot_rows
    
    % Read the displacements (SCC)
    tx_scc = tx_scc_spatial_mat(:, iteration_numbers(n));
    ty_scc = ty_scc_spatial_mat(:, iteration_numbers(n));
    
    % Read the displacements (RPC)
    tx_rpc = tx_rpc_spatial_mat(:, iteration_numbers(n));
    ty_rpc = ty_rpc_spatial_mat(:, iteration_numbers(n));
    
    % Read the displacements (APC)
    tx_apc = tx_apc_spatial_mat(:, iteration_numbers(n));
    ty_apc = ty_apc_spatial_mat(:, iteration_numbers(n));
    
    % Reshape them into arrays
    tx_mat_scc = reshape(tx_scc, ny, nx);
    ty_mat_scc = reshape(ty_scc, ny, nx);
    tx_mat_rpc = reshape(tx_rpc, ny, nx);
    ty_mat_rpc = reshape(ty_rpc, ny, nx);
    tx_mat_apc = reshape(tx_apc, ny, nx);
    ty_mat_apc = reshape(ty_apc, ny, nx);
    
    % Mean displacements along the horizontal direction
    tx_mean_scc = mean(tx_mat_scc, 2);
    tx_mean_rpc = mean(tx_mat_rpc, 2);
    tx_mean_apc = mean(tx_mat_apc, 2);
    
    % Standard deviations of displacements
    tx_std_scc = std(tx_mat_scc, [], 2);
    tx_std_rpc = std(tx_mat_rpc, [], 2);
    tx_std_apc = std(tx_mat_apc, [], 2);
    
    % % % %
    % Do the plotting
    % % %
    % Find the linear index of the plots
    ind_scc_quiver  = sub2ind([num_plot_cols, num_plot_rows], 1, n);
    ind_scc_profile = sub2ind([num_plot_cols, num_plot_rows], 2, n);
    
    ind_rpc_quiver  = sub2ind([num_plot_cols, num_plot_rows], 3, n);
    ind_rpc_profile = sub2ind([num_plot_cols, num_plot_rows], 4, n);
    
    ind_apc_quiver  = sub2ind([num_plot_cols, num_plot_rows], 5, n);
    ind_apc_profile = sub2ind([num_plot_cols, num_plot_rows], 6, n);
    
    
    % Figure out which y tick labels to use
    ytl_current = ytl_quiv;
    if n > 1
        ytl_current{end} = '';
    end
    
    % % % SCC PLOTS % % %
    
    % SCC quiver plot
    subplot(num_plot_rows, num_plot_cols, ind_scc_quiver);
    quiver(gx_mat(y_inds, x_inds),...
    gy_mat(y_inds, x_inds), ...
    Scale * tx_mat_scc(y_inds, x_inds), ...
    Scale * ty_mat_scc(y_inds, x_inds),...
    0, 'color', scc_color, 'linewidth', arrow_width);
    set(gca, 'FontSize', fSize);
    set(gca, 'ytick', []);
    box on
    ylim([1, image_height]);
    xlim([1, image_width]);
    set(gca, 'ytick', yt_quiv);
    set(gca, 'xtick', xt_quiv);
    if n == num_plot_rows
        set(gca, 'xticklabel', xtl_quiv);
    else
        set(gca, 'xticklabel', '');
    end
    set(gca, 'yticklabel', ytl_current);
    
    % Y label
    ylabel('$y / h$', 'interpreter', 'latex', 'fontsize', fSize_labels);
    
    % X label
    if n == num_plot_rows
        xlabel('$x / L$', 'interpreter', 'latex', 'FontSize', fSize_labels);
    end

    grid on
    set(gca, 'units', 'pixels');
    p = [x_left_quiver(1), y_bottom(n), axes_width, axes_height];
    set(gca, 'position', p);
    
    % SCC profile plot
    subplot(num_plot_rows, num_plot_cols, ind_scc_profile);
    shadedErrorBar(gy, tx_mean_scc, tx_std_scc, linespec_scc_profile);
    hold on
    plot(y_highres, u_exact, '--k', 'linewidth', lw);
    hold off
    set(gca, 'view', [90, 90]);
    box on;
    ylim([profile_y_min, profile_y_max]);
    xlim(image_height * [0, 1]);
    set(gca, 'ytick', yt_prof);
    set(gca, 'xtick', xt_prof);
    set(gca, 'xticklabel', '');
    set(gca, 'yticklabel', '');
    grid on
    set(gca, 'xdir', 'reverse')
    set(gca, 'FontSize', fSize);
    set(gca, 'units', 'pixels');
    p = [x_left_profile(1), y_bottom(n), profile_x_scale * axes_width, axes_height];
    set(gca, 'position', p);
    
    
     % % % RPC PLOTS % % %
     
    % RPC quiver plot
    subplot(num_plot_rows, num_plot_cols, ind_rpc_quiver);
    quiver(gx_mat(y_inds, x_inds),...
    gy_mat(y_inds, x_inds), ...
    Scale * tx_mat_rpc(y_inds, x_inds), ...
    Scale * ty_mat_rpc(y_inds, x_inds),...
    0, 'color', rpc_color, 'linewidth', arrow_width);
    set(gca, 'FontSize', fSize);
    set(gca, 'ytick', []);
    box on
    ylim([1, image_height]);
    xlim([1, image_width]);
    set(gca, 'ytick', yt_quiv);
    set(gca, 'xtick', xt_quiv);
    if n == num_plot_rows
        set(gca, 'xticklabel', xtl_quiv);
    else
        set(gca, 'xticklabel', '');
    end
    set(gca, 'yticklabel', '');
    grid on
    set(gca, 'units', 'pixels');
    p = [x_left_quiver(2), y_bottom(n), axes_width, axes_height];
    set(gca, 'position', p);
    % X label
    if n == num_plot_rows
        xlabel('$x / L$', 'interpreter', 'latex', 'FontSize', fSize_labels);
    end
    
    
    % RPC profile plot
    subplot(num_plot_rows, num_plot_cols, ind_rpc_profile);
    shadedErrorBar(gy, tx_mean_rpc, tx_std_rpc, linespec_rpc_profile);
    hold on
    plot(y_highres, u_exact, '--k', 'linewidth', lw);
    hold off
    set(gca, 'view', [90, 90]);
    box on;
    ylim([profile_y_min, profile_y_max]);
    xlim(image_height * [0, 1]);
    set(gca, 'ytick', yt_prof);
    set(gca, 'xtick', xt_prof);
    set(gca, 'xticklabel', '');
    set(gca, 'yticklabel', '');
    grid on
    set(gca, 'xdir', 'reverse')
    set(gca, 'FontSize', fSize);
    set(gca, 'units', 'pixels');
    p = [x_left_profile(2), y_bottom(n), profile_x_scale * axes_width, axes_height];
    set(gca, 'position', p);
    

    % % % APC PLOTS % % %    

    % APC quiver plot
    subplot(num_plot_rows, num_plot_cols, ind_apc_quiver);
    quiver(gx_mat(y_inds, x_inds),...
    gy_mat(y_inds, x_inds), ...
    Scale * tx_mat_apc(y_inds, x_inds), ...
    Scale * ty_mat_apc(y_inds, x_inds),...
    0, 'color', apc_color, 'linewidth', arrow_width);
    set(gca, 'FontSize', fSize);
    set(gca, 'ytick', []);
    box on
    ylim([1, image_height]);
    xlim([1, image_width]);
    set(gca, 'ytick', yt_quiv);
    set(gca, 'xtick', xt_quiv);
    if n == num_plot_rows
        set(gca, 'xticklabel', xtl_quiv);
    else
        set(gca, 'xticklabel', '');
    end
    set(gca, 'yticklabel', '');
    % X label
    if n == num_plot_rows
        xlabel('$x / L$', 'interpreter', 'latex', 'FontSize', fSize_labels);
    end
    grid on
    set(gca, 'units', 'pixels');
    p = [x_left_quiver(3), y_bottom(n), axes_width, axes_height];
    set(gca, 'position', p);
    
    % APC profile plot
    subplot(num_plot_rows, num_plot_cols, ind_apc_profile);
    shadedErrorBar(gy, tx_mean_apc, tx_std_apc, linespec_apc_profile);
    hold on
    plot(y_highres, u_exact, '--k', 'linewidth', lw);
    hold off
    set(gca, 'view', [90, 90]);
    box on;
    ylim([profile_y_min, profile_y_max]);
    xlim(image_height * [0, 1]);
    set(gca, 'ytick', yt_prof);
    set(gca, 'xtick', xt_prof);
    set(gca, 'xticklabel', '');
    set(gca, 'yticklabel', '');
    grid on
    set(gca, 'xdir', 'reverse')
    set(gca, 'FontSize', fSize);
    set(gca, 'units', 'pixels');
    p = [x_left_profile(3), y_bottom(n), profile_x_scale * axes_width, axes_height];
    set(gca, 'position', p);
    
  % Textbox positions
    textbox_x  =  x_left_profile(end) + 2 * axes_width + dx_textbox;
    textbox_y = y_bottom(n) + dy_textbox;
  % Textbox position
    textbox_pos = [textbox_x, textbox_y, 0, 0];
    textbox_str_current = textbox_str{n};
   
    % Make the annotation
    row_annotation = annotation('textbox', ...
        'string', textbox_str_current, ...
        'fontsize', fSize_textbox, ...
        'interpreter', 'latex', ...
        'Tag', 'text_box', 'units', 'pixels');
    
    % Shift the first label down a little bit
    if n == 1
        textbox_pos(2) = textbox_pos(2) - 25;
    end
    
    % Set the label position
    set(row_annotation, 'position', textbox_pos);
    
   
    
end

for p = 1 : 3
    % Textbox positions
    title_x = x_left_anchor + (p-1) * dx_pair + 1 * axes_width + titles_dx;
    title_y = y_bottom_anchor + 1 * axes_height + titles_dy;
  % Textbox position
    title_pos = [title_x, title_y, 0, 0];
    title_str_current = title_str{p};
   
    % Make the annotation
    title_annotation = annotation('textbox', ...
        'string', title_str_current, ...
        'fontsize', fSize_titles, ...
        'interpreter', 'latex', ...
        'Tag', 'text_box', 'units', 'pixels');
    set(title_annotation, 'position', title_pos);
    
end


set(gcf, 'color', 'white');

% tightfig;

% g = [-1989         711        1091         396];

% set(gcf, 'outerposition', g);


end