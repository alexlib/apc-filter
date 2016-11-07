function apc_quiver_plots_02(tx_apc, ty_apc, tx_rpc, ty_rpc, tx_scc, ty_scc, grid_x, grid_y, image_size, u_max, ...
     num_ens, diffusion_std_dev)

 % Particle size
particle_diameter = 3;
particle_std_dev = particle_diameter / 4;
 
nx = length(unique(grid_x(:)));
ny = length(unique(grid_y(:)));

gy = unique(grid_y);

image_height = image_size(1);
image_width = image_size(2);

Scale = 36;
lw = 2;
fSize = 16;
fSize_titles = 14;

profile_x_scale = 1;

Skip_x = 8;
Skip_y = 1;

gx_mat = reshape(grid_x, ny, nx);
gy_mat = reshape(grid_y, ny, nx);
tx_apc_mat = reshape(tx_apc, ny, nx);
ty_apc_mat = reshape(ty_apc, ny, nx);

tx_rpc_mat = reshape(tx_rpc, ny, nx);
ty_rpc_mat = reshape(ty_rpc, ny, nx);

tx_scc_mat = reshape(tx_scc, ny, nx);
ty_scc_mat = reshape(ty_scc, ny, nx);

x_inds = 1 : Skip_x : nx;
y_inds = 1 : Skip_y : ny;

% Average the velocity profiles
tx_mean_apc = mean(tx_apc_mat, 2);
tx_mean_rpc = mean(tx_rpc_mat, 2);
tx_mean_scc = mean(tx_scc_mat, 2);

% Standard deviations
tx_std_apc = std(tx_apc_mat, 0, 2);
tx_std_rpc = std(tx_rpc_mat, 0, 2);
tx_std_scc = std(tx_scc_mat, 0, 2);


% nticks_y = 4;
% xtl = num2cell(linspace(1, 10, nticks_x) / 10);
% xt = linspace(1, image_width, nticks_x);


profile_y_min = -5;
profile_y_max = 14;
n_ticks_profile = 6;

% Ticks for profiles
yt_prof = linspace(profile_y_min, profile_y_max, n_ticks_profile);

% Ticks for quiver

% Labels for quiver
n_ticks_y_quiv = 11;
yt_quiv = linspace(1, image_height, n_ticks_y_quiv);

ytl_quiv = num2cell(linspace(0, 1, n_ticks_y_quiv));
for n = 1 : n_ticks_y_quiv
   if mod(ytl_quiv{n}, 0.2)
        ytl_quiv{n} = '';
   end  
end



xtl_quiv = {'0', [], '0.5', [], '1'};

% Axes stuff
xt_prof = linspace(1 , image_height, 10);
xt_quiv = linspace(1, image_width, 5);
xtl_prof = num2cell(0.1 : 0.1 : 1);

x_min_fract = -0.0;

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

% Error bar plots
c_red = 1 / 255 * [178,34,34];
c_green = 	1 / 255 * [34, 139, 34] ;
c_gray = 0.5 * [1, 1, 1];
c_blue = 	1 / 255 * [30, 144, 255];

linespec_red = {'-', 'color', c_red};
linespec_green = {'-', 'color', c_green};
linespec_blue = {'-', 'color', c_blue};
linespec_gray = {'-', 'color', c_gray};


fig_pos_fract_y = 1.8;
% fig_pos_fract_x = 0.9;
fig_pos_x_shift = -0.05;

axes_width = 125;
axes_height = 250;
axes_y = 500;
x01 = 100;

dx_pair = axes_width * 2.2;
dx_sub = axes_width + 7;

% Title position
t_pos_fract = 2;

aw = 1;

subplot(1, 6, 1)
quiver(gx_mat(y_inds, x_inds),...
    gy_mat(y_inds, x_inds), ...
    Scale * tx_scc_mat(y_inds, x_inds), ...
    Scale * ty_scc_mat(y_inds, x_inds),...
    0, 'color', c_blue, 'linewidth', aw);
set(gca, 'units', 'pixels');
% axis image
title_str_01 = sprintf('$\\textrm{RPC}, \\, \\sigma_\\mathbf{d} / \\sigma_\\tau = %0.1f$', diffusion_std_dev / particle_std_dev);
t_01= title(title_str_01, 'interpreter', 'latex', 'FontSize', fSize_titles);


set(gca, 'FontSize', fSize);
set(gca, 'ytick', []);
box on
ylim([1, image_height]);
xlim([1, image_width]);
set(gca, 'ytick', yt_quiv);
set(gca, 'xtick', xt_quiv);
set(gca, 'xticklabel', xtl_quiv);
set(gca, 'yticklabel', ytl_quiv);
xlabel('$x / L$', 'interpreter', 'latex', 'FontSize', fSize);
ylabel('$y / h$', 'interpreter', 'latex', 'fontsize', fSize);
grid on
p = [x01, axes_y, axes_width, axes_height];
set(gca, 'position', p);
t_pos = get(t_01, 'Position');
t_pos(1) = t_pos_fract * t_pos(1);
set(t_01, 'Position', t_pos);


subplot(1, 6, 2);
% errorbar(gy, tx_mean_scc, tx_std_scc, '-', 'color', c_blue, 'linewidth', 2);
shadedErrorBar(gy, tx_mean_scc, tx_std_scc, linespec_blue);
hold on
plot(y_highres, u_exact, '--k');
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
p = [x01 + dx_sub, axes_y, profile_x_scale * axes_width, axes_height];
set(gca, 'position', p);


subplot(1, 6, 3)
quiver(gx_mat(y_inds, x_inds),...
    gy_mat(y_inds, x_inds), ...
    Scale * tx_rpc_mat(y_inds, x_inds), ...
    Scale * ty_rpc_mat(y_inds, x_inds),...
    0, 'color', c_red, 'linewidth', aw);
% axis image
title_str_02 = sprintf('$\\textrm{RPC}, \\, \\sigma_\\mathbf{d} / \\sigma_\\tau = %0.1f$', diffusion_std_dev / particle_std_dev);
t_02= title(title_str_02, 'interpreter', 'latex', 'FontSize', fSize_titles);



set(gca, 'FontSize', fSize);
set(gca, 'ytick', []);
box on
ylim([1, image_height]);
xlim([1, image_width]);
set(gca, 'ytick', yt_quiv);
set(gca, 'xtick', xt_quiv);
set(gca, 'xticklabel', xtl_quiv);
set(gca, 'yticklabel', '');
xlabel('$x / L$', 'interpreter', 'latex', 'FontSize', fSize);
grid on
set(gca, 'units', 'pixels');
t_pos = get(t_02, 'Position');
t_pos(1) = t_pos_fract * t_pos(1);
set(t_02, 'Position', t_pos);

x02 = x01 + dx_pair;

p = [x02, axes_y, axes_width, axes_height];
set(gca, 'position', p);

subplot(1, 6, 4);
shadedErrorBar(gy, tx_mean_rpc, tx_std_rpc, linespec_red);
hold on
plot(y_highres, u_exact, '--k');
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
p = [x02 + dx_sub, axes_y, profile_x_scale * axes_width, axes_height];
set(gca, 'position', p);


subplot(1, 6, 5)
quiver(gx_mat(y_inds, x_inds),...
    gy_mat(y_inds, x_inds), ...
    Scale * tx_apc_mat(y_inds, x_inds), ...
    Scale * ty_apc_mat(y_inds, x_inds),...
    0, 'color', c_gray, 'linewidth', aw);
% axis image


title_str_03 = sprintf('$\\textrm{APC}, \\, \\sigma_\\mathbf{d} / \\sigma_\\tau = %0.1f$', diffusion_std_dev / particle_std_dev);
t_03 = title(title_str_03, 'interpreter', 'latex', 'FontSize', fSize_titles);


set(gca, 'FontSize', fSize);
set(gca, 'ytick', []);
box on
ylim([1, image_height]);
xlim([1, image_width]);
set(gca, 'ytick', yt_quiv);
set(gca, 'xtick', xt_quiv);
set(gca, 'xticklabel', xtl_quiv);
set(gca, 'yticklabel', '');
xlabel('$x / L$', 'interpreter', 'latex', 'FontSize', fSize);
grid on
set(gca, 'units', 'pixels');
t_pos = get(t_03, 'Position');
t_pos(1) = t_pos_fract * t_pos(1);
set(t_03, 'Position', t_pos);

x03 = x02 + dx_pair;

p = [x03, axes_y, axes_width, axes_height];
set(gca, 'position', p);

subplot(1, 6, 6);
shadedErrorBar(gy, tx_mean_apc, tx_std_apc, linespec_gray);
hold on
plot(y_highres, u_exact, '--k');
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
p = [x03 + dx_sub, axes_y, profile_x_scale * axes_width, axes_height];
set(gca, 'position', p);


set(gcf, 'color', 'white');

tightfig;

g = [-1989         711        1091         396];

set(gcf, 'outerposition', g);


end