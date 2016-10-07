% image_dir = '~/Desktop/images/experimental';
% image_base_name = 'vortexring_d03_f60_t06_';
% num_digits = 6;
% image_ext = '.tif';
% start_image = 840;
% end_image = 859;
% skip_image = 1;
% c_step = 1;

addpaths('~/Desktop');

% image_dir = '/Users/matthewgiarra/Desktop/piv_test_images/poiseuille_diffusion_3.00/raw';
% image_base_name = 'poiseuille_diffusion_3.00_';
% num_digits = 6;
% image_ext = '.tiff';
% start_image = 1;
% end_image = 100;
% skip_image = 2;
% c_step = 1;

image_dir = '/Users/matthewgiarra/Desktop/Ball';
image_base_name = 'B';
num_digits = 3;
image_ext = '.bmp';
start_image = 1;
end_image = 50;
skip_image = 1;
c_step = 1;
trailer_a = 'a';
trailer_b = 'b';

% Region sizes
region_height = 64;
region_width  = 64;

% Window fraction
window_fraction = 0.4;

% Grid spacing
grid_spacing_y = 8;
grid_spacing_x = 8;

% Shuffle
shuffle_range = 0;
shuffle_step = 0;

% Region size vector
region_size = [region_height, region_width];

% LIst of image numbers
image_nums_01 = start_image : skip_image : end_image;
image_nums_02 = image_nums_01;

% Number of images
num_pairs = length(image_nums_01);

% Declare the image list
image_list_01 = {''};
image_list_02 = {''};

% Digit string
dig_str = ['%0' num2str(num_digits) 'd'];

% Form the image path lists
for k = 1 : num_pairs
   image_name_01 = [image_base_name num2str(image_nums_01(k), dig_str) trailer_a image_ext];
   image_name_02 = [image_base_name num2str(image_nums_02(k), dig_str) trailer_b image_ext];
    
   image_list_01{k} = fullfile(image_dir, image_name_01);
   image_list_02{k} = fullfile(image_dir, image_name_02);
    
end

% Load the first image and get its size
[image_height, image_width] = size(double(imread(image_list_01{1})));

% Grid the images
grid_spacing = [grid_spacing_y, grid_spacing_x];
grid_buffer_x = region_width/2  * [1, 1];
grid_buffer_y = region_height/2 * [1, 1];

% Grid the image
[grid_x, grid_y] = gridImage([image_height, image_width],...
    grid_spacing, grid_buffer_y, grid_buffer_x);

% Do the APC
[APC_STD_Y, APC_STD_X, disp_pdf_std_dev_y, disp_pdf_std_dev_x] = ...
    calculate_apc_filter_ensemble_no_shuffle(image_list_01, image_list_02, ...
    grid_y, grid_x, region_size,...
    window_fraction);

apc_std = sqrt(APC_STD_Y.^2 + APC_STD_X.^2);

xv = (1 : region_width) - fourier_zero(region_width);
yv = (1 : region_height) - fourier_zero(region_height);

[x, y] = meshgrid(xv, yv);

sx_apc_01 = APC_STD_X(1);
sy_apc_01 = APC_STD_X(1);


nx = length(unique(grid_x));
ny = length(unique(grid_y));

S = reshape(apc_std, [ny, nx]);


% Calculate the RPC filter size
rpc_filter = spectralEnergyFilter(region_height, region_width, 3);
apc_filt_rep = exp(-(x.^2) / (2 * sx_apc_01^2)) .* exp(-(y.^2) / (2 * sy_apc_01^2));

% Because what isn't these days
[~, rpc_std_y, rpc_std_x] = fit_gaussian_2D(rpc_filter);
 
rpc_dia = sqrt(rpc_std_y^2 + rpc_std_x^2);


% ens_lengths = [5, 10, 15, 20, 50, 250, 500];

ens_lengths = 5;

for e = 1 : length(ens_lengths)
    
    num_ens = ens_lengths(e);
    
  
% Now do the correlations with the calculated filters.
% This will be moved into its own function.
% Do the APC
[ty_apc, tx_apc] = ...
    apc_ensemble(image_list_01(1:num_ens), image_list_02(1:num_ens), ...
    grid_y, grid_x, region_size,...
    window_fraction, APC_STD_Y, APC_STD_X);

[ty_rpc, tx_rpc] = ...
     rpc_ensemble_spatial(image_list_01(1:num_ens), image_list_02(1:num_ens), ...
    grid_y, grid_x, region_size,...
    window_fraction, rpc_std_y, rpc_std_x);

[ty_scc, tx_scc] = ...
    scc_ensemble_spatial_full_image(image_list_01(1:num_ens), image_list_02(1:num_ens), ...
    grid_y, grid_x, region_size, window_fraction);

tx_true = 8;
ty_true = 0;

tx_err_apc = (tx_true - tx_apc);
ty_err_apc = (ty_true - ty_apc);
tx_err_rpc = (tx_true - tx_rpc);
ty_err_rpc = (ty_true - ty_rpc);

err_mag_apc = sqrt(tx_err_apc.^2 + ty_err_apc.^2);
err_mag_rpc = sqrt(tx_err_rpc.^2 + ty_err_rpc.^2);

% figure(1);
% plot(grid_x, err_mag_apc, 'ok', 'markerfacecolor', 'black');
% hold on
% plot(grid_x, err_mag_rpc, 'or', 'markerfacecolor', 'red');
% hold off
% axis square
% grid on;

Scale = 10;
lw = 2;
fSize = 16;

Skip_x = 16;
Skip_y_quiv = 4;
Skip_y_prof = 4; 

gx_mat = reshape(grid_x, ny, nx);
gy_mat = reshape(grid_y, ny, nx);
tx_apc_mat = reshape(tx_apc, ny, nx);
ty_apc_mat = reshape(ty_apc, ny, nx);

tx_rpc_mat = reshape(tx_rpc, ny, nx);
ty_rpc_mat = reshape(ty_rpc, ny, nx);

tx_scc_mat = reshape(tx_scc, ny, nx);
ty_scc_mat = reshape(ty_scc, ny, nx);

x_inds_quiv = 1 : Skip_x : nx;
y_inds_quiv = 1 : Skip_y_quiv : ny;

y_inds_prof = 1 : Skip_y_prof : ny;


% Average the velocity profiles
tx_mean_apc = mean(tx_apc_mat, 2);
tx_mean_rpc = mean(tx_rpc_mat, 2);
tx_mean_scc = mean(tx_scc_mat, 2);

% Standard deviations
tx_std_apc = std(tx_apc_mat, 0, 2);
tx_std_rpc = std(tx_rpc_mat, 0, 2);
tx_std_scc = std(tx_scc_mat, 0, 2);

gy = unique(grid_y);

nticks_x = 4;
nticks_y = 4;
xtl = num2cell(linspace(1, 10, nticks_x) / 10);
xt = linspace(1, image_width, nticks_x);

% Axes stuff
xt_prof = linspace(1 , image_height, 10);
xt_quiv = linspace(1, image_width, 10);
xtl_prof = num2cell(0.1 : 0.1 : 1);

for p = 1 : 2 : 9
    xtl_prof{p} = '';
end


% Error bar plots
figure(1);
c_red = 1 / 255 * [178,34,34];
c_green = 	1 / 255 * [34, 139, 34] ;
c_gray = 0.5 * [1, 1, 1];
c_blue = 	1 / 255 * [30, 144, 255];

linespec_red = {'-', 'color', c_red};
linespec_green = {'-', 'color', c_green};
linespec_blue = {'-', 'color', c_blue};
linespec_gray = {'-', 'color', c_gray};


vel_max = 14;
vel_min = 4;

fig_pos_fract_y = 2.5;
% fig_pos_fract_x = 0.9;
fig_pos_x_shift = -0.05;

subplot(2, 3, 1)
quiver(gx_mat(y_inds_quiv, x_inds_quiv),...
    gy_mat(y_inds_quiv, x_inds_quiv), ...
    Scale * tx_apc_mat(y_inds_quiv, x_inds_quiv), ...
    Scale * ty_apc_mat(y_inds_quiv, x_inds_quiv),...
    0, 'color', 'black', 'linewidth', lw);
axis image
title(sprintf('$\\textrm{APC, %d pairs}$', num_ens), ...
    'interpreter', 'latex', 'FontSize', fSize);
set(gca, 'FontSize', fSize);
box on
ylim([1, image_height]);
xlim([1, image_width]);
y_tick_label_current = get(gca, 'ytick');
set(gca, 'ytick', xt_prof);
set(gca, 'xtick', xt_quiv);
set(gca, 'xticklabel', xtl_prof);
set(gca, 'yticklabel', xtl_prof);
set(gca, 'ydir', 'reverse');
xlabel('$x / L$', 'interpreter', 'latex', 'FontSize', fSize);
ylabel('$y / h$', 'interpreter', 'latex', 'fontsize', fSize);
grid on

subplot(2, 3, 2)
quiver(gx_mat(y_inds_quiv, x_inds_quiv),...
    gy_mat(y_inds_quiv, x_inds_quiv), ...
    Scale * tx_rpc_mat(y_inds_quiv, x_inds_quiv), ...
    Scale * ty_rpc_mat(y_inds_quiv, x_inds_quiv),...
    0, 'color', c_red, 'linewidth', lw);
axis image
title(sprintf('$\\textrm{RPC, %d pairs}$', num_ens), ...
    'interpreter', 'latex', 'FontSize', fSize);
set(gca, 'FontSize', fSize);
set(gca, 'ytick', []);
box on
ylim([1, image_height]);
xlim([1, image_width]);
set(gca, 'ytick', xt_prof);
set(gca, 'xtick', xt_quiv);
set(gca, 'xticklabel', xtl_prof);
set(gca, 'yticklabel', '');
set(gca, 'ydir', 'reverse');
xlabel('$x / L$', 'interpreter', 'latex', 'FontSize', fSize);
p = get(gca, 'position');
p(1) =  p(1) + fig_pos_x_shift;
set(gca, 'position', p);
grid on

subplot(2, 3, 3)
quiver(gx_mat(y_inds_quiv, x_inds_quiv),...
    gy_mat(y_inds_quiv, x_inds_quiv), ...
    Scale * tx_scc_mat(y_inds_quiv, x_inds_quiv), ...
    Scale * ty_scc_mat(y_inds_quiv, x_inds_quiv),...
    0, 'color', c_blue, 'linewidth', lw);
axis image
title(sprintf('$\\textrm{SCC, %d pairs}$', num_ens), ...
    'interpreter', 'latex', 'FontSize', fSize);
set(gca, 'FontSize', fSize);
set(gca, 'ytick', []);
box on
ylim([1, image_height]);
xlim([1, image_width]);
set(gca, 'ytick', xt_prof);
set(gca, 'xtick', xt_quiv);
set(gca, 'xticklabel', xtl_prof);
set(gca, 'yticklabel', '');
set(gca, 'ydir', 'reverse');
xlabel('$x / L$', 'interpreter', 'latex', 'FontSize', fSize);
p = get(gca, 'position');
p(1) =  p(1) + 2 * fig_pos_x_shift;
set(gca, 'position', p);
grid on

subplot(2, 3, 4);
% errorbar(gy, tx_mean_apc, tx_std_apc, '-k', 'linewidth', 2);
shadedErrorBar(gy, tx_mean_apc, tx_std_apc, linespec_gray);
set(gca, 'view', [90, 90]);
axis square;
box on;
ylim([vel_min, vel_max]);
xlim([1, image_height]);
set(gca, 'ytick', 0 : 2 : vel_max);
ylabel('$\textrm{Axial velocity (pix / frame)}$',...
    'interpreter', 'latex', 'FontSize', fSize);
xlabel('$y / h$', 'interpreter', 'latex', 'FontSize', fSize);
set(gca, 'xtick', xt_prof);
set(gca, 'xticklabel', xtl_prof);
% set(gca, 'xdir', 'reverse')
set(gca, 'FontSize', fSize);
p = get(gca, 'position');
p(2) = fig_pos_fract_y * p(2);
set(gca, 'position', p);
grid on

subplot(2, 3, 5);
% errorbar(gy, tx_mean_rpc, tx_std_rpc, '-', 'color', c_red, 'linewidth', 2);
shadedErrorBar(gy, tx_mean_rpc, tx_std_rpc, linespec_red);
set(gca, 'view', [90, 90]);
axis square;
box on;
ylim([vel_min, vel_max]);
set(gca, 'ytick', 0 : 2 : vel_max);
set(gca, 'xtick', xt_prof);
xlim([1, image_height]);
ylabel('$\textrm{Axial velocity (pix / frame)}$',...
    'interpreter', 'latex', 'FontSize', fSize);
% set(gca, 'xdir', 'reverse')
set(gca, 'xticklabel', {''});
set(gca, 'FontSize', fSize);
p = get(gca, 'position');
p(2) = fig_pos_fract_y * p(2);
p(1) =  p(1) + fig_pos_x_shift;
set(gca, 'position', p);
grid on

subplot(2, 3, 6);
% errorbar(gy, tx_mean_scc, tx_std_scc, '-', 'color', c_gray, 'linewidth', 2);
shadedErrorBar(gy, tx_mean_scc, tx_std_scc, linespec_blue);
set(gca, 'view', [90, 90]);
axis square;
box on;
ylim([vel_min, vel_max]);
set(gca, 'ytick', 0 : 2 : vel_max);
set(gca, 'xtick', xt_prof);
xlim([1, image_height]);
ylabel('$\textrm{Axial velocity (pix / frame)}$',...
    'interpreter', 'latex', 'FontSize', fSize);
% set(gca, 'xdir', 'reverse')
set(gca, 'xticklabel', {''});
set(gca, 'FontSize', fSize);
p = get(gca, 'position');
p(1) =  p(1) + 2 * fig_pos_x_shift;
p(2) = fig_pos_fract_y * p(2);
set(gca, 'position', p);
grid on



set(gcf, 'color', 'white');

set(gcf, 'position', [-2115         524        1557         955]);


tightfig;
fig_save_name = sprintf('~/Desktop/figs/fig_ball_ens_%d.png', num_ens);


export_fig(fig_save_name, '-r300');


    
    
end

apc_std_x_mat = reshape(APC_STD_X ./ rpc_std_x, ny, nx);
apc_std_y_mat = reshape(APC_STD_Y ./ rpc_std_y, ny, nx);

figure; 
subplot(2, 1, 1);
imagesc(gx_mat(:), gy_mat(:), apc_std_x_mat);
hold on
quiver(gx_mat(y_inds_quiv, x_inds_quiv),...
    gy_mat(y_inds_quiv, x_inds_quiv), ...
    Scale * tx_apc_mat(y_inds_quiv, x_inds_quiv), ...
    Scale * ty_apc_mat(y_inds_quiv, x_inds_quiv),...
    0, 'color', 'black', 'linewidth', lw);
h = colorbar;
ylabel(h, '$\sigma_x / \sigma_{\textrm{RPC}}$', 'interpreter', 'latex', 'fontsize', 20);

hold off
axis image;
title({'PIV Challenge 2003 B, Filter from ensemble of 50 pairs', 'Vectors from ensemble of 5 pairs', 'APC diameter vs RPC diameter, x direction'}, 'FontSize', fSize);
set(gca, 'fontsize', fSize);
caxis([0.4, 0.9]);
xlim([1, 1.02 * image_width]);
ylim([1, image_height]);



subplot(2, 1, 2);
imagesc(gx_mat(:), gy_mat(:), apc_std_y_mat);
hold on
quiver(gx_mat(y_inds_quiv, x_inds_quiv),...
    gy_mat(y_inds_quiv, x_inds_quiv), ...
    Scale * tx_apc_mat(y_inds_quiv, x_inds_quiv), ...
    Scale * ty_apc_mat(y_inds_quiv, x_inds_quiv),...
    0, 'color', 'black', 'linewidth', lw);
h = colorbar;
ylabel(h, '$\sigma_y / \sigma_{\textrm{RPC}}$', 'interpreter', 'latex', 'fontsize', 20);

hold off
axis image;
title('APC diameter vs RPC diameter, y direction', 'FontSize', fSize);
caxis([0.4, 0.9]);
set(gcf, 'color', 'white');
set(gca, 'fontsize', fSize);
xlim([1, 1.02 * image_width]);
ylim([1, image_height]);



% 
% 
% 



