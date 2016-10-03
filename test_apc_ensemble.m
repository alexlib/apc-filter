% image_dir = '~/Desktop/images/experimental';
% image_base_name = 'vortexring_d03_f60_t06_';
% num_digits = 6;
% image_ext = '.tif';
% start_image = 840;
% end_image = 859;
% skip_image = 1;
% c_step = 1;

addpaths('~/Desktop');

image_dir = '/Users/matthewgiarra/Desktop/piv_test_images/synthetic/poiseuille';
image_base_name = 'poiseuille_';
num_digits = 6;
image_ext = '.tiff';
start_image = 1;
end_image = 110;
skip_image = 2;
c_step = 1;

% Region sizes
region_height = 64;
region_width  = 64;

% Window fraction
window_fraction = 0.5;

% Grid spacing
grid_spacing_y = 64;
grid_spacing_x = 64;

% Shuffle
shuffle_range = 0;
shuffle_step = 0;

% Region size vector
region_size = [region_height, region_width];

% LIst of image numbers
image_nums_01 = start_image : skip_image : end_image;
image_nums_02 = image_nums_01 + c_step;

% Number of images
num_images = length(image_nums_01);

% Declare the image list
image_list_01 = {''};
image_list_02 = {''};

% Digit string
dig_str = ['%0' num2str(num_digits) 'd'];

% Form the image path lists
for k = 1 : num_images
   image_name_01 = [image_base_name num2str(image_nums_01(k), dig_str) image_ext];
   image_name_02 = [image_base_name num2str(image_nums_02(k), dig_str) image_ext];
    
   image_list_01{k} = fullfile(image_dir, image_name_01);
   image_list_02{k} = fullfile(image_dir, image_name_02);
    
end

% Load the first image and get its size
[image_height, image_width] = size(double(imread(image_list_01{1})));

% Grid the images
grid_spacing = [grid_spacing_y, grid_spacing_x];
grid_buffer_y = grid_spacing_y * [1, 1];
grid_buffer_x = grid_spacing_x * [1, 1];
% Grid the image
[grid_x, grid_y] = gridImage([image_height, image_width],...
    grid_spacing, grid_buffer_y, grid_buffer_x);

% Do the APC
[APC_STD_Y, APC_STD_X, disp_pdf_std_dev_y, disp_pdf_std_dev_x] = ...
    calculate_apc_filter_ensemble(image_list_01, image_list_02, ...
    grid_y, grid_x, region_size,...
    window_fraction, shuffle_range, shuffle_step);

apc_std = sqrt(APC_STD_Y.^2 + APC_STD_X.^2);

nx = length(unique(grid_x));
ny = length(unique(grid_y));

S = reshape(apc_std, [ny, nx]);


% Calculate the RPC filter size
rpc_filter = spectralEnergyFilter(region_height, region_width, 3);

% Because what isn't these days
[~, rpc_std_y, rpc_std_x] = fit_gaussian_2D(rpc_filter);
 
rpc_dia = sqrt(rpc_std_y^2 + rpc_std_x^2);

plot(grid_x, apc_std, 'ok');
hold on
plot([0, max(grid_x(:))], rpc_dia * [1, 1], '--k', 'linewidth', 3);
hold off
axis square

num_ens = 100;

% Now do the correlations with the calculated filters.
% This will be moved into its own function.
% Do the APC
[ty_apc, tx_apc] = ...
    apc_ensemble(image_list_01(1:num_ens), image_list_02(1:num_ens), ...
    grid_y, grid_x, region_size,...
    window_fraction, APC_STD_Y, APC_STD_X);

[ty_rpc, tx_rpc] = ...
     apc_ensemble(image_list_01(1:num_ens), image_list_02(1:num_ens), ...
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

figure(1);
plot(grid_x, err_mag_apc, 'ok', 'markerfacecolor', 'black');
hold on
plot(grid_x, err_mag_rpc, 'or', 'markerfacecolor', 'red');
hold off
axis square
grid on;

Scale = 10;

figure(2);
subplot(1, 3, 1);
quiver(grid_x, grid_y, Scale * tx_apc, Scale * ty_apc, 0, 'black', 'linewidth', 2);
axis image

subplot(1, 3, 2);
quiver(grid_x, grid_y, Scale * tx_rpc, Scale * ty_rpc, 0, 'red', 'linewidth', 2);
axis image

subplot(1, 3, 3);
quiver(grid_x, grid_y, Scale * tx_scc, Scale * ty_scc, 0, 'blue', 'linewidth', 2);
axis image




% 
% 
% 



