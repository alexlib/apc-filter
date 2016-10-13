% image_dir = '~/Desktop/images/experimental';
% image_base_name = 'vortexring_d03_f60_t06_';
% num_digits = 6;
% image_ext = '.tif';
% start_image = 840;
% end_image = 859;
% skip_image = 1;
% c_step = 1;

addpaths('..');

diffusion_std = 5;

u_max = 10;

image_dir = sprintf('/Users/matthewgiarra/Desktop/piv_test_images/poiseuille_diffusion_%0.2f/raw', diffusion_std);
image_base_name = sprintf('poiseuille_diffusion_%0.2f_', diffusion_std);
num_digits = 6;
image_ext = '.tiff';
start_image = 1;

skip_image = 2;
c_step = 1;

rpc_diameter = 3.0;

ens_lengths = 10;

end_image = 2 * ens_lengths(1);
% end_image = 100;

% image_dir = '/Users/matthewgiarra/Desktop/Ball';
% image_base_name = 'poiseuille_diffusion_3.00_';
% num_digits = 6;
% image_ext = '.tiff';
% start_image = 1;
% end_image = 100;
% skip_image = 2;
% c_step = 1;
% Region sizes
region_height = 128;
region_width  = 128;

% Window fraction
window_fraction = 0.5;

% Grid spacing
grid_spacing_y = 64;
grid_spacing_x = 64;

% grid_buffer_y = 500 * [1, 1];
% grid_buffer_x = 630 * [1, 1];

grid_buffer_y = region_height/2 * [1, 1];
grid_buffer_x = region_width/2 * [1, 1];

% Shuffle
shuffle_range = [0, 0];
shuffle_step = [0, 0];

% Region size vector
region_size = [region_height, region_width];

% LIst of image numbers
image_nums_01 = start_image : skip_image : end_image;
image_nums_02 = image_nums_01 + c_step;

% Number of images
num_pairs = length(image_nums_01);

% Declare the image list
image_list_01 = {''};
image_list_02 = {''};

% Digit string
dig_str = ['%0' num2str(num_digits) 'd'];

% Form the image path lists
for k = 1 : num_pairs
   image_name_01 = [image_base_name num2str(image_nums_01(k), dig_str) image_ext];
   image_name_02 = [image_base_name num2str(image_nums_02(k), dig_str) image_ext];
    
   image_list_01{k} = fullfile(image_dir, image_name_01);
   image_list_02{k} = fullfile(image_dir, image_name_02);
    
end

% Load the first image and get its size
[image_height, image_width] = size(double(imread(image_list_01{1})));

% Grid the images
grid_spacing = [grid_spacing_y, grid_spacing_x];
% grid_buffer_y = region_height/2 * [1, 1];
% grid_buffer_x = region_width/2 * [1, 1];

% grid_buffer_x = image_width/2 * [1, 1];

% Grid the image
[grid_x, grid_y] = gridImage([image_height, image_width],...
    grid_spacing, grid_buffer_y, grid_buffer_x);

% Do the APC
[APC_STD_Y, APC_STD_X, disp_pdf_std_dev_y, disp_pdf_std_dev_x] = ...
    calculate_apc_filter_ensemble(image_list_01, image_list_02, ...
    grid_y, grid_x, region_size,...
    window_fraction, rpc_diameter, shuffle_range, shuffle_step);

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
rpc_filter = spectralEnergyFilter(region_height, region_width, rpc_diameter);
apc_filt_rep = exp(-(x.^2) / (2 * sx_apc_01^2)) .* exp(-(y.^2) / (2 * sy_apc_01^2));

% Because what isn't these days
[~, rpc_std_y, rpc_std_x] = fit_gaussian_2D(rpc_filter);
 
rpc_dia = sqrt(rpc_std_y^2 + rpc_std_x^2);

% Minimum diameters 
apc_std_x = min(rpc_std_x, APC_STD_X);
apc_std_y = min(rpc_std_y, APC_STD_Y);


for e = 1 : length(ens_lengths)
    
    num_ens = ens_lengths(e);
    
  
% Now do the correlations with the calculated filters.
% This will be moved into its own function.
% Do the APC
[ty_apc, tx_apc] = ...
    apc_ensemble(image_list_01(1:num_ens), image_list_02(1:num_ens), ...
    grid_y, grid_x, region_size,...
    window_fraction, apc_std_y, apc_std_x);

[ty_rpc, tx_rpc] = ...
     rpc_ensemble_spatial(image_list_01(1:num_ens), image_list_02(1:num_ens), ...
    grid_y, grid_x, region_size,...
    window_fraction, rpc_std_y, rpc_std_x);

[ty_scc, tx_scc] = ...
    scc_ensemble_spatial_full_image(image_list_01(1:num_ens), image_list_02(1:num_ens), ...
    grid_y, grid_x, region_size, window_fraction);

Scale = 10;
lw = 2;
fSize = 16;

Skip_x = 1;
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

% Make a plot
apc_quiver_plots_02(tx_apc, ty_apc, ...
    tx_rpc_mat, ty_rpc_mat, ...
    tx_scc_mat, ty_scc_mat, ...
    gx_mat, gy_mat, ...
    [image_height, image_width], ...
    u_max, ...
     num_ens, ...
     diffusion_std)

fig_save_name = sprintf('~/Desktop/piv_results/figs_dp_3/fig_diff_%0.2f_ens_%d.png', diffusion_std, num_ens);


export_fig(fig_save_name, '-r300');


    
    
end





% 
% 
% 



