% image_dir = '~/Desktop/images/experimental';
% image_base_name = 'vortexring_d03_f60_t06_';
% num_digits = 6;
% image_ext = '.tif';
% start_image = 840;
% end_image = 859;
% skip_image = 1;
% c_step = 1;

image_dir = '~/Desktop/images/synthetic';
image_base_name = 'lambvortex_h1024_w1024_';
num_digits = 6;
image_ext = '.tiff';
start_image = 1;
end_image = 19;
skip_image = 1;
c_step = 1;

% Region sizes
region_height = 64;
region_width  = 64;

% Window fraction
window_fraction = 0.4;

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
[APC_STD_Y, APC_STD_X] = ...
    calculate_apc_ensemble(image_list_01, image_list_02, ...
    grid_y, grid_x, region_size,...
    window_fraction, shuffle_range, shuffle_step);

apc_std = sqrt(APC_STD_Y.^2 + APC_STD_X.^2);

nx = length(unique(grid_x));
ny = length(unique(grid_y));

S = reshape(apc_std, [ny, nx]);










