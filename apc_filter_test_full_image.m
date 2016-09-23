addpaths;
% clear;


fSize_axes = 12;
fSize_title = 16;

% image_dir = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/fmc/analysis/data/synthetic/vortex/lambvortex_2014-05-22_centralDifference/c_0.0250/lambvortex_h1024_w1024_00001/raw';
% image_base_name = 'lambvortex_h1024_w1024_';
% num_digits = 6;
% image_ext = '.tiff';
% start_image = 1;
% end_image = 1;
% skip_image = 1;
% c_step = 1;

image_dir = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/fmc/analysis/data/experimental/vortex/vortexring_2013-11-12_d03_f60_t06/raw';
image_base_name = 'vortexring_d03_f60_t06_';
num_digits = 6;
image_ext = '.tif';
start_image = 840;
end_image = 840;
skip_image = 1;
c_step = 1;



% Image dimensions
region_height = 64;
region_width = 64;

% Grid spacing
grid_x = 64;
grid_y = 64;

% Grid buffer
grid_buffer_x = region_width * [1, 1];
grid_buffer_y = region_height * [1, 1];

% Window size
window_fraction = 0.4 * [1, 1];

% % Grid step
shuffle_range = 32;
shuffle_step = 16;

gx_range = shuffle_range;
gx_step = shuffle_step;

gy_range = shuffle_range;
gy_step = shuffle_step;

% Max number of reqions for the neq correlation
num_regions_neq = 1000;

% % % % %
% Center pixels
xc = (region_width  + 1) / 2 + 0.5 * (1 - mod(region_width,  2));
yc = (region_height + 1) / 2 + 0.5 * (1 - mod(region_height, 2));


% Window
g_win = gaussianWindowFilter([region_height, region_width], window_fraction, 'fraction');

% Allocate the particle shape
particle_shape_sum = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);

ncc_sum = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);

% List of image numbers
image_list_01 = start_image : skip_image : end_image;
image_list_02 = image_list_01 + c_step;

num_images = length(image_list_01);

% Load the first image to check its size
image_name_01 = sprintf('%s%06d%s', image_base_name, image_list_01(1), image_ext);
image_path_01 = fullfile(image_dir, image_name_01);
[image_height, image_width] = size(double(imread(image_path_01)));

[gx, gy] = gridImage([image_height, image_width], [grid_y, grid_x], grid_buffer_y, grid_buffer_x);

cc_sum = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);




for k = 1 : num_images;
    
    fprintf('On image %d of %d\n', k, num_images);
    
    image_name_01 = sprintf('%s%06d%s', image_base_name, image_list_01(k), image_ext);
    image_path_01 = fullfile(image_dir, image_name_01);
    
    image_name_02 = sprintf('%s%06d%s', image_base_name, image_list_02(k), image_ext);
    image_path_02 = fullfile(image_dir, image_name_02);
    
    image_01 = double(imread(image_path_01));
    image_02 = double(imread(image_path_02));
    
    SPECTRAL_FILTER = calculate_apc_filter_from_image_pair(image_01, image_02, ...
    gy, gx, region_height, region_width, window_fraction, shuffle_range, shuffle_step, num_regions_neq);

   
    
   
end









