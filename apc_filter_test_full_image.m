addpaths('..');
% clear;

fSize_axes = 12;
fSize_title = 16;

% Image directory
img_repo = '~/Desktop/piv_test_images/synthetic';

% Give this case a name
case_name = 'test_case';

% Image directory
image_dir = fullfile(img_repo, case_name);

% Make it 
if ~exist(image_dir, 'dir')
    mkdir(image_dir)
end;

num_digits = 6;
image_ext = '.tif';
start_image = 1;
end_image = 1;
skip_image = 1;
c_step = 1;

% % Image stuff
% Number of images
num_images = 1;

image_height = 1024;
image_width = 1280;

% % Region stuff
region_width = 64;
region_height = 64;

% Particle positions buffer
x_buffer = -100;
y_buffer = -100;
% Particle position max and min
x_min = 1 + x_buffer;
x_max = image_width  - x_buffer;
y_min = 1 + y_buffer;
y_max = image_height - y_buffer;

% Image noise
noise_mean_fract = 5E-2;
noise_std_fract  = 3E-4;

% noise_mean_fract = 1E-5;
% noise_std_fract = 3E-5;

% Particle stuff
d_mean = 1 * sqrt(8);
d_mean = 3;
d_std = 1;
particle_concentration = 3E-2;

% Displacement stuff
dx_mean = 5;
dy_mean = 0;
dx_rand_min = 0;
dx_rand_max = 4;
dy_rand_min = dx_rand_min;
dy_rand_max = dx_rand_max;

% APC stuff
shuffle_range_x = 64;
shuffle_step_x = 32;

shuffle_range_y = 64;
shuffle_step_y = 32;

shuffle_range = [shuffle_range_y, shuffle_range_x];
shuffle_step =  [shuffle_step_y,  shuffle_step_x];

% Grid stuff
% Grid spacing
grid_x = 128;
grid_y = 128;

% Grid buffers
grid_buffer_x = (region_width  + shuffle_range_x) * [1, 1];
grid_buffer_y = (region_height + shuffle_range_y) * [1, 1];

% Window stuff
window_fraction = 0.4 * [1, 1];

% % Do some werk
%
%
%
% 
% % Make the images
% Compute the total number of particles
aug_height = y_max - y_min + 1;
aug_width = x_max - y_min + 1;
num_particles = round(particle_concentration * aug_height * aug_width);

% Particle positions (image 1)
x_01 = (x_max - x_min) * rand(num_particles, 1) + x_min;
y_01 = (y_max - y_min) * rand(num_particles, 1) + y_min;

% Velocities
dx_rand = dx_rand_min + x_01 / image_width * dx_rand_max .* randn(num_particles, 1);
dy_rand = dy_rand_min + x_01 / image_width * dy_rand_max .* randn(num_particles, 1);


dx = - 3/2 * dx_mean * (y_01.^2 / image_height^2 - y_01/image_height) + dx_rand;

% Total displacement
% dx = dx_mean + dx_rand;
dy = dy_mean + dy_rand;

% Particle positions (image 2)
x_02 = x_01 + dx;
y_02 = y_01 + dy;

% Particle diameters
dp       = d_mean + d_std * randn(num_particles, 1);

% Particle intensities for the correlated images
particle_intensities     = d_mean ./ dp;

% Make the images
image_01_raw = generateParticleImage(image_height, image_width,...
    x_01, y_01, dp, particle_intensities);

% Generate the second image
image_02_raw = generateParticleImage(image_height, image_width,...
      x_02, y_02, dp, particle_intensities);

% Noise absolute values
noise_std = noise_std_fract * max(image_01_raw(:));
noise_mean = noise_mean_fract * max(image_01_raw(:));

% Noise matrices
noise_mat_full_01 = noise_mean + noise_std * ...
    randn(image_height, image_width);
noise_mat_full_02 = noise_mean + noise_std * ...
    randn(image_height, image_width);

% Add noise to the images.
image_01 = image_01_raw + noise_mat_full_01;
image_02 = image_02_raw + noise_mat_full_02;

% Create the image grid.
[gx, gy] = gridImage([image_height, image_width], [grid_y, grid_x], grid_buffer_y, grid_buffer_x);

gxm = reshape(gx, [length(unique(gy)), length(unique(gx))]);
gym = reshape(gy, [length(unique(gy)), length(unique(gx))]);

% Calculate the filters
[APC_STD_Y, APC_STD_X] = ...
    calculate_apc_filter_from_image_pair(image_01, image_02, ...
    gym, gxm, [region_height, region_width], window_fraction, shuffle_range, shuffle_step);

% imagesc(image_01);
% axis image;
% colormap gray;

nx = length(unique(gx(:)));
ny = length(unique(gy(:)));

A_X = reshape(APC_STD_X, [ny, nx]);
% A_Y = reshape(APC_STD_Y, size(gy));

% RPC filter
rpc_filt = spectralEnergyFilter(region_height, region_width, d_mean);

% Std dev of RPD filter
[~, rpc_std_y, rpc_std_x] = fit_gaussian_2D(rpc_filt);

plot(gx, A_X(:), '.k');
hold on;
plot([1, max(gx(:))], abs(rpc_std_x) * [1, 1], '--k');
hold off;






    

