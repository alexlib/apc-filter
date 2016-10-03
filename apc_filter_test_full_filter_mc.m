addpaths;
% clear;

do_ncc = 1;

if do_ncc
    clear
    do_ncc = 1;
end

fSize_axes = 9;
fSize_title = 12;

num_trials = 1;

% Number of corresponding regions
num_regions_eq = 20;
num_regions_neq = 1000;

cc_abs_mad = zeros(num_regions_eq, 1);

% 

% Image dimensions
region_height = 64;
region_width  = 64;
% % Grid step
gx_range = 0;
gx_step = 0;


% Center pixels
xc = fourier_zero(region_width);
yc = fourier_zero(region_height);

% Coordinates
x = (1 : region_width) - xc;
y = (1 : region_height) - yc;


% gx_range = 256;
% gx_step = 32;


gy_range = gx_range;
gy_step = gx_step;

% Random displacements
s_rand = 1;

%
%
%
%
% Bulk displacements (std dev)
sx_uniform_spread = region_width / 10;
sy_uniform_spread = 0;

%
%
x_dist = rectpuls(x, sx_uniform_spread);
y_dist = rectpuls(y, sy_uniform_spread);

dx_dist = y_dist' * x_dist;


% sx_lb =  -sx_uniform_spread/2;
% sx_ub =  -sx_uniform_spread/2;

sx_lb = 5;
sx_ub = 5;

% sy_lb = sx_lb;
% sy_ub = sx_ub;

sy_lb = 0;
sy_ub = 0;

sx_bulk_dist = (sx_ub - sx_lb) * rand(num_trials * num_regions_eq, 1) + sx_lb;
sy_bulk_dist = (sy_ub - sy_lb) * rand(num_trials * num_regions_eq, 1) + sy_lb;
% sx_bulk_dist = 5 * ones(num_trials, 1);

% sx_bulk_dist = linspace(1, 10, num_trials);

% Window size
window_fraction = 0.4 * [1, 1];

% Random displacements
sx_rand = s_rand;
sy_rand = s_rand;

% % Particle stuff

% Standard deviation of particle image diameters
d_std = 0;

% Mean particle diameter
d_mean = 3;
% d_mean = 1;
% Particle concentration in particles per pixel
particle_concentration = 6E-2;

% Image noise
noise_mean_fract = 1E-1;
noise_std_fract  = 5E-4;

% Particle positions buffer
x_buffer = -16;
y_buffer = -16;

% % % % %



% Window
g_win = gaussianWindowFilter([region_height, region_width], window_fraction, 'fraction');
g_win_ft = abs(fftshift(fft2(fftshift(g_win))));
% g_win = rect_lowpass_2D([region_height, region_width], [0.5, 0.5], 'fraction');

% Particle position max and min
x_min = 1 + x_buffer;
x_max = region_width  - x_buffer;
y_min = 1 + y_buffer;
y_max = region_width - y_buffer;

% Augmented size
aug_width = x_max - x_min + 1;
aug_height = y_max - y_min + 1;

% Compute the total number of particles
num_particles = round(particle_concentration * aug_height * aug_width);

% Allocate the particle shape
particle_shape_sum = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);

ncc_sum = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);

ncc_cur_max = zeros(num_regions_neq, 1);

% Make a test region to figure out the correct percent noise

max_val = 0;

for k = 1 : 10
% Particle diameters
dp_neq_01       = d_mean + d_std * randn(num_particles, 1);
% Particle intensities for the correlated images
particle_intensities_neq_01     = d_mean ./ dp_neq_01;

% Particle positions (image 1)
x_neq_01 = (x_max - x_min) * rand(num_particles, 1) + x_min;
y_neq_01 = (y_max - y_min) * rand(num_particles, 1) + y_min;

% % Generate the first image
region_neq_01 = generateParticleImage(region_height, region_width,...
    x_neq_01, y_neq_01, dp_neq_01, particle_intensities_neq_01);

max_val = max_val + max(region_neq_01(:));

end

% Average max value of the images before noise
max_val_mean = max_val / k;

% Standard deviation of the noise
noise_std = noise_std_fract * max_val_mean;

noise_mean = noise_mean_fract * max_val_mean;

% Allocate the CC sums
cc_abs_sum = zeros(region_height, region_width);
cc_abs_full_sum = zeros(region_height, region_width);

% Shake-the-box shfting coordinates
gx_shift = -1 * gx_range : gx_step : gx_range;
gy_shift = -1 * gx_range : gx_step : gx_range;

if isempty(gx_shift)
    gx_shift = 0;
end

if isempty(gy_shift)
    gy_shift = 0;
end

[GX, GY] = meshgrid(gx_shift, gy_shift);

for r = 1 : num_trials
    
    % Random mean displacements
    sx_bulk = sx_bulk_dist(r);
    sy_bulk = sy_bulk_dist(r);
    
     % Inform the user
    fprintf(1, 'Trial %d of %d, sx bulk = %0.2f\n', r, num_trials, sx_bulk);
    
    % Allocate running CC sum for the CCC
    cc_sum = zeros(region_height, region_width) + ...
        1i * zeros(region_height, region_width);
    
    ac_sum = zeros(region_height, region_width) + ...
        1i * zeros(region_height, region_width);
    
    % Allocate running CC sum for the CC minus the NCC
    cc_full_sum = zeros(region_height, region_width) + ...
        1i * zeros(region_height, region_width);
    
    
    cc_abs_cur_new = zeros(region_height, region_width);
    
    cc_test = zeros(region_height, region_width);

% Do the corresponding correlation
for k = 1 : num_regions_eq

    % Diffusion
    dx = sx_bulk + sx_rand * randn(num_particles, 1);
    dy = sy_bulk + sy_rand * randn(num_particles, 1);

     % Particle positions (image 1)
    x_01 = (x_max - x_min) * rand(num_particles, 1) + x_min;
    y_01 = (y_max - y_min) * rand(num_particles, 1) + y_min;
    
%     % Shearing
%     dx = sx_bulk + sx_uniform_spread * y_01 / region_height;
%     dy = sy_bulk + zeros(size(dx));
    
        
    % Uniform distribution
%     dx = (sx_uniform_spread) * rand(num_particles, 1) - sx_uniform_spread/2;
%     dy = (sy_uniform_spread) * rand(num_particles, 1) - sy_uniform_spread/2;
%     dy = zeros(size(dx));
    
    % Particle positions (image 2)
    % TLW
    x_02 = x_01 + dx;
    y_02 = y_01 + dy;

    % Particle diameters
    dp_eq       = d_mean + d_std * randn(num_particles, 1);

    % Particle intensities for the correlated images
    particle_intensities_eq     = d_mean ./ dp_eq;
    
    % Noise matrices
    noise_mat_full_01 = noise_mean + noise_std * ...
        randn(region_height + 2 * gx_range, region_width + 2 * gy_range);
    noise_mat_full_02 = noise_mean + noise_std * ...
        randn(region_height + 2 * gx_range, region_width + 2 * gy_range);
        
    % Shake-the-box
    for p = 1 : length(GX(:))
        
        x_eq_01 = x_01 + GX(p);
        x_eq_02 = x_02 + GX(p);
        
        y_eq_01 = y_01 + GY(p);
        y_eq_02 = y_02 + GY(p);
        
        % Noise coordinates
        noise_x_left = 1 + GX(p) + gx_range;
        noise_x_right = noise_x_left + region_width - 1;
        noise_y_top = 1 + GY(p) + gy_range;
        noise_y_bottom = noise_y_top + region_height - 1;

        % Extract noise
        noise_mat_eq_01 = noise_mat_full_01(...
            noise_y_top:noise_y_bottom, noise_x_left : noise_x_right);
        noise_mat_eq_02 = noise_mat_full_02(...
            noise_y_top:noise_y_bottom, noise_x_left : noise_x_right);

        % % Generate the images
        %
        % % Generate the first image
        region_eq_01 = (generateParticleImage(region_height, region_width,...
        x_eq_01, y_eq_01, dp_eq, particle_intensities_eq) + noise_mat_eq_01);
        %
        % Generate the second image
        region_eq_02 = (generateParticleImage(region_height, region_width,...
        x_eq_02, y_eq_02, dp_eq, particle_intensities_eq) + noise_mat_eq_02);
   
        % Transforms
        F1_eq = fft2(g_win .* (region_eq_01 - mean(region_eq_01(:))));
        F2_eq = fft2(g_win .* (region_eq_02 - mean(region_eq_02(:))));
    
        ac_01 = fftshift(F1_eq .* conj(F1_eq));
        ac_02 = fftshift(F2_eq .* conj(F2_eq));
        
        % Mean auto correlation
        ac_mean_cur = (ac_01 + ac_02) / 2;

        % Cross correlation
        cc_cur = fftshift(F1_eq .* conj(F2_eq));
        
        % Add the correlation to the sum
        cc_sum = cc_sum + cc_cur;
        
        % Auto correlation
        ac_sum = ac_sum + ac_mean_cur;
       
    end
    
end

% Add the sum of that trial to the sum-over-trials.
% This means we're taking the magnitude after summing.
cc_abs_sum = cc_abs_sum + (abs(cc_sum)).^2;

% keyboard;

end

% Square roots
cc_abs_sum_sqrt = sqrt(cc_abs_sum);
cc_abs_full_sum_sqrt = sqrt(cc_abs_full_sum);

surf(abs(cc_sum)); 
axis square;


% Fit the NCC
[A, sy, sx, YC, XC, B, ARRAY] = fit_gaussian_2D(cc_abs_sum_sqrt);

apc_filt = ARRAY - B;






