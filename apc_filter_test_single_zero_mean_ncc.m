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
s_rand = 0;

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


sx_lb =  -sx_uniform_spread/2;
sx_ub =  -sx_uniform_spread/2;

sy_lb = sx_lb;
sy_ub = sx_ub;

sy_lb = 0;
sy_ub = 0;

sx_bulk_dist = (sx_ub - sx_lb) * rand(num_trials * num_regions_eq, 1) + sx_lb;
sy_bulk_dist = (sy_ub - sy_lb) * rand(num_trials * num_regions_eq, 1) + sy_lb;
% sx_bulk_dist = 5 * ones(num_trials, 1);

% sx_bulk_dist = linspace(1, 10, num_trials);

% Window size
window_fraction = 1 * [1, 1];



% Random displacements
sx_rand = s_rand;
sy_rand = s_rand;

% % Particle stuff

% Standard deviation of particle image diameters
d_std = 0;

% Mean particle diameter
d_mean = 1 * sqrt(8);
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

for k = 1 : 1000
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

% Loop over the NCC
for k = 1 : num_regions_neq

    % Inform the user
    fprintf(1, 'NCC %d of %d\n', k, num_regions_neq);

    % Particle diameters
    dp_neq_01       = d_mean + d_std * randn(num_particles, 1);
    dp_neq_02       = d_mean + d_std * randn(num_particles, 1);

    % Particle intensities for the correlated images
    particle_intensities_neq_01     = d_mean ./ dp_neq_01;
    particle_intensities_neq_02     = d_mean ./ dp_neq_02;

    % Particle positions (image 1)
    x_neq_01 = (x_max - x_min) * rand(num_particles, 1) + x_min;
    y_neq_01 = (y_max - y_min) * rand(num_particles, 1) + y_min;

    % Particle positions (image 2)
    x_neq_02 = (x_max - x_min) * rand(num_particles, 1) + x_min;
    y_neq_02 = (y_max - y_min) * rand(num_particles, 1) + y_min;

    % Noise matrices
    noise_mat_neq_01      = noise_mean + noise_std * randn(region_height, region_width);
    noise_mat_neq_02      = noise_mean + noise_std * randn(region_height, region_width);

    % % Generate the images
    %
    % % Generate the first image
    region_neq_01 = (generateParticleImage(region_height, region_width,...
    x_neq_01, y_neq_01, dp_neq_01, particle_intensities_neq_01) + noise_mat_neq_01);
    %
    % Generate the second image
    region_neq_02 = (generateParticleImage(region_height, region_width,...
    x_neq_02, y_neq_02, dp_neq_02, particle_intensities_neq_02) + noise_mat_neq_02);

    % Transforms
    F1_neq = fft2(g_win .* (region_neq_01 - mean(region_neq_01(:))));
    F2_neq = fft2(g_win .* (region_neq_02 - mean(region_neq_02(:))));

    % Autocorrelations
    ac_01_neq = fftshift(F1_neq .* conj(F1_neq));
    ac_02_neq = fftshift(F2_neq .* conj(F2_neq));

    % Average auto correlation
    ac_mean = (ac_01_neq + ac_02_neq) / 2;

    % Cross correlation
    ncc_cur = fftshift(F1_neq .* conj(F2_neq));
    
    ncc_cur_max(k) = max(real(ncc_cur(:)));
    
    % Particle shape current
%     particle_shape_cur = ac_mean - ncc_cur;
    particle_shape_cur = abs(ac_mean);
    
    % Sum
    ncc_sum = ncc_sum + ncc_cur;
   
    % Sum particle shape
    particle_shape_sum = particle_shape_sum + particle_shape_cur;

end

% Normalize the particle so its maximum is one.
% Note that we don't want to force the minimum to zero
% Because that will cause the edges to blow up
% whenever we divide by it.
particle_shape_norm = real(particle_shape_sum ./ max(particle_shape_sum(:)));

% Fit a Gaussian to this;
% Fit the NCC
[A, ~, ~, ~, ~, B, particle_shape_fit] = fit_gaussian_2D(particle_shape_norm);

% % This line works!!!!
ncc_div = ncc_sum ./ particle_shape_sum;

% Normalize the NCC so that its maximum is one.
ncc_norm = ncc_div ./ max(ncc_div(:));

% Full NCC, which is not divided by the particle shape
ncc_full_norm = ncc_sum ./ max(real(ncc_sum(:)));

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
    
    % Allocate running CC sum for the CC minus the NCC
    cc_full_sum = zeros(region_height, region_width) + ...
        1i * zeros(region_height, region_width);
    
    
    cc_abs_cur_new = zeros(region_height, region_width);
    
    cc_test = zeros(region_height, region_width);

% Do the corresponding correlation
for k = 1 : num_regions_eq
    
%     dx = sx_bulk + sx_rand * randn(num_particles, 1);
%     dy = sy_bulk + sy_rand * randn(num_particles, 1);

     % Particle positions (image 1)
    x_01 = (x_max - x_min) * rand(num_particles, 1) + x_min;
    y_01 = (y_max - y_min) * rand(num_particles, 1) + y_min;
    
    % Shearing
    dx = sx_bulk + sx_uniform_spread * y_01 / region_height;
    dy = sy_bulk + zeros(size(dx));
    
        
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

        % Cross correlation
        cc_cur = fftshift(F1_eq .* conj(F2_eq));
        
        % Full CC sum
        cc_full_sum = cc_full_sum + cc_cur;

        % Divide by the particle shape
        cc_eq = cc_cur ./ particle_shape_norm;

        % Ensemble corresponding correlation
        % (divided by particle shape, minus NCC)
        cc_sum = cc_sum + cc_eq;
       
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


% Fit the NCC
[A, sy, sx, YC, XC, B, ARRAY] = fit_gaussian_2D(cc_abs_sum_sqrt);

% Fit the particle shape
[~, ~, ~, ~, ~, B_p, particle_shape_peak_fit] = fit_gaussian_2D(particle_shape_norm);

p_n = particle_shape_norm - B_p;

apc_filt = ARRAY - B;

filt_full = apc_filt ./ max(apc_filt(:)) .* p_n ./ max(p_n(:));

% Fit the NCC * Particle
% [A, sy, sx, YC, XC, B, ARRAY] = fit_gaussian_2D(cc_abs_full_sum_sqrt);

fprintf(1, 'stdx = %0.2f, stdy = %0.2f\n', sx, sy);

xv = (1 : region_width);
yv = (1 : region_height);

[X, Y] = meshgrid(xv, yv);

Z = (exp(-(X - xc).^2 / (2 * sx^2)) .* exp(-(Y - yc).^2 / (2 * sy^2)));

cc_abs_shift = abs(cc_abs_sum_sqrt - B);
cc_abs_full_shift = abs(cc_abs_full_sum_sqrt - B);

apc_filt_full = particle_shape_norm .* cc_abs_shift;


apc_filt_full_test = particle_shape_sum .* cc_abs_sum_sqrt;

apc_filt_norm_test = apc_filt_full_test ./ max(apc_filt_full_test(:));

apc_filt_sub = apc_filt_full - min(apc_filt_full(:));
apc_filt_norm = apc_filt_sub ./ max(apc_filt_sub(:));

[~, ~, ~, ~, ~, apc_fit_offset, apc_fit_full] = fit_gaussian_2D(apc_filt_norm);
apc_fit_sub = apc_fit_full - min(apc_fit_full(:));
apc_fit_norm = apc_fit_sub ./ max(apc_fit_sub(:));


cam_proj = 'orthographic';

v = [-34, 9];

lw = 1E-5;

figure(1);
% Particle shape
subplot(2, 2, 1);
set(gca, 'camerapositionmode', 'manual');
surf(particle_shape_norm ./ max(particle_shape_norm(:)), 'linewidth', lw);
camproj(cam_proj);
% surf(real(cc_sum) ./ max(real(cc_sum(:))));
xlim([1, region_width]);
ylim([1, region_height]);
zlim(1.1 * [0, 1]);
set(gca, 'view', v);
axis square;
axis vis3d
box on;
% set(gca, 'view', [0, 0]);
set(gca, 'FontSize', fSize_axes);
title(sprintf(...
    '$A \\left(k \\right) \\, , \\overline{d_p} = %0.1f$', d_mean),...
    'interpreter', 'latex', ...
    'fontSize', fSize_title);
set(gca, 'ztick', [0 : 0.25 : 1]);
set(gca, 'xtick', [16, 32, 48, 64]);
set(gca, 'ytick', [16, 32, 48, 64]);
set(gca, 'xticklabel', {''})
set(gca, 'yticklabel', {''})
set(gca, 'zticklabel', {''})
% xlabel('$k / \pi$', 'interpreter', 'latex', 'FontSize', fSize_axes);
% ylabel('$m / \pi$', 'interpreter', 'latex', 'FontSize', fSize_axes);


subplot(2, 2, 2);
% surf(cc_abs_shift ./ max(cc_abs_shift(:)), 'linewidth', lw);
surf(cc_abs_sum_sqrt ./ max(cc_abs_sum_sqrt(:)), 'linewidth', lw);
camproj(cam_proj);
xlim([1, region_width]);
ylim([1, region_height]);
zlim(1.1 * [0, 1]);
axis square;
box on;
set(gca, 'view', v);
axis square;
axis vis3d
set(gca, 'FontSize', fSize_axes);
title(sprintf(...
    '$P_\\mu \\left(k\\right) \\,\\, , \\sigma_{\\Delta x} = %0.1f$', s_rand),...
    'interpreter', 'latex', ...
    'fontSize', fSize_title);
set(gca, 'ztick', [0 : 0.25 : 1]);
set(gca, 'xtick', [16, 32, 48, 64]);
set(gca, 'ytick', [16, 32, 48, 64]);
set(gca, 'xticklabel', {''})
set(gca, 'yticklabel', {''})
set(gca, 'zticklabel', {''})
set(gca, 'view', [0, 0]);
% xlabel('$k / \pi$', 'interpreter', 'latex', 'FontSize', fSize_axes);
% ylabel('$m / \pi$', 'interpreter', 'latex', 'FontSize', fSize_axes);



subplot(2, 2, 3);
% surf(apc_filt_norm, 'linewidth', lw);
surf(apc_filt_norm_test, 'linewidth', lw);
xlim([1, region_width]);
ylim([1, region_height]);
zlim([0, 1.1]);
axis square;
box on;
set(gca, 'view', v);
axis square;
axis vis3d
set(gca, 'view', v);
set(gca, 'FontSize', fSize_axes);
title('$A\left(k\right) P_\mu \left(k\right)$', 'interpreter', 'latex', ...
    'FontSize', fSize_title);
set(gca, 'ztick', [0 : 0.25 : 1]);
set(gca, 'xtick', [16, 32, 48, 64]);
set(gca, 'ytick', [16, 32, 48, 64]);
set(gca, 'xticklabel', {''})
set(gca, 'yticklabel', {''})
set(gca, 'zticklabel', {''})
% xlabel('$k / \pi$', 'interpreter', 'latex', 'FontSize', fSize_axes);
% ylabel('$m / \pi$', 'interpreter', 'latex', 'FontSize', fSize_axes);



subplot(2, 2, 4);
surf(apc_fit_norm, 'linewidth', lw);
xlim([1, region_width]);
ylim([1, region_height]);
zlim([0, 1.1]);
axis square;
box on;
set(gca, 'view', v);
axis square;
axis vis3d
set(gca, 'view', v);
set(gca, 'FontSize', fSize_axes);
title('$\textrm{Gaussian Fit to } A\left(k\right) P_\mu \left(k\right)$', 'interpreter', 'latex', ...
    'FontSize', fSize_title);
set(gca, 'ztick', [0 : 0.25 : 1]);
set(gca, 'xtick', [16, 32, 48, 64]);
set(gca, 'ytick', [16, 32, 48, 64]);
set(gca, 'xticklabel', {''})
set(gca, 'yticklabel', {''})
set(gca, 'zticklabel', {''})
% xlabel('$k / \pi$', 'interpreter', 'latex', 'FontSize', fSize_axes);
% ylabel('$m / \pi$', 'interpreter', 'latex', 'FontSize', fSize_axes);

colormap parula;


% 
% figure(2);
% subplot(1, 3, 1);
% surf(apc_fit_norm, 'linewidth', lw);
% xlim([1, region_width]);
% ylim([1, region_height]);
% zlim([0, 1.1]);
% axis square;
% box on;
% set(gca, 'view', v);
% axis square;
% axis vis3d
% set(gca, 'view', v);
% set(gca, 'FontSize', fSize_axes);
% title(sprintf('$\\textrm{APC} \\left(k\\right) , \\,\\, \\overline{d_p} = %0.1f , \\, \\sigma_{\\Delta x} = %0.1f$', d_mean, s_rand), 'interpreter', 'latex', ...
%     'FontSize', ft);
% set(gca, 'ztick', [0 : 0.25 : 1]);
% set(gca, 'xtick', [16, 32, 48, 64]);
% set(gca, 'ytick', [16, 32, 48, 64]);
% set(gca, 'xticklabel', {''})
% set(gca, 'yticklabel', {''})
% set(gca, 'zticklabel', {''})
% 
% subplot(1, 3, 2);
% R = spectralEnergyFilter(region_height, region_width, d_mean);
% surf(R ./ max(R(:)), 'linewidth', lw);
% xlim([1, region_width]);
% ylim([1, region_height]);
% zlim([0, 1.1]);
% axis square;
% box on;
% set(gca, 'view', v);
% axis square;
% axis vis3d
% set(gca, 'view', v);
% set(gca, 'FontSize', fSize_axes);
% title(sprintf('$\\textrm{Analytical RPC Filter, } \\, \\, \\overline{d_p} = %0.1f$', d_mean), 'interpreter', 'latex', ...
%     'FontSize', ft);
% set(gca, 'ztick', [0 : 0.25 : 1]);
% set(gca, 'xtick', [16, 32, 48, 64]);
% set(gca, 'ytick', [16, 32, 48, 64]);
% set(gca, 'xticklabel', {''})
% set(gca, 'yticklabel', {''})
% set(gca, 'zticklabel', {''})
% 
% a = abs((R ./ max(R(:)) - apc_fit_norm));
% subplot(1, 3, 3);
% surf(a)
% xlim([1, region_width]);
% ylim([1, region_height]);
% % zlim([0, 1.1]);
% axis square;
% box on;
% set(gca, 'view', v);
% axis square;
% axis vis3d
% set(gca, 'view', v);
% set(gca, 'FontSize', fSize_axes);
% title('$ \vert \mathrm{APC}\left(k\right) - \mathrm{RPC} \left(k\right) \vert$', 'interpreter', 'latex', ...
%     'FontSize', ft);
% set(gca, 'xtick', [16, 32, 48, 64]);
% set(gca, 'ytick', [16, 32, 48, 64]);
% set(gca, 'xticklabel', {''})
% set(gca, 'yticklabel', {''})
% zl = zlim;
% 
% z_tick = linspace(0, max(zl), 5);
% set(gca, 'ztick', z_tick);
% zlim([0, 1.1 * max(z_tick)]);
% 
% % set(gca, 'zticklabel', {''})
% 
% 





plot_dir = '~/Desktop/apc_plots/quads';

file_name = sprintf('apc_full_filt_h%d_w%d_sr_%0.1f.png', ...
    region_height, region_width, s_rand);
file_path = fullfile(plot_dir, file_name);

% 
% do_print = input('Print? [y/N]\n', 's');
% switch lower(do_print)
%     case 'y'
%         print(1, '-dpng', '-r300', file_path);
%     otherwise
% end
% 
%         








