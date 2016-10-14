addpaths;
% clear;

do_ncc = 1;

if do_ncc
    clear
    do_ncc = 1;
end

fSize = 12;


num_trials = 1;

% Number of corresponding regions
num_regions_eq  = 1;
num_regions_neq = 1;

% 

% Image dimensions
region_height = 64;
region_width = 64;

% Grid step
gx_range = 0;
gx_step = 0;

gy_range = gx_range;
gy_step = gx_step;

% Random displacements
s_rand = 0;

sx_mean = 5;
sy_mean = 5;

sx_lb = 5;
sx_ub = 5;

sx_bulk_dist = (sx_ub - sx_lb) * rand(num_trials * num_regions_eq, 1) + sx_lb;
% sx_bulk_dist = 5 * ones(num_trials, 1);

% sx_bulk_dist = linspace(1, 10, num_trials);


% Window size
window_fraction = 0.5 * [1, 1];

% Bulk displacements (std dev)
sx_bulk_std = 0;
sy_bulk_std = 0;

% Random displacements
sx_rand = s_rand;
sy_rand = s_rand;

% % Particle stuff

% Standard deviation of particle image diameters
d_std = 0.0;

% Mean particle diameter
d_mean = 3;

% Particle concentration in particles per pixel
particle_concentration = 2E-2;

% Image noise
noise_std = 0E-2;
% noise_std = 0;

% Particle positions buffer
x_buffer = -100;
y_buffer = -100;

% % % % %

% Center pixels
xc = (region_width  + 1) / 2 + 0.5 * (1 - mod(region_width,  2));
yc = (region_height + 1) / 2 + 0.5 * (1 - mod(region_height, 2));


% Window
g_win = gaussianWindowFilter([region_height, region_width], window_fraction, 'fraction');
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

if do_ncc

% Allocate the particle shape
particle_shape_sum = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);

ncc_sum = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);

ac_mean_sum_neq = zeros(region_height, region_width) + ...
    1i * zeros(region_height, region_width);

ncc_cur_max = zeros(num_regions_neq, 1);

 
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
    noise_mat_neq_01      = noise_std * randn(region_height, region_width);
    noise_mat_neq_02      = noise_std * randn(region_height, region_width);

    % % Generate the images
    %
    % % Generate the first image
    region_neq_01 = generateParticleImage(region_height, region_width,...
    x_neq_01, y_neq_01, dp_neq_01, particle_intensities_neq_01) + noise_mat_neq_01;
    %
    % Generate the second image
    region_neq_02 = generateParticleImage(region_height, region_width,...
    x_neq_02, y_neq_02, dp_neq_02, particle_intensities_neq_02) + noise_mat_neq_02;

    region_neq_02 = circshift(region_neq_01, [0, sx_mean]);

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
    
    ac_mean_sum_neq = ac_mean_sum_neq + ac_mean;
    
    % Particle shape current
    particle_shape_cur = ac_mean - ncc_cur;
    
    % Sum
    ncc_sum = ncc_sum + ncc_cur;
   
    % Sum particle shape
    particle_shape_sum = particle_shape_sum + particle_shape_cur;

end

AN = (ac_mean_sum_neq - ncc_sum);
ANP = (ac_mean_sum_neq) - AN;

% Normalize the particle so its maximum is one.
% Note that we don't want to force the minimum to zero
% Because that will cause the edges to blow up
% whenever we divide by it.
particle_shape_norm = particle_shape_sum ./ max(particle_shape_sum(:));

% % This line works!!!!
ncc_div = ncc_sum ./ particle_shape_sum;

% Normalize the NCC so that its maximum is one.
ncc_norm = ncc_div ./ max(ncc_div(:));

% 
% [AMPLITUDE, STD_DEV_Y, STD_DEV_X, YC, XC, B, ncc_peak_fit] = fit_gaussian_2D(real(ncc_norm), 5);

% ncc_fit_sub = ncc_peak_fit - B;

% ncc_fit_norm = ncc_fit_sub ./ max(ncc_fit_sub(:));

% this worked when it was uncommented. Dealt with noise.
% Uncomment this.
ncc_fit_norm = ncc_norm;

% ncc_norm = (ncc_sum ./ particle_shape_sum) ./ num_regions_neq;

% % % % This line works
% ncc_norm = ncc_sum ./ max(ncc_sum(:));

end

cc_sum = zeros(region_height, region_width) + ...
    1i * zeros(region_height, region_width);

cc_abs_sum = zeros(region_height, region_width);
cc_cur_max = zeros(num_regions_eq, 1);
cc_abs_sum_02 = zeros(region_height, region_width);

% sx_bulk = sx_bulk_mean;
% sy_bulk = 0;

cc_peak_cols = (xc - 3) : (xc + 3);
cc_peak_rows = (yc - 3) : (yc + 3);

sy_bulk = 0;


gx_shift = -1 * gx_range : gx_step : gx_range;
gy_shift = -1 * gx_range : gx_step : gx_range;

[GX, GY] = meshgrid(gx_shift, gy_shift);

if isempty(gx_shift)
    GX = 0;
end

if isempty(gy_shift)
    GY = 0;
end


for r = 1 : num_trials
    
        % Random mean displacements
    sx_bulk = sx_bulk_dist(r);
    
    fprintf(1, 'Trial %d of %d, sx bulk = %0.2f\n', r, num_trials, sx_bulk);

    cc_sum = zeros(region_height, region_width) + ...
        1i * zeros(region_height, region_width);

    cc_sum_full = zeros(region_height, region_width) + ...
        1i * zeros(region_height, region_width);
    
    
    ac_sum_full = zeros(region_height, region_width) + ...
        1i * zeros(region_height, region_width);
    
    
% Do the corresponding correlation
for k = 1 : num_regions_eq
    
     % Particle positions (image 1)
    x_01 = (x_max - x_min) * rand(num_particles, 1) + x_min;
    y_01 = (y_max - y_min) * rand(num_particles, 1) + y_min;

    % Particle positions (image 2)
    x_02 = x_01 + sx_bulk + sx_rand * randn(num_particles, 1);
    y_02 = y_01 + sy_bulk + sy_rand * randn(num_particles, 1);
    
    % Particle diameters
    dp_eq       = d_mean + d_std * randn(num_particles, 1);

    % Particle intensities for the correlated images
    particle_intensities_eq     = d_mean ./ dp_eq;
    
    % Noise matrices
    noise_mat_full_01 = noise_std * ...
        randn(region_height + 2 * gx_range, region_width + 2 * gy_range);
    noise_mat_full_02 = noise_std * ...
        randn(region_height + 2 * gx_range, region_width + 2 * gy_range);
    
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
        region_eq_01 = generateParticleImage(region_height, region_width,...
        x_eq_01, y_eq_01, dp_eq, particle_intensities_eq) + noise_mat_eq_01;
        %
        % Generate the second image
        region_eq_02 = generateParticleImage(region_height, region_width,...
        x_eq_02, y_eq_02, dp_eq, particle_intensities_eq) + noise_mat_eq_02;

        % Transforms
        F1_eq = fftshift(fft2(g_win .* (region_eq_01 - mean(region_eq_01(:)))));
        F2_eq = fftshift(fft2(g_win .* (region_eq_02 - mean(region_eq_02(:)))));
        
        % Cross correlation
        cc_cur = (F1_eq .* conj(F2_eq));
        
        % Auto correlation
        ac_01 = F1_eq .* conj(F1_eq);
        ac_02 = F2_eq .* conj(F2_eq);
        ac_mean = (ac_01 + ac_02) / 2;
        
        
        % Full sums
        cc_sum_full = cc_sum_full + cc_cur;
        ac_sum_full = ac_sum_full + ac_mean;

        % Max value
        cc_div = cc_cur ./ particle_shape_norm;

        % Effective number of particles
        N = sqrt(max(abs(cc_div(:))));
        N = sqrt(max(abs(cc_cur(:))));

        % Real and imaginary parts of corrected CC
        cc_eq_real = real(cc_cur) ./ real(particle_shape_norm) - real((N^2 - 1 * N) * (ncc_fit_norm));
        cc_eq_imag = imag(cc_cur) ./ real(particle_shape_norm);

        % Complex corrected CC
        cc_eq = cc_eq_real + 1i * cc_eq_imag;

    %     cc_eq(cc_peak_rows, cc_peak_cols) = 0;
        % Ensemble corresponding correlation
        % This line works
        cc_sum = cc_sum + cc_eq;

        cc_abs_sum_02 = cc_abs_sum_02 + abs(cc_eq).^2;

    end
    
end

% Add the sum of that trial to the sum-over-trials.
% This means we're taking the magnitude after summing.
cc_abs_sum = cc_abs_sum + (abs(cc_sum)).^2;




end

% Square root
cc_abs_sum_sqrt = sqrt(cc_abs_sum);


% for k = 1 : 1 : size(g, 3);
%     
%     g_plot = g(:, :, k);
%     
%     surf(g_plot ./ max(g_plot(:))); 
%     drawnow; 
% end;

[A, sy, sx, B, ARRAY] = fit_gaussian_2D(cc_abs_sum_sqrt);

fprintf(1, 'stdx = %0.2f, stdy = %0.2f\n', sx, sy);

xv = (1 : region_width);
yv = (1 : region_height);

[X, Y] = meshgrid(xv, yv);

Z = (exp(-(X - xc).^2 / (2 * sx^2)) .* exp(-(Y - yc).^2 / (2 * sy^2)));

cc_abs_shift = abs(cc_abs_sum_sqrt - B);

cc_sum_real = real(cc_sum_full);
ac_sum_real = real(ac_sum_full);

cc_ifft = fftshift(abs(ifft2(fftshift(cc_sum_full))));

Neff = sqrt(max(ac_sum_real(:)));

% G = real(fftshift(fft2(fftshift(g_win))));
% G2 = G .* G;
% 
% Gn = G2 ./ max(G2(:));
% 
% G_full = Gn * Neff * (Neff - 1.5);

cc_sum_sub = cc_sum_full - Neff^2 / (Neff^2 - Neff) * ANP;

cc_sub_real = real(cc_sum_sub);


subplot(1, 3, 1);
imagesc(real(cc_sum_full) ./ max(real(cc_sum_full(:))));
% surf(real(cc_sum) ./ max(real(cc_sum(:))));
xlim([1, region_width]);
ylim([1, region_height]);
% zlim(1.1 * [0, 1]);
axis image;
% set(gca, 'view', [0, 0]);

subplot(1, 3, 2);
imagesc(real(cc_sub_real));
% surf(real(cc_sum) ./ max(real(cc_sum(:))));
xlim([1, region_width]);
ylim([1, region_height]);
% zlim(1.1 * [0, 1]);
axis image;
% set(gca, 'view', [0, 0]);

subplot(1, 3, 3);
imagesc(abs(cc_sum_sub));
% surf(real(cc_sum) ./ max(real(cc_sum(:))));
xlim([1, region_width]);
ylim([1, region_height]);
% zlim(1.1 * [0, 1]);
axis image;
% set(gca, 'view', [0, 0]);


% % Plot
% subplot(1, 2, 1);
% surf(real(cc_sum) ./ max(real(cc_sum(:))));
% set(gca, 'view', [0, 0]);
% axis square
% xlim([1, region_width]);
% ylim([1, region_height]);
% zlim(1.05 * [-1, 1]);
% title({'Ensemble CCC' , sprintf('$\\Delta s = %0.1f, \\sigma_s = %0.1f$', ...
%     sx_bulk_mean, sx_rand)}, 'interpreter', 'latex', 'fontsize', 16);
% 
% 
% % Plot
% subplot(1, 2, 2);
% surf(real(cc_abs_sum_sqrt) ./ max(cc_abs_sum_sqrt(:)));
% set(gca, 'view', [0, 0]);
% axis square
% xlim([1, region_width]);
% ylim([1, region_height]);
% zlim(1.05 * [0, 1]);
% title({'Ensemble CCC mag' , sprintf('$\\Delta s = %0.1f, \\sigma_s = %0.1f$', ...
%     sx_bulk_mean, sx_rand)}, 'interpreter', 'latex', 'fontsize', 16);
% 
% print(1, '-dpng', '-r200', sprintf('~/Desktop/figures_04/fig_sx_%0.2f.png', sx_bulk));

















