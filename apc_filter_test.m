% addpaths;
% clear;

clear


num_trials = 100;
% 

% Image dimensions
region_height = 64;
region_width = 64;

g = zeros(region_height, region_width, num_trials);


sx_lb = -10;
sx_ub = 10;

sx_bulk_dist = (sx_ub - sx_lb) * rand(num_trials, 1) + sx_lb;
sx_bulk_dist = 5 * ones(num_trials, 1);

% sx_bulk_dist = linspace(1, 10, num_trials);


fSize = 12;

do_ncc = 1;

% Number of corresponding regions
num_regions_eq = 10;
num_regions_neq = 1000;


% Window size
window_fraction = 0.5 * [1, 1];


% Bulk displacements (std dev)
sx_bulk_std = 0;
sy_bulk_std = 0;

% Random displacements
s_rand = 1;
sx_rand = s_rand;
sy_rand = s_rand;

% % Particle stuff

% Standard deviation of particle image diameters
d_std = 0.5;

% Mean particle diameter
d_mean = 1 * sqrt(8) ;

% Particle concentration in particles per pixel
particle_concentration = 1E-2;

% Image noise
noise_std = 1E-2;

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


% Allocate the particle shape
particle_shape_sum = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);

ncc_sum = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);

ncc_cur_max = zeros(num_regions_neq, 1);

 
% Loop over the NCC
for k = 1 : num_regions_neq

    % Inform the user
%     fprintf(1, 'NCC %d of %d\n', k, num_regions_neq);

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

    % Transforms
    F1_neq = fft2(g_win .* region_neq_01);
    F2_neq = fft2(g_win .* region_neq_02);

    % Autocorrelations
    ac_01_neq = fftshift(F1_neq .* conj(F1_neq));
    ac_02_neq = fftshift(F2_neq .* conj(F2_neq));

    % Average auto correlation
    ac_mean = (ac_01_neq + ac_02_neq) / 2;

    % Cross correlation
    ncc_cur = fftshift(F1_neq .* conj(F2_neq));
    
    ncc_cur_max(k) = max(real(ncc_cur(:)));
    
    % Particle shape current
    particle_shape_cur = ac_mean - ncc_cur;
    
    % Sum
    ncc_sum = ncc_sum + ncc_cur;
  
    
    % Sum particle shape
    particle_shape_sum = particle_shape_sum + particle_shape_cur;
    
%     % Plot
%     surf(abs(particle_shape_sum) ./ max(abs(particle_shape_sum(:))));
%     xlim([1, region_width]);
%     ylim([1, region_height]);
%     
%     pause(0.1);



end



% Normalize the particle shape by the number of regions
particle_shape_norm = particle_shape_sum / (num_regions_neq - num_regions_eq);

% This should be the average of:
% ncc_norm = A * N^2 * Pu
% ncc_norm =  ncc_sum / num_regions_neq;
ncc_norm = ncc_sum ./ max(ncc_sum(:));

cc_sum = zeros(region_height, region_width) + ...
    1i * zeros(region_height, region_width);

cc_abs_sum = zeros(region_height, region_width);
cc_cur_max = zeros(num_regions_eq, 1);


% sx_bulk = sx_bulk_mean;
% sy_bulk = 0;



for r = 1 : num_trials
    
    fprintf(1, 'Trial %d of %d\n', r, num_trials);

cc_sum = zeros(region_height, region_width) + ...
    1i * zeros(region_height, region_width);

% Random mean displacements
    sx_bulk = sx_bulk_dist(r);
    sy_bulk = 0;

% Do the corresponding correlation
for k = 1 : num_regions_eq

    % Inform the user
%         fprintf(1, 'CCC %d of %d\n', k, num_regions_eq);

    % Particle diameters
    dp_eq       = d_mean + d_std * randn(num_particles, 1);

    % Particle intensities for the correlated images
    particle_intensities_eq     = d_mean ./ dp_eq;

    % Particle positions (image 1)
    x_eq_01 = (x_max - x_min) * rand(num_particles, 1) + x_min;
    y_eq_01 = (y_max - y_min) * rand(num_particles, 1) + y_min;

    % Particle positions (image 2)
    x_eq_02 = x_eq_01 + sx_bulk + sx_rand * randn(num_particles, 1);
    y_eq_02 = y_eq_01 + sy_bulk + sy_rand * randn(num_particles, 1);

    % Noise matrices
    noise_mat_eq_01      = noise_std * randn(region_height, region_width);
    noise_mat_eq_02      = noise_std * randn(region_height, region_width);

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
    F1_eq = fft2(g_win .* region_eq_01);
    F2_eq = fft2(g_win .* region_eq_02);

    % Cross correlation
    cc_cur = fftshift(F1_eq .* conj(F2_eq));
    
    N = sqrt(max(real(cc_cur(:))));
    

    cc_cur_max(k) = max(real(cc_cur(:)));

    % Ensemble corresponding correlation
    cc_sum = cc_sum + cc_cur ./ particle_shape_norm - (N^2 - N) * (ncc_norm ./ particle_shape_norm);

    cc_abs_sum = cc_abs_sum + (abs(cc_sum)).^2;

end
    
    


cc_abs_sum_sqrt = sqrt(cc_abs_sum);

g(:, :, r) = cc_abs_sum;

end

for k = 1 : size(g, 3);
    
    g_plot = g(:, :, k);
    
    surf(g_plot ./ max(g_plot(:))); 
    drawnow; 
end;

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

















