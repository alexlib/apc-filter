% addpaths;


fSize_axes = 20;
fSize_title = 30;
fSize_labels = 30;
fSize_legend = 16;

sx_mean = 8;
sy_mean = 0;
sz_mean = 0;

diffusion_list = [0, 2, 5];

% diffusion_list = 0;

num_diffusions = length(diffusion_list);

num_trials = 1;

% Number of corresponding regions
num_regions_eq = 1000;

% Image dimensions
region_height = 128;
region_width  = 128;

% % Grid step
gx_range = 0;
gx_step = 0;

% Center pixels
xc = fourier_zero(region_width);
yc = fourier_zero(region_height);

% Coordinates
xv = (1 : region_width) - xc;
yv = (1 : region_height) - yc;

% Grid
[x, y] = meshgrid(xv, yv);


% Window size
window_fraction = 0.5 * [1, 1];


% % Particle stuff

% Standard deviation of particle image diameters
d_std = 0;

% Mean particle diameter
d_mean = 3;
% d_mean = 1;
% Particle concentration in particles per pixel
particle_concentration = 2E-2;

% Image noise
noise_mean_fract = 0E-2;
noise_std_fract  = 5E-2;
% noise_mean_fract = 0;
% noise_std_fract  = 0;

% Particle positions buffer
x_buffer = -16;
y_buffer = -16;


% % % % %

% RPC filter
rpc_filter = spectralEnergyFilter(region_height, region_width, d_mean);

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

% Magnification in microns per pixel
mag_um_pix = 7.5;

% Light sheet thickness in microns = 4 * std dev
sheet_thickness_microns = 0.5E3;

% Light sheet thickness in pixels
sheet_thickness_pixels = sheet_thickness_microns / mag_um_pix;

% Light sheet standard deviation in pixels
sheet_std_dev_pix = sheet_thickness_pixels / 4;

% Constrain the z coordinates
z_min = -1 * sheet_thickness_pixels / 2;
z_max =      sheet_thickness_pixels / 2;

% Line width
lw = 3;

xt = linspace(1, region_width, 11);
xtl = cell(11, 1);
xtl{1} = '-1';
xtl{end} = '1';
xtl{6} = '0';

ax_fract = 1;
dy_subplot = 20;
dx_subplot = 40;
height_fract = 2;

close all;
f = figure;
f_pos = [-1244         497         972         694];
set(f, 'position', f_pos);


for n = 1 : num_diffusions
    
    P1 = sub2ind([num_diffusions, 2], n, 1);
    P2 = sub2ind([num_diffusions, 2], n, 2);
    
    ax{2 * n - 1} = subplot(2, num_diffusions, P1, 'parent', f);
    set(gca, 'units', 'pixels');
    
    ax{2 * n} = subplot(2, num_diffusions, P2, 'parent', f);
    set(gca, 'units', 'pixels');

end

% Allocate running CC sum for the CCC
cc_full_sum = zeros(region_height, region_width, num_diffusions) + ...
    1i * zeros(region_height, region_width, num_diffusions);

cc_rows = zeros(num_regions_eq, region_width, num_diffusions) + ...
    1i * zeros(num_regions_eq, region_width, num_diffusions);


parfor n = 1 : num_diffusions
 

% Random displacements
s_rand = diffusion_list(n);
sx_rand = s_rand;
sy_rand = s_rand;
sz_rand = s_rand;


% Allocate running CC sum for the CC minus the NCC
cc_sum = zeros(region_height, region_width) + ...
    1i * zeros(region_height, region_width);
  
% Do the corresponding correlation
for k = 1 : num_regions_eq
    
    fprintf(1, 'On %d of %d...\n', k, num_regions_eq);
    
    dx = sx_mean + sx_rand * randn(num_particles, 1);
    dy = sy_mean + sy_rand * randn(num_particles, 1);
    dz = sz_mean + sz_rand * randn(num_particles, 1);

    % Particle positions (image 1)
    x_01 = (x_max - x_min) * rand(num_particles, 1) + x_min;
    y_01 = (y_max - y_min) * rand(num_particles, 1) + y_min;
    z_01 = (z_max - z_min) * rand(num_particles, 1) + z_min; 

    % Particle positions (image 2)
    % TLW
    x_02 = x_01 + dx;
    y_02 = y_01 + dy;
    z_02 = z_01 + dz;

    % Particle diameters
    dp       = d_mean + d_std * randn(num_particles, 1);

    % Particle intensities for the correlated images
    particle_intensities_01 = exp(-z_01.^2 / (2 * sheet_std_dev_pix^2));
    particle_intensities_02 = exp(-z_02.^2 / (2 * sheet_std_dev_pix^2));
    

    % % Generate the images
    %
    % % Generate the first image
    region_eq_clean_01 = generateParticleImage(region_height, region_width,...
    x_01, y_01, dp, particle_intensities_01);
    %
    % Generate the second image
    region_eq_clean_02 = generateParticleImage(region_height, region_width,...
    x_02, y_02, dp, particle_intensities_02);

    % Noise levels
    max_val = max([region_eq_clean_01(:); region_eq_clean_02(:)]);
    
    % Noise levels
    noise_mean = max_val * noise_mean_fract * 0.99;
    noise_std_dev = max_val * noise_std_fract * 0.99;

    % Noise levels
    noise_mat_01 = noise_mean + noise_std_dev * randn(region_height, region_width);
    noise_mat_02 = noise_mean + noise_std_dev * randn(region_height, region_width);
    
    % Add the noise to the regions
    region_eq_01 = region_eq_clean_01 + noise_mat_01;
    region_eq_02 = region_eq_clean_02 + noise_mat_02;

    % Transforms
    F1_eq = fft2(g_win .* (region_eq_01 - mean(region_eq_01(:))));
    F2_eq = fft2(g_win .* (region_eq_02 - mean(region_eq_02(:))));
    
    % Cross correlation
    cc_cur = fftshift(F1_eq .* conj(F2_eq));

    % Full CC sum
    cc_rows(k, :, n) = cc_cur(yc + 2, :);

end

end


cc_rows_sum = squeeze(sum(cc_rows, 1));
cc_rows_mean = cc_rows_sum ./ num_regions_eq;

ph_rows = angle(cc_rows);



cc_abs = abs(cc_rows_sum);
cc_rows_std_dev = squeeze(std(ph_rows, [], 1)) ./ cc_abs;

z = zeros(region_width, 1);
xv = (x(1, :))';

for n = 1 : num_diffusions
   
cc_std_dev = cc_rows_std_dev(:, n);
cc_sum_vect = cc_rows_sum(:, n);
cc_sum_vect_real = real(cc_sum_vect);
cc_sum_abs_vect = abs(cc_sum_vect);
cc_sum_phase = real(phaseOnlyFilter(cc_sum_vect));

cc_sum_abs_norm = cc_sum_abs_vect ./ max(cc_sum_abs_vect);

% subplot(2, num_diffusions, P1, 'parent', f);
axes(ax{2 * n - 1})
plot(cc_sum_vect_real ./ max(cc_sum_vect_real), '-k', 'linewidth', lw);
hold on;
plot(cc_sum_abs_norm, '--r', 'linewidth', lw)
set(gca, 'xticklabel', '');
set(gca, 'yticklabel', '');
set(gca, 'xtick', xt);
set(gca, 'ytick', linspace(-1, 1, 5))
box on;
grid on;
title_str = sprintf('$\\sigma_\\mathbf{d} = %d$', diffusion_list(n));
title(title_str, 'fontsize', fSize_title, 'interpreter', 'latex');
   

if n == 1
    p_anchor = get(gca, 'position');
    p_left = p_anchor(1);
    p_bottom = p_anchor(2);
    p_height = p_anchor(4);
    p_width = p_height * height_fract;    
end


p = p_anchor;
p(1) = p_left + (p_width + dx_subplot) * (n - 1);
p(3) = p_width;
p(4) = p_height;
set(gca, 'position', p);

ylim(1.05 * [-1, 1]);
xlim([1, region_width]);

if n == 1
ylabel('$ \mathcal{R} \{ R\left(k\right) \}$',...
    'interpreter', 'latex', 'fontsize', fSize_labels);

h_leg = legend('$\mathcal{R}\{ R \left(k\right) \}$',...
    '$\left| \mathcal{R}\{ R \left(k\right)\} \right| $');
set(h_leg, 'interpreter', 'latex')
set(h_leg, 'fontsize', fSize_legend);

end

% subplot(2, num_diffusions, P2, 'parent', f);
axes(ax{2 * n});
plot(cc_sum_phase, '-k', 'linewidth', lw);

p = get(gca, 'position');
p(1) = p_left + (p_width + dx_subplot) * (n - 1);
p(2) = p_bottom - p_height - dy_subplot;
p(3) = p_width;
p(4) = p_height;
set(gca, 'position', p);
xlim([1, region_width]);
ylim([-1, 1]);
set(gca, 'yticklabel', '');
set(gca, 'ytick', linspace(-1, 1, 10))
set(gca, 'xtick', xt);
set(gca, 'xticklabel', xtl);
box on;
grid on;
set(gca, 'fontsize', fSize_axes);
xlabel('$k / \pi$', ...
    'interpreter', 'latex', 'fontsize', fSize_labels);

if n == 1
    ylabel('$\mathcal{R} \{ \phi \left(k\right) \} $',...
    'interpreter', 'latex', 'fontsize', fSize_labels);
end


    
end


set(gcf, 'color', 'white');
tightfig;








% surf(cc_mag_norm);

% % This is for making the SCC ensemble
% % versus phase angle plots
% subplot(2, 1, 1);
% % surf(cc_ifft);
% surf(scc_full_sum);
% axis square
% axis vis3d
% axis off
% xlim([1, region_width]);
% ylim([1, region_height]);
% set(gca, 'view', [0, 0]);
% set(gca, 'units', 'normalized')
% p = get(gca, 'position');
% p(2) = 0.78 * p(2);
% p(4) = 1.8 * p(4);
% set(gca, 'position', p);
% % %
% subplot(2, 1, 2);
% imagesc(phase_angle);
% axis image;
% axis off
% p = get(gca,'position');
% p(2) = 1.75 * p(2);
% set(gca, 'position', p);
% % % % 
% % % % plot_path = sprintf('~/Desktop/corr_figs/fig_01.%02d.eps', s_rand);
% % % % print(1, '-depsc', plot_path);
% 







