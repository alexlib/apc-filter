% addpaths;


fSize_axes = 9;
fSize_title = 12;

sx_mean = 5;
sy_mean = 0;
sz_mean = 0;

% Window size
window_fraction = 1 * [1, 1];


num_trials = 1;

% Number of corresponding regions
num_regions_eq = 100;
num_regions_neq = 1;
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
xv = (1 : region_width) - xc;
yv = (1 : region_height) - yc;

% Grid
[x, y] = meshgrid(xv, yv);





% % Particle stuff

% Standard deviation of particle image diameters
d_std = 0;

% Mean particle diameter
d_mean = 3;

% Diffusion ratio
R = 2;

% Diffusion list
diffusion_list = R .* (d_mean / 4);

% Number of diffusions
num_diffusions = length(diffusion_list);


% Range of displacements
sx_range = R * d_mean;
sy_range = R * d_mean;

% d_mean = 1;
% Particle concentration in particles per pixel
particle_concentration = 2E-2;

% Image noise
noise_mean_fract = 10E-2;
noise_std_fract  = 5.0E-2;
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
lw = 1E-5;

for n = 1 : num_diffusions
    
% Random displacements
s_rand = diffusion_list(n);
sx_rand = s_rand;
sy_rand = s_rand;
sz_rand = s_rand;


% Allocate running CC sum for the CCC
cc_sum = zeros(region_height, region_width) + ...
    1i * zeros(region_height, region_width);

% Allocate running CC sum for the CC minus the NCC
cc_full_sum = zeros(region_height, region_width) + ...
    1i * zeros(region_height, region_width);

scc_full_sum = zeros(region_height, region_width);

rpc_full_sum = zeros(region_height, region_width);


cc_abs_cur_new = zeros(region_height, region_width);

cc_test = zeros(region_height, region_width);

 scc_max = zeros(num_regions_eq, 1);
  
% Do the corresponding correlation
for k = 1 : num_regions_eq
    
    
    % Uncomment this for random displacements
    dx = sx_mean + sx_rand * randn(num_particles, 1);
    dy = sy_mean + sy_rand * randn(num_particles, 1);
    dz = sz_mean + sz_rand * randn(num_particles, 1);

    % Particle positions (image 1)
    x_01 = (x_max - x_min) * rand(num_particles, 1) + x_min;
    y_01 = (y_max - y_min) * rand(num_particles, 1) + y_min;
    z_01 = (z_max - z_min) * rand(num_particles, 1) + z_min; 
    
    % Uncomment this for shearing
%     dx = sx_range * (y_01 - yc) /region_height + sx_mean;
%     dy = sy_range * (x_01 - xc) /region_width + sy_mean;
%     dz = zeros(num_particles, 1);

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

    scc_cur = fftshift(abs(ifft2(fftshift(cc_cur))));
    
    scc_max(k) = max(scc_cur(:));

    rpc_cur = fftshift(abs(ifft2(fftshift(phaseOnlyFilter(cc_cur) .* rpc_filter))));

    scc_full_sum = scc_full_sum + scc_cur;

    rpc_full_sum = rpc_full_sum + rpc_cur;

    % Full CC sum
    cc_full_sum = cc_full_sum + cc_cur;

end


cc_mag = abs(cc_full_sum);
cc_mag_norm = cc_mag ./ max(cc_mag(:));


% Fit the APC filter
[A, sy, sx, B, ARRAY] = fit_gaussian_2D(cc_mag_norm);

% Create the filter

% Filter the APC
apc_filter = exp(-x.^2 / (2 * sx^2)) .* exp(-y.^2 / (2 * sy^2));

apc_plane_spect = apc_filter .* phaseOnlyFilter(cc_full_sum);

rpc_plane_spectral_ensemble = rpc_filter .* phaseOnlyFilter(cc_full_sum);

% Spatial APC plane
apc_plane_spatial = fftshift(abs(ifft2(fftshift(apc_plane_spect))));

rpc_plane_spectral_ensemble_spatial = fftshift(abs(ifft2(fftshift(rpc_plane_spectral_ensemble))));

scc_plane_spectral_ensemble_spatial = fftshift(abs(ifft2(fftshift(cc_full_sum))));


P1 = sub2ind([4, num_diffusions], 1, n);
P2 = sub2ind([4, num_diffusions], 2, n);
P3 = sub2ind([4, num_diffusions], 3, n);
P4 = sub2ind([4, num_diffusions], 4, n);

subplot(num_diffusions, 4, P1);
surf(scc_plane_spectral_ensemble_spatial ./ max(scc_plane_spectral_ensemble_spatial(:)), 'linewidth', lw);
axis square
xlim([1, region_width]);
ylim([1, region_height]);

if n == 1
    p_origin = get(gca, 'position');
    p_origin_left  = p_origin(1);
    p_origin_bottom = p_origin(2);
    p_origin_width = p_origin(3);
    p_origin_height = p_origin(4);
else
  p = get(gca, 'position');
  p(2) =  p_origin_bottom - (n - 1) * p_origin_height;
  set(gca, 'position', p);  
end
axis off
axis vis3d

subplot(num_diffusions, 4, P2 );
surf(rpc_plane_spectral_ensemble_spatial ./ max(rpc_plane_spectral_ensemble_spatial(:)), 'linewidth', lw);
axis square
xlim([1, region_width]);
ylim([1, region_height]);

if n > 1
    p = get(gca, 'position');
    p(2) = p_origin_bottom - (n - 1) * p_origin_height;
    set(gca, 'position', p); 
end

p = get(gca, 'position');
p(1) = p_origin_left + p_origin_width;
set(gca, 'position', p);

axis off
axis vis3d

subplot(num_diffusions, 4, P3);
surf(apc_plane_spatial ./ max(apc_plane_spatial(:)), 'linewidth', lw);
axis square
xlim([1, region_width]);
ylim([1, region_height]);
if n > 1
    p = get(gca, 'position');
    p(2) = p_origin_bottom - (n - 1) * p_origin_height;
    set(gca, 'position', p); 
end

p = get(gca, 'position');
p(1) = p_origin_left + 2 * p_origin_width;
set(gca, 'position', p);

axis off
axis vis3d

subplot(num_diffusions, 4, P4);
surf(real(cc_full_sum) ./ max(real(cc_full_sum(:))), 'linewidth', lw);
axis square
xlim([1, region_width]);
ylim([1, region_height]);
if n > 1
    p = get(gca, 'position');
    p(2) = p_origin_bottom - (n - 1) * p_origin_height;
    set(gca, 'position', p); 
end

p = get(gca, 'position');
p(1) = p_origin_left + 3 * p_origin_width;
set(gca, 'position', p);

axis off
axis vis3d
set(gca, 'view', [0, 0]);



end

s_max = 62;

scc_ens_norm = scc_full_sum ./ num_regions_eq;

cur_max = max(scc_ens_norm(:) / s_max);

fprintf(1, 'Peak height: %0.2f\n', cur_max);

scc_plane = scc_cur;

% subplot(2, 1, 1);
figure(1);
subplot(2, 1, 1);
mesh(scc_full_sum ./ num_regions_eq, 'edgecolor', 'black', 'linewidth', 1E-5);
axis square
xlim([1, region_width]);
ylim([1, region_height]);
zlim([0, s_max]);
axis off;
set(gca, 'view', [-35, 15]);

scc_angle = angle(cc_full_sum);
scc_real = real(cc_full_sum);
scc_abs = abs(cc_full_sum);
scc_phase = real(phaseOnlyFilter(cc_full_sum));



plot_row = 4;

% figure(2);
% imagesc(scc_real);
% axis image;


subplot(2, 1, 2);
imagesc(scc_phase);
axis image;
axis off
colormap winter
caxis(0.6 * [-1, 1]);

figure(2);
plot(scc_phase(yc + plot_row, :), '-k', 'linewidth', 2);
axis square;
axis off
colormap winter
caxis(0.6 * [-1, 1]);

% hold on;
% plot(scc_real(yc + plot_row, :), '--k', 'linewidth', 2);
% hold off

% subplot(2, 1, 2);
% plot((abs((cc_cur(yc + 2, :)))), '-k', 'linewidth', 2);
% axis square;
% xlim([1, region_width]);





% p = phaseOnlyFilter(cc_full_sum(yc, :));
% plot(real(p), '-k', 'linewidth', 2);
% pbaspect([3, 1, 1]);


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







