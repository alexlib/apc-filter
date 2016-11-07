% addpaths;


% Number of corresponding regions
num_regions_eq = 200;
num_regions_neq = 1;

% Image dimensions
region_height = 128;
region_width  = 128;

% Fraction of diffusion std dev 
% in the Y direction compared to X
sy_rand_fract = 1;

% List of diffusion ratios
R_list = [0, 3, 5];
R_list = 0.6;

% MAD thresholds
mad_x_thresh = 0.1;
mad_y_thresh = 0.1;

fSize_axes = 9;
fSize_title = 12;
fSize_textbox = 16;

% Number of plot columns
n_plot_cols = 3;

units_str = 'normalized';


p_scale = 1;

dx_plot = p_scale * 0.4;
dx_plot_fract = p_scale * 0.75;
dy_plot_fract = p_scale * 0.55;
subplot_width = p_scale * 0.25;
subplot_height = p_scale * 0.25;
p_origin_left = 0.00;
p_origin_bottom = 0.7073;


dy_plot_mag = 0;

dx_plot_2D = 0.06;
dy_plot_2D = 0.1;

line_plot_fract = 0.4;


% Textbox spacing
dy_textbox = 0.19;
dx_textbox = -0.015;

% axis view vector
ax_view = [16.8000; 4.8000];
ax_view = [45, 10];

% Z limits 
zl = 1.1 * [-1, 1];

% Figure position
f_pos  = [ -0.5844    0.1458    0.4742    0.6882];

% Line width for mesh plots
lw_mesh = 1E-5;
lw_line = 3;

sx_mean = 5;
sy_mean = 0;
sz_mean = 0;



% Window size
window_fraction = 0.5 * [1, 1];

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

s_range_list = R_list * d_mean;
diffusion_list = R_list .* d_mean / 4;
num_diffusions = length(diffusion_list);


R = 0;

% Range of displacements
sx_range = R * d_mean;
sy_range = R * d_mean;

% d_mean = 1;
% Particle concentration in particles per pixel
particle_concentration = 2E-2;

% Image noise
noise_mean_fract = 0E-2;
noise_std_fract  = 10E-2;
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

figure(1);
% set(gcf, 'position', f_pos);
delete(findall(gcf,'Tag','text_box'))


for n = 1 : num_diffusions

accept_this = false;



while accept_this == false

% if n == num_diffusions
%     keyboard
% end
    
% Random displacements
s_rand = diffusion_list(n);
sx_rand = s_rand;
sy_rand = s_rand * sy_rand_fract;
% sy_rand = 0;
sz_rand = s_rand;

sx_range = s_range_list(n);
sy_range = s_range_list(n);

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
    
%      % Uncomment this for shearing
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

% Plot indices
P1 = sub2ind([n_plot_cols, num_diffusions], 1, n);
P2 = sub2ind([n_plot_cols, num_diffusions], 2, n);
P3 = sub2ind([n_plot_cols, num_diffusions], 3, n);


p_colors = get_plot_colors(3);

% Fit the Gaussian;
[~, ~, ~, ~, gauss_fit] = fit_gaussian_2D(abs(cc_full_sum) ./ max(abs(cc_full_sum(:))));

set(gcf, 'units', units_str);
set(gcf, 'color', 'white');
set(gcf, 'position', f_pos);

% % % % PLOT 1
% ax1 = subplot(num_diffusions, n_plot_cols, P1);
% if exist('ax1', 'var')
%     cla(ax1);
% end
ax1 = subplot(num_diffusions, n_plot_cols, P1);
set(gca, 'units', units_str);
mesh(real(cc_full_sum) ./ max(abs(real(cc_full_sum(:)))), ...
    'edgecolor', 'black', 'linewidth', lw_mesh);

axis square
xlim([1, region_width]);
ylim([1, region_height]);
zlim(zl);

p = get(gca, 'position');
p(1) = p_origin_left;
p(2) = p_origin_bottom - (n - 1) * subplot_height * dy_plot_fract;
p(3) = subplot_width;
p(4) = subplot_height; 
set(gca, 'position', p); 
axis off
axis vis3d
set(gca, 'view', ax_view);
 
% % % % % PLOT 2
% subplot(num_diffusions, n_plot_cols, P2 );
% set(gca, 'units', units_str);
% mesh(imag(cc_full_sum) ./ max(abs(imag(cc_full_sum(:)))),...
%     'edgecolor', 'black', 'linewidth', lw_mesh);
% axis square
% xlim([1, region_width]);
% ylim([1, region_height]);
% zlim(zl);
% p = get(gca, 'position');
% p(1) = p_origin_left + 1 * subplot_width * dx_plot_fract;
% p(2) = p_origin_bottom - (n - 1) * subplot_height * dy_plot_fract;
% p(3) = subplot_width;
% p(4) = subplot_height;
% set(gca, 'position', p);
% axis off
% axis vis3d
% set(gca, 'view', ax_view);

% % % % PLOT 3
ax2 = subplot(num_diffusions, n_plot_cols, P2);
set(gca, 'units', units_str);
mesh(abs(cc_full_sum) ./ max(abs(cc_full_sum(:))), ...
    'edgecolor', 'black', 'linewidth', lw_mesh);
axis square
xlim([1, region_width]);
ylim([1, region_height]);
zlim(zl);
p = get(gca, 'position');
p(1) = p_origin_left + 1 * subplot_width * dx_plot_fract;
p(2) = p_origin_bottom - (n - 1) * subplot_height * dy_plot_fract;
p(3) = subplot_width;
p(4) = subplot_height;
set(gca, 'position', p); 
axis off
axis vis3d
set(gca, 'view', ax_view);

% Plot 3
ax3 = subplot(num_diffusions, n_plot_cols, P3);
set(gca, 'units', units_str);

x = (1 : region_width) - fourier_zero(region_width);
k = x / region_width;
sx = d_mean / 4;

% This is the FT of the particle shape 
F = sqrt(pi) * sx * sqrt(2) * exp(-2  * sx^2 * pi^2 * k.^2);

% This is the FT of the auto correlation
A = F.* F;

% This is the FT of the displacement PDF
P_x = exp(-2 * pi^2 * k.^2 * sx_rand^2);
P_y = exp(-2 * pi^2 * k.^2 * sy_rand^2);

% This is the total SNR envelope
E_x = A .* P_x;
E_y = A .* P_y;

E_norm_x = E_x ./ max(E_x(:));
E_norm_y = E_y ./ max(E_y(:));

E_slice_yproj = gauss_fit(yc, :);
E_slice_norm_yproj = E_slice_yproj ./ max(E_slice_yproj(:));

E_slice_xproj = gauss_fit(:, xc);
E_slice_norm_xproj = E_slice_xproj ./ max(E_slice_xproj(:));

mad_y = mean(abs(E_norm_y' - E_slice_xproj));
mad_x = mean(abs(E_norm_x - E_slice_yproj));

mesh(gauss_fit ./ max(gauss_fit(:)), ...
    'edgecolor', 'black', 'linewidth', lw_mesh);

hold on;
n_grid_points = 7;
xpv_yproj = (1 : region_width);
ypv_yproj = region_height * ones(1, region_width);
gxv_yproj = linspace(1, region_width, n_grid_points);
gzv_yproj = linspace(0.01, 1, n_grid_points);
[gx_yproj, gz_yproj] = meshgrid(gxv_yproj, gzv_yproj);
gy_yproj = region_height * ones(size(gz_yproj));
mesh(gx_yproj, gy_yproj, gz_yproj, 'edgecolor', 'black', 'linewidth', lw_mesh, 'facealpha', 0);
box on;
plot3(xpv_yproj, ypv_yproj, E_norm_x, '-', 'linewidth', lw_line, 'color', p_colors(2, :));
plot3(xpv_yproj, ypv_yproj - 0.2, E_slice_norm_yproj, '--', 'linewidth', lw_line, 'color', p_colors(3, :));


 
xpv_xproj = (1 : region_height);
ypv_xproj = ones(1, region_height);

gyv_xproj = linspace(1, region_height, n_grid_points);
gzv_xproj = linspace(0, 1, n_grid_points);

[gy_xproj, gz_xproj] = meshgrid(gyv_xproj, gzv_xproj);

gx_xproj = ones(size(gz_xproj));
mesh(gx_xproj, gy_xproj, gz_xproj, 'edgecolor', 'black', 'linewidth', lw_mesh, 'facealpha', 0);
plot3(ypv_xproj, xpv_xproj, E_norm_y, '-', 'linewidth', lw_line, 'color', p_colors(2, :));
plot3(ypv_xproj + 0.2, xpv_xproj, E_slice_norm_xproj, '--', 'linewidth', lw_line, 'color', p_colors(3, :));

axis square
xlim([1, region_width]);
ylim([1, region_height]);
zlim(zl);
p = get(gca, 'position');
p(1) = p_origin_left + 2 * subplot_width * dx_plot_fract;
p(2) = p_origin_bottom - (n - 1) * subplot_height * dy_plot_fract + dy_plot_mag;
p(3) = subplot_width;
p(4) = subplot_height;
set(gca, 'position', p); 
axis off
axis vis3d
set(gca, 'view', ax_view);

hold off;

fprintf(1, 'MAD y: %0.3f\tMAD x: %0.3f\n', mad_y, mad_x);

drawnow;

accept = min([mad_y < mad_y_thresh, mad_x < mad_x_thresh]);

if accept;
    accept_this = true;
else
    cla(ax1);
    cla(ax2);
    cla(ax3);
end

end

end
 
clear label_strings;
label_strings{1} = '$\mathcal{R} \langle \tilde{R} \left( \mathbf{k} \right) \rangle$';
% label_strings{2} = '$\mathcal{I} \langle \tilde{R} \left( \mathbf{k} \right) \rangle$';
label_strings{2} = '$\left| \langle \tilde{R} \left( \mathbf{k} \right) \rangle \right|$';
label_strings{3} = '$\textrm{Gaussian fit}$';

for p = 1 : 3
    
% Textbox positions
textbox_x  = p_origin_left + (p-1) * subplot_width * 1 * dx_plot_fract - (p == 3) * 0.01 + 0.1;
textbox_y = p_origin_bottom + subplot_height * dy_plot_fract + 0.10;
textbox_width = 0.5 * subplot_width;
textbox_height = 0.25 * subplot_height;

% Textbox position
textbox_pos = [textbox_x, textbox_y, 0, 0];

textbox_str = label_strings{p};

% Make the annotation
annotation('textbox', 'position', textbox_pos, ...
    'string', textbox_str, ...
    'fontsize', fSize_textbox, ...
    'interpreter', 'latex', ...
    'Tag', 'text_box');


end

textbox_sep_y = -0.02;

% Diffusion labels
for p = 1 : num_diffusions
    % Textbox positions
textbox_x = p_origin_left - dx_textbox;
textbox_y_xdir = p_origin_bottom - (p - 1) * subplot_height * dy_plot_fract + dy_textbox;
textbox_y_ydir = p_origin_bottom - (p - 1) * subplot_height * dy_plot_fract + dy_textbox + textbox_sep_y;

textbox_width = 0.5 * subplot_width;
textbox_height = 0.25 * subplot_height;

% Textbox position
textbox_pos_xdir = [textbox_x, textbox_y_xdir, 0, 0];
textbox_pos_ydir = [textbox_x, textbox_y_ydir, 0, 0];

textbox_str_xdir = sprintf('$\\sigma_{d_x} / \\sigma_\\tau = %0.2f$', diffusion_list(p) / (d_mean / 4));
textbox_str_ydir = sprintf('$\\sigma_{d_y} / \\sigma_\\tau = %0.2f$', diffusion_list(p) / (d_mean / 4) * sy_rand_fract);

% Make the annotation
annotation('textbox', 'position', textbox_pos_xdir, ...
    'string', textbox_str_xdir, ...
    'fontsize', fSize_textbox, ...
    'interpreter', 'latex', ...
    'Tag', 'text_box');

% Make the annotation
annotation('textbox', 'position', textbox_pos_ydir, ...
    'string', textbox_str_ydir, ...
    'fontsize', fSize_textbox, ...
    'interpreter', 'latex', ...
    'Tag', 'text_box');
    
end

% Make the legend by hand. 
% Like our forefathers used to.

dx_legend = 0.05;
dy_legend = 0.09;
leg_line_length = subplot_width / 7;
leg_line_spacing_x = 0.08;

text_spacing_x = 0;
text_spacing_y = 0.0125;
c_blue = p_colors(2, :);
c_black = p_colors(3, :);

delete(findall(gcf,'Tag','legend_line_01'))
delete(findall(gcf,'Tag','legend_line_02'))
delete(findall(gcf,'Tag','legend_text_01'))
delete(findall(gcf,'Tag','legend_text_02'))

h1 = annotation('line', 'Tag', 'legend_line_01');
h1.Color = c_black;
h1.LineWidth = lw_line;
h1.LineStyle = '--';

h1_left = p_origin_left + 2 * subplot_width * dx_plot_fract + dx_legend;
h1_bottom = p_origin_bottom - (num_diffusions - 1) * subplot_height * dy_plot_fract + dy_plot_mag + dy_legend;
h1_width = leg_line_length;
h1_height = 0;
h1.Position= [h1_left, h1_bottom, h1_width, h1_height] ;


% Legend line 2
h2 = annotation('line', 'Tag', 'legend_line_01');
h2.Color = c_blue;
h2.LineWidth = lw_line;
h2.LineStyle = '-';

h2_left = h1_left + leg_line_spacing_x;
h2_bottom = h1_bottom;
h2_width = h1_width;
h2_height = h1_height;
h2.Position= [h2_left, h2_bottom, h2_width, h2_height] ;


t1 = annotation('textbox', 'position', textbox_pos, ...
    'string', 'Fit', ...
    'fontsize', fSize_textbox, ...
    'interpreter', 'latex', ...
    'Tag', 'legend_text_01', 'linestyle', 'none');

t1_pos = h1.Position;
t1_pos(1) = t1_pos(1) + leg_line_length + text_spacing_x;
t1_pos(2) = t1_pos(2) + text_spacing_y;
t1.Position = t1_pos;

t2 = annotation('textbox', 'position', textbox_pos, ...
    'string', 'Theory', ...
    'fontsize', fSize_textbox, ...
    'interpreter', 'latex', ...
    'Tag', 'legend_text_02', 'linestyle', 'none');

t2_pos = h2.Position;
t2_pos(1) = h2_left + leg_line_length + text_spacing_x;
t2_pos(2) = t1_pos(2);
t2.Position = t2_pos;


figure(2);
surf(real(cc_full_sum) ./ max(abs(real(cc_full_sum(:)))), ...
   'linewidth', lw_mesh);
axis off
axis vis3d
set(gca , 'view', [0, 0]);

figure(3);
plot(imag(cc_full_sum(yc, :)), '-k', 'linewidth', lw_line);
xlim([1, region_width]);
axis square
axis off

figure(4);
plot(imag(phaseOnlyFilter(cc_full_sum(yc, :))), '-k', 'linewidth', lw_line);
xlim([1, region_width]);
axis square;
axis off

figure(5);
plot(abs(cc_full_sum(yc+1, :)), '-k', 'linewidth', lw_line);
xlim([1, region_width]);
axis square;
axis off

figure(6);
plot(gauss_fit(yc, :), '-k', 'linewidth', lw_line);
xlim([1, region_width]);
axis square;
axis off;


g_spatial = fftshift(abs(ifft2(fftshift(phaseOnlyFilter(cc_full_sum) .* gauss_fit))));

figure(7);
plot(g_spatial(yc, :), '-k', 'linewidth', lw_line);
xlim([1, region_width]);
axis square;
axis off;

filtered_phase = imag(phaseOnlyFilter(cc_full_sum)) .* gauss_fit;


figure(8);
plot(filtered_phase(yc, :), '-k', 'linewidth', lw_line);
xlim([1, region_width]);
axis square;
axis off;









