addpath('~/Desktop/codes/tightfig');
addpath('~/Desktop/codes/export_fig');

load('~/Desktop/apc_test_s_mean_16_s_rand_3_single_region.mat');


rpc_filter = spectralEnergyFilter(region_width, region_height, 3);
cc_rpc_spect = phaseOnlyFilter(cc_full_sum) .* rpc_filter;
g3 = fftshift(abs(ifft2(fftshift(cc_rpc_spect))));

fSize_titles = 32;
fSize_axislabels = 16;
fSize_axes = 16;

xv = linspace(0, 1, region_width);
yv = linspace(0, 1, region_height);

[x, y] = meshgrid(xv, yv);

plot_dx = -0.05;

figure(2); 
subplot(1, 3, 1); 
surf(x, y, g ./ max(g(:)), 'linewidth', 1E-5); 
axis square; 
xlim([0, 1]);
ylim([0, 1]);
zlim([0, 1]);
box on
set(gca, 'view', [0, 0]);
% set(gca, 'fontsize', fSize_axes);
title({'$\textrm{Unfiltered cross correlation}$'}, ...
    'interpreter', 'latex', ...
    'fontsize',fSize_titles);
set(gca, 'zticklabel', {''});
set(gca, 'xticklabel', {''});
set(gca, 'ytick', 0 : 0.2 : 1);

subplot(1, 3, 2); 
surf(x, y, g3 ./ max(g3(:)), 'linewidth', 1E-5); 
axis square; 
xlim([0, 1]);
ylim([0, 1]);
zlim([0, 1]);
box on;
set(gca, 'view', [0, 0]);
p = get(gca, 'position');
p(1) = p(1) + plot_dx;
set(gca, 'position', p);
set(gca, 'zticklabel', {''});
set(gca, 'xticklabel', {''});
title({'$\textrm{Phase correlation with}$', '$\textrm{particle-size filter (RPC)}$'}, ...
    'interpreter', 'latex', ...
    'fontsize',fSize_titles);

subplot(1, 3, 3); 
surf(x, y, g2 ./ max(g2(:)), 'linewidth', 1E-5); 
axis square; 
xlim([0, 1]);
ylim([0, 1]);
zlim([0, 1]);
box on;
set(gca, 'view', [0, 0]);
p = get(gca, 'position');
p(1) = p(1) + 2 *  plot_dx;
set(gca, 'position', p);
set(gca, 'zticklabel', {''});
set(gca, 'xticklabel', {''});
title({'$\textrm{Phase correlation with exact filter}$'}, ...
    'interpreter', 'latex', ...
    'fontsize',fSize_titles);

set(gcf, 'color', 'white');

tightfig;
set(gcf, 'position', [  -1644         369        1519         710]);

export_fig('~/Desktop/correlation_plane_subplots_single_region.png', '-r300');




