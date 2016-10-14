% 

addpath ~/Documents/School/VT/Research/Aether/Codes/apc-filter/
addpaths('..');

fSize_axes = 24;
fSize_labels = 30;
fsize_title = 18;

lw = 3;

% Size of the domain
region_width = 512;

% Window fraction
win_fract = 0.5;

% This is the location of
% the centroid pixel.
xc = fourier_zero(region_width);

% This is the domain
x = (1 : region_width) - xc;

k = x / region_width;

% Particle diameter
% dp = 1 * sqrt(8);
dp = 3;
% dp = 5;

% dp = 1E-3;

% This is the displacent PDF
% dx_range = 10;
dx_range = region_width / 80;

% Mean displacement
dx_mean = dx_range / 2;
dx_mean = 8;

% Ratio between spread and diameter
R = dx_range / dp;

% Particle diameter is defined as four times the
% standard deviation of the Gaussian particle shape
sx = dp / 4;

% This is the Gaussian that
% defines the particle shape.
f = exp(-x.^2 / (2 * sx^2));

% This is the Fourier transform
% of the particle shape, 
% evaluated numerically
F_c = real(fftshift(fft(fftshift(f))));

% This is the FT of the particle shape 
F = sqrt(pi) * sx * sqrt(2) * exp(-2  * sx^2 * pi^2 * k.^2);

% This is the FT of the auto correlation
A = F.* F;

% This is the auto correlation of the particle,
% computed numerically
A_c = abs(F_c .* conj(F_c));

% This is the Prana RPC filter
A_rpc = spectralEnergyFilter(1, 128, dp);

% This is the FT of the displacement distribution
% (Normal distribution)
P = (sinc(R * dp * k));

% Correlation SNR evaluated analytically
W = P .* A;

% Correlation SNR evaluated numerically
W_c = (P .* A_c);

W_c_abs = abs(W_c);
W_abs = abs(W);

% Sinusoidal displacement correlation
pc = exp(-1i * 2 * pi * dx_mean * k);

% Displacement correlation windowed by
% the magnitude.
pc_full = real(W_c .* pc);

% This is the FT of the displacement PDF

% Define the 'gray' color
c_gray = 0.7 * [1, 1, 1];

figure(1);
% plot(k, A_rpc, '-k', 'linewidth', 5);
% subplot(1, 2, 1)
plot(k, pc_full, '-k', 'linewidth', lw, 'color', c_gray);
hold on
plot(k, A , '--k', 'linewidth', lw);
plot(k, (P) , '--', 'linewidth', lw, 'color', c_gray);
plot(k, W_abs , '-', 'linewidth', lw, 'color', 'black');
hold off

% Options
axis square;
grid on
box on
set(gca, 'fontsize', fSize_axes);
% Label the axes
xlabel('$k / \pi$', 'interpreter', 'latex', 'fontsize', fSize_labels)
ylabel('$E / E_o$',...
    'interpreter', 'latex', 'fontsize', fSize_labels);

% Axes limits
xlim(0.5* [-1, 1])
ylim([-2.6, 3.6]);

% Format the axes
set(gca, 'xtick', 0.5 * [-1, 0, 1]);
set(gca, 'xticklabel', {'-1', '0', '1'});


% Legend entries
legend_entries{1} = '$\mathcal{R} \{ \tilde{R}_D \left(k\right) \}$';
legend_entries{2} = '$A \left(k\right)$';
legend_entries{3} = '$P_{\mathbf{\overline{d}}}\left(k\right)$';
legend_entries{4} = '$E \left(k\right)$';

% Make the legend
h = legend(legend_entries);

% Format the legend
set(h, 'interpreter', 'latex',...
    'location', 'southeast', ...
    'fontsize', 16);

% Title
title(sprintf(...
    'Spectral correlation, shear gradient 0.1 pix/pix, $d_p = %0.1f$ pix', dp),...
    'interpreter', 'latex', ...
    'fontsize', fsize_title);

% Figure color (need to set
% this to use export_fig
set(gcf, 'color', 'white');

tightfig

% % This is the command used to print to PNG
export_fig ~/Desktop/apc_1D_example_fig.png -r300 -opengl

% This command prints to eps





