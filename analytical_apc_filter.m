% 
addpaths('~/Desktop');

% Size of the domain
region_width = 128;

% This is the location of
% the centroid pixel.
xc = fourier_zero(region_width);

% This is the domain
x = (1 : region_width) - xc;

k = x / region_width;

% Particle diameter
dp = sqrt(8);

% Particle diameter is defined as four times the
% standard deviation of the Gaussian particle shape
sx = dp / 4;

% This is the Gaussian that
% defines the particle shape.
f = exp(-x.^2 / (2 * sx^2));

% This is the FT of the particle shape 
F = sqrt(2 * pi) * sx * exp(-2  * sx^2 * pi^2 * k.^2);

% This is the FT of the auto correlation
A = F.* F;

% This is the Prana RPC filter
A_rpc = spectralEnergyFilter(1, 128, dp);

% This is the displacent PDF
dx_range = 8;

pdf_amp = 1 / dx_range;

N = 10;

P = N * pdf_amp * abs(sinc(dx_range  *  k));

W = P .* A;

% This is the FT of the displacement PDF
figure(2);
% plot(k, A_rpc, '-k', 'linewidth', 5);

plot(k, A, '-k', 'linewidth', 3);
hold on
plot(k, P, '-', 'linewidth', 3, 'color', 0.5 * [1, 1, 1]);
plot(k, W, ':', 'linewidth', 3, 'color', 0.5 * [1, 1, 1]);
hold off
xlim(0.5* [-1, 1])
axis square;

h = legend('$A_c \left(k\right)$', '$P_{\mathbf{\overline{d}}}\left(k\right)$', '$\left| E \left(k\right) \right|$');
set(h, 'interpreter', 'latex');
set(h, 'fontsize', 16);









