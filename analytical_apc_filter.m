% 
addpaths('~/Desktop');

fSize_axes = 16;
fSize_labels = 30;

lw = 3;

% Size of the domain
region_width = 128;

% This is the location of
% the centroid pixel.
xc = fourier_zero(region_width);

% This is the domain
x = (1 : region_width) - xc;

k = x / region_width;

% Particle diameter
dp = 1 * sqrt(8);
% dp = 5;

% dp = 1E-3;

% This is the displacent PDF
% dx_range = 10;
dx_range = 8;

% Particle diameter is defined as four times the
% standard deviation of the Gaussian particle shape
sx = dp / 4;

% This is the Gaussian that
% defines the particle shape.
f = exp(-x.^2 / (2 * sx^2));
% f = exp(-x.^2 / (2 * sx^2)) ./ sx;

F_c = real(fftshift(fft(fftshift(f))));

% This is the FT of the particle shape 
F = sqrt(pi) * sx * sqrt(2) * exp(-2  * sx^2 * pi^2 * k.^2);

% figure(3); 
% plot(F, '-k', 'linewidth', 5);
% hold on;
% plot(F_c, '-r', 'linewidth', 2);
% plot(f, '-', 'color', 0.5 * [1, 1, 1]);
% hold off;
% axis square
% h = legend('Exact', 'Comp', 'p(x)');
% set(h, 'fontsize', 16);

% F = exp(-2  * sx^2 * pi^2 * k.^2);

% This is the FT of the auto correlation
A = F.* F;

% This is the auto correlation of the particle,
% computed numerically
A_c = abs(F_c .* conj(F_c));

% This is the Prana RPC filter
A_rpc = spectralEnergyFilter(1, 128, dp);


% P = N * pdf_amp * abs(sinc(dx_range  *  k));
P = abs(sinc(dx_range * k));

W = P .* A;

W_c = P .* A_c;

% This is the FT of the displacement PDF
figure(2);
% plot(k, A_rpc, '-k', 'linewidth', 5);
subplot(1, 2, 1)
plot(k, A, '--k', 'linewidth', lw);
hold on
plot(k, P, '--', 'linewidth', lw, 'color', 0.5 * [1, 1, 1]);
plot(k, W, '-', 'linewidth', 0.6 * lw, 'color', 'black');
hold off
xlim(0.5* [-1, 1])
axis square;
grid on
set(gca, 'fontsize', fSize_axes);
xlabel('$\omega / \pi$', 'interpreter', 'latex', 'fontsize', fSize_labels)
ylabel('$\frac{E }{NE_o}$',...
    'interpreter', 'latex', 'fontsize', fSize_labels);
h = legend('$A_c \left(k\right)$', '$P_{\mathbf{\overline{d}}}\left(k\right)$', '$\left| E \left(k\right) \right|$');
set(h, 'interpreter', 'latex');
set(h, 'fontsize', 16);

title_str = sprintf('$d_p = %0.1f, \\sigma_d = %0.1f$', dp, dx_range);

title('Exact', 'FontSize', fSize_labels);

subplot(1, 2, 2)
plot(k, A_c, '--k', 'linewidth', lw);
hold on
plot(k, P, '--', 'linewidth', lw, 'color', 0.5 * [1, 1, 1]);
plot(k, W_c, '-', 'linewidth', 0.6 * lw, 'color', 'black');
hold off
xlim(0.5* [-1, 1])
axis square;
grid on
set(gca, 'fontsize', fSize_axes);
xlabel('$\omega / \pi$', 'interpreter', 'latex', 'fontsize', fSize_labels)
ylabel('$\frac{E }{NE_o}$',...
    'interpreter', 'latex', 'fontsize', fSize_labels);
h = legend('$A_c \left(k\right)$', '$P_{\mathbf{\overline{d}}}\left(k\right)$', '$\left| E \left(k\right) \right|$');
set(h, 'interpreter', 'latex');
set(h, 'fontsize', 16);
title('Numerical', 'FontSize', fSize_labels);











