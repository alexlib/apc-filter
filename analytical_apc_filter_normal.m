% 
clear;

addpaths('~/Desktop');



fSize_axes = 16;
fSize_labels = 30;

lw = 3;

% Size of the domain
region_width = 64;

% This is the location of
% the centroid pixel.
xc = fourier_zero(region_width);

% This is the domain
x = (1 : region_width) - xc;

% Wave number
k = x / region_width;

% Particle diameter
dp = 1 * sqrt(8);

% Particle diameter is defined as four times the
% standard deviation of the Gaussian particle shape
sx = dp / 4;

% This is the std dev of the
% displacent PDF
dx_std = 1;

% Ratio of the dx std dev to the particle diameter
R = dx_std / dp;

% This is the Gaussian that
% defines the particle shape.
f = exp(-x.^2 / (2 * sx^2));

% This is the FT of the particle shape,
% evaluated analytically
F = sqrt(2 * pi) * dp/4 * exp(-2 * pi^2 * (dp/4)^2 * k.^2);

% This is the Fourier transform
% of the particle shape, 
% evaluated numerically
F_c = real(fftshift(fft(fftshift(f))));

% This is the FT of the auto correlation
% A = F.* F;

A = pi / 8 * dp^2 * exp(- 1/4 * pi^2 * dp^2 * k.^2);

% This is the auto correlation of the particle,
% computed numerically
A_c = abs(F_c .* conj(F_c));

% This is the Prana RPC filter
A_rpc = spectralEnergyFilter(1, 128, dp);

% P = N * pdf_amp * abs(sinc(dx_range  *  k));
% This is the FT of the displacement distribution
% (Normal distribution)
P = exp(-2 * pi^2 * dx_std^2 * k.^2);

% Correlation SNR evaluated analytically
W = P .* A;

% Correlation SNR evaluated numerically
W_c = P .* A_c;

% Exact solution for the correlation SNR
% using the ratio b/t the dx std and the particle diameter
% E = pi/8 * dp^2 * exp(-pi^2 * k.^2 * (1/4 * dp^2 + 2 * dx_std^2));
E = pi/8 * dp^2 * exp(-1/4 * pi^2 * dp^2 * k.^2 * (1 + 8 * R^2));

% This is the FT of the displacement PDF
figure(1);
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

title_str = sprintf('$d_p = %0.1f, \\sigma_d = %0.1f$', dp, dx_std);

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


figure(2)
r_range = 0 : 0.5 : 1.5;

c = linspace(0, 0.75, length(r_range));

for n = 1 : length(r_range)
    
    R = r_range(n);
    
    E = pi/8 * dp^2 * exp(-1/4 * pi^2 * dp^2 * k.^2 * (1 + 8 * R^2));
    
    legend_entries{n} = sprintf('$R = %0.1f$', R);
    
    plot(k, E./max(E(:)), '-', 'color', c(n) * [1, 1, 1], 'linewidth', lw);

    hold on;
   
    
end


hold off;
axis square;
grid on
set(gca, 'fontsize', fSize_axes);
xlabel('$\omega / \pi$', 'interpreter', 'latex', 'fontsize', fSize_labels)
ylabel('$\mathrm{SNR}\left(k\right)$',...
    'interpreter', 'latex', 'fontsize', fSize_labels);
h = legend(legend_entries);
set(h, 'interpreter', 'latex', 'fontsize', 16);




