% 
addpaths('~/Desktop');

fSize_axes = 16;
fSize_labels = 30;

lw = 2;

% Size of the domain
region_width = 128;

% Window fraction
win_fract = 0.5;

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
dx_range = region_width / 20;

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
P = abs(sinc(R * dp * k));

% Correlation SNR evaluated analytically
W = P .* A;

% Correlation SNR evaluated numerically
W_c = P .* A_c;

% This is the FT of the displacement PDF
figure(figure('rend','painters','pos',[10 10 900 600]));
% plot(k, A_rpc, '-k', 'linewidth', 5);
% subplot(1, 2, 1)
plot(k, A , '--k', 'linewidth', lw);
hold on
plot(k, P , '--', 'linewidth', lw, 'color', 0.5 * [1, 1, 1]);
plot(k, W , '-', 'linewidth', lw, 'color', 'black');
hold off
xlim(0.5* [-1, 1])
axis square;
grid on
set(gca, 'fontsize', fSize_axes);
xlabel('$k / \pi$', 'interpreter', 'latex', 'fontsize', fSize_labels)
ylabel('$\frac{E }{NE_o}$',...
    'interpreter', 'latex', 'fontsize', fSize_labels);
h = legend('$A_c \left(k\right)$', '$\left| P_{\mathbf{\overline{d}}}\left(k\right) \right|$', '$E \left(k\right)$');
set(h, 'interpreter', 'latex');
set(h, 'fontsize', 16);
set(gca, 'xtick', 0.5 * [-1, 0, 1]);
set(gca, 'xticklabel', {'-1', '0', '1'});
box on


% title_str = sprintf('$d_p = %0.1f, \\sigma_d = %0.1f$', dp, dx_range);
% 
% title('Exact', 'FontSize', fSize_labels);
% 
% subplot(1, 2, 2)
% 
% % r_range = linspace(0, 5, 5);
% ppp = [0, 0.1, 0.5];
% dx_range = region_width * ppp;
% R = dx_range ./ dp;
% c = linspace(0, 0.75, length(ppp));
% clear 'legend_entries';
% for n = 1 : length(ppp)
%     
%     E = abs(sinc(R(n) * dp * k)) * pi / 8 * dp^2 .* ...
%         exp(-1/4 * pi^2 * dp^2 * k.^2);
%     
%     legend_entries{n} = sprintf('$R = %0.1f$', ppp(n));
%     
%     plot(k, E./max(E(:)), '-', 'color', c(n) * [1, 1, 1], 'linewidth', lw);
% 
%     hold on;
%    
%     
% end
% 
% 
% hold off;
% axis square;
% grid on
% set(gca, 'fontsize', fSize_axes);
% xlabel('$\omega / \pi$', 'interpreter', 'latex', 'fontsize', fSize_labels)
% ylabel('$\mathrm{SNR}\left(k\right)$',...
%     'interpreter', 'latex', 'fontsize', fSize_labels);
% h = legend(legend_entries);
% set(h, 'interpreter', 'latex', 'fontsize', 16);
% set(gca, 'xtick', 0.5 * [-1, 0, 1]);
% set(gca, 'xticklabel', {'-1', '0', '1'});
% 


% figure(2); 
% plot(k, P, '--', 'linewidth', lw, 'color', 0.5 * [1, 1, 1]);
% axis square;
% 







