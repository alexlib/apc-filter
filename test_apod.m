

lw = 3;



region_width = 128;
xc = fourier_zero(region_width);

win_fract = 0.5;

sx = findGaussianWidth(region_width, 1, win_fract * region_width, 1);

dx_range = 8;

sx_p = region_width / ( sx * dx_range); 

x = (1 : region_width) - xc;

xp = x * dx_range / region_width;

% Wave number
k = x / region_width;

g = exp(-x.^2 / (2 * sx^2));

gp = exp(-x.^2 / (2 * sx_p^2));

g2 = g.*g;

g_dist = exp(-(x * region_width / dx_range ).^2 / (2 * sx^2));

F_dist = sqrt(pi) * dx_range * sx / region_width * exp(-k.^2 * pi^2 * ...
    dx_range^2 * sx^2 / region_width^2);


plot(x, g, '-k', 'linewidth', lw);
hold on;
plot(x, gp, '-r', 'linewidth', lw);
hold off