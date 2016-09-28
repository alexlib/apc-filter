function xc = fourier_zero(region_width)

% This is the location of the DC component
% of the Fourier spectrium along
% a dimension containing a discrete
% number of points equal to region_width
xc = (region_width  + 1) / 2 + ...
    0.5 * (1 - mod(region_width,  2));

end