function [TRANSLATION_Y, TRANSLATION_X, SPATIAL_PLANE, ...
    CORR_HEIGHT, CORR_DIAMETER] = ...
    scc_ensemble_spatial(REGION_MATRIX_01, REGION_MATRIX_02, ...
    PEAK_FIT_METHOD);

if nargin < 3
    PEAK_FIT_METHOD = 1;
end

% Calculate size of interrogation regions (homogeneous) (pixels)
[region_height, region_width, num_regions] = size(REGION_MATRIX_01);

% Initialize the cross correlation
SPATIAL_PLANE = zeros(region_height, region_width);

% Loop over all the regions
for k = 1 : num_regions
    
    % Extract the regions
   region_01 = REGION_MATRIX_01(:, :, k);
   region_02 = REGION_MATRIX_02(:, :, k);
   
   % Fourier transforms
   f1 = fft2(region_01);
   f2 = fft2(region_02);
   
   % Add the complex correlation to the ensemble
   cc_spectral = f1 .* conj(f2);
   
   SPATIAL_PLANE = SPATIAL_PLANE + fftshift(abs(ifft(cc_spectral)));
    
end

% Subpixel fit
% Prana subpixel implmentation of the sub-pixel fit
[TRANSLATION_X, TRANSLATION_Y, CORR_HEIGHT, CORR_DIAMETER] = ...
    subpixel(SPATIAL_PLANE, ...
    region_width, region_height, ones(size(SPATIAL_PLANE)), PEAK_FIT_METHOD, ...
    0, sqrt(8) * [1, 1]);


end