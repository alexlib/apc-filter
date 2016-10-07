function [TY, TX] = ...
    rpc_ensemble_spatial(image_list_01, image_list_02, ...
    grid_y, grid_x, region_size, window_fraction, apc_std_y, apc_std_x)
% APC_STD_Y, APC_STD_X, APC_FILTER] = ...
%     calculate_apc_filter_from_image_pair(image_01, image_02, ...
%     grid_y, grid_x, region_size, window_fraction, shuffle_range, shuffle_step)
% 
% CALCULATE_APC_FILTER_FROM_IMAGE_PAIR calculates the standard deviations
% of Gaussian-shaped spectral filters from a pair of images
%
% INPUTS
%     image_01: [M x N] array containing the image data
%         the first image in the pair.
% 
%     image_02: [M x N] array containing the image data
%         the second image in the pair.
% 
%     grid_y: Array or vector containing the pixel
%         row coordinates of the interrogation regions
%         at which the filters will be calculated.
% 
%     grid_x: Array or vector containing the pixel
%         column coordinates of the interrogation regions
%         at which the filters will be calculated.
% 
%     region_size: [2 x 1] vector specifying the
%         height and width of the interrogation regions
%         in pixels.
%         region_size = [region_size_rows, region_size_columns];
% 
%     window_fraction: [2 x 1] vector specifying the 
%         fractional height and width of the interrogation
%         regions' apodization window. Each element is between 0 and 1.
%         0 specifies the image is completely apodized in that direction
%         and 1 specifies the image is un-apodized. 
%         window_fraction = [window_fraction_rows, window_fraction_cols];
% 
%     shuffle_range: [2 x 1] vector specifying the row and column 
%         range across which each interrogation is shifted for calculating
%         the filter. This is the range of the shake-the-box analog for this
%         method. 
%         shuffle_range = [shuffle_range_rows, shuffle_rows_cols];
% 
%     shuffle_step: [2 x 1] vector specifying the row and column 
%         spacing between shifted locations of regions for calculating
%         the filter. This is the spacing of the of the shake-the-box
%         analog for this method. 
%         shuffle_step = [shuffle_step_rows, shuffle_step_cols];
% 
% OUTPUTS
%     APC_STD_Y = Vector containing the standard deviations
%         of the APC filters in the row direction,
%         located at each grid point.
%         The order of this list of values corresponds 
%         to the order of the grid points located
%         at [grid_y(:), grid_x(:)].
%
%     APC_STD_X = Vector containing the standard deviations
%         of the APC filters in the column direction,
%         located at each grid point.
%         The order of this list of values corresponds 
%         to the order of the grid points located
%         at [grid_y(:), grid_x(:)].
%
%

% This is the number of images that will be correlated.
num_images = length(image_list_01);



% % Extract grid coordinates
gy = grid_y(:);
gx = grid_x(:);

% Number of regions
num_regions = length(gy);

if length(apc_std_y) == 1
    apc_std_y = apc_std_y .* ones(num_regions, 1);
end

if length(apc_std_x) == 1
    apc_std_x = apc_std_x .* ones(num_regions,1);
end

% Size of the region
region_height = region_size(1);
region_width = region_size(2);

% Make the Gaussian window
% and save it for later.
g_win = gaussianWindowFilter(...
    [region_height, region_width], ...
    window_fraction, 'fraction');

% Allocate array to hold all of the ensemble
% subregions. 
% This is a complex array. 
% This could be large in memory, be careful.
spatial_correlation_array = ...
    zeros(region_height, region_width, num_regions);
    
% Center of fourier domain
xc = fourier_zero(region_width);
yc = fourier_zero(region_height);
    
% Grid of region coordinates
X = (1 : region_width)  - xc;
Y = (1 : region_height) - yc;

% Pixel coordinates of each region
[x, y] = meshgrid(X, Y);

% Allocate the translation vectors
TX = zeros(num_regions, 1);
TY = zeros(num_regions, 1);

% Loop over all the images.
for p = 1 : num_images
    
     % Inform the user.
    fprintf(1, 'On image %d of %d\n', p, num_images);

    % Load the images from disk.
    image_01 = double(imread(image_list_01{p}));
    image_02 = double(imread(image_list_02{p}));

% Loop over all the interrogation regions.
    for k = 1 : num_regions
        
        % Create the Gaussian fitlers
        apc_filt = exp(-(x.^2)./ (2 * apc_std_x(k)^2)) .* ...
            exp(-(y.^2)./ (2 * apc_std_y(k)^2));

        % Extract the subregions.
        region_01 = extractSubRegions(image_01,...
            [region_height, region_width], gx(k), gy(k));
        
        region_02 = extractSubRegions(image_02,...
            [region_height, region_width], gx(k), gy(k));

        % Transforms
        F1_eq = fft2(g_win .* (region_01 - mean(region_01(:))));
        F2_eq = fft2(g_win .* (region_02 - mean(region_02(:))));

        % Cross correlation
        complex_correlation_current = fftshift(F1_eq .* conj(F2_eq));
        
        complex_correlation_filt_current = phaseOnlyFilter(complex_correlation_current) .* apc_filt;
        
        spatial_corr_current = fftshift(abs(ifft2(fftshift(complex_correlation_filt_current))));
        
        spatial_correlation_array(:, :, k) = ...
            spatial_correlation_array(:, :, k) + spatial_corr_current;

    end % End (for k = 1 : num_regions)

end % End (for p = 1 : num_images)


% Do the Gaussian fitting
for k = 1 : num_regions
    
    % Particle diameter
    dp_x = 4 * abs(pi^2 / (sqrt(2) * apc_std_x(k)));
    dp_y = 4 * abs(pi^2 / (sqrt(2) * apc_std_y(k)));
    
    % Extract the correlation
    cc = spatial_correlation_array(:, :, k);
 
    % Subpixel fit
    [tx, ty ] = subpixel(cc,...
        region_width, region_height, ones(size(cc)), ...
        1, 0, [dp_x, dp_y]);
    
    TX(k) = -1 * tx;
    TY(k) = -1 * ty;
   
end


end




