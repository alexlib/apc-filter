function [APC_STD_Y, APC_STD_X] = ...
    calculate_apc_filter_from_image_pair(image_01, image_02, ...
    grid_y, grid_x, region_size, window_fraction, shuffle_range, shuffle_step)
% APC_STD_Y, APC_STD_X, APC_FILTER] = ...
%     calculate_apc_filter_from_image_pair(image_01, image_02, ...
%     grid_y, grid_x, region_size, window_fraction, shuffle_range, shuffle_step)
% 
% CALCULATE_APC_FILTER_FROM_IMAGE_PAIR calculates the standard deviations
% of Gaussian-shaped spectral filters from a pair of
% images at a list of grid points. 
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


% Load images if they're specified as paths
if ischar(image_01);
    image_01 = double(imread(image_01));
end
if ischar(image_02);
    image_02 = double(imread(image_02));
end

% % Extract grid coordinates
gy = grid_y(:);
gx = grid_x(:);

% Number of regions
num_regions = length(gy);

% Size of the region
region_height = region_size(1);
region_width = region_size(2);

% Make the Gaussian window
% and save it for later.
g_win = gaussianWindowFilter(...
    [region_height, region_width], ...
    window_fraction, 'fraction');

% Default to a single shuffle range.
if length(shuffle_range) == 1
    shuffle_range_y = shuffle_range;
    shuffle_range_x = shuffle_range;
else
    shuffle_range_y = shuffle_range(1);
    shuffle_range_x = shuffle_range(2);
end

% Default to a single shuffle step.
if length(shuffle_step) == 1
    shuffle_step_y = shuffle_step;
    shuffle_step_x = shuffle_step;
else
    shuffle_step_y = shuffle_step(1);
    shuffle_step_x = shuffle_step(2);
end

% Shake-the-box shfting coordinates
gx_shift = -1 * shuffle_range_x : shuffle_step_x : shuffle_range_x;
gy_shift = -1 * shuffle_range_y : shuffle_step_y : shuffle_range_y;

% Default to not shifting
% the regions (columns)
if isempty(gx_shift)
    gx_shift = 0;
end

% Default to not shifting
% the regions (rows)
if isempty(gy_shift)
    gy_shift = 0;
end

% Make a grid of the shuffle points
% These are the locations to which
% each grid point will be shifted 
% ("shuffled"), and the origin
% of these shifts is specified
% as the grid point itself.
[shuffle_X, shuffle_Y] = meshgrid(gx_shift, gy_shift);   
   
% Allocate vectors to hold the
% paramters of the Gaussian
% functions that will be fit
% to the magnitudes of
% the correlations.
APC_STD_X = zeros(num_regions, 1);
APC_STD_Y = zeros(num_regions, 1);

% Do the corresponding correlation
for k = 1 : num_regions
    
    % Inform the user.
    fprintf(1, 'On region %d of %d\n', k, num_regions);
    
    % Allocate running correlation sum.
    % This is re-zeroed for each region
    % because it is the magnitude of 
    % this complex array that is taken
    % to be the APC filter, and a 
    % different filter is calculated for each region!
    cc_sum = zeros(region_height, region_width) + ...
        1i * zeros(region_height, region_width);
   
    % Shake-the-box
    for p = 1 : length(shuffle_X(:))
        
        % Shift the grid (this is the shuffle)
        % The Kobayashi Shake of PIV
        % The 
        shifted_grid_x = gx(k) + shuffle_X(p);
        shifted_grid_y = gy(k) + shuffle_Y(p);
        
        % Extract the subregions.
        region_01 = extractSubRegions(image_01,...
            [region_height, region_width], shifted_grid_x, shifted_grid_y);
        region_02 = extractSubRegions(image_02,...
            [region_height, region_width], shifted_grid_x, shifted_grid_y);

        % Transforms
        F1_eq = fft2(g_win .* (region_01 - mean(region_01(:))));
        F2_eq = fft2(g_win .* (region_02 - mean(region_02(:))));

        % Cross correlation
        cc_cur = fftshift(F1_eq .* conj(F2_eq));

        % Ensemble correlation
        cc_sum = cc_sum + cc_cur;
        
    end
    
    % Take the magnitude of
    % the summed correlation. 
    cc_mag = abs(cc_sum);
    
    % Fit a Gaussian to this
    % since, let's be honest, 
    % it's probably Gaussian...
    % Because what isn't these days
    [~, APC_STD_Y(k), APC_STD_X(k)] =...
    fit_gaussian_2D(cc_mag);

    subplot(1, 2, 1);
    mesh(cc_mag ./ max(cc_mag(:)), 'edgecolor', 'black');
    xlim([1, region_width]);
    ylim([1, region_height]);
    zlim([0, 1.1]);
    axis square;
    
    subplot(1, 2, 2);
    imagesc(image_01); axis image;
    hold on;
    plot(gx, gy, '.r');
    plot(gx(k), gy(k), 'oy', 'markerfacecolor', 'yellow');
    hold off
    
    drawnow;
    


    
end


end




