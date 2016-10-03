function [APC_STD_Y, APC_STD_X, disp_pdf_std_dev_y, disp_pdf_std_dev_x] = ...
    calculate_apc_filter_ensemble(image_list_01, image_list_02, ...
    grid_y, grid_x, region_size, window_fraction, shuffle_range, shuffle_step)
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

% Number of shuffle locations
num_shuffles = length(shuffle_X(:));
   
% Allocate vectors to hold the
% paramters of the Gaussian
% functions that will be fit
% to the magnitudes of
% the correlations.
APC_STD_X = zeros(num_regions, 1);
APC_STD_Y = zeros(num_regions, 1);

% Allocate arrays to hold the
% displacement PDF vectors.
disp_pdf_std_dev_y = zeros(num_regions, 1);
disp_pdf_std_dev_x = zeros(num_regions, 1);


% Allocate array to hold all of the ensemble subregions. 
% This is a complex array. 
% This could be large in memory, be careful.
spectral_correlation_array = ...
    zeros(region_height, region_width, num_regions) + ...
        1i * zeros(region_height, region_width, num_regions);
    
% ALlocate array to hold the autocorrelations.
auto_correlation_array = ...
    zeros(region_height, region_width, num_regions);

% Loop over all the images.
for p = 1 : num_images
    
     % Inform the user.
    fprintf(1, 'On image %d of %d\n', p, num_images);

    % Load the images from disk.
    image_01 = double(imread(image_list_01{p}));
    image_02 = double(imread(image_list_02{p}));

% Loop over all the interrogation regions.
    for k = 1 : num_regions

        % Allocate running correlation sum.
        % This is re-zeroed for each region
        % because it is the magnitude of 
        % this complex array that is taken
        % to be the APC filter, and a 
        % different filter is calculated for each region!
        complex_correlation_sum = zeros(region_height, region_width) + ...
            1i * zeros(region_height, region_width);
        
        % Re zero the auto correlation sum 
        auto_corr_sum = zeros(region_height, region_width);

        % Shake-the-box
        for g = 1 : num_shuffles

            % Shift the grid (this is the shuffle)
            % The Kobayashi Shake of PIV
            % The 
            shifted_grid_x = gx(k) + shuffle_X(g);
            shifted_grid_y = gy(k) + shuffle_Y(g);

            % Extract the subregions.
            region_01 = extractSubRegions(image_01,...
                [region_height, region_width], shifted_grid_x, shifted_grid_y);
            region_02 = extractSubRegions(image_02,...
                [region_height, region_width], shifted_grid_x, shifted_grid_y);

            % Transforms
            F1 = fftshift(fft2(g_win .* (region_01 - mean(region_01(:)))));
            F2 = fftshift(fft2(g_win .* (region_02 - mean(region_02(:)))));

            % Cross correlation
            complex_cross_correlation_current = F1 .* conj(F2);
            
            % Auto correlations
            auto_corr_01 = F1 .* conj(F1);
            auto_corr_02 = F2 .* conj(F2);
            
            % Mean auto correlatoin
            auto_corr_mean =  (auto_corr_01 + auto_corr_02) / 2;
            
            % Add the mean auto correlation 
            % to the running sum of auto correlations
            % The auto correlation captures the shape
            % of the correlating pattern (i.e. 
            % the shape of the particle)
            % and dividing the cross correlation
            % by this array gives the FT
            % of the PDF of displacements of the 
            % correlating tracer particles.
            auto_corr_sum = auto_corr_sum + auto_corr_mean;

            % Ensemble correlation
            complex_correlation_sum = ...
                complex_correlation_sum + complex_cross_correlation_current;
            
        end % End (for g = 1 : num_shuffles)

        % Add the correlation to the big array
        % containing all of the ensemble correlation planes.
        spectral_correlation_array(:, :, k) = ...
            spectral_correlation_array(:, :, k) + complex_correlation_sum;
        
        % Add the auto correlation to the big array
        % containing all of the autocorrelations
        auto_correlation_array(:, :, k) = ...
            auto_correlation_array(:, :, k) + auto_corr_sum;
        
    end % End (for k = 1 : num_regions)

end % End (for p = 1 : num_images)

% Do the Gaussian fitting
parfor k = 1 : num_regions
    
    % Inform the user
    fprintf(1, 'Fitting region %d of %d\n', k, num_regions);
    
    % Extract data from their big arrays
    spectral_corr = spectral_correlation_array(:, :, k);
    auto_corr = auto_correlation_array(:, :, k);
    
    % This is the cross correlation divided by the auto correlation
    cc_div = spectral_corr ./ ...
        auto_corr;
 
    % Fit a Gaussian function to the remaining magnitude
    % which should represent the Fourer transform
    % of the PDF of displacements, centered about the
    % mean displacement.
    [~, ft_pdf_std_dev_y, ft_pdf_std_dev_x] =...
        fit_gaussian_2D(abs(cc_div));
    
    % Convert the PDF standard deviations in to pixels per frame.
    disp_pdf_std_dev_y(k) = pi^2 / (ft_pdf_std_dev_y);
    disp_pdf_std_dev_x(k) = pi^2 / (ft_pdf_std_dev_x);
    
    % Fit a Gaussian function to the magnitude
    % of the complex correlation, 
    % which should represent the SNR versus wavenumber.
    [~, APC_STD_Y(k), APC_STD_X(k)] =...
        fit_gaussian_2D(abs(spectral_corr));

   
end


end




