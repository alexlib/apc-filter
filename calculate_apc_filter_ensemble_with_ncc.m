function [APC_STD_Y, APC_STD_X, DX_STD_DEV_Y, DX_STD_DEV_X] = ...
    calculate_apc_filter_ensemble_with_ncc(image_list_01, image_list_02, ...
    grid_y, grid_x, region_size, window_fraction, rpc_diameter, shuffle_range, shuffle_step, num_ncc)
% APC_STD_Y, APC_STD_X, APC_FILTER] = ...
%     calculate_apc_filter_from_image_pair(image_01, image_02, ...
%     grid_y, grid_x, region_size, window_fraction, shuffle_range, shuffle_step)
% 
% CALCULATE_APC_FILTER_FROM_IMAGE_PAIR calculates the standard deviations
% of Gaussian-shaped spectral filters from a pair of images
%
% INPUTS
%     image_list_01: This is a cell array containing the paths
%       to images the filesystem. Each element
%       in the cell array is a string that specifies
%       the path to an image. The images refered to
%       in image_list_01 specify the first image
%       in each pair of images. See below for example usage.
% 
%     image_02: Same as image_list_02, but specifying
%       the paths to the second image in each pair.
% 
%     grid_y: 2D Array or 1D vector containing the pixel
%         row coordinates of the interrogation regions
%         at which the filters will be calculated.
%         Note: Grids can be created using the function gridImage()
% 
%     grid_x: 2D Array or 1D vector containing the pixel
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
%         The window will be created using the function
%         gaussianWindowFilter().
%         For example, to create a 50% window (i.e., a 
%         64x64 window on a 128x128 region), 
%         then set window_fraction = 0.5 * [1, 1];
% 
%     shuffle_range: [2 x 1] vector specifying the row and column 
%         range across which each interrogation is shifted for calculating
%         the filter. This is the range of the shake-the-box analog for this
%         method. 
%         shuffle_range = [shuffle_range_rows, shuffle_rows_cols];
%         If no value is input then this defaults to zero.
% 
%     shuffle_step: [2 x 1] vector specifying the row and column 
%         spacing between shifted locations of regions for calculating
%         the filter. This is the spacing of the of the shake-the-box
%         analog for this method. 
%         shuffle_step = [shuffle_step_rows, shuffle_step_cols];
%         % If no value is input then this defaults to zero.
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
%     DX_STD_DEV_Y = Vector containing the standard deviations
%       of the distribution of vertical displacements in
%       the interrogation regions centered about each grid point.
%       This code assumes that distributions are Gaussian.
%
%     DX_STD_DEV_X = Vector containing the standard deviations
%       of the distribution of horizontal displacements in
%       the interrogation regions centered about each grid point.
%       This code assumes that distributions are Gaussian.
% 
% SEE ALSO
%     gridImage; subpixel, fit_gaussian_2D, extractSubRegions
% 

% Default to not shuffling
% the interrogation regions.
if nargin < 9
    shuffle_step = 0;
end

% Default to not shuffling
% the interrogation regions.
if nargin < 8
    shuffle_range = 0;
end;

% Default to rpc diameter of 3 pixels
if nargin < 7
    rpc_diameter = 3;
end

% Convert the RPC diameter to the std dev of the filter
rpc_std_dev = 8 * pi^2 / rpc_diameter;

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
DX_STD_DEV_Y = zeros(num_regions, 1);
DX_STD_DEV_X = zeros(num_regions, 1);

% Figure out lots of possible
% combinations of pairs
% for the NCC
img_num_ncc_01 = ceil((num_images - 1) * rand(2 * num_ncc-1, 1));
img_num_ncc_02 = ceil((num_images - 1) * rand(2 * num_ncc-1, 1));

% index of the image numbers
img_index = (sub2ind(...
    [num_images, num_images], img_num_ncc_01, img_num_ncc_02));

% Number of NCC to run
num_ncc_total = min(length(img_index), num_ncc);

% Which ones to use
img_index_ncc = img_index(1 : num_ncc_total);

% Number of NCC to run
num_ncc_total = length(img_index_ncc);

% Allocate array to hold all of the ensemble subregions. 
% This is a complex array. 
% This could be large in memory, be careful.
spectral_correlation_array = ...
    zeros(region_height, region_width, num_regions) + ...
        1i * zeros(region_height, region_width, num_regions);
    
% ALlocate array to hold the autocorrelations.
auto_correlation_array = ...
    zeros(region_height, region_width, num_regions);

% NCC array
ncc_array = zeros(region_height, region_width, num_regions) + ...
        1i * zeros(region_height, region_width, num_regions);
    
% Do the NCC
for p = 1 : num_ncc_total
    
     % Inform the user.
    fprintf(1, 'On NCC %d of %d\n', p, num_ncc_total);

    % Numbers of the images
    [img_num_01, img_num_02] = ind2sub(...
        [num_images, num_images], img_index_ncc(p));
    
    % Load the images from disk.
    image_01 = double(imread(image_list_01{img_num_01}));
    image_02 = double(imread(image_list_02{img_num_02}));

% Loop over all the interrogation regions.
    for k = 1 : num_regions

        % Extract the subregions.
        region_01 = extractSubRegions(image_01,...
            [region_height, region_width], gx(k), gy(k));
        region_02 = extractSubRegions(image_02,...
            [region_height, region_width], gx(k), gy(k));
        
        % Transforms
        F1 = fftshift(fft2(g_win .* (region_01 - mean(region_01(:)))));
        F2 = fftshift(fft2(g_win .* (region_02 - mean(region_02(:)))));
                
        % Cross correlation
        ncc_current = F1 .* conj(F2);
        
        % This is the NCC
        ncc_array(:, :, k) = ...
            ncc_array(:, :, k) + ncc_current;
        
    end % End (for k = 1 : num_regions)

end % End (for p = 1 : num_images)

% Divide the NCC by the number of images
ncc_array_div = ncc_array ./ num_ncc_total;
    

% Loop over all the images.
for p = 1 : num_images
    
     % Inform the user.
    fprintf(1, 'On image %d of %d\n', p, num_images);

    % Load the images from disk.
    image_01 = double(imread(image_list_01{p}));
    image_02 = double(imread(image_list_02{p}));

% Loop over all the interrogation regions.
    for k = 1 : num_regions

        % Extract the subregions.
        region_01 = extractSubRegions(image_01,...
            [region_height, region_width], gx(k), gy(k));
        region_02 = extractSubRegions(image_02,...
            [region_height, region_width], gx(k), gy(k));
        
        % Transforms
        F1 = fftshift(fft2(g_win .* (region_01 - mean(region_01(:)))));
        F2 = fftshift(fft2(g_win .* (region_02 - mean(region_02(:)))));
                
        % Cross correlation
        complex_cross_correlation_current = F1 .* conj(F2);
        
        % Cross correlation minus the NCC
%         ccc = complex_cross_correlation_current - real(ncc_array_div);
        
        % This is the spectral correlation
        spectral_correlation_array(:, :, k) = ...
            spectral_correlation_array(:, :, k) + complex_cross_correlation_current;

       ca = real(spectral_correlation_array(:, :, k));
       surf(ca ./ max(ca(:)));
       axis square;
       xlim([1 region_width]);
       ylim([1, region_height]);
       zlim([-1, 1]);
       drawnow;
        
    end % End (for k = 1 : num_regions)

end % End (for p = 1 : num_images)

% Normalize the spectral corr
spc_norm = spectral_correlation_array ./ num_images;

% Subtract the NCC
for p = 1 : num_regions
    ncc_cur = ncc_array_div(:, :, k);
    ncc_cur_norm = ncc_cur ./ max(real(ncc_cur(:)));
    cc_cur = spc_norm(:, :, k);
    cc_cur_norm = cc_cur ./ max(cc_cur(:));
    
    N = sqrt(max(cc_cur(:)));
    
    ncc_scaled = ncc_cur_norm * (N^2 - N);
    
    ncc_ratio = max(real(ncc_cur(:))) / max(real(cc_cur(:)));
     
    
end


% Do the Gaussian fitting
for k = 1 : num_regions
    
    % Inform the user
    fprintf(1, 'Fitting region %d of %d\n', k, num_regions);
    
    % Extract data from their big arrays
    spectral_corr = spectral_correlation_array(:, :, k);
    auto_corr = auto_correlation_array(:, :, k);
    
    % Fit a Gaussian function to the magnitude
    % of the complex correlation, 
    % which should represent the SNR versus wavenumber.
    [~, sy, sx, ~ , A] =...
        fit_gaussian_2D(abs(spectral_corr));
    
    surf(A);
    
    % The fit can crap out and come back with
    % a standard deviation of less than 1. This is nonphysical
    % and can be used as a flag.
    if sx <= 1
        sx = rpc_std_dev;
    end
    if sy <= 1
        sy = rpc_std_dev;
    end
    
    % Take the APC diameter as the minimum
    % between the RPC equivalent std dev
    % and the standard deviation diameter.
    APC_STD_Y(k) = min(rpc_std_dev, sy);
    APC_STD_X(k) = min(rpc_std_dev, sx);
    
    % This is the cross correlation divided by the auto correlation
    cc_div = spectral_corr ./ auto_corr;
 
    % Fit a Gaussian function to the remaining magnitude
    % which should represent the Fourer transform
    % of the PDF of displacements, centered about the
    % mean displacement.
    [~, ft_pdf_std_dev_y, ft_pdf_std_dev_x] =...
        fit_gaussian_2D(abs(cc_div));
    
    % Convert the PDF standard deviations in to pixels per frame.
    DX_STD_DEV_Y(k) = pi^2 / (ft_pdf_std_dev_y);
    DX_STD_DEV_X(k) = pi^2 / (ft_pdf_std_dev_x);
    
end

end




