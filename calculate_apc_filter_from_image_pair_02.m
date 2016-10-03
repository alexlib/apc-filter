function [APC_STD_Y, APC_STD_X] = ...
    calculate_apc_filter_from_image_pair_02(image_01, image_02, ...
    grid_y, grid_x, region_size, window_fraction, shuffle_range, shuffle_step)

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

% Shake-the-box shfting coordinates
gx_shift = -1 * shuffle_range : shuffle_step : shuffle_range;
gy_shift = -1 * shuffle_range : shuffle_step : shuffle_range;

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
        region_eq_01 = extractSubRegions(image_01,...
            [region_height, region_width], shifted_grid_x, shifted_grid_y);
        region_eq_02 = extractSubRegions(image_02,...
            [region_height, region_width], shifted_grid_x, shifted_grid_y);

        % Transforms
        F1_eq = fft2(g_win .* (region_eq_01 - mean(region_eq_01(:))));
        F2_eq = fft2(g_win .* (region_eq_02 - mean(region_eq_02(:))));

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
    
end


end




