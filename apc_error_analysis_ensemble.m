function [ty_apc, tx_apc, ...
    ty_rpc, tx_rpc, ...
    ty_scc, tx_scc, ...
    apc_std_y, apc_std_x] = apc_error_analysis_ensemble(...
    image_list_01, image_list_02, grid_y, grid_x, region_size, ...
    window_fraction, rpc_diameter)

% This is the number of images that will be correlated.
num_images = length(image_list_01);

% Static particle diameter
dp_static = rpc_diameter / 4;

% % Extract grid coordinates
gy = grid_y(:);
gx = grid_x(:);

% Number of grid points in each direction
nx = length(unique(gx));
ny = length(unique(gy));

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

% Make the RPC filter
rpc_filter = spectralEnergyFilter(...
    region_height, region_width, rpc_diameter);

% Allocate the running ensemble correlation
spectral_correlation_array = ...
    zeros(region_height, region_width, num_regions) + ...
    1i * zeros(region_height, region_width, num_regions);

% SCC ensemble
scc_ensemble = zeros(region_height, region_width, num_regions);

% RPC ensemble
rpc_ensemble = zeros(region_height, region_width, num_regions);

% Standard deviations of the APC filters
apc_std_x = zeros(num_regions, num_images);
apc_std_y = zeros(num_regions, num_images);

% Allocate APC displacements
ty_apc = zeros(num_regions, num_images);
tx_apc = zeros(num_regions, num_images);

% Allocate RPC displacements
ty_rpc = zeros(num_regions, num_images);
tx_rpc = zeros(num_regions, num_images);

% Allocate SCC displacements
ty_scc = zeros(num_regions, num_images);
tx_scc = zeros(num_regions, num_images);

% Make the interrogation region coordinate vectors
xv2 = ((1 : region_width) - fourier_zero(region_width)) .^2;
yv2 = ((1 : region_height) - fourier_zero(region_height)).^2;

% Make the interrogation region coordinate arrays
[x2, y2] = meshgrid(xv2, yv2);

% Subpixel weighting matrix
subpix_weights = ones(region_height, region_width);

% Timer
t = tic;
% Loop over the images
for p = 1 : num_images
         
    % Inform the user.
    fprintf(1, 'On image %d of %d\n', p, num_images);

    % Load the images from disk.
    image_01 = double(imread(image_list_01{p}));
    image_02 = double(imread(image_list_02{p}));
    
    % Loop over all the interrogation regions
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
        spectral_cc_current = F1 .* conj(F2);
        
        % Spectral correlation full
        spectral_correlation_array(:, :, k) = ...
            spectral_correlation_array(:, :, k) + spectral_cc_current;
        
        % Inverse FT of the cross correlation. This is the SCC plane.
        scc_ensemble(:, :, k) = scc_ensemble(:, :, k) + ...
            fftshift(abs(ifft2(fftshift(spectral_cc_current))));
        
        % Inverse FT of the RPC. 
        rpc_ensemble(:, :, k) = rpc_ensemble(:, :, k) + ...
            fftshift(abs(ifft2(fftshift(...
            phaseOnlyFilter(spectral_cc_current) .* rpc_filter))));
         
    end % End (for k = 1 : num_regions)

    % Do the peak fitting and APC filtering
    % This is separated from the first loop
    % so that it can be run in parallel,
    % and for some reason running the xa
    % first "region" loop in parallel
    % is slower than running it serially.
    for k = 1 : num_regions
        % Inform the user
        fprintf(1, 'On image %d, APC fit %d of %d\n', p, k, num_regions);

        % Extract the data we need
        %
        % Spatial planes
        scc_spatial = scc_ensemble(:, :, k);
        rpc_spatial = rpc_ensemble(:, :, k);

        % Spectral Planes
        cc_spectral = spectral_correlation_array(:, :, k);

        % First calculate the APC filters
        % Fit a Gaussian function to the magnitude
        % of the complex correlation, 
        % which should represent the SNR versus wavenumber.
        [~, sy_apc, sx_apc] =...
            fit_gaussian_2D(abs(cc_spectral ./ max(abs(cc_spectral(:)))));
        
        % APC filter
        apc_filt = exp(-x2 ./ (2 * sx_apc^2) - y2 ./ (2 * sy_apc^2));

        % Equivalent particle diameters
        dp_equiv_x = 8 * pi^2 / sx_apc;
        dp_equiv_y = 8 * pi^2 / sy_apc;

        % Save APC standard deviations to output variables
        apc_std_x(k, p) = sx_apc;
        apc_std_y(k, p) = sy_apc;

        % Filter the spectral correlation and invert it
        apc_spatial = fftshift(abs(ifft2(fftshift(...
            phaseOnlyFilter(cc_spectral) .* apc_filt))));

        % APC Subpixel fit
        [tx_apc(k, p), ty_apc(k, p)] = subpixel(...
            apc_spatial, region_width, region_height, subpix_weights,...
            1, 0, [dp_equiv_x, dp_equiv_y]);

        % RPC Subpixel fit
        [tx_rpc(k, p), ty_rpc(k, p)] = subpixel(...
            rpc_spatial, region_width, region_height, subpix_weights,...
            1, 0, dp_static * [1, 1]);

        % SCC Subpixel fit
        [tx_scc(k, p), ty_scc(k, p)] = subpixel(...
            scc_spatial, region_width, region_height, subpix_weights,...
            1, 0, dp_static * [1, 1]);

    end % End (parfor k = 1 : num_regions);    

end % End (for p = 1 : num_images)

% End timer
tf = toc(t);

% Print elapsed time
fprintf(1, 'Elapsed time: %d seconds for %d images.\n', tf, num_images);
fprintf(1, '%0.1f seconds per image.\n\n', tf / num_images);


end































