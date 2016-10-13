function [ty_apc, tx_apc, ...
    ty_rpc, tx_rpc, ...
    ty_scc, tx_scc, ...
    apc_std_y, apc_std_x] = apc_error_analysis_ensemble_from_planes(...
    planes_path_list, vector_save_path_list, rpc_diameter)

% This is the number of images that were correlated
num_pairs = length(planes_path_list);

% Load the first file to count some things
load(planes_path_list{1});

% Size and number of planes
[region_height, region_width, regions_per_pair] = size(spectral_correlation_array);

% Default to RPC diameter of 3
if nargin < 3
    rpc_diameter = 3;
end

% Make the RPC filter
rpc_filter = spectralEnergyFilter(...
    region_height, region_width, rpc_diameter);

% Fit to the RPC filter
rpc_std_dev = 8 * pi^2 / rpc_diameter;

% SCC ensemble
scc_ensemble = zeros(region_height, region_width, regions_per_pair);

% RPC ensemble
rpc_ensemble = zeros(region_height, region_width, regions_per_pair);

% Complex ensemble
complex_cc_ensemble = ...
    zeros(region_height, region_width, regions_per_pair) ...
    + 1i *  zeros(region_height, region_width, regions_per_pair);

% Standard deviations of the APC filters
apc_std_x = zeros(regions_per_pair, num_pairs);
apc_std_y = zeros(regions_per_pair, num_pairs);

% Allocate APC displacements
ty_apc = zeros(regions_per_pair, num_pairs);
tx_apc = zeros(regions_per_pair, num_pairs);

% Allocate RPC displacements
ty_rpc = zeros(regions_per_pair, num_pairs);
tx_rpc = zeros(regions_per_pair, num_pairs);

% Allocate SCC displacements
ty_scc = zeros(regions_per_pair, num_pairs);
tx_scc = zeros(regions_per_pair, num_pairs);

% Make the interrogation region coordinate vectors
xv2 = ((1 : region_width) - fourier_zero(region_width)).^2;
yv2 = ((1 : region_height) - fourier_zero(region_height)).^2;

% Make the interrogation region coordinate arrays
[x2, y2] = meshgrid(xv2, yv2);

% Subpixel weighting matrix
subpix_weights = ones(region_height, region_width);

% Timer
t = tic;
% Loop over the images
for p = 1 : num_pairs
         
    % Inform the user.
    fprintf(1, 'On image %d of %d\n', p, num_pairs);

    % Load the planes data
    load(planes_path_list{p});
    
    % Add the spectral correlation to the ensemble
    complex_cc_ensemble = complex_cc_ensemble + spectral_correlation_array;

    % Loop over all the interrogation regions
    for k = 1 : regions_per_pair
       
        % Inverse FT of the cross correlation. This is the SCC plane.
        scc_ensemble(:, :, k) = scc_ensemble(:, :, k) + ...
            fftshift(abs(ifft2(fftshift(spectral_correlation_array(:, :, k)))));
        
        % Inverse FT of the RPC. 
        rpc_ensemble(:, :, k) = rpc_ensemble(:, :, k) + ...
            fftshift(abs(ifft2(fftshift(...
            phaseOnlyFilter(spectral_correlation_array(:, :, k)) .* rpc_filter))));
         
    end % End (for k = 1 : num_regions)

    % Do the peak fitting and APC filtering
    % This is separated from the first loop
    % so that it can be run in parallel,
    % and for some reason running the xa
    % first "region" loop in parallel
    % is slower than running it serially.  
    parfor k = 1 : regions_per_pair
        
        % Inform the user
        fprintf(1, 'Image %d of %d, APC fit %d of %d\n', p, num_pairs, k, regions_per_pair);

        % Extract the data we need
        %
        % Spatial planes
        scc_spatial = scc_ensemble(:, :, k);
        rpc_spatial = rpc_ensemble(:, :, k);

        % Spectral Planes
        cc_spectral = complex_cc_ensemble(:, :, k);
        
        % Fit to the correlation magnitude
        [~, sy_fit, sx_fit] =...
            fit_gaussian_2D(abs(cc_spectral));
        
        % The fit can crap out and come back with
        % a standard deviation of less than 1. This is nonphysical
        % and can be used as a flag.
        if sx_fit <= 1
            sx_fit = rpc_std_dev;
        end
        if sy_fit <= 1
            sy_fit = rpc_std_dev;
        end
        
        % Pick the minimum between the calculated
        % APC diameter and the analytical RPC diameter.
        % In other words, set the RPC diameter 
        % as the upper bound to the filter size.
        sx_apc = min(sx_fit, rpc_std_dev);
        sy_apc = min(sy_fit, rpc_std_dev);

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
            1, 0, rpc_diameter * [1, 1]);

        % SCC Subpixel fit
        [tx_scc(k, p), ty_scc(k, p)] = subpixel(...
            scc_spatial, region_width, region_height, subpix_weights,...
            1, 0, rpc_diameter * [1, 1]);

    end % End (parfor k = 1 : num_regions);
    
    % Extract the data for this pair
    tx_pair_scc = tx_scc(:, p);
    ty_pair_scc = ty_scc(:, p);
    
    tx_pair_rpc = tx_rpc(:, p);
    ty_pair_rpc = ty_rpc(:, p);
    
    tx_pair_apc = tx_apc(:, p);
    ty_pair_apc = ty_apc(:, p);
    
    apc_std_x_pair = apc_std_x(:, p);
    apc_std_y_pair = apc_std_y(:, p);
    
    % Save the data for this pair
    save(vector_save_path_list{p}, ...
        'tx_pair_scc', 'ty_pair_scc', ...
        'tx_pair_rpc', 'ty_pair_rpc', ...
        'tx_pair_apc', 'ty_pair_apc', ...
        'apc_std_x_pair', 'apc_std_y_pair', ...
        'gx', 'gy');
    
%     % Make a plot
%     apc_quiver_plots(tx_apc(:, p), ty_apc(:, p), ...
%         tx_rpc(:, p), ty_rpc(:, p), ...
%         tx_scc(:, p), ty_scc(:, p), ...
%         gx, gy, [2048, 2048], 10);
%     drawnow;

end % End (for p = 1 : num_images)

% End timer
tf = toc(t);

% Print elapsed time
fprintf(1, 'Elapsed time: %d seconds for %d images.\n', tf, regions_per_pair);
fprintf(1, '%0.1f seconds per image.\n\n', tf / regions_per_pair);


end































