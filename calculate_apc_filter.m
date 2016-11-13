function [APC_STD_Y, APC_STD_X, APC_FILTER] = ...
    calculate_apc_filter(spectral_correlation_complex, rpc_diameter)

    % Size of the region
    [region_height, region_width] = size(spectral_correlation_complex);

    % Standard deviations of the RPC filter.
    rpc_std_dev_x = sqrt(2) / (pi * rpc_diameter) * region_width;
    rpc_std_dev_y = sqrt(2) / (pi * rpc_diameter) * region_height;

    % Fit a Gaussian function to the magnitude
    % of the complex correlation, 
    % which should represent the SNR versus wavenumber.
    [~, sy, sx] =...
        fit_gaussian_2D(abs(spectral_correlation_complex) ...
        ./ max(abs(spectral_correlation_complex(:))));
    
    % The fit can crap out and come back with
    % a standard deviation of less than 1. This is nonphysical
    % and can be used as a flag.
    if sx <= 1
        sx = rpc_std_dev_x;
    end
    if sy <= 1
        sy = rpc_std_dev_y;
    end
    
    % Take the APC diameter as the minimum
    % between the RPC equivalent std dev
    % and the standard deviation diameter.
    APC_STD_Y = min(rpc_std_dev_y, sy);
    APC_STD_X = min(rpc_std_dev_x, sx);
    
    % % Calculate the filter
    %
    % Make coordinate vectors
    xv = (1 : region_width) - fourier_zero(region_width);
    yv = (1 : region_height) - fourier_zero(region_height);
    
    % Make coordinate arrays
    [X, Y] = meshgrid(xv, yv);
    
    % Calculate filter
    APC_FILTER = exp(-X.^2 / (2 * APC_STD_X^2) - Y.^2 / (2 * APC_STD_Y^2));
    
end