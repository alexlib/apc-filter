function [spectral_correlation_array, spatial_correlation_array, gx, gy] = apc_save_planes(...
    image_01, image_02, grid_y, grid_x, region_size, ...
    window_fraction, save_path)

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

% Allocate the running ensemble correlation
spectral_correlation_array = ...
    zeros(region_height, region_width, num_regions) + ...
    1i * zeros(region_height, region_width, num_regions);

% Read the image if it's passed as a string
% i.e. if it's a path
if ischar(image_01)
    image_01 = double(imread(image_01));
else
    image_01 = double(image_01);
end

% Read the image if it's passed as a string
% i.e. if it's a path
if ischar(image_02)
    image_02 = double(imread(image_02));
else
    image_02 = double(image_02);
end

% Extract the region matrices from image 01
region_mat_01 = extractSubRegions(image_01,...
    [region_height, region_width], gx, gy);

% Extract the region matrices from image 02
region_mat_02 = extractSubRegions(image_02,...
    [region_height, region_width], gx, gy);

t1 = tic;
% Loop over all the interrogation regions
parfor k = 1 : num_regions
    
    % Extract the subregions.
    region_01 = region_mat_01(:, :, k);
    region_02 = region_mat_02(:, :, k);

    % Transforms
    F1 = fftshift(fft2(g_win .* (region_01 - mean(region_01(:)))));
    F2 = fftshift(fft2(g_win .* (region_02 - mean(region_02(:)))));

    % Spectral cross correlation
    spectral_cc_current = F1 .* conj(F2);

    % Save data to the output arrays
    spectral_correlation_array(:, :, k) = spectral_cc_current;

end % End (for k = 1 : num_regions)

% Save all the data
t2 = double(toc(t1));
regions_per_sec = num_regions / t2;
fprintf('Correlated %0d regions in %0.2f sec.\n', num_regions, t2);
fprintf(1, '%0.0f regions per second.\n', regions_per_sec);

% % % % % % % SAVING % % % % % % 
t1 = tic;
save(save_path, ...
    'gx', 'gy', ...
    'spectral_correlation_array', ...
    '-v6');
% % % % % % % % % % % % % % % % 

t2 = double(toc(t1));

file_info = dir(save_path);
file_size_MB = file_info.bytes / 10^6;
MB_per_sec = file_size_MB / t2;
fprintf(1, 'Saved results to: %s\n', save_path);
fprintf(1, 'Save time: %0.2f sec\n', t2);
fprintf(1, '%0.2f MB per sec\n', MB_per_sec);


end































