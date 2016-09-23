function SPECTRAL_FILTER = calculate_apc_filter_from_image_pair(image_01, image_02, ...
    grid_y, grid_x, region_height, region_width, window_fraction, shuffle_range, shuffle_step, num_regions_neq);

% Image size
[image_height, image_width, num_channels] = size(image_01);

% Grid the image. This is temporary line.


% % Extract grid coordinates
gy = grid_y(:);
gx = grid_x(:);

% Number of regions
num_regions_eq = length(gy);

% Make the Gaussian window
g_win = gaussianWindowFilter([region_height, region_width], window_fraction, 'fraction');


% Calculate the NCC
num_regions = length(gy);

% Augment the list by a lot so we can get rid of the equal-region pairs
% and still have plenty left over
% List of the first random correlation regions
region_numbers_neq_aug_01 = round((num_regions - 1) * rand(2 * num_regions_neq, 1) + 1);

% List of the second random correlation regions
region_numbers_neq_aug_02 = round((num_regions - 1) * rand(2 * num_regions_neq, 1) + 1);

% Find the equivalent indices
equiv_inds = find(region_numbers_neq_aug_01 == region_numbers_neq_aug_02);

% Pop the equivalent indices
region_numbers_neq_aug_01(equiv_inds) = [];
region_numbers_neq_aug_02(equiv_inds) = [];

% Take the first num_ens_neq of regions
region_numbers_neq_01 = region_numbers_neq_aug_01(1 : num_regions_neq);
region_numbers_neq_02 = region_numbers_neq_aug_02(1 : num_regions_neq);

% % % % % shuffle_range_x = shuffle_range;
% % % % % shuffle_step_x = shuffle_step;
% % % % % 
% % % % % shuffle_range_y = shuffle_range;
% % % % % shuffle_step_y = shuffle_step;
% % % % % 
% % % % % % Center pixels
% % % % % xc = (region_width  + 1) / 2 + 0.5 * (1 - mod(region_width,  2));
% % % % % yc = (region_height + 1) / 2 + 0.5 * (1 - mod(region_height, 2));

% Allocate the NCC sum
ncc_sum = zeros(region_height, region_width) ...
    + 1i * zeros(region_height, region_width);

% Allocate the particle shape sum
particle_shape_sum = zeros(region_height, region_width) ...
    + 1i * zeros(region_height, region_width);

for k = 1 : num_regions_neq

    % Inform the user
    fprintf(1, 'NCC %d of %d\n', k, num_regions_neq);

    % % Extract the region from the first image
    region_neq_01 = extractSubRegions(image_01, ...
        [region_height, region_width], gx(region_numbers_neq_01(k)), gy(region_numbers_neq_01(k)));

    % Extract the region from the second image
    region_neq_02 = extractSubRegions(image_02, ...
        [region_height, region_width], gx(region_numbers_neq_02(k)), gy(region_numbers_neq_02(k)));

    % Transforms
    F1_neq = fft2(g_win .* (region_neq_01 - mean(region_neq_01(:))));
    F2_neq = fft2(g_win .* (region_neq_02 - mean(region_neq_02(:))));

    % Autocorrelations
    ac_01_neq = fftshift(F1_neq .* conj(F1_neq));
    ac_02_neq = fftshift(F2_neq .* conj(F2_neq));

    % Average auto correlationd
    ac_mean = (ac_01_neq + ac_02_neq) / 2;

    % Cross correlation
    ncc_cur = fftshift(F1_neq .* conj(F2_neq));
    
    % Particle shape current
    particle_shape_cur = abs(ac_mean);
   
    % Sum
    ncc_sum = ncc_sum + ncc_cur;
   
    % Sum particle shape
    particle_shape_sum = particle_shape_sum + particle_shape_cur;


end

% Normalize the particle so its maximum is one.
% Note that we don't want to force the minimum to zero
% Because that will cause the edges to blow up
% whenever we divide by it.
particle_shape_norm = real(particle_shape_sum ./ max(particle_shape_sum(:)));

% Allocate the CC sums
cc_abs_sum = zeros(region_height, region_width);
cc_abs_full_sum = zeros(region_height, region_width);

% Shake-the-box shfting coordinates
gx_shift = -1 * shuffle_range : shuffle_step : shuffle_range;
gy_shift = -1 * shuffle_range : shuffle_step : shuffle_range;

if isempty(gx_shift)
    gx_shift = 0;
end

if isempty(gy_shift)
    gy_shift = 0;
end

[GX, GY] = meshgrid(gx_shift, gy_shift);   
    

cc_sum_full = zeros(region_height, region_width) + ...
    1i * zeros(region_height, region_width);

% Do the corresponding correlation
for k = 1 : num_regions_eq
    
    fprintf(1, 'CCC region %d of %d\n', k, num_regions_eq);
    
    % Allocate running CC sum for the CCC
    cc_sum = zeros(region_height, region_width) + ...
    1i * zeros(region_height, region_width);
   
    % Shake-the-box
    for p = 1 : length(GX(:))
        
        shifted_grid_x = gx(k) + GX(p);
        shifted_grid_y = gy(k) + GY(p);
        
        region_eq_01 = extractSubRegions(image_01,...
            [region_height, region_width], shifted_grid_x, shifted_grid_y);
        region_eq_02 = extractSubRegions(image_02,...
            [region_height, region_width], shifted_grid_x, shifted_grid_y);

        % Transforms
        F1_eq = fft2(g_win .* (region_eq_01 - mean(region_eq_01(:))));
        F2_eq = fft2(g_win .* (region_eq_02 - mean(region_eq_02(:))));

        % Cross correlation
        cc_cur = fftshift(F1_eq .* conj(F2_eq));

        % Divide by the particle shape
        cc_eq = cc_cur ./ particle_shape_norm;

        % Ensemble corresponding correlation
        % (divided by particle shape, minus NCC)
        cc_sum = cc_sum + cc_eq;
        
        cc_sum_full = cc_sum_full + cc_cur;
       
    end
    
    % Add the sum of that trial to the sum-over-trials.
    % This means we're taking the magnitude after summing.
    cc_abs_sum = cc_abs_sum + (abs(cc_sum)).^2;
    
end

% Square root of the abs sum
cc_abs_sum_sqrt = abs(sqrt(cc_abs_sum));

% Subtract the minimum
cc_abs_sum_sub = cc_abs_sum_sqrt - min(cc_abs_sum_sqrt(:));

% Normalize
cc_abs_sum_norm = cc_abs_sum_sub ./ max(cc_abs_sum_sub(:));


% Full filter;
SPECTRAL_FILTER = cc_abs_sum_norm .* particle_shape_norm;


end




