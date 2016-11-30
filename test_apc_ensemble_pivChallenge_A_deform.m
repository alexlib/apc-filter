% image_dir = '~/Desktop/images/experimental';
% image_base_name = 'vortexring_d03_f60_t06_';
% num_digits = 6;
% image_ext = '.tif';
% start_image = 840;
% end_image = 859;
% skip_image = 1;
% c_step = 1;

addpaths('..');

lw = 2 ;
fSize = 16;


% image_dir = '/Users/matthewgiarra/Desktop/piv_test_images/poiseuille_diffusion_3.00/raw';
% image_base_name = 'poiseuille_diffusion_3.00_';
% num_digits = 6;
% image_ext = '.tiff';
% start_image = 1;
% end_image = 100;
% skip_image = 2;
% c_step = 1;

% image_dir = fullfile(get_image_repo, 'Aall', 'raw');
% image_base_name = 'A';
% num_digits = 3;
% image_ext = '.tif';
% start_image = 1;
% end_image = 100;
% skip_image = 1;
% c_step = 0;
% trailer_a = 'a';
% trailer_b = 'b';


% image_dir = fullfile(get_image_repo, 'Ball', 'raw');
% image_base_name = 'B';
% num_digits = 3;
% image_ext = '.bmp';
% start_image = 1;
% end_image = 100;
% skip_image = 1;
% c_step = 0;
% trailer_a = 'a';
% trailer_b = 'b';

% image_dir = fullfile(get_image_repo, 'grasshopper', 'grasshopper_3', 'mng_2-072-B', 'raw');
% image_base_name = 'mng_2-072-B';
% num_digits = 6;
% image_ext = '.tiff';
% start_image = 1;
% end_image = 100;
% skip_image = 1;imgAmask3.tif
% c_step = 3;
% trailer_a = '';
% trailer_b = '';

% image_dir = fullfile(get_image_repo, 'piv_challenge', '2014', 'piv_challenge_2014_A', 'raw');
image_dir = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/piv_test_images/pivchallenge/2014/A/images/proc/ghost';
image_base_name = 'A_deghost_';

image_dir = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/piv_test_images/pivchallenge/2014/A/images/proc/mean_sub';
image_base_name = 'A_mean_sub_';

num_digits = 5;
image_ext = '.tif';
start_image = 300;
skip_image = 1;
c_step = 0;
trailer_a = '_a';
trailer_b = '_b';

% Locate the grid mask.
mask_dir = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/piv_test_images/pivchallenge/2014/A/images/masks';
mask_name = 'imgAmask3.tif';
mask_path = fullfile(mask_dir, mask_name);

% image_dir = fullfile(get_image_repo, 'confocal', '100nm_resonant_30uL_90deg_1450um', 'tif');
% image_base_name = '100nm_resonant_30uL_90deg_1450um_';
% num_digits = 6;
% image_ext = '.tif';
% start_image = 1;
% end_image = 1999;
% skip_image = 1;
% c_step = 1;
% trailer_a = '';
% trailer_b = '';

% mask_path = '~/Desktop/A_mask.tif';
% grid_mask = double(imread(mask_path));

% Number of images to use to calculate the APC filter
num_images_filter = 50;

% Number of images to process after calculating the filter.
num_images_correlate = 10;

% Number of Passes
num_passes_spec = 1;

% RPC diameter
rpc_diameter = 8;

% Smoothing kernel diameter
smoothing_kernel_diameter = 7;

% Smoothing kernel standard deviation
smoothing_kernel_std = 1;

% Deform method
deform_method = 'interp2';

% Faux job file
Parameters.Processing(1).Region.Height = 64;
Parameters.Processing(1).Region.Width = 128;
Parameters.Processing(1).Grid.Spacing.Y = 64;
Parameters.Processing(1).Grid.Spacing.X = 64;
Parameters.Processing(1).Grid.Shift.Y = -16;
Parameters.Processing(1).Grid.Shift.X = 0;
Parameters.Processing(1).Window.Fraction.Image_01 = [0.25, 0.5];
Parameters.Processing(1).Window.Fraction.Image_02 = [0.25, 1];
% Parameters.Processing(1).Window.Fraction.Image_02 = [0.25, 0.5];

% Copy the first processing pass to the second one.
Parameters.Processing(2) = Parameters.Processing(1);

% % Create a third processing pass.
% Parameters.Processing(3).Region.Height = 32;
% Parameters.Processing(3).Region.Width = 128;
% Parameters.Processing(3).Grid.Spacing.Y = 24;
% Parameters.Processing(3).Grid.Spacing.X = 24;
% Parameters.Processing(3).Grid.Shift.Y = -16;
% Parameters.Processing(3).Grid.Shift.X = 0;
% Parameters.Processing(3).Window.Fraction = 0.5;

% Number of iterations
num_iterations = min(num_passes_spec, length(Parameters.Processing));

% Grid mask
grid_mask = double(imread(mask_path));

% Final image for filtering
end_image_filter = start_image + num_images_filter - 1;

% Final image
end_image_correlate = start_image + num_images_correlate - 1;

% List of image numbers for the ensemble
image_nums_correlate_01 = start_image : skip_image : end_image_correlate;
image_nums_correlate_02 = image_nums_correlate_01 + c_step;

% List of image numbers for the filter calculation
image_nums_filter_01 = start_image : skip_image : end_image_filter;
image_nums_filter_02 = image_nums_filter_01 + c_step;

% Number of images for filtering
num_pairs_filter = length(image_nums_filter_01);

% Number of images for correlating
num_pairs_correlate = length(image_nums_correlate_01);

% Declare the image list
image_list_01 = {''};
image_list_02 = {''};

% Digit string
dig_str = ['%0' num2str(num_digits) 'd'];

% Form the image path lists for the filter calculation
for k = 1 : num_pairs_filter
    
   % Names of the images
   image_name_01 = [image_base_name num2str(image_nums_filter_01(k), dig_str) trailer_a image_ext];
   image_name_02 = [image_base_name num2str(image_nums_filter_02(k), dig_str) trailer_b image_ext];
    
   % List of paths to the images
   image_list_filter_01{k} = fullfile(image_dir, image_name_01);
   image_list_filter_02{k} = fullfile(image_dir, image_name_02);
    
end

% Paths to the images to correlate
for k = 1 : num_pairs_correlate
    
   % Names of the images
   image_name_01 = [image_base_name num2str(image_nums_correlate_01(k), dig_str) trailer_a image_ext];
   image_name_02 = [image_base_name num2str(image_nums_correlate_02(k), dig_str) trailer_b image_ext];
    
   % List of paths to the images
   image_list_correlate_01{k} = fullfile(image_dir, image_name_01);
   image_list_correlate_02{k} = fullfile(image_dir, image_name_02);
    
end

% Load the first image and get its size
[image_height, image_width] = size(double(imread(image_list_filter_01{1})));

for p = 1 : num_iterations
    
    % Region sizes
    region_height = Parameters.Processing(p).Region.Height;
    region_width  = Parameters.Processing(p).Region.Width;

    % Window fraction
    window_fraction_01 = Parameters.Processing(p).Window.Fraction.Image_01;
    window_fraction_02 = Parameters.Processing(p).Window.Fraction.Image_02;

    % Grid spacing
    grid_spacing_y = Parameters.Processing(p).Grid.Spacing.Y;
    grid_spacing_x = Parameters.Processing(p).Grid.Spacing.X;

    % Grid shift
    grid_shift_y = Parameters.Processing(p).Grid.Shift.Y;
    grid_shift_x = Parameters.Processing(p).Grid.Shift.X;

    % Grid the images
    grid_spacing = [grid_spacing_y, grid_spacing_x];
    grid_buffer_x = region_width/2  * [1, 1];
    grid_buffer_y = region_height/2 * [1, 1];

    % Calculate the RPC filter size
    rpc_filter = spectralEnergyFilter(...
        region_height, region_width, rpc_diameter);

    % Spatial window for the first image
    g_win_01 = gaussianWindowFilter(...
        [region_height, region_width], window_fraction_01, 'fraction');
    
    g_win_02 = gaussianWindowFilter(...
        [region_height, region_width], window_fraction_02, 'fraction');

    % Grid the image
    % DO not put the mask in here
    % This is the fixed grid.
    [grid_x, grid_y] = gridImage([image_height, image_width],...
        grid_spacing, grid_buffer_y, grid_buffer_x, grid_shift_y, grid_shift_x);
    
    % Temporary grid points for debugging
%     grid_x = 1856;
%     grid_y = 440;
    
    % Save the grid to the cell structure
    gx{p} = grid_x;
    gy{p} = grid_y;

    % Number of x and y grid points
    nx = length(unique(grid_x(:)));
    ny = length(unique(grid_y(:)));

    % Number of regions
    num_regions = length(grid_x(:));

    % Allocate arrays for the ensemble
    % correlation displacements
    % for all of the regions
    % even the masked ones
    % This is so the grid is compatible
    % with UOD, deform etc.
    tx_scc{p} = nan(num_regions, num_images_correlate);
    ty_scc{p} = nan(num_regions, num_images_correlate);
    tx_rpc{p} = nan(num_regions, num_images_correlate);
    ty_rpc{p} = nan(num_regions, num_images_correlate);
    tx_apc{p} = nan(num_regions, num_images_correlate);
    ty_apc{p} = nan(num_regions, num_images_correlate);

    % Validated
    tx_val_scc{p} = nan(num_regions, num_images_correlate);
    ty_val_scc{p} = nan(num_regions, num_images_correlate);
    tx_val_rpc{p} = nan(num_regions, num_images_correlate);
    ty_val_rpc{p} = nan(num_regions, num_images_correlate);
    tx_val_apc{p} = nan(num_regions, num_images_correlate);
    ty_val_apc{p} = nan(num_regions, num_images_correlate);

    % Smoothed vector fields
    tx_smoothed_scc{p} = zeros(num_regions, num_images_correlate);
    ty_smoothed_scc{p} = zeros(num_regions, num_images_correlate);
    tx_smoothed_rpc{p} = zeros(num_regions, num_images_correlate);
    ty_smoothed_rpc{p} = zeros(num_regions, num_images_correlate);
    tx_smoothed_apc{p} = zeros(num_regions, num_images_correlate);
    ty_smoothed_apc{p} = zeros(num_regions, num_images_correlate);

    % Flags for outliers
    is_outlier_scc{p} = zeros(num_regions, num_images_correlate);
    is_outlier_rpc{p} = zeros(num_regions, num_images_correlate);
    is_outlier_apc{p} = zeros(num_regions, num_images_correlate);

    % % Make coordinates for the filters
    %
    xv = (1 : region_width) - fourier_zero(region_width);
    yv = (1 : region_height) - fourier_zero(region_height);    
    % Region coordinates
    [x, y] = meshgrid(xv, yv);

    % Subpixel weighting matrix
    subpix_weights = ones(region_height, region_width);

    % Image indices specified by the grid
    grid_inds = sub2ind([image_height, image_width], grid_y, grid_x);

    % Indices to correlate
    correlate_inds = find(grid_mask(grid_inds));

    % Grid points to correlate
    grid_x_correlate = grid_x(correlate_inds);
    grid_y_correlate = grid_y(correlate_inds);

    % Number of regions to correlate
    num_regions_correlate = length(grid_x_correlate);

    % Allocate arrays for the ensemble
    % correlation displacements
    % for only the correlated regions
    tx_valid_scc = zeros(num_regions_correlate, 1);
    ty_valid_scc = zeros(num_regions_correlate, 1);
    tx_valid_rpc = zeros(num_regions_correlate, 1);
    ty_valid_rpc = zeros(num_regions_correlate, 1);
    tx_valid_apc = zeros(num_regions_correlate, 1);
    ty_valid_apc = zeros(num_regions_correlate, 1);

    % Allocate the spatial correlation plane arrays (real valued)
    scc_ens = zeros(region_height, region_width, num_regions_correlate);
    rpc_ens = zeros(region_height, region_width, num_regions_correlate);
    
    % Allocate the spectral correlation plane arrays (complex valued)
    spc_ens = zeros(region_height, region_width, num_regions_correlate) + ...
        1i * zeros(region_height, region_width, num_regions_correlate);
    
    % Deform parameters
    if p > 1
        % On subsequent passes, make a deform structure to
        % give to the APC filter code
        deform_params.Grid.X = gx{p-1};
        deform_params.Grid.Y = gy{p-1};
        deform_params.Displacement.Y = ty_smoothed_apc{p-1}(:, end);
        deform_params.Displacement.X = tx_smoothed_apc{p-1}(:, end);
        deform_params.Method = deform_method;
    else
        % On the first pass, make an empty deform parameters structure.
        deform_params = [];
    end

    % Calculate the APC filter
    [APC_STD_Y{p}, APC_STD_X{p}] = ...
        calculate_apc_filter_ensemble_deform(image_list_filter_01, image_list_filter_02, ...
        grid_y, grid_x,...
        [region_height, region_width], window_fraction_01, window_fraction_02, rpc_diameter, grid_mask, deform_params);

    % Extract the filters to correlate
    apc_std_x_correlate = APC_STD_X{p}(correlate_inds);
    apc_std_y_correlate = APC_STD_Y{p}(correlate_inds);

    % Loop over the images.
    for n = 1 : num_images_correlate

        % Inform the user.
        fprintf(1, 'On image %d of %d\n', n, num_images_correlate);
       
        % Load the images
        image_raw_01 = double(imread(image_list_correlate_01{n}));
        image_raw_02 = double(imread(image_list_correlate_02{n}));
        
        % Deform the image if requested.
        if p > 1  
            % Source velocity
            dx_source = tx_smoothed_apc{p-1}(:, end);
            dy_source = ty_smoothed_apc{p-1}(:, end);
            gx_source = gx{p-1};
            gy_source = gy{p-1};
            
            % Deform the images
            [image_01, image_02] = deform_image_pair(...
            image_raw_01, image_raw_02, ...
            gx_source, gy_source,...
            dx_source, dy_source);
        else
            % If deform wasn't specified
            % then use the raw images. 
            image_01 = image_raw_01;
            image_02 = image_raw_02;
        end
        
        % Extract the regions
        region_mat_01 = extractSubRegions(image_01, [region_height, region_width], grid_x_correlate, grid_y_correlate);
        region_mat_02 = extractSubRegions(image_02, [region_height, region_width], grid_x_correlate, grid_y_correlate);

        % Loop over the regions
        for k = 1 : num_regions_correlate
            
            % Determine which overall region is being correlated
            global_ind = correlate_inds(k);

            % APC standard deviations
            sx_apc = apc_std_x_correlate(k);
            sy_apc = apc_std_y_correlate(k);

            % Equivalent particle diameters
            dp_equiv_apc_x = equiv_particle_diameter(sx_apc, region_width);
            dp_equiv_apc_y = equiv_particle_diameter(sy_apc, region_height);

            % Extract the regions
            region_01 = double(region_mat_01(:, :, k));
            region_02 = double(region_mat_02(:, :, k));

            % Make the spectral filter
            apc_filter = exp(-x.^2 / (2 * sx_apc^2) - y.^2 / (2 * sy_apc^2));

            % Fourier transforms
            F1 = fftshift(fft2(g_win_01 .* (region_01 - mean(region_01(:)))));
            F2 = fftshift(fft2(g_win_02 .* (region_02 - mean(region_02(:)))));

            % Correlation
            cc_spect = F1 .* conj(F2);

            % SCC plane
            scc_plane = abs(fftshift(ifft2(fftshift(cc_spect))));

            % RPC plane
            rpc_plane = abs(fftshift(ifft2(fftshift(phaseOnlyFilter(cc_spect) .* rpc_filter))));

            % Ensemble planes
            scc_ens(:, :, k) = scc_ens(:, :, k) + scc_plane;
            rpc_ens(:, :, k) = rpc_ens(:, :, k) + rpc_plane;
            spc_ens(:, :, k) = spc_ens(:, :, k) + cc_spect;

            % APC plane
            apc_plane = abs(fftshift(ifft2(fftshift(phaseOnlyFilter(spc_ens(:, :, k)) .* apc_filter))));

            % Calculate the SCC displacement
            % from the ensemble SCC plane
            % % % % REPLACE WITH CORRECT INDEX
            [tx_valid_scc(k), ty_valid_scc(k)] = subpixel(scc_ens(:, :, k),...
            region_width, region_height, subpix_weights, ...
                1, 0, rpc_diameter * [1, 1]);

            % Calculate the RPC displacement
            % from the ensemble RPC plane
            [tx_valid_rpc(k), ty_valid_rpc(k)] = subpixel(rpc_ens(:, :, k),...
                region_width, region_height, subpix_weights, ...
                1, 0, rpc_diameter * [1, 1]);

            % Calculate the APC displacement
            % from the ensemble APC plane
            [tx_valid_apc(k), ty_valid_apc(k)] = subpixel(apc_plane,...
                region_width, region_height, subpix_weights, ...
                1, 0, [dp_equiv_apc_x, dp_equiv_apc_y]);        
        end

        % Put the unmasked displacement estimates
        % into the overall displacement array
        tx_scc{p}(correlate_inds, n) = tx_valid_scc;
        ty_scc{p}(correlate_inds, n) = ty_valid_scc;
        tx_rpc{p}(correlate_inds, n) = tx_valid_rpc;
        ty_rpc{p}(correlate_inds, n) = ty_valid_rpc;
        tx_apc{p}(correlate_inds, n) = tx_valid_apc;
        ty_apc{p}(correlate_inds, n) = ty_valid_apc;

        % Validate the vector fields
        [tx_val_scc{p}(:, n), ty_val_scc{p}(:, n), is_outlier_scc{p}(:, n)] = validateField_prana(grid_x, grid_y, tx_scc{p}(:, n), ty_scc{p}(:, n));
        [tx_val_rpc{p}(:, n), ty_val_rpc{p}(:, n), is_outlier_rpc{p}(:, n)] = validateField_prana(grid_x, grid_y, tx_rpc{p}(:, n), ty_rpc{p}(:, n));
        [tx_val_apc{p}(:, n), ty_val_apc{p}(:, n), is_outlier_apc{p}(:, n)] = validateField_prana(grid_x, grid_y, tx_apc{p}(:, n), ty_apc{p}(:, n));

        % Smooth the vector fields
        tx_smoothed_scc{p}(:, n) = smoothField(tx_val_scc{p}(:, n), smoothing_kernel_diameter, smoothing_kernel_std);
        ty_smoothed_scc{p}(:, n) = smoothField(ty_val_scc{p}(:, n), smoothing_kernel_diameter, smoothing_kernel_std);
        tx_smoothed_rpc{p}(:, n) = smoothField(tx_val_rpc{p}(:, n), smoothing_kernel_diameter, smoothing_kernel_std);
        ty_smoothed_rpc{p}(:, n) = smoothField(ty_val_rpc{p}(:, n), smoothing_kernel_diameter, smoothing_kernel_std);
        tx_smoothed_apc{p}(:, n) = smoothField(tx_val_apc{p}(:, n), smoothing_kernel_diameter, smoothing_kernel_std);
        ty_smoothed_apc{p}(:, n) = smoothField(ty_val_apc{p}(:, n), smoothing_kernel_diameter, smoothing_kernel_std);

    end
end

for p = 1 : num_iterations
   
    nx = length(unique(gx{p}));
    ny = length(unique(gy{p}));
    
    gx_mat = reshape(gx{p}, [ny, nx]);
    gy_mat = reshape(gy{p}, [ny, nx]);
    
    tx_mat_scc = reshape(tx_scc{1}(:, end), [ny, nx]);
    ty_mat_scc = reshape(ty_scc{1}(:, end), [ny, nx]);
    
    tx_mat_rpc = reshape(tx_rpc{1}(:, end), [ny, nx]);
    ty_mat_rpc = reshape(ty_rpc{1}(:, end), [ny, nx]);
    
    tx_mat_apc = reshape(tx_apc{1}(:, end), [ny, nx]);
    ty_mat_apc = reshape(ty_apc{1}(:, end), [ny, nx]);
    
end

% Mask points
[mask_y, mask_x] = get_mask_outline(grid_mask);

dx_plot = 5;

gv = [-15, 16];

Skip_x = 1;
Skip_y = 1;
Scale = 4;

lw_mesh = 1;

scc_plot = scc_ens(1 : Skip : end, 1 : Skip : end);
rpc_plot = rpc_ens(1 : Skip : end, 1 : Skip : end);
apc_plot = apc_plane(1 : Skip : end, 1 : Skip : end);

% First and last rows to plot
min_row = min(gy_mat(:));
max_row = max(gy_mat(:));
min_col = min(gx_mat(:));
max_col = max(gx_mat(:));

subplot(3, 1, 1);
imagesc(gx_mat(:), gy_mat(:), tx_mat_scc);
axis image;
hold on;
quiver(gx_mat(1 : Skip_y : end, 1 : Skip_x : end), ...
       gy_mat(1 : Skip_y : end, 1 : Skip_x : end), ...
       Scale * tx_mat_scc(1 : Skip_y : end, 1 : Skip_x : end), ...
       Scale * ty_mat_scc(1 : Skip_y : end, 1 : Skip_x : end), ...
       0, 'black', 'linewidth', 1.5);
plot(mask_x, mask_y, 'ok');   
hold off
xlim([min_col, max_col]);
yl = ylim;

% First and last rows to plot
min_row = min(gy_mat(:));
max_row = max(gy_mat(:));
min_col = min(gx_mat(:));
max_col = max(gx_mat(:));

subplot(3, 1, 2);
imagesc(gx_mat(:), gy_mat(:), tx_mat_rpc);
axis image;
hold on;
quiver(gx_mat(1 : Skip_y : end, 1 : Skip_x : end), ...
       gy_mat(1 : Skip_y : end, 1 : Skip_x : end), ...
       Scale * tx_mat_rpc(1 : Skip_y : end, 1 : Skip_x : end), ...
       Scale * ty_mat_rpc(1 : Skip_y : end, 1 : Skip_x : end), ...
       0, 'black', 'linewidth', 1.5);
plot(mask_x, mask_y, 'ok');   
hold off
xlim([min_col, max_col]);
ylim(yl);

subplot(3, 1, 3);
% First and last rows to plot
min_row = min(gy_mat(:));
max_row = max(gy_mat(:));
min_col = min(gx_mat(:));
max_col = max(gx_mat(:));

% subplot(3, 1, 1);
imagesc(gx_mat(:), gy_mat(:), tx_mat_apc);
axis image;
hold on;
quiver(gx_mat(1 : Skip_y : end, 1 : Skip_x : end), ...
       gy_mat(1 : Skip_y : end, 1 : Skip_x : end), ...
       Scale * tx_mat_apc(1 : Skip_y : end, 1 : Skip_x : end), ...
       Scale * ty_mat_apc(1 : Skip_y : end, 1 : Skip_x : end), ...
       0, 'black', 'linewidth', 1.5);
plot(mask_x, mask_y, 'ok');   
hold off
xlim([min_col, max_col]);
ylim(yl);








% subplot(3, 1, 1);
% % surf(scc_plot ./ max(scc_plot(:)), 'linewidth', 1E-5);
% mesh(scc_plot ./ max(scc_plot(:)), 'linewidth', lw_mesh, 'edgecolor', 'black');
% axis square;
% xlim([1, region_width / Skip]);
% ylim([1, region_height / Skip]);
% zlim([0, 1.1]);
% axis off;
% % set(gca, 'units', 'pixels');
% p1 = get(gca, 'position');
% set(gca, 'view', gv);
% 
% subplot(3, 1, 2);
% % surf(rpc_plot ./ max(rpc_plot(:)), 'linewidth', 1E-5);
% mesh(rpc_plot ./ max(rpc_plot(:)), 'linewidth', lw_mesh, 'edgecolor', 'black');
% axis square;
% xlim([1, region_width / Skip]);
% ylim([1, region_height / Skip]);
% zlim([0, 1.1]);
% axis off;
% % set(gca, 'units', 'pixels');
% p2 = get(gca, 'position');
% p2(2) = p1(2) - p1(4) + dx_plot;
% % set(gca, 'position', p2);
% set(gca, 'view', gv);
% 
% 
% subplot(3, 1, 3);
% % surf(apc_plot ./ max(apc_plot(:)), 'linewidth', 1E-5);
% mesh(apc_plot ./ max(apc_plot(:)), 'linewidth', lw_mesh, 'edgecolor', 'black');
% axis square;
% xlim([1, region_width / Skip]);
% ylim([1, region_height / Skip]);
% zlim([0, 1.1]);
% axis off;
% % set(gca, 'units', 'pixels');
% p3 = get(gca, 'position');
% p3(2) = p2(2) - p2(4) + dx_plot;
% % set(gca, 'position', p3);
% set(gca, 'view', gv);
% set(gcf, 'color', 'white');


% =======

















