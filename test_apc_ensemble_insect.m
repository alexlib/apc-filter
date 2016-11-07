% image_dir = '~/Desktop/images/experimental';
% image_base_name = 'vortexring_d03_f60_t06_';
% num_digits = 6;
% image_ext = '.tif';
% start_image = 840;
% end_image = 859;
% skip_image = 1;
% c_step = 1;

addpaths('..');

% Plotting garbage
% Skip vectors
skip_x = 1;
skip_y = 1;
Scale = 8;

lw = 2 ;


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
% c_step = 1;
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

image_dir = fullfile(get_image_repo, 'grasshopper', 'grasshopper_3', 'mng_2-072-B', 'raw');
image_base_name = 'mng_2-072-B';
num_digits = 6;
image_ext = '.tiff';
start_image = 1;
end_image = 100;
skip_image = 1;
c_step = 4;
trailer_a = '';
trailer_b = '';

% White field image directory
white_field_dir = fullfile(image_dir, '..');
white_field_name = 'mng_2-072-E_white_field_ave.tiff';
white_field_path = fullfile(white_field_dir, white_field_name);

% Load the white field image
white_field_image = double(imread(white_field_path));


% Number of images to process after calculating the filter.
num_images_correlate = 1;

% Ensemble length
ens_lengths = 100;

% Region sizes
region_height = 128;
region_width  = 128;

% Window fraction
window_fraction = 0.5;

% Grid spacing
grid_spacing_y = 64;
grid_spacing_x = 64;

% Shuffle
shuffle_range = [0, 0];
shuffle_step = [0, 0];

% RPC diameter
rpc_diameter = 3;

% Region size vector
region_size = [region_height, region_width];

% LIst of image numbers
image_nums_01 = start_image : skip_image : end_image;
image_nums_02 = image_nums_01 + c_step;

% Number of images
num_pairs = length(image_nums_01);

% Declare the image list
image_list_01 = {''};
image_list_02 = {''};

% Digit string
dig_str = ['%0' num2str(num_digits) 'd'];

% Form the image path lists
for k = 1 : num_pairs
   image_name_01 = [image_base_name num2str(image_nums_01(k), dig_str) trailer_a image_ext];
   image_name_02 = [image_base_name num2str(image_nums_02(k), dig_str) trailer_b image_ext];
    
   image_list_01{k} = fullfile(image_dir, image_name_01);
   image_list_02{k} = fullfile(image_dir, image_name_02);
    
end

% Load the first image and get its size
[image_height, image_width] = size(double(imread(image_list_01{1})));

% Grid the images
grid_spacing = [grid_spacing_y, grid_spacing_x];
grid_buffer_x = region_width/2  * [1, 1];
grid_buffer_y = region_height/2 * [1, 1];

% Grid the image
[grid_x, grid_y] = gridImage([image_height, image_width],...
    grid_spacing, grid_buffer_y, grid_buffer_x);

grid_x = 900; grid_y = 345;

% Calculate the APC filter
[APC_STD_Y, APC_STD_X, disp_pdf_std_dev_y, disp_pdf_std_dev_x] = ...
    calculate_apc_filter_ensemble(image_list_01, image_list_02, ...
    grid_y, grid_x, region_size,...
    window_fraction, rpc_diameter);

% Vectors for region coordinates
xv = (1 : region_width) - fourier_zero(region_width);
yv = (1 : region_height) - fourier_zero(region_height);

% Region coordinates
[x, y] = meshgrid(xv, yv);

% Calculate the RPC filter size
rpc_filter = spectralEnergyFilter(region_height, region_width, rpc_diameter);

% Spatial window
g_win = gaussianWindowFilter([region_height, region_width], window_fraction, 'fraction');

% Subpixel weighting matrix
subpix_weights = ones(region_height, region_width);

% Number of regions
num_regions = length(grid_x(:));

% Number of grid points
ny = size(unique(grid_y));
nx = size(unique(grid_x));

% Grid in matrix format
gx_mat = reshape(grid_x, [ny, nx]);
gy_mat = reshape(grid_y, [ny, nx]);

% Allocate translations
tx_scc = zeros(num_regions, num_images_correlate);
ty_scc = zeros(num_regions, num_images_correlate);
tx_rpc = zeros(num_regions, num_images_correlate);
ty_rpc = zeros(num_regions, num_images_correlate);
tx_apc = zeros(num_regions, num_images_correlate);
ty_apc = zeros(num_regions, num_images_correlate);

% Loop over the images.
for n = 1 : num_images_correlate
    
    image_01 = double(imread(image_list_01{n})) ./ white_field_image;
    image_02 = double(imread(image_list_02{n})) ./ white_field_image;
    
    image_01(isinf(image_01)) = 0;
    image_02(isinf(image_02)) = 0;
    
    % Extract the regions
    region_mat_01 = extractSubRegions(image_01, [region_height, region_width], grid_x, grid_y);
    region_mat_02 = extractSubRegions(image_02, [region_height, region_width], grid_x, grid_y);
    
    % Loop over the regions
    for k = 1 : num_regions
        
        % Inform user
        fprintf(1, 'On region %d of %d\n', k, num_regions);
        
        % APC standard deviations
        sx_apc = APC_STD_X(k);
        sy_apc = APC_STD_Y(k);
        
        % Equivalent particle diameters
        dp_equiv_apc_x = equiv_particle_diameter(sx_apc, region_width);
        dp_equiv_apc_y = equiv_particle_diameter(sy_apc, region_height);
        
        % Extract the regions
        region_01 = double(region_mat_01(:, :, k));
        region_02 = double(region_mat_02(:, :, k));
        
        % Make the spectral filter
        apc_filter = exp(-x.^2 / (2 * sx_apc^2) - y.^2 / (2 * sy_apc^2));
    
        % Fourier transforms
        F1 = fftshift(fft2(g_win .* (region_01 - mean(region_01(:)))));
        F2 = fftshift(fft2(g_win .* (region_02 - mean(region_02(:)))));
   
        % Correlation
        cc_spect = F1 .* conj(F2);
      
        % SCC plane
        scc_plane = abs(fftshift(ifft2(fftshift(cc_spect))));
        
        % RPC plane
        rpc_plane = abs(fftshift(ifft2(fftshift(phaseOnlyFilter(cc_spect) .* rpc_filter))));
        
        % APC plane
        apc_plane = abs(fftshift(ifft2(fftshift(phaseOnlyFilter(cc_spect) .* apc_filter))));
                
        [tx_scc(k, n), ty_scc(k, n)] = subpixel(scc_plane,...
        region_width, region_height, subpix_weights, ...
            3, 0, rpc_diameter * [1, 1]);

        [tx_rpc(k, n), ty_rpc(k, n)] = subpixel(rpc_plane,...
            region_width, region_height, subpix_weights, ...
            3, 0, rpc_diameter * [1, 1]);
          
        [tx_apc(k, n), ty_apc(k, n)] = subpixel(apc_plane,...
            region_width, region_height, subpix_weights, ...
            3, 0, rpc_diameter * [1, 1]); 
        
    end
    
    tx_scc_mat = reshape(tx_scc(:, n), [ny, nx]);
    ty_scc_mat = reshape(ty_scc(:, n), [ny, nx]);
    
    tx_rpc_mat = reshape(tx_rpc(:, n), [ny, nx]);
    ty_rpc_mat = reshape(ty_rpc(:, n), [ny, nx]);
    
    tx_apc_mat = reshape(tx_apc(:, n), [ny, nx]);
    ty_apc_mat = reshape(ty_apc(:, n), [ny, nx]);
    


    subplot(3, 1, 1);
    quiver(gx_mat(1 : skip_y : end, 1 : skip_x : end), ...
        gy_mat(1 : skip_y : end, 1 : skip_x : end), ...
        Scale * tx_scc_mat(1 : skip_y : end, 1 : skip_x : end), ...
        Scale * ty_scc_mat(1 : skip_y : end, 1 : skip_x : end), ...
        0, 'b', 'linewidth', lw);
    axis image
    xlim([1, image_width]);
    ylim([1, image_height]);
    set(gca, 'ydir', 'reverse');

    subplot(3, 1, 2);
    quiver(gx_mat(1 : skip_y : end, 1 : skip_x : end), ...
        gy_mat(1 : skip_y : end, 1 : skip_x : end), ...
        Scale * tx_rpc_mat(1 : skip_y : end, 1 : skip_x : end), ...
        Scale * ty_rpc_mat(1 : skip_y : end, 1 : skip_x : end), ...
        0, 'r', 'linewidth', lw);
    axis image
    xlim([1, image_width]);
    ylim([1, image_height]);
    set(gca, 'ydir', 'reverse');

    subplot(3, 1, 3);
    quiver(gx_mat(1 : skip_y : end, 1 : skip_x : end), ...
        gy_mat(1 : skip_y : end, 1 : skip_x : end), ...
        Scale * tx_apc_mat(1 : skip_y : end, 1 : skip_x : end), ...
        Scale * ty_apc_mat(1 : skip_y : end, 1 : skip_x : end), ...
        0, 'k', 'linewidth', lw);
    axis image
    xlim([1, image_width]);
    ylim([1, image_height]);
    set(gca, 'ydir', 'reverse');

    drawnow;
    
end

% Shape some grids
nx = length(unique(grid_x));
ny = length(unique(grid_y));

apc_std = sqrt(APC_STD_Y.^2 + APC_STD_X.^2);

S = reshape(apc_std, [ny, nx]);
apc_x_mat = reshape(APC_STD_X, [ny, nx]); 

figure(2); imagesc(grid_x, grid_y, apc_x_mat); axis image;


