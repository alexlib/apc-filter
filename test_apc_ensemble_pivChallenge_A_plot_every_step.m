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
% skip_image = 1;
% c_step = 3;
% trailer_a = '';
% trailer_b = '';

% image_dir = fullfile(get_image_repo, 'piv_challenge', '2014', 'piv_challenge_2014_A', 'raw');
image_dir = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/piv_test_images/piv_challenge/2014/piv_challenge_2014_A/proc/min_sub';
image_base_name = 'A_';
num_digits = 5;
image_ext = '.tif';
start_image = 300;
end_image = 309;
skip_image = 1;
c_step = 0;
trailer_a = '_a';
trailer_b = '_b';
mask_path = fullfile(image_dir, '..', '..', 'mask', 'A_mask.tif');

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


% Number of images to process after calculating the filter.
num_images_correlate = 10;

% Region sizes
region_height = 48;
region_width  = 192;

% Window fraction
window_fraction = 0.5;

% Grid spacing
grid_spacing_y = 24;
grid_spacing_x = 48;

% Shuffle
shuffle_range = [0, 0];
shuffle_step = [0, 0];

% RPC diameter
rpc_diameter = 8;

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

% Grid buffer
grid_buffer_x = 0 * [1, 1];
grid_buffer_y = 0 * [1, 1];

% Grid shift
grid_shift_x = 0;
grid_shift_y = -16;

% Grid mask
grid_mask = double(imread(mask_path));

% Grid the image
[gx, gy] = gridImage([image_height, image_width],...
    grid_spacing, grid_buffer_y, grid_buffer_x, grid_mask);

grid_x = gx + grid_shift_x;
grid_y = gy + grid_shift_y;

figure(1);
image_01 = imread(image_list_01{1});
imagesc(image_01); axis image; colormap gray;
caxis([2^7, 2^8]);
hold on;
plot(grid_x, grid_y, 'oy', 'markerfacecolor', 'y');
hold off
drawnow;

% % Number of grid points in each direction
% nx = length(unique(grid_x));
% ny = length(unique(grid_y));
% 


% grid_x = 1465; grid_y = 995;

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
% 
% % Number of grid points
% ny = size(unique(grid_y));
% nx = size(unique(grid_x));
% 
% % Grid in matrix format
% gx_mat = reshape(grid_x, [ny, nx]);
% gy_mat = reshape(grid_y, [ny, nx]);

% Allocate translations
tx_scc = zeros(num_regions, num_images_correlate);
ty_scc = zeros(num_regions, num_images_correlate);
tx_rpc = zeros(num_regions, num_images_correlate);
ty_rpc = zeros(num_regions, num_images_correlate);
tx_apc = zeros(num_regions, num_images_correlate);
ty_apc = zeros(num_regions, num_images_correlate);

% Allocate plane arrays
scc_ens = zeros(region_height, region_width, num_regions);
rpc_ens = zeros(region_height, region_width, num_regions);
spc_ens = zeros(region_height, region_width, num_regions) + ...
    1i * zeros(region_height, region_width, num_regions);

tx_ens_scc = zeros(num_regions, 1);
ty_ens_scc = zeros(num_regions, 1);
tx_ens_rpc = zeros(num_regions, 1);
ty_ens_rpc = zeros(num_regions, 1);
tx_ens_apc = zeros(num_regions, 1);
ty_ens_apc = zeros(num_regions, 1);


% Loop over the images.
for n = 1 : num_images_correlate
    
    fprintf(1, 'On image %d of %d\n', n, num_images_correlate);
    
    % Extract the regions
    region_mat_01 = extractSubRegions(image_list_01{n}, [region_height, region_width], grid_x, grid_y);
    region_mat_02 = extractSubRegions(image_list_02{n}, [region_height, region_width], grid_x, grid_y);
    
    % Loop over the regions
    parfor k = 1 : num_regions
        
        % Inform user
        
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
        
        
        
        % Ensemble planes
        scc_ens(:, :, k) = scc_ens(:, :, k) + scc_plane;
        rpc_ens(:, :, k) = rpc_ens(:, :, k) + rpc_plane;
        spc_ens(:, :, k) = spc_ens(:, :, k) + cc_spect;
        
        
        
        % APC plane
        apc_plane = abs(fftshift(ifft2(fftshift(phaseOnlyFilter(cc_spect) .* apc_filter))));
        
        
                
        [tx_scc(k, n), ty_scc(k, n)] = subpixel(scc_plane,...
        region_width, region_height, subpix_weights, ...
            1, 0, rpc_diameter * [1, 1]);
 
        [tx_rpc(k, n), ty_rpc(k, n)] = subpixel(rpc_plane,...
            region_width, region_height, subpix_weights, ...
            1, 0, rpc_diameter * [1, 1]);
%           
        [tx_apc(k, n), ty_apc(k, n)] = subpixel(apc_plane,...
            region_width, region_height, subpix_weights, ...
            1, 0, [dp_equiv_apc_x, dp_equiv_apc_y]); 
        

        
    end
    


end



% Do the APC fits and calculate the displacements
parfor k = 1 : num_regions
    [tx_ens_scc(k), ty_ens_scc(k)] = subpixel(scc_ens(:, :, k),...
        region_width, region_height, subpix_weights, ...
            1, 0, rpc_diameter * [1, 1]);
        
    [tx_ens_rpc(k), ty_ens_rpc(k)] = subpixel(rpc_ens(:, :, k),...
        region_width, region_height, subpix_weights, ...
            1, 0, rpc_diameter * [1, 1]);
        
    % APC diameters
    sx_apc = APC_STD_X(k);
    sy_apc = APC_STD_Y(k);
    dp_equiv_apc_x = equiv_particle_diameter(sx_apc, region_width);
    dp_equiv_apc_y = equiv_particle_diameter(sy_apc, region_height);
    
    
     % Make the spectral filter
    apc_filter = exp(-x.^2 / (2 * sx_apc^2) - y.^2 / (2 * sy_apc^2));    
    % Filter and invert the ensemble spectral plane
    apc_plane = abs(fftshift(ifft2(fftshift(phaseOnlyFilter(spc_ens(:, :, k)) .* apc_filter))));
    
    [tx_ens_apc(k), ty_ens_apc(k)] = subpixel(apc_plane,...
            region_width, region_height, subpix_weights, ...
            1, 0, [dp_equiv_apc_x, dp_equiv_apc_y]); 
    
end

% % Shape some grids
% nx = length(unique(grid_x));
% ny = length(unique(grid_y));
% 
apc_std = sqrt(APC_STD_Y.^2 + APC_STD_X.^2);

% S = reshape(apc_std, [ny, nx]);
% apc_x_mat = reshape(APC_STD_X, [ny, nx]); 
% apc_y_mat = reshape(APC_STD_Y, [ny, nx]);

% tx_ens_scc_mat = reshape(tx_ens_scc, [ny, nx]);
% ty_ens_scc_mat = reshape(ty_ens_scc, [ny, nx]);
% tx_ens_rpc_mat = reshape(tx_ens_rpc, [ny, nx]);
% ty_ens_rpc_mat = reshape(ty_ens_rpc, [ny, nx]);
% tx_ens_apc_mat = reshape(tx_ens_apc, [ny, nx]);
% ty_ens_apc_mat = reshape(ty_ens_apc, [ny, nx]);

% Plotting garbage
% Skip vectors
Skip = 2;
Scale = 4;
c_max = 45;

figure(2);
subplot(3, 2, 1);
quiver(grid_x(1 : Skip : end), grid_y(1 : Skip : end), ...
    Scale * tx_ens_scc(1 : Skip : end), ...
    Scale * ty_ens_scc(1 : Skip : end), ...
    0, 'b', 'linewidth', lw);
axis image
xlim([1, image_width]);
ylim([1, image_height]);
set(gca, 'ydir', 'reverse');
title(sprintf('$\\textrm{SCC}, \\, %d \\,\\, \\textrm{pairs}$', num_images_correlate), 'FontSize', fSize, 'interpreter', 'latex');
set(gca, 'FontSize', fSize)
box on


subplot(3, 2, 3);

quiver(grid_x(1 : Skip : end), grid_y(1 : Skip : end), ...
    Scale * tx_ens_rpc(1 : Skip : end), ...
    Scale * ty_ens_rpc(1 : Skip : end), ...
    0, 'r', 'linewidth', lw);
axis image
xlim([1, image_width]);
ylim([1, image_height]);
set(gca, 'ydir', 'reverse');
title(sprintf('$\\textrm{RPC}, d_\\mathrm{RPC} = %0.1f,\\, %d  \\,\\, \\textrm{pairs}$', rpc_diameter, num_images_correlate), 'FontSize', fSize, 'interpreter', 'latex')
set(gca, 'FontSize', fSize)
box on


subplot(3, 2, 5 ,'FontSize', fSize);
quiver(grid_x(1 : Skip : end), grid_y(1 : Skip : end), ...
    Scale * tx_ens_apc(1 : Skip : end), ...
    Scale * ty_ens_apc(1 : Skip : end), ...
    0, 'k', 'linewidth', lw);
axis image
xlim([1, image_width]);
ylim([1, image_height]);
set(gca, 'ydir', 'reverse');
drawnow;
title(sprintf('$\\textrm{APC}, \\, %d \\,\\, \\textrm{pairs}$', num_images_correlate), 'FontSize', fSize, 'interpreter', 'latex')
set(gca, 'FontSize', fSize)
box on



% Full grids
% Grid the image
[gx_full_full, gy_full_full] = gridImage([image_height, image_width],...
    grid_spacing, grid_buffer_y, grid_buffer_x);

gx_full = gx_full_full + grid_shift_x;
gy_full = gy_full_full + grid_shift_y;

tx_scc_full = zeros(size(gx_full));
tx_rpc_full = zeros(size(gx_full));
tx_apc_full = zeros(size(gx_full));

for k = 1 : length(grid_x)
    
    ind = find(gx_full == grid_x(k) & gy_full == grid_y(k));
    
    tx_scc_full(ind) = tx_ens_scc(k);
    tx_rpc_full(ind) = tx_ens_rpc(k);
    tx_apc_full(ind) = tx_ens_apc(k);
end

nx = length(unique(gx_full));
ny = length(unique(gy_full));
tx_scc_full_mat = reshape(tx_scc_full, [ny, nx]);
tx_rpc_full_mat = reshape(tx_rpc_full, [ny, nx]);
tx_apc_full_mat = reshape(tx_apc_full, [ny, nx]);


subplot(3, 2, 2);
imagesc(gx_full, gy_full, ...
    tx_scc_full_mat);
axis image
xlim([1, image_width]);
ylim([1, image_height]);
set(gca, 'ydir', 'reverse');
title(sprintf('$\\Delta x \\left( \\textrm{SCC} \\right), \\, \\, %d \\,\\, \\textrm{pairs}$', num_images_correlate), 'FontSize', fSize, 'interpreter', 'latex');
set(gca, 'FontSize', fSize)
caxis([0, c_max]); 
c = get(gca, 'position');
h = colorbar;
ylabel(h, '$\Delta x \left( \textrm{pixels} \right)$', 'interpreter', 'latex', 'fontsize', fSize);
set(gca, 'position', c);



subplot(3, 2, 4);
imagesc(gx_full, gy_full, ...
    tx_rpc_full_mat);
axis image
xlim([1, image_width]);
ylim([1, image_height]);
set(gca, 'ydir', 'reverse');
title(sprintf('$\\Delta x \\, \\, \\left(\\textrm{RPC}, d_\\mathrm{RPC} = %0.1f \\right),\\, %d  \\,\\, \\textrm{pairs}$', rpc_diameter, num_images_correlate), 'FontSize', fSize, 'interpreter', 'latex')
set(gca, 'FontSize', fSize)
caxis([0, c_max]); 
c = get(gca, 'position');
h = colorbar;
ylabel(h, '$\Delta x \left( \textrm{pixels} \right)$', 'interpreter', 'latex', 'fontsize', fSize);
set(gca, 'position', c);



subplot(3, 2, 6 ,'FontSize', fSize);
imagesc(gx_full, gy_full, ...
    tx_apc_full_mat);
axis image
xlim([1, image_width]);
ylim([1, image_height]);
set(gca, 'ydir', 'reverse');
drawnow;
title(sprintf('$\\Delta x \\,\\, \\left(\\textrm{APC}\\right), \\, %d \\,\\, \\textrm{pairs}$', num_images_correlate), 'FontSize', fSize, 'interpreter', 'latex')
set(gca, 'FontSize', fSize)
caxis([0, c_max]); 
c = get(gca, 'position');
h = colorbar;
ylabel(h, '$\Delta x \left( \textrm{pixels} \right)$', 'interpreter', 'latex', 'fontsize', fSize);
set(gca, 'position', c);

set(gcf, 'color', 'white');

set(gcf, 'outerposition', [ -1696         214        1293        1054]);


