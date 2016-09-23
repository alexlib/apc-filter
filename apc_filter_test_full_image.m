addpaths;
% clear;


fSize_axes = 12;
fSize_title = 16;


image_dir = '~/Desktop/piv_test_images/tiffs/0.5_ul_min';
image_base_name = 'img_';
num_digits = 2;
image_ext = '.tiff';
start_image = 1;
end_image = 1;
skip_image = 1;
c_step = 1;

JobList = MonteCarloImageGenerationJobFile_micro;
[imageMatrix1, imageMatrix2] = generateMonteCarloImageSet_micro(JobList);


% image_dir = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/fmc/analysis/data/synthetic/vortex/lambvortex_2014-05-22_centralDifference/c_0.0250/lambvortex_h1024_w1024_00001/raw';
% image_base_name = 'lambvortex_h1024_w1024_';
% num_digits = 6;
% image_ext = '.tiff';
% start_image = 1;
% end_image = 1;
% skip_image = 1;
% c_step = 1;

% image_dir = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/fmc/analysis/data/experimental/vortex/vortexring_2013-11-12_d03_f60_t06/raw';
% image_base_name = 'vortexring_d03_f60_t06_';
% num_digits = 6;
% image_ext = '.tif';
% start_image = 840;
% end_image = 840;
% skip_image = 1;
% c_step = 1;



% Image dimensions
region_height = 64;
region_width = 64;

% Grid spacing
grid_x = 64;
grid_y = 64;

% Grid buffer
grid_buffer_x = region_width * [1, 1];
grid_buffer_y = region_height * [1, 1];

% Window size
window_fraction = 0.4 * [1, 1];

% % Grid step
shuffle_range = 32;
shuffle_step = 16;

gx_range = shuffle_range;
gx_step = shuffle_step;

gy_range = shuffle_range;
gy_step = shuffle_step;

% Max number of reqions for the neq correlation
num_regions_neq = 1000;

% % % % %
% Center pixels
xc = (region_width  + 1) / 2 + 0.5 * (1 - mod(region_width,  2));
yc = (region_height + 1) / 2 + 0.5 * (1 - mod(region_height, 2));

% Window
g_win = gaussianWindowFilter([region_height, region_width], window_fraction, 'fraction');

% Allocate the particle shape
particle_shape_sum = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);

ncc_sum = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);
% 
% % List of image numbers
% image_list_01 = start_image : skip_image : end_image;
% image_list_02 = image_list_01 + c_step;
% 
% % Number of images
% num_images = length(image_list_01);
% 
% % String Formatting 
% num_str = ['%0' num2str(num_digits) 'd'];
% 
% % Image name string
% image_name_01 = [image_base_name num2str(image_list_01(1), num_str) image_ext];
% image_path_01 = fullfile(image_dir, image_name_01);
% [image_height, image_width] = size(double(imread(image_path_01)));

[image_height, image_width, num_images] = size(imageMatrix1);

% Create the image grid.
[gx, gy] = gridImage([image_height, image_width], [grid_y, grid_x], grid_buffer_y, grid_buffer_x);

% Number of regions
num_regions = numel(gx);

% Allocate velocity matrices
v = zeros(num_regions, num_images);
u = zeros(num_regions, num_images);

% Basic RPC filter
rpc_filter = spectralEnergyFilter(region_height, region_width, sqrt(8));

for k = 1 : num_images;
    
    fprintf('On image %d of %d\n', k, num_images);
    
    % Image 1 name and path
%     image_name_01 = [image_base_name num2str(image_list_01(k), num_str) image_ext];
%     image_path_01 = fullfile(image_dir, image_name_01);
%     
%     % Image 2 name and path
%     image_name_02 = [image_base_name num2str(image_list_02(k), num_str) image_ext];
%     image_path_02 = fullfile(image_dir, image_name_02);
%     
%     % Read the images
%     image_01 = double(imread(image_path_01));
%     image_02 = double(imread(image_path_02));
%     
    image_01 = double(imageMatrix1(:, :, k));
    image_02 = double(imageMatrix2(:, :, k));
    
    % Create the spectral filter.
    
    [SPECTRAL_FILTER, PARTICLE_SHAPE, FILTER_FIT] =...
        calculate_apc_filter_from_image_pair(image_01, image_02, ...
    gy, gx, region_height, region_width, window_fraction, ...
    shuffle_range, shuffle_step, num_regions_neq);

    % The image is already gridded. Just do the correlations now.
    %
    % Extract regions for image 1
    region_matrix_01 = extractSubRegions(image_01, ...
        [region_height, region_width], gx, gy);
    % Extract regions for image 2
    region_matrix_02 = extractSubRegions(image_02, ...
        [region_height, region_width], gx, gy);
    
    % Do all the mean subtraction
    fprintf(1, 'Doing mean subtraction\n');
    for p = 1 : num_regions
        region_matrix_01(:, :, p) = g_win .* (region_matrix_01(:, :, p) - mean(mean(region_matrix_01(:, :, p))));
        region_matrix_02(:, :, p) = g_win .* (region_matrix_02(:, :, p) - mean(mean(region_matrix_02(:, :, p))));
    end
    
    % Allocate RPC planes
    apc_plane_list = zeros(region_height, region_width, num_regions);
    scc_ensemble_plane = zeros(region_height, region_width);
    
    % Loop over all the regions
    for p = 1 : num_regions
        
        % Velocity measurements
        fprintf('On PIV region %d of %d\n', p, num_regions);
        
        % Extract regions
        region_01 = region_matrix_01(:, :, p);
        region_02 = region_matrix_02(:, :, p);
        
        % Do the displacement estimate
        [v(p, k), u(p, k), apc_plane] = RPC(region_01, region_02, FILTER_FIT);
        
        [~, ~, scc_plane] = SCC(region_01, region_02, 1);
        
        scc_ensemble_plane = scc_ensemble_plane + scc_plane;
        
        % Save the RPC Plane to the matrix
        apc_plane_list(:, :, p) = apc_plane;
        
        
    end
    
    % Now do the ensemble with APC to compare.
    fprintf(1, 'On ensemble correlation...\n');
    [v_ensemble_apc, u_ensemble_apc, apc_ensemble_plane] = ...
    rpc_ensemble(region_matrix_01, region_matrix_02, ...
    FILTER_FIT);
    
 %   Now do the ensemble with RPC to compare.
    fprintf(1, 'On ensemble correlation...\n');
    [v_ensemble_rpc, u_ensemble_rpc, rpc_ensemble_plane] = ...
    rpc_ensemble(region_matrix_01, region_matrix_02, ...
    rpc_filter);

end

fSize = 16;

subplot(1, 4, 1);
surf(rpc_ensemble_plane ./ max(rpc_ensemble_plane(:)));
title('RPC Ensemble', 'fontsize', fSize);
xlim([1, region_width]);
ylim([1, region_height]);
zlim([0, 1]);
set(gca, 'fontsize', fSize);
axis square;

subplot(1, 4, 2);
surf(imrotate(scc_ensemble_plane, 180) ./ max(scc_ensemble_plane(:)));
title('SCC Ensemble', 'fontsize', fSize);
xlim([1, region_width]);
ylim([1, region_height]);
zlim([0, 1]);
set(gca, 'fontsize', fSize);
axis square;

subplot(1, 4, 3);
surf(apc_ensemble_plane ./ max(apc_ensemble_plane(:)));
title('APC Ensemble', 'fontsize', fSize);
xlim([1, region_width]);
ylim([1, region_height]);
zlim([0, 1]);
set(gca, 'fontsize', fSize);
axis square;

subplot(1, 4, 4);
surf(FILTER_FIT);
title('APC Filter', 'fontsize', fSize);
xlim([1, region_width]);
ylim([1, region_height]);
zlim([0, 1]);
set(gca, 'fontsize', fSize);
axis square;










