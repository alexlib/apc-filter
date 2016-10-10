
addpaths('..');

% Path to the solution file
solution_file_path = JobFile.Solution.Path;

image_dir = JobFile.Images.Directory;
image_base_name = JobFile.Images.BaseName;
img_ext = JobFile.Images.Extension;
num_digits = JobFile.Images.NumDigits;
start_image = JobFile.Images.Start;
end_image = JobFile.Images.End;
skip_image = JobFile.Images.Skip;
img_trailer_01 = JobFile.Images.Trailer_01;
img_trailer_02 = JobFile.Images.Trailer_02;


c_step = JobFile.Parameters.Processing.CorrelationStep;
region_width = JobFile.Parameters.Processing.RegionWidth;
region_height = JobFile.Parameters.Processing.RegionHeight;
window_fraction = JobFile.Parameters.Processing.WindowFraction;
shuffle_range_x = JobFile.Parameters.Processing.Shuffle.Range.X;
shuffle_range_y = JobFile.Parameters.Processing.Shuffle.Range.Y;
shuffle_step_x = JobFile.Parameters.Processing.Shuffle.Step.X;
shuffle_step_y = JobFile.Parameters.Processing.Shuffle.Step.Y;

grid_spacing_x = JobFile.Parameters.Processing.Grid.Spacing.X;
grid_spacing_y = JobFile.Parameters.Processing.Grid.Spacing.Y;

grid_buffer_x = JobFile.Parameters.Processing.Grid.Buffer.X;
grid_buffer_y = JobFile.Parameters.Processing.Grid.Buffer.Y;

% Load the Exact solution file
solution_file = load(solution_file_path);

% True average particle diameter in pixels.
rpc_diameter = solution_file.JobFile.Parameters.Particles.Diameter.Mean;

% True mean horizontal displacement
% of the Poiseuille displacement profile
tx_mean_true = solution_file.JobFile.Parameters.Displacement.Mean.X;

% Max velocity
tx_max_true = 3 / 2 * tx_mean_true;

% Diffusion velocity std dev
diffusion_std_dev = solution_file.JobFile.Parameters.Displacement.Rand.X;

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

% Remove any dots in the extension
img_ext = strrep(img_ext, '.', '');

% Form the image path lists
for k = 1 : num_pairs
    
   % Image names
   image_name_01 = [image_base_name num2str(image_nums_01(k), dig_str) img_trailer_01 '.' img_ext];
   image_name_02 = [image_base_name num2str(image_nums_02(k), dig_str) img_trailer_02 '.' img_ext];
    
   % List of image paths
   image_list_01{k} = fullfile(image_dir, image_name_01);
   image_list_02{k} = fullfile(image_dir, image_name_02);
    
end

% Load the first image and get its size
[image_height, image_width] = size(double(imread(image_list_01{1})));

% Save image size to vector
image_size = [image_height, image_width];

% Grid the images
grid_spacing = [grid_spacing_y, grid_spacing_x];
% grid_buffer_y = region_height/2 * [1, 1];
% grid_buffer_x = region_width/2 * [1, 1];

% Grid the image
[grid_x, grid_y] = gridImage([image_height, image_width],...
    grid_spacing, grid_buffer_y, grid_buffer_x);

% Number of grid points in each direction
nx = length(unique(grid_x(:)));
ny = length(unique(grid_y(:)));

% Do the correlation error analysis
[ty_apc, tx_apc, ...
    ty_rpc, tx_rpc, ...
    ty_scc, tx_scc, ...
    apc_std_y, apc_std_x] = apc_error_analysis_ensemble(...
    image_list_01, image_list_02, grid_y, grid_x, region_size, ...
    window_fraction, rpc_diameter);

gx_mat = reshape(grid_x, [ny, nx]);
gy_mat = reshape(grid_y, [ny, nx]);

tx_apc_mat = -1 * reshape(tx_apc, [ny, nx, num_pairs]);
ty_apc_mat = -1 * reshape(ty_apc, [ny, nx, num_pairs]);

tx_rpc_mat = -1 * reshape(tx_rpc, [ny, nx, num_pairs]);
ty_rpc_mat = -1 * reshape(ty_rpc, [ny, nx, num_pairs]);

tx_scc_mat = -1 * reshape(tx_scc, [ny, nx, num_pairs]);
ty_scc_mat = -1 * reshape(ty_scc, [ny, nx, num_pairs]);

apc_std_x_mat =  reshape(apc_std_x, [ny, nx, num_pairs]);
apc_std_y_mat =  reshape(apc_std_y, [ny, nx, num_pairs]);

apc_std_x_vect = reshape(apc_std_x_mat(:, :, end), [ny * nx, 1]);

figure(2);
apc_quiver_plots(...
    tx_apc_mat(:, :, end), ty_apc_mat(:, :, end), ...
    tx_rpc_mat(:, :, end), ty_rpc_mat(:, :, end), ...
    tx_scc_mat(:, :, end), ty_scc_mat(:, :, end), ...
    grid_x, grid_y, ...
    image_size, ...
    tx_max_true);

% Results file name
results_file_name = sprintf('poiseuille_piv_h%d_w%d_diff_std_%0.2f.mat', ...
    image_height, image_width,  diffusion_std_dev);

% Results file directory
results_dir = fullfile(image_dir, '..', 'results');

% Make the results dir
if ~exist(results_dir, 'dir')
    mkdir(results_dir)
end

% Results path
results_path = fullfile(results_dir, results_file_name);

% Save results
save(results_path, ...
    'solution_file', 'JobFile', ...
    'tx_apc_mat', 'ty_apc_mat', ...
    'tx_rpc_mat', 'ty_rpc_mat', ...
    'tx_scc_mat', 'ty_scc_mat', ...
    'apc_std_x_mat', 'apc_std_y_mat');

% Save results

% for p = 1 : num_pairs
%     
%     sx = apc_std_x(p);
%     sy = apc_std_y(p);
% 
%     apc_filt = exp(-X.^2 / (2 * sx^2))  .* exp(-Y.^2 / (2 * sy^2));
%     
%     surf(apc_filt, 'linewidth', 1E-5);
%     axis square
%     xlim([1, region_width]);
%     ylim([1, region_height]);
%     zlim([0, 1]);
%     
%     pause(0.01);
% 
% end




