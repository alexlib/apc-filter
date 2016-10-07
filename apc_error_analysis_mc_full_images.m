
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

% Load the solution file
solution_file = load(solution_file_path);

% True average particle diameter in pixels.
rpc_diameter = solution_file.JobFile.Parameters.Particles.Diameter.Mean;

% True mean horizontal displacement
% of the Poiseuille displacement profile
tx_mean_true = solution_file.JobFile.Parameters.Displacement.Mean.X;

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

% Grid the images
grid_spacing = [grid_spacing_y, grid_spacing_x];
grid_buffer_y = region_height/2 * [1, 1];
grid_buffer_x = region_width/2 * [1, 1];

% grid_buffer_x = image_width/2 * [1, 1];

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







