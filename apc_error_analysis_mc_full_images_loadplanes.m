
function apc_error_analysis_mc_full_images_loadplanes(JOBLIST)

% Count the number of jobs
num_jobs = length(JOBLIST);

% Add some paths
addpaths('..');

% Loop over the jobs
for n = 1 : num_jobs
    JobFile = JOBLIST(n);
    
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
    
    % Correlation planes file directory
    planes_dir = fullfile(image_dir, '..', 'planes');
    
    % Vector save directory
    vector_dir = fullfile(image_dir, '..', 'vect');

    c_step = JobFile.Parameters.Processing.CorrelationStep;
    region_width = JobFile.Parameters.Processing.RegionWidth;
    region_height = JobFile.Parameters.Processing.RegionHeight;
    window_fraction = JobFile.Parameters.Processing.WindowFraction;

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

    % List of image numbers
    image_nums_01 = start_image : skip_image : end_image;
    image_nums_02 = image_nums_01 + c_step;

    % Number of images
    num_pairs = length(image_nums_01);
    
    % Declare the results list
    plane_save_paths = {''};
    vector_save_paths = {''};
    
    % Make the vector dir if it doesn't exist
    if ~exist(vector_dir, 'dir')
        mkdir(vector_dir);
    end

    % Form the image path lists
    for k = 1 : num_pairs
        
        % Planes file name
        planes_file_name = sprintf('poiseuille_piv_h%d_w%d_diff_std_%0.2f_%06d_%06d.mat', ...
            region_height, region_width, diffusion_std_dev, image_nums_01(k), image_nums_02(k));
        
        % Vector file name
        vector_file_name = sprintf('poiseuille_vect_h%d_w%d_diff_std_%0.2f_%06d_%06d.mat', ...
            region_height, region_width, diffusion_std_dev, image_nums_01(k), image_nums_02(k));

        % Planes paths (these will be read / loaded)
        plane_save_paths{k} = fullfile(planes_dir, planes_file_name);
        
        % Vector save paths (these will be written / saved)
        vector_save_paths{k} = fullfile(vector_dir, vector_file_name);

    end
    
    % Do the displacement measurements
    apc_error_analysis_ensemble_from_planes(...
        plane_save_paths, rpc_diameter, ...
        vector_save_paths);
  
% 
%     % Load the first image and get its size
%     [image_height, image_width] = size(double(imread(image_list_01{1})));
% 
%     % Save image size to vector
%     image_size = [image_height, image_width];
% 
% 
% 
%     % Number of grid points in each direction
%     nx = length(unique(grid_x(:)));
%     ny = length(unique(grid_y(:)));
% 
%     % Do the correlation error analysis
%     [ty_apc, tx_apc, ...
%         ty_rpc, tx_rpc, ...
%         ty_scc, tx_scc, ...
%         apc_std_y, apc_std_x] = apc_error_analysis_ensemble(...
%         image_list_01, image_list_02, grid_y, grid_x, region_size, ...
%         window_fraction, rpc_diameter);
% 
%     gx_mat = reshape(grid_x, [ny, nx]);
%     gy_mat = reshape(grid_y, [ny, nx]);
% 
%     tx_apc_mat = -1 * reshape(tx_apc, [ny, nx, num_pairs]);
%     ty_apc_mat = -1 * reshape(ty_apc, [ny, nx, num_pairs]);
% 
%     tx_rpc_mat = -1 * reshape(tx_rpc, [ny, nx, num_pairs]);
%     ty_rpc_mat = -1 * reshape(ty_rpc, [ny, nx, num_pairs]);
% 
%     tx_scc_mat = -1 * reshape(tx_scc, [ny, nx, num_pairs]);
%     ty_scc_mat = -1 * reshape(ty_scc, [ny, nx, num_pairs]);
% 
%     apc_std_x_mat =  reshape(apc_std_x, [ny, nx, num_pairs]);
%     apc_std_y_mat =  reshape(apc_std_y, [ny, nx, num_pairs]);
% 
%     apc_std_x_vect = reshape(apc_std_x_mat(:, :, end), [ny * nx, 1]);
% 
% 
%     % Save results
%     save(results_path, ...
%         'solution_file', 'JobFile', ...
%         'tx_apc_mat', 'ty_apc_mat', ...
%         'tx_rpc_mat', 'ty_rpc_mat', ...
%         'tx_scc_mat', 'ty_scc_mat', ...
%         'apc_std_x_mat', 'apc_std_y_mat', ...
%         'gx_mat', 'gy_mat');
% 
%     fprintf(1, 'Saved results to %s\n', results_path);

end

end


