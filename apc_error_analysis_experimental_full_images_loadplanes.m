
function apc_error_analysis_experimental_full_images_loadplanes(JOBLIST)

% Count the number of jobs
num_jobs = length(JOBLIST);

% Add some paths
addpaths('..');

% Loop over the jobs
for n = 1 : num_jobs
    JobFile = JOBLIST(n);
    
    calculate_displacements = JobFile.JobOptions.CalculateDisplacements;
    calculate_error = JobFile.JobOptions.CalculateError;
    
    
    image_dir = JobFile.Images.Directory;
    image_base_name = JobFile.Images.BaseName;
  
  
    start_image = JobFile.Images.Start;
    end_image = JobFile.Images.End;
    skip_image = JobFile.Images.Skip;
    
    % Correlation planes file directory
    planes_dir = fullfile(image_dir, '..', 'planes');
    
    % Vector save directory
    vector_dir = fullfile(image_dir, '..', 'vect');

    c_step = JobFile.Parameters.Processing.CorrelationStep;
    region_width = JobFile.Parameters.Processing.RegionWidth;
    region_height = JobFile.Parameters.Processing.RegionHeight;
    
    % RPC Diameter
    rpc_diameter = JobFile.Parameters.Processing.RpcDiameter;
    
    % Load filter fit?
    load_filter_fit = JobFile.JobOptions.LoadFilters;

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
        
        % Paths to the save planes
        planes_file_name = sprintf('%s_h%d_w%d_%06d_%06d.mat', ...
            image_base_name, region_height, region_width, image_nums_01(k), image_nums_02(k));

        % Results path
        plane_save_paths{k} = fullfile(planes_dir, planes_file_name);
        
        % Vector file name
        vector_file_name = sprintf('%s_vect_h%d_w%d_%06d_%06d.mat', ...
            image_base_name, region_height, region_width, image_nums_01(k), image_nums_02(k));
        
        % Planes paths (these will be read / loaded)
        plane_save_paths{k} = fullfile(planes_dir, planes_file_name);
        
        % Vector save paths (these will be written / saved)
        vector_save_paths{k} = fullfile(vector_dir, vector_file_name);

    end
    
    if calculate_displacements
        % Do the displacement measurements
        apc_error_analysis_ensemble_from_planes(...
            plane_save_paths, ...
            vector_save_paths, rpc_diameter, load_filter_fit);
    end
  
    % Calculate errors?
    if calculate_error
          
    % Image size
    image_height = JobFile.Images.Height;
    image_width = JobFile.Images.Width;
    
    % Image size
    image_size = [image_height, image_width];
        
        apc_calculate_error_from_vectors_poiseuille(...
            vector_save_paths, image_size, dx_max);
    end
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


