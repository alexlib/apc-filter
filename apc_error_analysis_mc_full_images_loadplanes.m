
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
    load_filters = JobFile.JobOptions.LoadFilters;
    
    calculate_error = JobFile.JobOptions.CalculateError;
    calculate_displacements = JobFile.JobOptions.CalculateDisplacements;
    
    % Correlation planes file directory
    planes_dir = fullfile(image_dir, '..', 'planes');
    
    % Vector save directory
    vector_dir = fullfile(image_dir, '..', 'vect');
    
    % Results save directory
    error_dir = fullfile(image_dir, '..', 'error');

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
    
    % Size of the images
    image_height = solution_file.JobFile.Parameters.Image.Height;
    image_width = solution_file.JobFile.Parameters.Image.Width;
    image_size = [image_height, image_width];

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
    
    % Make the directory for the results if it doesn't exist
    if ~exist(error_dir, 'dir')
        mkdir(error_dir);
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
    
    % Calculate displacements?
    if calculate_displacements
        % Do the displacement measurements
        apc_error_analysis_ensemble_from_planes(...
        plane_save_paths, ...
        vector_save_paths, rpc_diameter, load_filters);
    
    end
    
    % Results file path
    error_save_name = sprintf(...
        'poiseuille_error_h%d_w%d_diff_std_%0.2f_%06d_%06d.mat', ...
            region_height, region_width, ...
            diffusion_std_dev, ...
            start_image, end_image);
        
        % results file path
        results_save_path = fullfile(error_dir, error_save_name);
    
        % Calculate errors?
    if calculate_error
        apc_calculate_error_from_vectors_poiseuille(...
            vector_save_paths, ...
            image_size, ...
            tx_max_true, ...
            diffusion_std_dev, ...
            results_save_path);
    end
    
    results_save_paths{n} = results_save_path;
    
end

% Make a plot
plot_apc_error_analysis_results_with_spectral(results_save_paths);

end


