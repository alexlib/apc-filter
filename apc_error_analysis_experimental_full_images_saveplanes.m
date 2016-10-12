
function apc_error_analysis_experimental_full_images_saveplanes(JOBLIST)

% Count the number of jobs
num_jobs = length(JOBLIST);

% Add some paths
addpaths('..');

% Loop over the jobs
for n = 1 : num_jobs
    JobFile = JOBLIST(n);
   
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

    grid_spacing_x = JobFile.Parameters.Processing.Grid.Spacing.X;
    grid_spacing_y = JobFile.Parameters.Processing.Grid.Spacing.Y;

    grid_buffer_x = JobFile.Parameters.Processing.Grid.Buffer.X;
    grid_buffer_y = JobFile.Parameters.Processing.Grid.Buffer.Y;

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
    
    % Declare the results list
    plane_save_paths = {''};

    % Digit string
    dig_str = ['%0' num2str(num_digits) 'd'];

    % Remove any dots in the extension
    img_ext = strrep(img_ext, '.', '');
    
    % Correlation planes output file directory
    planes_dir = fullfile(image_dir, '..', 'planes');

    % Make the results dir
    if ~exist(planes_dir, 'dir')
        mkdir(planes_dir)
    end

    % Form the image path lists
    for k = 1 : num_pairs

        % Image names
        image_name_01 = [image_base_name num2str(image_nums_01(k), dig_str) img_trailer_01 '.' img_ext];
        image_name_02 = [image_base_name num2str(image_nums_02(k), dig_str) img_trailer_02 '.' img_ext];

        % List of image paths
        image_list_01{k} = fullfile(image_dir, image_name_01);
        image_list_02{k} = fullfile(image_dir, image_name_02);

        % Paths to the save planes
        planes_file_name = sprintf('%s_h%d_w%d_%06d_%06d.mat', ...
        image_base_name, region_height, region_width, image_nums_01(k), image_nums_02(k));

        % Results path
        plane_save_paths{k} = fullfile(planes_dir, planes_file_name);

    end
   
    % Load the first image to measure its size for gridding
    [image_height, image_width, ~] = size(imread(image_list_01{1}));
    
    % Grid the images
    grid_spacing = [grid_spacing_y, grid_spacing_x];

    % Grid the image
    [grid_x, grid_y] = gridImage([image_height, image_width],...
        grid_spacing, grid_buffer_y, grid_buffer_x);
  
    % Do the correlations and save the planes.
    for k = 1 : num_pairs      
        t1 = tic;

        fprintf(1, '\nOn pair %d of %d...\n', k, num_pairs);

        plane_save_path = plane_save_paths{k};
        image_path_01 = image_list_01{k};
        image_path_02 = image_list_02{k};

        % Do the correlations and save the planes
        apc_save_planes(...
        image_path_01, image_path_02, grid_y, grid_x, region_size, ...
        window_fraction, plane_save_path);     

        % End pair timer
        t2 = toc(t1);

        fprintf(1, 'Total image pair time: %0.2f sec.\n', t2);
    end
   
    
    

end

end


