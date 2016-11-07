

animal_number = 3;
case_name = 'mng_2-072-B';

% Image extension
image_ext = '.tiff';

% White field name
white_field_name = 'mng_2-072-E_white_field_ave.tiff';

% First and last images
start_image = 0;
end_image = 1000;
skip_image = 1;

% Get the image repository
image_repo = get_image_repo();

% Image input directory
input_dir_raw = fullfile(image_repo, ...
    'grasshopper', ...
    sprintf('grasshopper_%d', animal_number),...
    case_name, 'raw');

% Output input directory
output_dir_div = fullfile(input_dir_raw, '..', 'div', 'output');

% Make it
if ~exist(output_dir_div, 'dir');
    mkdir(output_dir_div);
end

% Name of the average white field image
white_field_path = fullfile(output_dir_div, '..', white_field_name);

% Load the white field image;
white_field = double(imread(white_field_path));

% Image numbers of the images to be divided
image_nums = start_image : skip_image : end_image;

% Number of images
num_images = length(image_nums);

% Create the image names
for k = 1 : num_images
   input_image_name  = sprintf('%s%06d%s', case_name, image_nums(k), image_ext); 
   output_image_name = sprintf('%s_div_%06d%s', case_name, image_nums(k), image_ext);
   
   % Image paths
   input_image_path{k} = fullfile(input_dir_raw, input_image_name);
   output_image_path{k} = fullfile(output_dir_div, output_image_name);
   
   
end

% Process the images
for k = 1 : num_images
   
    % Load the input image
    input_image = double(imread(input_image_path{k}));
    
    % Divide it by the white field
    output_image = input_image ./ white_field;
    
    % Fix the infinities
    output_image(isinf(output_image)) = 0;
    
    % Display the image
    imagesc(output_image);
    colormap gray
    caxis([0, 2]);
    axis image;
    pause(1/15);
    
end




















