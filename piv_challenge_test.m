img_repo = get_image_repo();

% Image dir
im_dir = fullfile(img_repo, 'Ball', 'raw');

# Image base name
im_base = 'B';

% Number of digits
n_digits = 3;


% Extension
img_ext = '.bmp';

img_trailer_01 = 'a';
img_trailer_02 = 'b';

% Image numbers
start_pair = 1;
end_pair = 1;

img_nums = start_pair : end_pair;
num_images = length(img_nums);

% Name of the first image
img_name_01 = sprintf('im_base %03d img_trailer_01 img_ext', img_nums(1));

% Image path
image_path_01 = fullfile(im_dir, img_name_01);

% Load the image
img_01 = double(imread(image_path_01));

% Get the image size
[img_height, img_width] = size(img_01);

% Grid the image
[gx, gy] = gridImage([img_height, img_width], 











