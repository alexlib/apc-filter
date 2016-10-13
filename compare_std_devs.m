% compare_std_devs

pair_01 = 1;
pair_02 = 10;

ca = [20, 35];

vect_dir = '/Users/matthewgiarra/Desktop/piv_test_images/Ball/vect/';

name_01 = sprintf('B_vect_h128_w128_%06d_%06d.mat', pair_01, pair_01);
name_02 = sprintf('B_vect_h128_w128_%06d_%06d.mat', pair_02, pair_02);

path_01 = fullfile(vect_dir, name_01);
path_02 = fullfile(vect_dir, name_02);

% Load them
load(path_01);
nx = length(unique(gx));
ny = length(unique(gy));
sx_01 = reshape(apc_std_x_pair, [ny, nx]);
sy_01 = reshape(apc_std_y_pair, [ny, nx]);
apc_std_dev_01 = sqrt(sx_01.^2 + sy_01.^2);

load(path_02);
sx_02 = reshape(apc_std_x_pair, [ny, nx]);
sy_02 = reshape(apc_std_y_pair, [ny, nx]);
apc_std_dev_02 = sqrt(sx_02.^2 + sy_02.^2);

subplot(2, 1, 2);
imagesc(gx, gy, apc_std_dev_02);
% ca = caxis;
caxis(ca);
axis image;

subplot(2, 1, 1);
imagesc(gx, gy, apc_std_dev_01);
caxis(ca);
axis image;

