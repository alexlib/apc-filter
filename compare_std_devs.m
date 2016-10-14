% compare_std_devs

pair_01_list = [1 : 100];
pair_02 = 100;

lw = 1;

rpc_std = 8 * pi^2 / 3;

Skip_x = 6;
Skip_y = 1;
Scale = 6;

ca = [0.32, 0.6];

num_pairs = length(pair_01_list);

% ca = [20, 35];

vect_dir = '/Users/matthewgiarra/Desktop/piv_test_images/Ball/vect/';

name_02 = sprintf('B_vect_h128_w128_%06d_%06d.mat', pair_02, pair_02);
path_02 = fullfile(vect_dir, name_02);
load(path_02);
nx = length(unique(gx));
ny = length(unique(gy));
sx_02 = reshape(apc_std_x_pair ./ rpc_std, [ny, nx]);
sy_02 = reshape(apc_std_y_pair ./ rpc_std, [ny, nx]);
gx_mat = reshape(gx, [ny, nx]);
gy_mat = reshape(gy, [ny, nx]);
apc_std_dev_02 = sqrt(sx_02.^2 + sy_02.^2);

subplot(2, 1, 2);
imagesc(gx, gy, sx_02);
caxis(ca);
% ca = caxis
axis image;



std_frac_mean = mean(sx_02, 2);

gyu = unique(gy);

yl = [5, 13];

for k = 1 : num_pairs
    
name_01 = sprintf('B_vect_h128_w128_%06d_%06d.mat', pair_01_list(k), pair_01_list(k));
path_01 = fullfile(vect_dir, name_01);

% Load them
load(path_01);

sx_01 = reshape(apc_std_x_pair ./ rpc_std, [ny, nx]);
sy_01 = reshape(apc_std_y_pair ./ rpc_std, [ny, nx]);
apc_std_dev_01 = sqrt(sx_01.^2 + sy_01.^2);
apc_std_dev_01 = sx_01;

tx_apc_mat = reshape(tx_pair_apc, [ny, nx]);
ty_apc_mat = reshape(ty_pair_apc, [ny, nx]);

tx_rpc_mat = reshape(tx_pair_rpc, [ny, nx]);
tx_scc_mat = reshape(tx_pair_scc, [ny, nx]);

tx_scc_mean = mean(tx_scc_mat, 2);
tx_scc_std  = 2 * std(tx_scc_mat, [], 2);

tx_rpc_mean = mean(tx_rpc_mat, 2);
tx_rpc_std = 2 * std(tx_rpc_mat, [], 2);

tx_apc_mean = mean(tx_apc_mat, 2);
tx_apc_std  = 2 * std(tx_apc_mat, [], 2);

subplot(2, 3, 4:6);
imagesc(gx, gy, apc_std_dev_01);
caxis(ca);
axis image;
hold on;
quiver(...
    gx_mat(1 : Skip_y : end, 1 : Skip_x : end), ...
    gy_mat(1 : Skip_y : end, 1 : Skip_x : end), ...
    Scale * tx_apc_mat(1 : Skip_y : end, 1 : Skip_x : end),...
    Scale * ty_apc_mat(1 : Skip_y : end, 1 : Skip_x : end),...
    0, 'black', 'linewidth', lw);
hold off

subplot(2, 3, 1);
shadedErrorBar(gyu, tx_scc_mean, tx_scc_std, 'blue');
set(gca, 'view', [90, 90]);
axis square;
title('SCC');
grid on;
box on;
ylim(yl)

subplot(2, 3, 2);
shadedErrorBar(gyu, tx_rpc_mean, tx_rpc_std, 'red');
set(gca, 'view', [90, 90]);
axis square;
title('SCC');
grid on;
box on;
ylim(yl)

subplot(2, 3, 3);
shadedErrorBar(gyu, tx_apc_mean, tx_apc_std, 'black');
set(gca, 'view', [90, 90]);
axis square;
title('SCC');
grid on;
box on;
ylim(yl)

set(gcf, 'color', 'white')

pause(0.01);
    
end







