% test_quiver_plot
% load('/Users/matthewgiarra/Desktop/piv_test_images/poiseuille_diffusion_3.00/vect/poiseuille_vect_h128_w128_diff_std_3.00_000099_000100.mat');

num_pairs = 50;
diffusion_std_dev = 5;

img_01 = 2 * num_pairs  - 1;
img_02 = img_01 + 1;

data_repo = '/Users/matthewgiarra/Desktop/piv_test_images';
vect_dir = fullfile(data_repo, sprintf('poiseuille_diffusion_%0.2f', diffusion_std_dev), 'vect');

vect_name = sprintf('poiseuille_vect_h128_w128_diff_std_%0.2f_%06d_%06d.mat', diffusion_std_dev, img_01, img_02);
vect_path = fullfile(vect_dir, vect_name);

load(vect_path);


image_size = 2048 * [1, 1];
u_max = 10;

tx_apc = tx_pair_apc;
ty_apc = ty_pair_apc;

tx_rpc = tx_pair_rpc;
ty_rpc = ty_pair_rpc;

tx_scc = tx_pair_scc;
ty_scc = ty_pair_scc;

apc_quiver_plots_02(tx_apc, ty_apc, tx_rpc, ty_rpc, tx_scc, ty_scc, gx, gy, image_size, u_max, num_pairs, diffusion_std_dev)
