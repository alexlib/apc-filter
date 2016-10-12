% test_quiver_plot
load('/Users/matthewgiarra/Desktop/piv_test_images/poiseuille_diffusion_3.00/vect/poiseuille_vect_h128_w128_diff_std_3.00_000099_000100.mat');

image_size = 2048 * [1, 1];
u_max = 10;

tx_apc = tx_pair_apc;
ty_apc = ty_pair_apc;

tx_rpc = tx_pair_rpc;
ty_rpc = ty_pair_rpc;

tx_scc = tx_pair_scc;
ty_scc = ty_pair_scc;


apc_quiver_plots_02(tx_apc, ty_apc, tx_rpc, ty_rpc, tx_scc, ty_scc, gx, gy, image_size, u_max)
