old_results_path = '~/Desktop/debug/old_code_vars.mat';

new_results_path = '~/Desktop/debug/new_code_vars.mat';

new_mats = load(new_results_path);
old_mats = load(old_results_path);

cc_array_new = new_mats.spectral_correlation_array;

cc_array_old = old_mats.spectral_correlation_array;

std_x_old = old_mats.apc_std_x;
std_y_old = old_mats.apc_std_y;

std_x_new = new_mats.apc_std_x(:, end);
std_y_new = new_mats.apc_std_y(:, end);

num_regions = size(cc_array_new, 3);

% for k = 1 : num_regions
%    
%     
%     subplot(1,2,1)
%     imagesc(real(cc_array_old(:, :, k)));
%     axis image;
%     title('old');
%     
%     subplot(1,2,2)
%     imagesc(real(cc_array_new(:, :, k)));
%     axis image;
%     title('new');
%     
%     pause
%     
%     
% end