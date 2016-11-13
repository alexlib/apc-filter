function [convergence_iteration, tx_diff_mean_smoothed, ty_diff_mean_smoothed, pair_nums] = ...
    calculate_vector_convergence(tx, ty, convergence_criterion, kernel_length)

% This is the convolution kernel for the moving mean
kernel = 1 / kernel_length * ones(kernel_length, 1);

% This is the difference in error between passes
tx_diff = abs(diff(tx, 1, 2));
ty_diff = abs(diff(ty, 1, 2));

% Average change in estimates
tx_diff_mean = mean(tx_diff, 1);
ty_diff_mean = mean(ty_diff, 1);

% Smoothed change in estimates
tx_diff_mean_smoothed = conv(tx_diff_mean, kernel, 'valid');
ty_diff_mean_smoothed = conv(ty_diff_mean, kernel, 'valid');

% Number of pairs correlated
num_pairs = size(tx, 2);

% This is the anchor point of the convolution kernel
kernel_len_half = (kernel_length  + 1) / 2 + ...
    0.5 * (1 - mod(kernel_length,  2));

% This vector holds the image pair numbers.
pair_nums = kernel_len_half : (num_pairs - kernel_len_half);

% Find where convergence is reached
convergence_iteration = pair_nums(find(...
    tx_diff_mean_smoothed < convergence_criterion ...
    & ty_diff_mean_smoothed < convergence_criterion, 1));

end