function apc_calculate_error_from_vectors_poiseuille(...
    vector_save_paths, ...
       image_size, ...
       dx_max,...
       diffusion_std,...
       results_save_path)

% This is the number of images that were correlated
num_pairs = length(vector_save_paths);

% Image height
image_height = image_size(1);

% Timer
t = tic;

% Load the first vector path to get the grid
load(vector_save_paths{1});

% Calculate the true solution at the grid points
% Center of the image in the height direction
yc = image_height / 2;
    
% Radial coordinate in the height direction
r = abs((gy - yc) / (image_height/2));

% Number of points
ny = length(unique(gy));
nx = length(unique(gx));
regions_per_pair = ny * nx;

% True mean displacement along the velocity profile
tx_true = - dx_max * (r.^2 - 1);
ty_true = zeros(size(tx_true));

% Vectors to hold the mean errors
mean_err_scc_spatial = zeros(num_pairs, 1);
mean_err_rpc_spatial = zeros(num_pairs, 1);
mean_err_scc_spectral = zeros(num_pairs, 1);
mean_err_rpc_spectral = zeros(num_pairs, 1);
mean_err_apc = zeros(num_pairs, 1);

% Vectors to hold the std devs of error
std_err_scc_spatial = zeros(num_pairs, 1);
std_err_rpc_spatial = zeros(num_pairs, 1);
std_err_scc_spectral = zeros(num_pairs, 1);
std_err_rpc_spectral = zeros(num_pairs, 1);
std_err_apc = zeros(num_pairs, 1);

% Allocate arrays to hold the filters
apc_sx_mean = zeros(regions_per_pair, 1);
apc_sy_mean = zeros(regions_per_pair, 1);
apc_sx_std = zeros(regions_per_pair, 1);
apc_sy_std = zeros(regions_per_pair, 1);

% Loop over the images
for p = 1 : num_pairs
        
    % Inform the user.
    fprintf(1, 'On image %d of %d\n', p, num_pairs);

    % Load the planes data
    load(vector_save_paths{p});
    
    % Velocity components (spatial ensembles)
    tx_scc_spatial = 1 * tx_pair_scc_spatial;
    ty_scc_spatial = 1 * ty_pair_scc_spatial;
    tx_rpc_spatial = 1 * tx_pair_rpc_spatial;
    ty_rpc_spatial = 1 * ty_pair_rpc_spatial;
    
    % Velocity components (spectral ensembles)
    tx_scc_spectral = 1 * tx_pair_scc_spectral;
    ty_scc_spectral = 1 * ty_pair_scc_spectral;
    tx_rpc_spectral = 1 * tx_pair_rpc_spectral;
    ty_rpc_spectral = 1 * ty_pair_rpc_spectral;
    
    % Velocity components (APC spectral ensembles)
    tx_apc = 1 * tx_pair_apc;
    ty_apc = 1 * ty_pair_apc;
    
    % APC filter standard deviations
    apc_sx = apc_std_x_pair;
    apc_sy = apc_std_y_pair;
    
    % Absolute value of the error
    tx_err_scc_spatial = (tx_true - tx_scc_spatial);
    ty_err_scc_spatial = (ty_true - ty_scc_spatial);
    tx_err_rpc_spatial = (tx_true - tx_rpc_spatial);
    ty_err_rpc_spatial = (ty_true - ty_rpc_spatial);
    
    % Spectral correlations
    tx_err_scc_spectral = (tx_true - tx_scc_spectral);
    ty_err_scc_spectral = (ty_true - ty_scc_spectral);
    tx_err_rpc_spectral = (tx_true - tx_rpc_spectral);
    ty_err_rpc_spectral = (ty_true - ty_rpc_spectral);
    
    % APC error
    tx_err_apc = (tx_true - tx_apc);
    ty_err_apc = (ty_true - ty_apc);
    
    % Error magnitudes
    err_mag_scc_spatial = sqrt(tx_err_scc_spatial.^2 + ty_err_scc_spatial.^2);
    err_mag_rpc_spatial = sqrt(tx_err_rpc_spatial.^2 + ty_err_rpc_spatial.^2);
    err_mag_scc_spectral = sqrt(tx_err_scc_spectral.^2 + ty_err_scc_spectral.^2);
    err_mag_rpc_spectral = sqrt(tx_err_rpc_spectral.^2 + ty_err_rpc_spectral.^2);
    err_mag_apc = sqrt(tx_err_apc.^2 + ty_err_apc.^2);
    
    % Average the errors
    mean_err_scc_spatial(p) = mean(err_mag_scc_spatial);
    mean_err_rpc_spatial(p) = mean(err_mag_rpc_spatial);
    mean_err_scc_spectral(p) = mean(err_mag_scc_spectral);
    mean_err_rpc_spectral(p) = mean(err_mag_rpc_spectral);
    mean_err_apc(p) = mean(err_mag_apc);
    
    % Standard deviations of error
    std_err_scc_spatial(p) = std(err_mag_scc_spatial);
    std_err_rpc_spatial(p) = std(err_mag_rpc_spatial);
    std_err_scc_spectral(p) = std(err_mag_scc_spectral);
    std_err_rpc_spectral(p) = std(err_mag_rpc_spectral);
    std_err_apc(p) = std(err_mag_apc);
    
    % Mean and std devs of filter.
    apc_sx_mean(p) = mean(apc_sx);
    apc_sy_mean(p) = mean(apc_sy);
    apc_sx_std(p) = std(apc_sx);
    apc_sy_std(p) = std(apc_sy);
    
end % End (for p = 1 : num_images)

% End timer
tf = toc(t);

% Save the results
save(results_save_path, ...
    'mean_err_scc_spatial', 'std_err_scc_spatial', ...
    'mean_err_rpc_spatial', 'std_err_rpc_spatial', ...
    'mean_err_scc_spectral', 'std_err_scc_spectral', ...
    'mean_err_rpc_spectral', 'std_err_rpc_spectral', ...
    'mean_err_apc', 'std_err_apc', ...
    'apc_sx_mean', 'apc_sx_std', ...
    'apc_sy_mean', 'apc_sy_std', ...
    'diffusion_std');

% Print elapsed time
fprintf(1, 'Elapsed time: %d seconds for %d image pairs.\n', tf, num_pairs);
fprintf(1, '%0.1f seconds per image.\n\n', tf / num_pairs);




end































