function apc_calculate_error_from_vectors_poiseuille(vector_save_paths, ...
       image_size, dx_max)

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
mean_err_scc = zeros(num_pairs, 1);
mean_err_rpc = zeros(num_pairs, 1);
mean_err_apc = zeros(num_pairs, 1);

% Allocate arrays to hold the filters
apc_sx = zeros(regions_per_pair, num_pairs);
apc_sy = zeros(regions_per_pair, num_pairs);

% Loop over the images
for p = 1 : num_pairs
        
    % Inform the user.
    fprintf(1, 'On image %d of %d\n', p, num_pairs);

    % Load the planes data
    load(vector_save_paths{p});
    
    % Velocity components
    tx_scc = 1 * tx_pair_scc;
    ty_scc = 1 * ty_pair_scc;
    tx_rpc = 1 * tx_pair_rpc;
    ty_rpc = 1 * ty_pair_rpc;
    tx_apc = 1 * tx_pair_apc;
    ty_apc = 1 * ty_pair_apc;
    
    apc_sx(:, p) = apc_std_x_pair;
    apc_sy(:, p) = apc_std_y_pair;
    
    % Absolute value of the error
    tx_err_scc = (tx_true - tx_scc);
    ty_err_scc = (ty_true - ty_scc);

    tx_err_rpc = (tx_true - tx_rpc);
    ty_err_rpc = (ty_true - ty_rpc);

    tx_err_apc = (tx_true - tx_apc);
    ty_err_apc = (ty_true - ty_apc);
    
    % Error magnitudes
    err_mag_scc = sqrt(tx_err_scc.^2 + ty_err_scc.^2);
    err_mag_rpc = sqrt(tx_err_rpc.^2 + ty_err_rpc.^2);
    err_mag_apc  = sqrt(tx_err_apc.^2 + ty_err_apc.^2);
    
    % Average the errors
    mean_err_scc(p) = mean(err_mag_scc);
    mean_err_rpc(p) = mean(err_mag_rpc);
    mean_err_apc(p) = mean(err_mag_apc);

end % End (for p = 1 : num_images)

% End timer
tf = toc(t);

apc_sx_mean = mean(apc_sx, 1);
apc_sy_mean = mean(apc_sy, 1);

apc_sx_std = std(apc_sx, [], 1);
apc_sy_std = std(apc_sy, [], 1);

% Error bar plots
c_red = 1 / 255 * [178,34,34];
c_green = 	1 / 255 * [34, 139, 34] ;
c_gray = 0.5 * [1, 1, 1];
c_blue = 	1 / 255 * [30, 144, 255];

lw = 2;

fSize = 16;

figure;
semilogy(mean_err_scc, ':', 'color', c_blue, 'linewidth', lw);
hold on
plot(mean_err_rpc, ':', 'color', c_red, 'linewidth', lw);
plot(mean_err_apc, ':', 'color', c_gray, 'linewidth', lw);
hold off
% ylim([0, 0.2]);
xlim([0, 1000]);
axis square
% set(gca, 'ytick', 0 : 0.5 : max(ylim));
% set(gca, 'ytick', linspace(0, max(ylim), 5));

grid on
box on
xlabel('Number of ensemble pairs', ...
    'FontSize', fSize, 'interpreter', 'latex')
ylabel('Average translation error magnitude (pixels)',...
    'FontSize', fSize, 'interpreter', 'latex');
set(gca, 'fontsize', fSize);

h = legend('SCC', 'RPC', 'APC');
set(h, 'FontSize', fSize)
set(h, 'interpreter', 'latex');

% title('$\textrm{Poiseuille flow,} \, \sigma_{\Delta x} = 1.0 \,\, \textrm{pix / frame}$', 'interpreter', 'latex');

set(gcf, 'color', 'white');


% figure;
% hold on;
% g1 = cdfplot(err_mag_scc);
% g2 = cdfplot(err_mag_rpc);
% g3 = cdfplot(err_mag_apc);
% hold off;
% 
% % Set line widths
% set(g1, 'linewidth', lw);
% set(g2, 'linewidth', lw);
% set(g3, 'linewidth', lw);
% 
% % Set colors
% set(g1, 'color', c_blue);
% set(g2, 'color', c_red);
% set(g3, 'color', c_gray);
% 
% % Format the axes
% axis square;
% box on;
% set(gca, 'fontsize', fSize);
% 
% % Labels
% xlabel('$\textrm{Translation error magnitude (pixels)}$',...
%     'fontsize', fSize, ...
%     'interpreter', 'latex');
% 
% ylabel('$\textrm{Cumulative probability}$',...
%     'fontsize', fSize, ...
%     'interpreter', 'latex');
% 
% title('');
% 
% % Colors
% set(gcf, 'color', 'white');
% tightfig;
% 
% Print elapsed time
fprintf(1, 'Elapsed time: %d seconds for %d image pairs.\n', tf, num_pairs);
fprintf(1, '%0.1f seconds per image.\n\n', tf / num_pairs);




end































