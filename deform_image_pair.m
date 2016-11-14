function [IMAGE_OUT_01, IMAGE_OUT_02] = ...
    deform_image_pair(IMAGE_01, IMAGE_02, X, Y, U, V, METHOD);

% Default to not specifying
% an interpolation method
% which means use Matlab interp2
if nargin < 7
    METHOD = '';
end
%
% % All this commented out section is for 
% later on modifying the grid to be rectangular
% so that the deform code can accept non-rectangular grids.
%
%
% % Make the grid rectangular.
% x_min = min(X(:));
% x_max = max(X(:));
% y_min = min(Y(:));
% y_max = max(Y(:));
% dx = max(diff(unique(X(:))));
% dy = max(diff(unique(Y(:))));
% 
% % Grid vectors
% xv_new = x_min : dx : x_max;
% yv_new = y_min : dy : y_max;
% 
% % Size of the whole array
% nx = length(xv_new);
% ny = length(yv_new);
% 
% % Number of grid measurements
% num_measurements = length(X(:));
% 
% % Grid mesh (whole thing)
% [x_rect, y_rect] = meshgrid(xv_new, yv_new);
% 
% % Make arrays for rectangular
% % grids of measurements
% u_rect = zeros(ny, nx);
% v_rect = zeros(ny, nx);
% 
% % Allocate the indices
% inds = zeros(num_measurements, 1);
% 
% % Find grid indices of measured values
% for k = 1 : num_measurements
%     inds(k) = find(x_rect == X(k) & y_rect == Y(k));   
% end
% 
% % Find indices where the velocity was measured
% measured_inds = find(x_rect == X & y_rect == Y);

% The first image gets deformed "forwards"
% which means the input velocity 
% should be in the opposite direction
% as the measured velocity
u_source_01 = -1 * U / 2;
v_source_01 = -1 * V / 2;

% The second image gets deformed "backwards"
% which means the input velocity 
% should be in the direction
% of the measured velocity
u_source_02 = 1 * U / 2;
v_source_02 = 1 * V / 2;

% Deform the images.
IMAGE_OUT_01 = deformImage(...
    IMAGE_01, X, Y, u_source_01, v_source_01, METHOD);

IMAGE_OUT_02 = deformImage(...
    IMAGE_02, X, Y, u_source_02, v_source_02, METHOD);


end



