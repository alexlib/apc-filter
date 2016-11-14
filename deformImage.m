function IMAGEOUT = deformImage(IMAGEIN, X, Y, U, V, METHOD)

% Default to no method sepecified
if nargin < 6
    METHOD = '';
end

% Measure the number of grid points in each direction
% The input has to be on a rectangular grid
% one way or another.
% Make the calling code deal with that part.
nx = length(unique(X));
ny = length(unique(Y));

% Image dimensions
[imageHeight, imageWidth] = size(IMAGEIN); 

% Create the pixel coordinates.
[xi_integer, yi_integer] = meshgrid(1:imageWidth, 1:imageHeight);

% Shift the pixel coordinates by 0.5 pixels
XI = xi_integer - 0.5;
YI = yi_integer - 0.5;

% Coordinates as matrices
y_array = reshape(Y, [ny, nx]);
x_array = reshape(X, [ny, nx]);
u_array = reshape(U, [ny, nx]);
v_array = reshape(V, [ny, nx]);

% Create interpolation structures for the velocity field.
interpolant_tx = griddedInterpolant(y_array, x_array, u_array, 'spline', 'linear');
interpolant_ty = griddedInterpolant(y_array, x_array, v_array, 'spline', 'linear');

% This is the velocity field upsampled to every pixel.
UI = interpolant_tx(YI, XI);
VI = interpolant_ty(YI, XI);

% These are the coordinates at which to resample image 2.
XD = XI + UI;
YD = YI + VI;

switch lower(METHOD)
    case 'blackman'
        % Resample the images using sinc-blacman interpolation (slower)   
        IMAGEOUT = whittaker_blackman(IMAGEIN, XD + 0.5, YD + 0.5, 8, true);
    otherwise
        % Resample the image using matlab interp2
        IMAGEOUT = interp2(IMAGEIN, XD + 0.5, YD + 0.5, 'spline', 0);
end

end