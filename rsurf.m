function rsurf(X, Y, Z)

[height, width] = size(X);

if nargin == 1
    
    Z = X;
   
    [X, Y] = meshgrid( 1:size(Z, 2), 1 : size(Z, 1));
    
end

surf(X, Y, real(Z));

xlim([1, width]);
ylim([1, height]);
axis square;
axis vis3d;

end