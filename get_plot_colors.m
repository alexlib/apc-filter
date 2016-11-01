function rgb_list =  get_plot_colors(N)

if nargin == 0
    num_colors = 3;
else
    num_colors = N;
end

if num_colors == 5
    rgb_list = zeros(5, 3);
    rgb_list(1, :) = [215, 25, 28];
    rgb_list(2, :) = [253, 174, 97];
    rgb_list(3, :) = [0, 0, 0];
    rgb_list(4, :) = [171, 221, 164];
    rgb_list(5, :) = [43, 131, 186];

else
    rgb_list = zeros(3, 3);
    rgb_list(1, :) = [215, 25, 28];
    rgb_list(2, :) = [30, 180, 255];
    rgb_list(3, :) = [0, 0, 0];
end

rgb_list = rgb_list ./ 255;



end