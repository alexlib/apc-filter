function [y, x] = get_mask_outline(grid_mask);

% Measure the size of the image
[image_height, image_width] = size(grid_mask);

% Do edge detection on the binary mask
mask_edges = edge(grid_mask, 'sobel');

% Indices of non-zero points (these should be the mask edges)
non_zero_inds = find(mask_edges);

% Convert to [x, y] coordinates
[y, x] = ind2sub([image_height, image_width], non_zero_inds);

end