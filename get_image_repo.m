function img_repo = get_image_repo();

% Decide where the image repo is
if ismac
    img_repo = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/piv_test_images';
elseif islinux
    img_repo = '/home/shannon/c/aether/piv_test_images';
end


end