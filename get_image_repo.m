function img_repo = get_image_repo();

% Decide where the image repo is
if ismac
    img_repo = '~/Desktop/piv_test_images';
%     img_repo = '/Volumes/aether_c/piv_test_images';
elseif islinux
    img_repo = '/home/shannon/c/aether/piv_test_images';
end


end