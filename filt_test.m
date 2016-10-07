apc_filter_test_single_zero_mean_ncc;

load('~/Desktop/apc_filt_s_rand_3.mat');

g = fftshift(abs(ifft2(fftshift(cc_full_sum))));

cc_apc_spect = phaseOnlyFilter(cc_full_sum) .* apc_filt; 

g2 = fftshift(abs(ifft2(fftshift(cc_apc_spect))));

rpc_filter = spectralEnergyFilter(region_width, region_height, 3);
cc_rpc_spect = phaseOnlyFilter(cc_full_sum) .* rpc_filter;
g3 = fftshift(abs(ifft2(fftshift(cc_rpc_spect))));



figure(2); 
subplot(1, 3, 1); 
surf(g); 
axis square; 
set(gca, 'view', [0, 0]); 

subplot(1, 3, 2); 
surf(g3); 
axis square; 
set(gca, 'view', [0, 0]);

subplot(1, 3, 3); 
surf(g2); 
axis square; 
set(gca, 'view', [0, 0]);