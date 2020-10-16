%% paired_t_roi_snr
clear all; 
load /users/PAS1342/osu9912/fmri/isss_multi/isss_multiband_122119/batch/roi/mask_lng_noi_042020/hybrid_isss/mean_SNR_unique.mat

x = mean_SNR(1, :)'; % hybrid
y = mean_SNR(2, :)'; % isss

[h1, p1] = ttest(x, y)

load /users/PAS1342/osu9912/fmri/isss_multi/isss_multiband_122119/batch/roi/mask_lng_noi_042020/hybrid_multiband/mean_SNR_unique.mat
x = mean_SNR(1, :)'; % hybrid
y = mean_SNR(2, :)'; % multiband

[h2, p2] = ttest(x, y)
