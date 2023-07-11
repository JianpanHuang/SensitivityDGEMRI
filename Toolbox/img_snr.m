function snr = img_snr(img, mask_sig, mask_noi)
    if nargin<3
        mask_noi = mask_sig;
    end
    mean_val = mean2(img(mask_sig == 1));
    std_val = std2(img(mask_noi == 1));
    snr = 20*log10(mean_val/std_val);
end