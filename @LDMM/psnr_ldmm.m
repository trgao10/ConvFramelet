function y = psnr_ldmm(f, g, MAX_PIXEL_VAL)

y = -20*log10(norm(f-g,'fro')/sqrt(numel(f))/MAX_PIXEL_VAL);

end