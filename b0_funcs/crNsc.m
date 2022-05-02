function [imout] = crNsc(im_in,ccrop)

im_in = im_in(:,ccrop+1:end-ccrop);
im_in_sc = max(abs(im_in(:)));
imout = im_in ./(im_in_sc + eps);

end

