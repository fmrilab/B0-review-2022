%% Support functions
function Imask = circmask(I)
npix = size(I,1);
[xmask,ymask] = meshgrid(1:npix, 1:npix);
bw_circ = (xmask - ((npix/2)+0.5)).^2 + (ymask - ((npix/2)+0.5)).^2 <= (npix/2)^2;
Idims = size(I);
Imask = I.*repmat(bw_circ,[1 1 Idims(3:end)]);
end