function [adj_dcf_im] = adj_dcf(kspace,nshot,kdata,mask,FOV,dcf_plot)

% reshape multishot spiral and get intial DCF from voronoi method
kx = kspace(:,1); ky = kspace(:,2);
npts = numel(kx)/nshot;
dcf0 = voronoidens(kx(:),ky(:));

%%%% crop NaNs in voronoi
cr_extra = .2;
dcf_crp = round( (1+cr_extra) * sum(isnan(dcf0))/nshot );

% put dcf back into shot dimensions and crop
dcf_resh = reshape(dcf0,[numel(dcf0)/nshot, nshot]);
dcf = dcf_resh(1:end-dcf_crp,:);

% crop kspace locations at edge of each spiral shot
kspace_tmp = reshape(kspace,[npts, nshot, 2]);
kspace_cr = kspace_tmp(1:end-dcf_crp,:,:);
kx_dcf = kspace_cr(:,:,1);
ky_dcf = kspace_cr(:,:,2);

% resize data and crop data points with bad DCF from voronoi
kdata_tmp = reshape(kdata,[npts, nshot]);
kdata_cr = kdata_tmp(1:end-dcf_crp,:);

if nargin > 5
    if dcf_plot
        figure(11); scatter(kx(:),ky(:),'.'); hold on; scatter(kx_dcf(:),ky_dcf(:),'.')
    end
end

% make new imaging operator with cropped kspace
A = Gmri([kx_dcf(:),ky_dcf(:)], mask, 'fov', FOV);

% scale data using the DCF
kdata_cr_DCF = kdata_cr(:).*dcf(:);

% get new image from adjoint+DCF operation, and reshape
adj_dcf_im = reshape(A'*kdata_cr_DCF(:),size(mask));


end

