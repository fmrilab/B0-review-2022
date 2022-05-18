function plot_multi(fignum,subnum,im1,title_str,gt_im,scale_ims)
% plot multiple images one by one

% scale_ims == 2 will plot the error image

% Melissa W. Haskell, UMich, 2022

%%
figure(fignum); 
subplot(subnum); 
if isreal(im1), imagesc(im1)
else, imagesc(abs(im1));  end
title(title_str)
axis image; colormap gray; colorbar; axis off

if nargin > 4
    if nargin > 5
        if scale_ims == 1 || scale_ims == 2
            im1 = im1 ./ (max(abs(im1(:))) + eps);
            gt_im = gt_im ./ (max(abs(gt_im(:))) + eps);
        end
    end
    rmse = 100*norm(im1(:)-gt_im(:))/norm(gt_im(:));
    title(sprintf('%s - rmse %2.1f',title_str,rmse))
    if scale_ims == 2
        figure(fignum); 
        subplot(subnum); 
        imagesc(abs(im1-gt_im))
        title(title_str)
        axis image; colormap gray; colorbar; axis off
    end
end

end

