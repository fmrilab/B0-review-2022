% run_sims.m
%
% Main script to run the simulations in Haskell, Nielsen, and Noll,
% NMR in Biomedicine 2022 off-resonance correction review.
%
%
% Requires Michigan Image Reconstruction Toolbox for MATLAB:
%       https://github.com/JeffFessler/mirt
%
%
% Melissa Haskell, University of Michigan
% Spring 2022


%%
clear; close all
addpath('./b0_funcs')
addpath('./mat_files')


%% Initialize variables needed for simulation

im_spi = 30; im_epi = 40;

% image parameters
n = 180;        % image size
FOV = 24;      % field of view (cm)

% additional variables for Gmri object
L = 8;
ov = 8;

%% Load image and bmap

load('t1_image.mat'); 


%% Resize and prep for simulation

xtrue = imresize(t1_im,[n n]);
mask = true(n);
mask_ov = true(ov*n);
xtrue_ov = imresize(t1_im,ov*[n n]);

b0_scale = 1;
bmap = b0_scale*imresize(bmap,[n n]);
bmap_ov = b0_scale*imresize(bmap,ov.*[n n]);

% create zmap (in i*radians/s instead of Hz)
zmap = 1i*2*pi*bmap;
zmap_ov = 1i*2*pi*bmap_ov;


%% Load and plot k-space data trajectories

load('spiral_traj.mat')
ksp_spi = ktraj; ti_spi = t_s; clear ktraj; 

load('epi_traj.mat')
ksp_epi = kspace; ti_epi = t_s; clear kspace; 
ksp_epi = [ksp_epi(end:-1:1,2),ksp_epi(end:-1:1,1)]; % for data axis & direction flipping

figure(1); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
ti_epi_ms = ti_epi * 1000;
ti_spi_ms = ti_spi * 1000;
zm = 0.25; % zoom in crop for ktraj data points
sz1 = 8; sz2 = 6; % scatter plot dot sizes
subplot(231); plot(ti_epi_ms); title('epi time (ms)'); axis square
ylabel('time (ms)'); xlabel('data point index')
subplot(234); plot(ti_spi_ms); title('spiral time (ms)'); axis square
ylabel('time (ms)'); xlabel('data point index')
subplot(235); scatter(ksp_spi(:,1),ksp_spi(:,2),sz1,ti_spi_ms,'filled'); axis image; colorbar
ylabel('kx (1/cm)'); xlabel('ky (1/cm)'); title('spiral traj, colorbar in ms')
subplot(232); scatter(ksp_epi(:,1),ksp_epi(:,2),sz1,ti_epi_ms,'filled'); axis image; colorbar
ylabel('kx (1/cm)'); xlabel('ky (1/cm)'); title('epi traj, colorbar in ms')
 
subplot(236); scatter(ksp_spi(:,1),ksp_spi(:,2),sz2,ti_spi_ms,'filled'); 
axis image; colorbar; xlim([-zm zm]); ylim([-zm zm])
ylabel('kx (1/cm)'); xlabel('ky (1/cm)'); title('spiral traj (zoomed), colorbar in ms')
subplot(233); scatter(ksp_epi(:,1),ksp_epi(:,2),sz2,ti_epi_ms,'filled'); axis image; colorbar
axis image; colorbar; xlim([-zm zm]); ylim([-zm zm])
ylabel('kx (1/cm)'); xlabel('ky (1/cm)'); title('epi traj (zoomed), colorbar in ms')


%% Create imaging forward models with and without B0 effects

A_spi = Gmri(ksp_spi, mask,'fov', FOV);
A_spi_ov = Gmri(ksp_spi, mask_ov,'fov', FOV);
A_wB0_spi = Gmri(ksp_spi, mask,'fov', FOV, 'zmap', zmap, 'ti', ti_spi, 'L', L);
A_wB0_spi_ov = Gmri(ksp_spi, mask_ov,'fov', FOV, 'zmap', zmap_ov, 'ti', ti_spi, 'L', L);

A_epi = Gmri(ksp_epi, mask,'fov', FOV);
A_wB0_epi = Gmri(ksp_epi, mask,'fov', FOV, 'zmap', zmap, 'ti', ti_epi, 'L', L);


%% Simulate data acquired trajectories
kdata_spi = A_spi *  xtrue(:);
kdata_spi_ov = A_spi_ov *  xtrue_ov(:);
kdata_wB0_spi = A_wB0_spi *  xtrue(:);
kdata_wB0_spi_ov = A_wB0_spi_ov *  xtrue_ov(:);

kdata_epi = A_epi *  xtrue(:);
kdata_wB0_epi = A_wB0_epi *  xtrue(:);


%% Recon using density compensated adjoint spiral, adjoint if epi

x_recon_spi = circmask(adj_dcf(ksp_spi,nshot,kdata_spi,mask,FOV));
x_recon_spi_ov = circmask(adj_dcf(ksp_spi,nshot,kdata_spi_ov,mask,FOV));
x_recon_wB0art_spi = circmask(adj_dcf(ksp_spi,nshot,kdata_wB0_spi,mask,FOV));
x_recon_wB0art_spi_ov = circmask(adj_dcf(ksp_spi,nshot,kdata_wB0_spi_ov,mask,FOV));

% direct "inversion" with adjoint, works ok for EPI bc roughly a FFT
x_recon_epi = reshape(A_epi'*kdata_epi,[n n])./numel(xtrue);
x_recon_wB0art_3T_epi = reshape(A_epi'*kdata_wB0_epi,[n n])./numel(xtrue);


%% images used in paper (Figure 6)
cr = 20; cind = cr+1:n-cr;
row1 = cat(2,crNsc(xtrue,cr),crNsc(x_recon_wB0art_3T_epi,cr),...
    crNsc(x_recon_wB0art_spi_ov,cr));

figure(100); hold off
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
imagesc(abs(cat(1,row1)),[0 .7]); colormap gray; axis image; axis off
set(gcf,'color','w')

load('xtrue_mask.mat')
bounds = bwboundaries(xtrue_mask); outline=bounds{1};

lw = 1;
hold on; 
plot(outline(:,2)-cr,outline(:,1),'g','LineWidth',lw)
plot(outline(:,2)-(3*cr)+n,outline(:,1),'g','LineWidth',lw)
plot(outline(:,2)-(5*cr)+(2*n),outline(:,1),'g','LineWidth',lw)



