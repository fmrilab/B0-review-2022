% run_fig6.m
%
% Script to run Figure 6 in Haskell, Nielsen, and Noll,
% "Off-resonance artifact correction for magnetic resonance
% imaging: a review", NMR in Biomedicine 2022.
%
%
% Requires Michigan Image Reconstruction Toolbox for MATLAB:
%       https://github.com/JeffFessler/mirt
%
%
% Melissa Haskell, University of Michigan
% Jon-Fredrik Nielsen, University of Michigan
% 2022
%


%%
clear; close all
addpath('./b0_funcs')
addpath('./mat_files')


%% Initialize variables needed for simulation

% image parameters
n = 180;        % image size
FOV = 24;      % field of view (cm)

% additional variables for Gmri object
L = 39;
ov = 8;

%% Load image and simulate linear B0 map

load('t1_image.mat','t1_im');
n1 = 76;
bmap = 70*ndgrid(linspace(-1,1,n1), linspace(-1,1,n1));


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

% cartesian kspace
dk = 1/FOV;  % kspace sample spacing
kmax = n*dk/2 - dk/2;
[ksp_epi_x, ksp_epi_y] = ndgrid(linspace(-kmax, kmax, n), linspace(-kmax, kmax, n));
ksp_epi = [ksp_epi_y(:) ksp_epi_x(:)];
ti_epi = (1 + ksp_epi_y(:)/kmax)/2 * 100e-3;

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
ylabel('ky (1/cm)'); xlabel('kx (1/cm)'); title('spiral traj, colorbar in ms')
subplot(232); scatter(ksp_epi(:,2),ksp_epi(:,1),sz1,ti_epi_ms,'filled'); axis image; colorbar
ylabel('ky (1/cm)'); xlabel('kx (1/cm)'); title('epi traj, colorbar in ms')

subplot(236); scatter(ksp_spi(:,1),ksp_spi(:,2),sz2,ti_spi_ms,'filled');
axis image; colorbar; xlim([-zm zm]); ylim([-zm zm])
ylabel('ky (1/cm)'); xlabel('kx (1/cm)'); title('spiral traj (zoomed), colorbar in ms')
subplot(233); scatter(ksp_epi(:,2),ksp_epi(:,1),sz2,ti_epi_ms,'filled'); axis image; colorbar
axis image; colorbar; xlim([-zm zm]); ylim([-zm zm])
ylabel('ky (1/cm)'); xlabel('kx (1/cm)'); title('epi traj (zoomed), colorbar in ms')


%% Create imaging forward models with and without B0 effects

A_spi = Gmri(ksp_spi, mask,'fov', FOV);
A_spi_ov = Gmri(ksp_spi, mask_ov,'fov', FOV);
A_wB0_spi = Gmri(ksp_spi, mask,'fov', FOV, 'zmap', zmap, 'ti', ti_spi, 'L', L);
A_wB0_spi_ov = Gmri(ksp_spi, mask_ov,'fov', FOV, 'zmap', zmap_ov, 'ti', ti_spi, 'L', L);

A_epi = Gmri(ksp_epi, mask,'fov', FOV);
N = [n n];
nufft_args = {N, 6*ones(size(N)), 2*N, N/2};
A_wB0_epi = Gmri(ksp_epi, mask,'fov', FOV, 'zmap', zmap, 'ti', ti_epi, 'L', L); %, 'nufft', nufft_args);


%% Simulate data acquired with both trajectories
kdata_spi = A_spi *  xtrue(:);
kdata_spi_ov = A_spi_ov *  xtrue_ov(:);
kdata_wB0_spi = A_wB0_spi *  xtrue(:);
kdata_wB0_spi_ov = A_wB0_spi_ov *  xtrue_ov(:);

kdata_epi = A_epi *  xtrue(:);
kdata_wB0_epi = A_wB0_epi *  xtrue(:);


%% Recon using density compensated adjoint if spiral, iFFT if epi

x_recon_spi = circmask(adj_dcf(ksp_spi,nshot,kdata_spi,mask,FOV));
x_recon_spi_ov = circmask(adj_dcf(ksp_spi,nshot,kdata_spi_ov,mask,FOV));
x_recon_wB0art_spi = circmask(adj_dcf(ksp_spi,nshot,kdata_wB0_spi,mask,FOV));
x_recon_wB0art_spi_ov = circmask(adj_dcf(ksp_spi,nshot,kdata_wB0_spi_ov,mask,FOV));

% kspace is now cartesian so do ift
x_recon_wB0art_3T_epi = fftshift(ifftn(fftshift(reshape(kdata_wB0_epi, [n n]))))';


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



