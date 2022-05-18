% mbir.m
%
% Performs a model-based image reconstruction (MBIR) that includes the
% effects of B0 off-resonance in the forward model using the Michigan Image 
% Reconstruction Toolbox (MIRT) for MATLAB.
%
% This script was written to accompany the NMR in Biomedicine 2022 
% off-resonance correction review by Haskell, Nielsen, and Noll.
%
%
% Requires Michigan Image Reconstruction Toolbox for MATLAB:
%       https://github.com/JeffFessler/mirt
%
%
% Melissa W. Haskell, University of Michigan
% Spring 2022


clear; close all
addpath('./b0_funcs')
addpath('./mat_files')


%% Initialize variables needed for simulation

% image parameters
n = 180;        % image size
FOV = 24;       % field of view (cm)

% additional variables for Gmri object (more info below)
L = 8;

% MBIR parameters
W = 1; % data weighting matrix (can be scalar if all values the same)
C = 0; % penalty 'differencing matrix' (0 for unregularized)
niter = 20; 

%% Load image and its B0 fieldmap, then load kspace trajectory

load('t1_image.mat'); 
load('spiral_traj.mat')


%% Resize image and prep for simulation

xtrue = imresize(t1_im,[n n]);
mask = true(n);

bmap = imresize(bmap,[n n]);

% create zmap (in i*radians/s instead of Hz)
zmap = 1i*2*pi*bmap;


%% Create imaging forward model with B0 effects

% no B0 effect forward model for comparison
A = Gmri(ktraj, mask,'fov', FOV);

% Forward model including B0 effects, which is done by adding:
%  1. zmap - fieldmap plus relaxation map, but relaxation ignored here
%  2. ti - time of each data sample, in units inverse to the zmap, here 
%          the fieldmap is in Hz, so t is in seconds (s)
%  3. L - number of basis functions in the complex exponential. This is the
%         same as L in Eqns. 17 & 23 in the review paper
A_wB0 = Gmri(ktraj, mask,'fov', FOV, 'zmap', zmap, 'ti', t_s, 'L', L);


%% Simulate data by applying the forward model to the ground truth image
kdata = A *  xtrue(:);
kdata_wB0 = A_wB0 *  xtrue(:);


%% Model-based image recon using qpwls_pcg1 from MIRT
% run recons and crop out edges from spiral undersampling aliasing using
% circmask

x_mbir_no_B0 = circmask(qpwls_pcg1_wrap([n n],A,W,kdata,C,niter)); % (not plotted)
x_mbir_B0_corrupted = circmask(qpwls_pcg1_wrap([n n],A,W,kdata_wB0,C,niter));
x_mbir_B0_corrected = circmask(qpwls_pcg1_wrap([n n],A_wB0,W,kdata_wB0,C,niter));

%% Plot images together

mbir_fig = 10; 
plot_multi(mbir_fig,131,xtrue,'ground truth')
plot_multi(mbir_fig,132,x_mbir_B0_corrupted,sprintf('B0 corrupted - %d iters',niter),xtrue,1)
plot_multi(mbir_fig,133,x_mbir_B0_corrected,sprintf('B0 corrected - %d iters',niter),xtrue,1)







