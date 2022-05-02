% Figure 5 from B0 review


clear; close all

%% Initialize parameters

% plotting parameters
lw1 = 3; % line width
lw2 = 2;

% Gaussian parameters for simulated fieldmap
g1 = .16; g2 = .5; 
max_b0 = 400; % max B0 in Hz
b0_bump_size = .5;   % size of gaussian relative to full image
b0_bump_start = .5; % where to start gaussian within the image

% Input parameters in standard physical units
pos_x = -10:.01:10;   % position in cm
pos_y = pos_x;
npos = numel(pos_x);
Gx = 1;   % gradient strength in mT/m
Gy = 0;   % gradient strength in mT/m
TE = 1;    % echo time in ms


%% Simulate 2d field 

% convert units to Hz
ideal_field_2d = repmat(Gx.*(pos_x/100),[npos 1]) + repmat(Gy.*(pos_y'/100),[1 npos]);  % ideal field in mT
ideal_field_2d = ideal_field_2d/1000; % ideal field in T
ideal_freq_2d = 42.577 * ideal_field_2d; % ideal freq in MHz
ideal_freq_2d = ideal_freq_2d * 1e6; % ideal freq in Hz

% calculate ideal phase assuming perfect rect gradient pulse
ideal_phase_2d = mod((ideal_freq_2d * 2*pi * (TE/1000))+pi,2*pi)-pi;  % ideal encoded phase in radians

% create fieldmap from params
x = linspace(0,1,round(b0_bump_size*npos));
y = max_b0*gaussmf(x,[g1 g2]);
y_2d= (y'*y)./max(y(:));
df_start = round(b0_bump_start*npos);
df_2d = zeros(npos,npos);
df_2d(df_start:df_start+numel(y)-1,round(end/2-numel(y)/2):round(end/2+numel(y)/2-1)) = y_2d;
df_2d = flipud(df_2d);


%% calculate actual phase with off-resonance

act_freq_2d = ideal_freq_2d + df_2d; % (again assumes perfect rect gradient)
act_phase_2d = mod((act_freq_2d * 2*pi * (TE/1000))+pi,2*pi)-pi;


%% Plot comparison of ideal and actual phase (Figure 5)

load('axial_mask.mat'); 
roi2d = double(roi2d); roi2d = imresize(roi2d,[npos npos]);
roi2d = logical(roi2d);
figure('units','normalized','outerposition',[0 0 1 1])
lpos = npos-df_start - round(numel(y)/2);

ax2 = subplot(221); imagesc(df_2d); colorbar; axis image; title('fieldmap (Hz)'); axis off
line([0,npos],[lpos,lpos],'Color','red','LineStyle','--','LineWidth',lw1)

subplot(222);hold on;
plot(pos_x,ideal_freq_2d(lpos,:),'--','LineWidth',lw2,'Color','green'); 
plot(pos_x,df_2d(lpos,:),'--','Color','red','LineWidth',lw2); 
plot(pos_x,act_freq_2d(lpos,:),'-','Color','blue','LineWidth',lw2) 
legend('ideal gradient field','delta f','actual encoding field')
xlabel('position (cm)'); ylabel('Frequency (Hz)'); 
axis square
title(sprintf('Gx = %s mT/m, Gy = %s mT/m, TE = %s ms',num2str(Gx),num2str(Gy),num2str(TE)))

ax4 = subplot(223); imagesc(roi2d.*ideal_phase_2d)
colorbar; axis image; title('ideal phase encoding in radians (masked)'); axis off
ax6 = subplot(224); imagesc(roi2d.*act_phase_2d); colorbar; axis image; title('actual phase encoding in radians (masked)'); axis off
colormap(ax2,gray)
colormap(ax4,hsv)
colormap(ax6,hsv)

set(gcf,'color','w');


%% Add fieldmaps to plot (not shown in paper)

figure('units','normalized','outerposition',[0 0 1 1])

ax1 = subplot(232); imagesc(roi2d.*ideal_freq_2d); colorbar; axis image; axis off; title('ideal encoding field in Hz (masked)')
line([0,npos],[lpos,lpos],'Color','green','LineStyle','--','LineWidth',lw1)
ax2 = subplot(231); imagesc(df_2d); colorbar; axis image; title('fieldmap (Hz)'); axis off
line([0,npos],[lpos,lpos],'Color','red','LineStyle','--','LineWidth',lw1)
ax3 = subplot(233); imagesc(roi2d.*act_freq_2d); colorbar; axis image; title('actual encoding field (in Hz masked)'); axis off
line([0,npos],[lpos,lpos],'Color','blue','LineStyle','--','LineWidth',lw1)

subplot(234);hold on;
plot(pos_x,ideal_freq_2d(lpos,:),'--','LineWidth',lw2,'Color','green'); 
plot(pos_x,df_2d(lpos,:),'--','Color','red','LineWidth',lw2); 
plot(pos_x,act_freq_2d(lpos,:),'-','Color','blue','LineWidth',lw2) 
legend('ideal gradient field','delta f','actual encoding field')
xlabel('position (cm)'); ylabel('Frequency (Hz)'); 
axis square
title(sprintf('Gx = %s mT/m, Gy = %s mT/m, TE = %s ms',num2str(Gx),num2str(Gy),num2str(TE)))

ax4 = subplot(235); imagesc(roi2d.*ideal_phase_2d)
colorbar; axis image; title('ideal phase encoding in radians (masked)'); axis off
ax6 = subplot(236); imagesc(roi2d.*act_phase_2d); colorbar; axis image; title('actual phase encoding in radians (masked)'); axis off
colormap(ax1,parula)
colormap(ax2,gray)
colormap(ax3,parula)
colormap(ax4,hsv)
colormap(ax6,hsv)

set(gcf,'color','w');

