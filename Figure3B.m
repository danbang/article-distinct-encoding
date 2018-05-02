% Bang & Fleming (2018) Distinct encoding of decision confidence in human
% medial prefrontal cortex
%
% Visualises ROI GLM1 contrast estimates
%
% ROIs: pgACC (leave-one-out), MT+ (localiser), IPS (leave-one-out),
% ventral striatum (anatomical), pre-SMA (leave-one-out), rlPFC (anatomical)
%
% Reproduces panel B in Figure 3
%
% Dan Bang danbang.db@gmail.com 2018

%% -----------------------------------------------------------------------
%% PREPARATION

% fresh memory
clear;

% Subjects
subjects = [1:9 11:13 15:20 22:35];

% Paths [change 'repoBase' according to local setup]
fs = filesep;
repoBase = [getDropbox(1),fs,'Ego',fs,'Matlab',fs,'ucl',fs,'sensory_vs_decision',fs,'Repository'];
dataScanDir = [repoBase,fs,'Data',fs,'Scan',fs,'Contrast'];

% ROIs
ROIs= {'pgacc','V5','IPS','vstriatum','preSMA','rlPFC'};
ROIs_for_plot= {'pgACC','MT+','IPS','ventral striatum','pre-SMA','rlPFC'};

%% -----------------------------------------------------------------------
%% ANALYSIS

% Loop through ROIs
for i_roi = 1:length(ROIs)
% Loop through subjects
for i_sbj = 1:length(subjects)
    % Load subject data
    file= [dataScanDir,fs,ROIs{i_roi},'_s',num2str(subjects(i_sbj)),'.mat'];
    load(file);
    % Log subject data
    c_data(i_sbj,:)= Cestimate;
end
% Plot ROI data
colorz=[255 51 51;
        255 0 255;
        0 255 0]./255;
figure('color',[1 1 1]);
plot([0 4.5],[0 0],'k-','LineWidth',4); hold on;
for i_con= 1:3;
    y= c_data(:,i_con);
    n = numel(y);
    mu = mean(y);
    sem = std(y)./sqrt(n);
    bar(i_con,mu,'FaceColor',colorz(i_con,:),'FaceAlpha',.7,'LineWidth',4); hold on;
    plot([i_con i_con],[mu-sem mu+sem],'k-','LineWidth',4); hold on;
end
title([ROIs_for_plot{i_roi}],'Fontsize',46,'FontWeight','normal')
ylim([-3.5 3.5]); 
set(gca,'YTick',-3:1:3);
xlim([0 4]);
set(gca,'Xtick',[1:3],'XTickLabel',{'C','D','CxD'});
set(gca,'FontSize',50,'LineWidth',4);
box('off')
end