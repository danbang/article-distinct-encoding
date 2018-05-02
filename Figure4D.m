% Bang & Fleming (2018) Distinct encoding of decision confidence in human
% medial prefrontal cortex
%
% Individual-difference analysis of relationship between pgACC single-trial 
% activity estimates and confidence estimates
%
% Reproduces panel D in Figure 4
%
% Dan Bang danbang.db@gmail.com 2018

%% -----------------------------------------------------------------------
%% PREPARATION

% fresh memory
clear;

% Subjects
subjects = [1:9 11:13 15:20 22:35];

% ROIs
ROIs= {'pgacc'};
ROIs_for_plot= {'pgACC'};

% Paths [change 'repoBase' according to local setup]
fs = filesep;
repoBase = [getDropbox(1),fs,'Ego',fs,'Matlab',fs,'ucl',fs,'sensory_vs_decision',fs,'Repository'];
dataBehaviourDir = [repoBase,fs,'Data',fs,'Behaviour',fs,'Scan'];
dataScanDir = [repoBase,fs,'Data',fs,'Scan',fs,'Activity'];
dataModelDir = [repoBase,fs,'Model'];

% Add custom scripts
addpath('Custom');

% Load out-of-sample predictions
load([dataModelDir,fs,'conmodel_outOFsample_predictions.mat']);

% Load exclusions
load([dataScanDir,fs,'Tactivity_exclusions_RT1.mat']);

%% -----------------------------------------------------------------------
%% ANALYSIS

% Loop through ROIs
for i_roi = 1:length(ROIs)
% Loop through subjects
for i_sbj = 1:length(subjects)
    
    % Load ROI data
    file= [dataScanDir,fs,ROIs{i_roi},'_s',num2str(subjects(i_sbj)),'_RT1.mat'];
    load(file);    
    roi_ts = Tactivity;
    
    % load behavioural data
    
    % loop through scan runs
    for i_blk = 1:5     
        
    % load file
    file = [dataBehaviourDir,fs,'s',num2str(subjects(i_sbj)),'_task_b',num2str(i_blk),'.mat'];
    load(file);
    
    % get data field names
    fn = fieldnames(data);
    
    % if first block, then initialise temporary storage structure
    if i_blk == 1; for i_field = 1:length(fn); eval(['tmp.',fn{i_field},'=[];']); end; end
   
    % add data to temporary storage structure
    for i_field = 1:length(fn); eval(['tmp.',fn{i_field},'=[tmp.',fn{i_field},' data.',fn{i_field},'];']); end
    
    end
    
    % Re-assign
    data = tmp;
    
    % Include trials based on deviation from grand mean: RT
    rt1= log(data.rt1./1000);
    centre= mean(rt1);
    stdval= std(rt1)*2.5;
    include= (rt1>(centre-stdval))&(rt1<(centre+stdval));
    
    % Include trials based on deviation from grand mean: single-trial beta
    include2= exclusion.subject{i_sbj}==0;
    include= include & include2;
    
    % z-scoring
    roi_ts(include)= zscore(roi_ts(include)); 
    data.con(data.concat&include)=zscore(data.con(data.concat&include));

        
    %% NEURAL INTERACTION: DIFFERENCE OF DIFFERENCES
    interactionROI{i_roi}(i_sbj,1)= (mean(roi_ts(data.cohcat==2&data.bodcat==2&include))-...
                                    mean(roi_ts(data.cohcat==2&data.bodcat==1&include)))-...
                                    (mean(roi_ts(data.cohcat==1&data.bodcat==2&include))-...
                                    mean(roi_ts(data.cohcat==1&data.bodcat==1&include)));
                         
    %% NEURAL INTERACTION: DIFFERENCE OF DIFFERENCES
    interactionCON{i_roi}(i_sbj,1)= (mean(data.con(data.cohcat==2&data.bodcat==2&include&data.concat))-...
                             mean(data.con(data.cohcat==2&data.bodcat==1&include&data.concat)))-...
                             (mean(data.con(data.cohcat==1&data.bodcat==2&include&data.concat))-...
                             mean(data.con(data.cohcat==1&data.bodcat==1&include&data.concat)));
                         
    %% NEURAL INTERACTION: DIFFERENCE OF DIFFERENCES
    interactionACC{i_roi}(i_sbj,1)= (mean(data.acc(data.cohcat==2&data.bodcat==2&include&data.concat))-...
                             mean(data.acc(data.cohcat==2&data.bodcat==1&include&data.concat)))-...
                             (mean(data.acc(data.cohcat==1&data.bodcat==2&include&data.concat))-...
                             mean(data.acc(data.cohcat==1&data.bodcat==1&include&data.concat)));
    
    end
end

%% -----------------------------------------------------------------------
%% FIGURES

%% Loop through ROIs
for i_roi= 1:length(ROIs);
figure('color',[1 1 1]);
plot([1000 1001],[1000 1001],'-','color',[255/255 51/255 51/255],'LineWidth',4); hold on;
plot([-2 +2],[0 0],'k--','LineWidth',2); hold on;
plot([0 0],[-2 +2],'k--','LineWidth',2); hold on;
X= interactionROI{i_roi}; Y=interactionCON{i_roi}/3;
plot(X,Y,'ko','MarkerSize',8,'MarkerFaceColor','k'); hold on;
[B,STATS] = robustfit(X,Y);
plot(X,B(1)+B(2)*X,'-','color',[255/255 51/255 51/255],'LineWidth',4); hold on;
set(gca,'XTick',[-.4:.2:.4]); set(gca,'YTick',[-.4:.2:.4]);
xlim([-0.5 +0.5]);ylim([-.5 +.5]);
set(gca,'FontSize',24,'LineWidth',2);
title(ROIs_for_plot{i_roi},'FontSize',28,'FontWeight','normal')
legend(['p = ', num2str(STATS.p(end))],'Location','NorthWest');
legend('boxoff');
xlabel('neural interaction [a.u.]','FontSize',28,'FontWeight','normal');
ylabel('conf. interaction [a.u.]','FontSize',28,'FontWeight','normal');
end