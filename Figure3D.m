% Bang & Fleming (2018) Distinct encoding of decision confidence in human
% medial prefrontal cortex
%
% Reward-magnitude analysis of ROI activity timecourses
%
% Loads permutation-testing statistics for indicating significance (see folder
% 'Permutation' for permutation procedure)
%
% ROIs: pgACC (leave-one-out), MT+ (localiser), IPS (leave-one-out),
% ventral striatum (anatomical), pre-SMA (leave-one-out), rlPFC (anatomical)
%
% Reproduces panel D in Figure 3
%
% Dan Bang danbang.db@gmail.com 2018

%% -----------------------------------------------------------------------
%% PREPARATION

% fresh memory
clear;

% Subjects
subjects = [1:9 11:13 15:20 22:35];

% ROIs
ROIs= {'pgacc','V5','IPS','vstriatum','preSMA','rlPFC'};
ROIs_for_plot= {'pgACC','MT+','IPS','ventral striatum','pre-SMA','rlPFC'};

% Paths [change 'repoBase' according to local setup]
fs = filesep;
repoBase = [getDropbox(1),fs,'Ego',fs,'Matlab',fs,'ucl',fs,'sensory_vs_decision',fs,'Repository'];
dataBehaviourDir = [repoBase,fs,'Data',fs,'Behaviour',fs,'Scan'];
dataScanDir = [repoBase,fs,'Data',fs,'Scan',fs,'Timecourse'];
dataPermutationDir = [repoBase,fs,'Permutation'];

% Add custom scripts
addpath('Custom');

%% -----------------------------------------------------------------------
%% ANALYSIS

% Loop through ROIs
for i_roi = 1:length(ROIs)
% Loop through subjects
for i_sbj = 1:length(subjects)
    
    % Load ROI data
    file= [dataScanDir,fs,ROIs{i_roi},'_s',num2str(subjects(i_sbj)),'_scale.mat'];
    load(file);    
    roi_ts = Tcourse;
    
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
    
    % Include trials based on deviation from grand mean
    rt1= log(data.rt1./1000);
    centre= mean(rt1);
    stdval= std(rt1)*2.5;
    include= (rt1>(centre-stdval))&(rt1<(centre+stdval));
    
    % Include trials where confidence estimates were required
    for i= 1:size(data.concat,1); if data.concat(i)==0; include(i)=0; end; end;
    
    % Reshape ROI time series to all trials
    zeropad= NaN(560,85); j=0;
    for i= 1:size(zeropad,1); if data.concat(i); j=j+1; zeropad(i,:)= roi_ts(j,:); end; end
    roi_ts= zeropad;
    
    % Include trials where final time-point estimate is ~NaN
    for i= 1:size(roi_ts,1); if isnan(roi_ts(i,end)); include(i)=0; end; end;
    
    % UP-SAMPLED GLM
    roi_Zts = zscore(roi_ts(include,:));
    rwd= zscore(data.rwdcat(include)-1.5);
    t= 0;
    for j= 1:size(roi_ts,2)
        t= t+1;
        x= rwd';
        y= roi_Zts(:,j);
        beta= glmfit(x,y);
        beta_ts{i_roi}.rwd(i_sbj,t)= beta(end);
    end

    end
end

%% -----------------------------------------------------------------------
%% FIGURES

%% Load permutation output
load([dataPermutationDir,fs,'timecourse_permutation_reward.mat']);
pmt_coh=permutation.emp_rwd_pval;

%% Loop through ROIs
max_t = 85;
srate = .144;
for i_roi= 1:length(ROIs);
figure('color',[1 1 1]);
plot([0 max_t+20],[0 0],'k-','LineWidth',4); hold on
plot([2/srate 2/srate],[-1 +1],'k--','LineWidth',4); hold on
fillsteplotk(beta_ts{i_roi}.rwd,4); hold on
for t=1:max_t
    if pmt_coh{i_roi}(t); plot(t,-.06,'ko','MarkerSize',12,'MarkerFaceColor','k'); end
end
ylim([-.07 .07]); 
xlim([1 max_t]);
xlabel('time [seconds]','FontSize',32,'FontWeight','bold');
title(ROIs_for_plot{i_roi},'FontSize',32,'FontWeight','normal')
ylabel('beta [a.u.]','FontSize',32,'FontWeight','bold');
set(gca,'YTick',[-.06:.03:.06]);
set(gca,'XTick',1:14:max_t-1)
set(gca,'XTickLabel',{'0','2','4','6','8','10'})
box('off')
set(gca,'FontSize',40,'LineWidth',4);
end