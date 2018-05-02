% Bang & Fleming (2018) Distinct encoding of decision confidence in human
% medial prefrontal cortex
%
% Confidence analysis of pgACC activity timecourses
%
% Loads predictions from confidence model (see folder 'Model' for
% generation of model predictions)
%
% Loads permutation-testing statistics for indicating significance (see folder
% 'Permutation' for permutation procedure)
%
% Reproduces panel C in Figure 4
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
dataScanDir = [repoBase,fs,'Data',fs,'Scan',fs,'Timecourse'];
dataModelDir = [repoBase,fs,'Model'];
dataPermutationDir = [repoBase,fs,'Permutation'];

% Add custom scripts
addpath('Custom');

% Load out-of-sample predictions
load([dataModelDir,fs,'conmodel_outOFsample_predictions.mat']);

%% -----------------------------------------------------------------------
%% ANALYSIS

% Loop through ROIs
for i_roi = 1:length(ROIs)
% Loop through subjects
for i_sbj = 1:length(subjects)
    
    % Load ROI data
    file= [dataScanDir,fs,ROIs{i_roi},'_s',num2str(subjects(i_sbj)),'_stimulus.mat'];
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
    rt1= (data.rt1./1000);
    centre= mean(rt1);
    stdval= std(rt1)*2.5;
    include= (rt1>(centre-stdval))&(rt1<(centre+stdval));
    
    % Include trials where final time-point estimate is ~NaN
    for i= 1:size(roi_ts,1); if isnan(roi_ts(i,end)); include(i)=0; end; end;
    
    % UP-SAMPLED GLM: MODEL
    include= include;
    roi_Zts = zscore(roi_ts(include,:));
    con= conmodel.Yint{i_sbj}(include);
    t= 0;
    for j= 1:size(roi_ts,2)
        t= t+1;
        x= con';
        y= roi_Zts(:,j);
        beta= glmfit(x,y);
        beta_ts{i_roi}.conModel(i_sbj,t)= beta(end);
    end
    
    % UP-SAMPLED GLM: EMPIRICAL
    include= include & data.concat==1;
    roi_Zts = zscore(roi_ts(include,:));
    con= data.con(include).*10-4;
    t= 0;
    for j= 1:size(roi_ts,2)
        t= t+1;
        x= con';
        y= roi_Zts(:,j);
        beta= glmfit(x,y);
        beta_ts{i_roi}.conEmpirical(i_sbj,t)= beta(end);
    end

    end
end

%% -----------------------------------------------------------------------
%% FIGURES

%% Load permutation output
load([dataPermutationDir,fs,'timecourse_permutation_confidence.mat']);
pmt_conModel=permutation.emp_conModel_pval;
pmt_conEmpirical=permutation.emp_conEmpirical_pval;

%% Loop through ROIs
max_t = 85;
srate = .144;
for i_roi= 1:length(ROIs);
figure('color',[1 1 1]);
plot([0 max_t+20],[0 0],'k-','LineWidth',4); hold on
plot([1/srate 1/srate],[-1 +1],'k--','LineWidth',4); hold on
plot([2/srate 2/srate],[-1 +1],'k--','LineWidth',4); hold on
fillsteplotk(beta_ts{i_roi}.conModel,4); hold on
for t=1:max_t
    if pmt_conModel{i_roi}(t); plot(t,-.155,'o','MarkerSize',10,'Color','k','MarkerFaceColor','k'); end
end
ylim([-.17 .17]); 
xlim([1 max_t]);
xlabel('time [seconds]','FontSize',28,'FontWeight','bold');
title(ROIs_for_plot{i_roi},'FontSize',28,'FontWeight','normal')
ylabel('model: beta [a.u.]','FontSize',28,'FontWeight','bold');
set(gca,'YTick',[-.16:.04:.16]);
set(gca,'XTick',1:14:max_t-1)
set(gca,'XTickLabel',{'0','2','4','6','8','10'})
box('off')
set(gca,'FontSize',24,'LineWidth',4);
end

%% Loop through ROIs
max_t = 85;
srate = .144;
for i_roi= 1:length(ROIs);
figure('color',[1 1 1]);
plot([0 max_t+20],[0 0],'k-','LineWidth',4); hold on
plot([1/srate 1/srate],[-1 +1],'k--','LineWidth',4); hold on
plot([2/srate 2/srate],[-1 +1],'k--','LineWidth',4); hold on
fillsteplotk(beta_ts{i_roi}.conEmpirical,4); hold on
for t=1:max_t
    if pmt_conEmpirical{i_roi}(t); plot(t,-.155,'o','MarkerSize',10,'Color','k','MarkerFaceColor','k'); end
end
ylim([-.17 .17]); 
xlim([1 max_t]);
xlabel('time [seconds]','FontSize',28,'FontWeight','normal');
title(ROIs_for_plot{i_roi},'FontSize',28,'FontWeight','normal')
ylabel('empirical: beta [a.u.]','FontSize',28,'FontWeight','normal');
set(gca,'YTick',[-.16:.04:.16]);
set(gca,'XTick',1:14:max_t-1)
set(gca,'XTickLabel',{'0','2','4','6','8','10'})
box('off')
set(gca,'FontSize',24,'LineWidth',4);
end