% Bang & Fleming (2018) Distinct encoding of decision confidence in human
% medial prefrontal cortex
%
% Confidence analysis of pgACC single-trial activity estimates
%
% Loads predictions from confidence model (see folder 'Model' for
% generation of model predictions)
%
% Statistical tests can be reproduced by applying one-sample t-tests to
% regression betas: e.g. [H,P,CI,STATS]= ttest(beta_conModel{i_roi},0)
%
% Reproduces panel B in Figure 4
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
load([dataScanDir,fs,'Tactivity_exclusions_1sec.mat']);

%% -----------------------------------------------------------------------
%% ANALYSIS

% Loop through ROIs
for i_roi = 1:length(ROIs)
% Loop through subjects
for i_sbj = 1:length(subjects)
    
    % Load ROI data
    file= [dataScanDir,fs,ROIs{i_roi},'_s',num2str(subjects(i_sbj)),'_1sec.mat'];
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
    
    %% parse
    concat = data.concat(include);
    conEmpirical = data.con(include);
    conModel = conmodel.Yint{i_sbj}(include)/10+.4;
    roi= roi_ts(include);
    
    %% binned ROI activity
    n_bins= 6;
    thresh1= quantile(roi(concat==1),n_bins-1); % confidence trials
    thresh2= quantile(roi,n_bins-1); % all trials
    clear roi_emp_cat roi_mod_cat;
    for i= 1:sum(include); roi_emp_cat(i)=sum(roi(i)>thresh1)+1; end;
    for i= 1:sum(include); roi_mod_cat(i)=sum(roi(i)>thresh2)+1; end;

    %% binned ROI activity
    for i= 1:length(thresh1)+1; roi_emp_con{i_roi}(i_sbj,i)= mean(conEmpirical(roi_emp_cat==i&concat==1)); end;
    for i= 1:length(thresh2)+1; roi_mod_con{i_roi}(i_sbj,i)= mean(conModel(roi_mod_cat==i)); end;
    
    %% trial-by-trial regression

    X= roi(concat==1)';
    Y= conEmpirical(concat==1)';
    [B,DEV,STATS]= glmfit(X,Y);
    beta_conEmpirical{i_roi}(i_sbj,1)= B(end);
    
    %% trial-by-trial regression
    X= roi';
    Y= conModel';
    [B,DEV,STATS]= glmfit(X,Y);
    beta_conModel{i_roi}(i_sbj,1)= B(end);

    end
end

%% -----------------------------------------------------------------------
%% FIGURES

%% Loop through ROIs
for i_roi= 1:length(ROIs);
figure('color',[1 1 1]);
plot([0 4],[0 0],'k-','LineWidth',4); hold on;
c_mean= nanmean(roi_emp_con{i_roi});
c_std= nanstd(roi_emp_con{i_roi})/sqrt(32);
for i= 1:numel(c_mean);
    plot([i i],[c_mean(i)-c_std(i) c_mean(i)+c_std(i)],'k-','LineWidth',4); hold on;
end
plot(c_mean,'ko','MarkerFaceColor',[.5 .5 .5],'MarkerSize',20,'LineWidth',4);
set(gca,'FontSize',24,'LineWidth',4);
set(gca,'Xtick',[1:n_bins]);
set(gca,'Ytick',[.76:.02:.82]);
xlim([0 n_bins+1]);
ylim([.76 .82]);
xlabel('pgACC beta bin','FontSize',28,'FontWeight','normal');
title(ROIs_for_plot{i_roi},'FontSize',28,'FontWeight','normal')
ylabel('empirical confidence','FontSize',28,'FontWeight','normal');
box('off')  
end

%% Loop through ROIs
for i_roi= 1:length(ROIs);
figure('color',[1 1 1]);
plot([0 7],[0 0],'k-','LineWidth',4); hold on;
c_mean= nanmean(roi_mod_con{i_roi});
c_std= nanstd(roi_mod_con{i_roi})/sqrt(32);
for i= 1:numel(c_mean);
    plot([i i],[c_mean(i)-c_std(i) c_mean(i)+c_std(i)],'k-','LineWidth',4); hold on;
end
plot(c_mean,'ko','MarkerFaceColor',[.5 .5 .5],'MarkerSize',20,'LineWidth',4);
set(gca,'FontSize',24,'LineWidth',4);
set(gca,'Xtick',[1:n_bins]);
set(gca,'Ytick',[.76:.02:.82]);
xlim([0 n_bins+1]);
ylim([.76 .82]);
xlabel('pgACC beta bin','FontSize',28,'FontWeight','normal');
title(ROIs_for_plot{i_roi},'FontSize',28,'FontWeight','normal')
ylabel('model confidence','FontSize',28,'FontWeight','normal');
box('off')  
end
