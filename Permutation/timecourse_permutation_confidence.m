% Bang & Fleming (2018) Distinct encoding of decision confidence in human
% medial prefrontal cortex
%
% Permutation test for confidence analysis of ROI activity timecourses
%
% Loads predictions from confidence model (see folder 'Model' for
% generation of model predictions)
%
% ROIs: pgACC (leave-one-out), MT+ (localiser), IPS (leave-one-out),
% ventral striatum (anatomical), pre-SMA (leave-one-out), rlPFC (anatomical)
%
% Dan Bang danbang.db@gmail.com 2018

%% -----------------------------------------------------------------------
%% PREPARATION

% Fresh memory
clear; close all;

% Number of permutation
n_sim = 1e3+1;

% Subjects
subjects = [1:9 11:13 15:20 22:35];

% Paths [change 'repoBase' according to local setup]
fs = filesep;
repoBase = [getDropbox(1),fs,'Ego',fs,'Matlab',fs,'ucl',fs,'sensory_vs_decision',fs,'Repository'];
dataBehaviourDir = [repoBase,fs,'Data',fs,'Behaviour',fs,'Scan'];
dataScanDir = [repoBase,fs,'Data',fs,'Scan',fs,'Timecourse'];
dataModelDir = [repoBase,fs,'Model'];

% ROIs
ROIs= {'pgacc','V5','IPS','vstriatum','preSMA','rlPFC'};
ROIs_for_plot= {'pgACC','MT+','IPS','ventral striatum','pre-SMA','rlPFC'};
permutation.ROI_order= ROIs;

% Load out-of-sample predictions
load([dataModelDir,fs,'conmodel_outOFsample_predictions.mat']);

%% -----------------------------------------------------------------------
%% PERMUTATION PROCEDURE

%% loop through ROIs
for i_roi = 1:length(ROIs)

%% customer service
fprintf(['===== Performing permutation for ',ROIs_for_plot{i_roi},' ===== \n']);        

%% loop through permutations
for i_sim = 1:n_sim

%% loop through subjects
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
    rt1= log(data.rt1./1000);
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
    rindx = randperm(size(roi_Zts,1)); % scramble trials (not time points within trial)
    for j= 1:size(roi_ts,2)
        t= t+1;
        x= con';
        y= roi_Zts(rindx,j);
        beta= glmfit(x,y);
        permutation.sim_conModel_beta{i_roi}(i_sbj,t)= beta(end);
    end
    
    % UP-SAMPLED GLM: EMPIRICAL
    include= include & data.concat==1;
    roi_Zts = zscore(roi_ts(include,:));
    con= data.con(include).*10-4;
    t= 0;
    rindx = randperm(size(roi_Zts,1)); % scramble trials (not time points within trial)
    for j= 1:size(roi_ts,2)
        t= t+1;
        x= con';
        y= roi_Zts(rindx,j);
        beta= glmfit(x,y);
        permutation.sim_conEmpirical_beta{i_roi}(i_sbj,t)= beta(end);
    end
    
end

%% statistical test against zero
[h p ci stats] = ttest(permutation.sim_conModel_beta{i_roi});
permutation.sim_conModel_tval{i_roi}(i_sim,:) = stats.tstat;
[h p ci stats] = ttest(permutation.sim_conEmpirical_beta{i_roi});
permutation.sim_conEmpirical_tval{i_roi}(i_sim,:) = stats.tstat;

%% mean permuted beta
permutation.sim_conModel_beta_mean{i_roi}(i_sim,:) = mean(permutation.sim_conModel_beta{i_roi});
permutation.sim_conEmpirical_beta_mean{i_roi}(i_sim,:) = mean(permutation.sim_conEmpirical_beta{i_roi});

end

end

%% -----------------------------------------------------------------------
%% COMPARE AGAINST EMPIRICAL RESULTS

%% loop through ROIs
for i_roi = 1:length(ROIs)

%% loop through subjects
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
    rt1= log(data.rt1./1000);
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
        permutation.emp_conModel_beta{i_roi}(i_sbj,t)= beta(end);
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
        permutation.emp_conEmpirical_beta{i_roi}(i_sbj,t)= beta(end);
    end
        
end

%% statistical tests against zero
[h p ci stats] = ttest(permutation.emp_conModel_beta{i_roi});
permutation.emp_conModel_tval{i_roi} = stats.tstat;
[h p ci stats] = ttest(permutation.emp_conEmpirical_beta{i_roi});
permutation.emp_conEmpirical_tval{i_roi} = stats.tstat;

%% statistical testing against permutation -- thresholds
permutation.sim_conModel_LB{i_roi}=quantile(permutation.sim_conModel_tval{i_roi},.025);
permutation.sim_conModel_UB{i_roi}=quantile(permutation.sim_conModel_tval{i_roi},.975);
permutation.sim_conEmpirical_LB{i_roi}=quantile(permutation.sim_conEmpirical_tval{i_roi},.025);
permutation.sim_conEmpirical_UB{i_roi}=quantile(permutation.sim_conEmpirical_tval{i_roi},.975);

%% statistical testing against permutation -- comparison
permutation.emp_conModel_pval{i_roi}=(permutation.emp_conModel_tval{i_roi}<permutation.sim_conModel_LB{i_roi})|(permutation.emp_conModel_tval{i_roi}>permutation.sim_conModel_UB{i_roi});
permutation.emp_conEmpirical_pval{i_roi}=(permutation.emp_conEmpirical_tval{i_roi}<permutation.sim_conEmpirical_LB{i_roi})|(permutation.emp_conEmpirical_tval{i_roi}>permutation.sim_conEmpirical_UB{i_roi});

%% statistical testing against permutation -- comparison
permutation.emp_conModel_bdev{i_roi}=mean(permutation.emp_conModel_beta{i_roi})-mean(permutation.sim_conModel_beta_mean{i_roi});
permutation.emp_conEmpirical_bdev{i_roi}=mean(permutation.emp_conEmpirical_beta{i_roi})-mean(permutation.sim_conEmpirical_beta_mean{i_roi});

end

%% Save output
save('timecourse_permutation_confidence','permutation');