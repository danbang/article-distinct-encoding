% Bang & Fleming (2018) Distinct encoding of decision confidence in human
% medial prefrontal cortex
%
% Generates out-of-sample predictions about confidence in scan session using 
% multinomial ordinal regression model fitted to data from pre-scan session
%
% Visualises empirically observed and predicted summary statistics (note
% that visualisation for scan session only uses confidence trials)
%
% Dan Bang danbang.db@gmail.com 2018

%% -----------------------------------------------------------------------
%% PREPARATION

% Fresh memory
clear; close all;

% Subjects
subjects = [1:9 11:13 15:20 22:35];

% Paths [change 'repoBase' according to local setup]
fs = filesep;
repoBase = [getDropbox(1),fs,'Ego',fs,'Matlab',fs,'ucl',fs,'sensory_vs_decision',fs,'Repository'];
dataPrescanDir = [repoBase,fs,'Data',fs,'Behaviour',fs,'Prescan'];
dataScanDir = [repoBase,fs,'Data',fs,'Behaviour',fs,'Scan'];
dataHDDMDir = [repoBase,fs,'HDDM'];

% Add custom scripts
addpath([repoBase,fs,'Custom']);

%% -----------------------------------------------------------------------
%% FIT TO PRE-SCAN DATA

%% loop through subjects
for i_sbj = 1:length(subjects)
        
    %% HOUSE KEEEPING
    
    % Load file
    file = [dataPrescanDir,fs,'s',num2str(subjects(i_sbj)),'_task.mat'];
    load(file);
    
    % Include trials based on deviation from grand mean
    rt1 = log(data.rt1./1000);
    centre = mean(rt1);
    stdval = std(rt1)*2.5;
    include = (rt1>(centre-stdval))&(rt1<(centre+stdval));

    % Distance categories
    sindx = 1:length(data.acc);
    c_del = abs(data.deltaz(sindx));
    u_del = unique(c_del);
    for t = 1:length(sindx); data.deltacat(sindx(t)) = sum(c_del(t)>=u_del); end
        
    % regression model
    % variables
    c_con = (data.con(include).*10-4);
    c_coh = data.cohz(include);
    c_del = abs(data.deltaz(include))./360;
    c_rt1 = log(data.rt1(include)./1000);
    c_acc = data.acc(include)-.5;
    % predictor matrix
    X = [c_rt1;
         c_acc;
         c_coh;
         c_del;
         c_coh.*c_del;
         c_coh.*c_acc;
         c_del.*c_acc;
         c_coh.*c_del.*c_acc;
         ]';
    % outcome
    Y = (7-c_con)';
    % run regression
    [B,~,STATS] = mnrfit(X,Y,'model','ordinal','link','logit');
    % save betas
    fit.betas{i_sbj} = B;
    % generate prediction
    Yhat          = mnrval(fit.betas{i_sbj},X,'model','ordinal','link','logit');
    Yint          = sum(repmat([1:6],size(Yhat,1),1).*Yhat,2);
    % zero-pad
    zeropad= NaN(1,length(include)); j=0;
    for i= 1:length(zeropad); if include(i); j=j+1; zeropad(i)= Yint(j); end; end
    Yint= zeropad;
    
    % MODEL PREDICTION
    % wrong trials
    sindx = data.acc==0 & include;
    % load relevant data
    c_con = (7-Yint(sindx))/10+.4;
    c_coh = data.cohz(sindx);
    c_del = abs(data.deltaz(sindx));
    u_coh = unique(c_coh);
    u_del = unique(c_del);    
    % loop through conditions
    j = 1; for l = 1:numel(u_del); prescan.sim_con1w(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
    j = 2; for l = 1:numel(u_del); prescan.sim_con2w(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
    % corret trials
    sindx = data.acc==1 & include;
    % load relevant data
    c_con = (7-Yint(sindx))/10+.4;
    c_coh = data.cohz(sindx);
    c_del = abs(data.deltaz(sindx));
    u_coh = unique(c_coh);
    u_del = unique(c_del);    
    % loop through conditions
    j = 1; for l = 1:numel(u_del); prescan.sim_con1c(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
    j = 2; for l = 1:numel(u_del); prescan.sim_con2c(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
   
    % EMPIRICAL OBSERVATION
    % wrong trials
    sindx = data.acc==0 & include;
    % load relevant data
    c_con = data.con(sindx);
    c_coh = data.cohz(sindx);
    c_del = abs(data.deltaz(sindx));
    u_coh = unique(c_coh);
    u_del = unique(c_del);    
    % loop through conditions
    j = 1; for l = 1:numel(u_del); prescan.emp_con1w(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
    j = 2; for l = 1:numel(u_del); prescan.emp_con2w(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
    % corret trials
    sindx = data.acc==1 & include;
    % load relevant data
    c_con = data.con(sindx);
    c_coh = data.cohz(sindx);
    c_del = abs(data.deltaz(sindx));
    u_coh = unique(c_coh);
    u_del = unique(c_del);    
    % loop through conditions
    j = 1; for l = 1:numel(u_del); prescan.emp_con1c(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
    j = 2; for l = 1:numel(u_del); prescan.emp_con2c(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
    
end

%% PLOT DATA
figure('color',[1 1 1]);
subplot(2,2,1)
fillsteplotblue(prescan.emp_con1c,2,'-'); hold on
fillsteplotred(prescan.emp_con2c,2,'-'); hold on
fillsteplotblue(prescan.emp_con1w,2,'--'); hold on
fillsteplotnanred(prescan.emp_con2w,2,'--'); hold on
ylim([.5 1]); xlim([0 5]);
set(gca,'YTick',.5:.1:1);
set(gca,'XTiCk',1:4);
set(gca,'XTickLabel',{'1','2','3','4'});
set(gca,'FontSize',16,'LineWidth',2)
legend('boxoff')
ylabel('confidence','FontSize',20,'FontWeight','normal')
xlabel('distance [4=easiest]','FontSize',20,'FontWeight','normal')
title('pre-scan: empirical','FontWeight','normal');
box('off')
subplot(2,2,2)
fillsteplotblue(prescan.sim_con1c,2,'-'); hold on
fillsteplotred(prescan.sim_con2c,2,'-'); hold on
fillsteplotblue(prescan.sim_con1w,2,'--'); hold on
fillsteplotnanred(prescan.sim_con2w,2,'--'); hold on
ylim([.5 1]); xlim([0 5]);
set(gca,'YTick',.5:.1:1);
set(gca,'XTiCk',1:4);
set(gca,'XTickLabel',{'1','2','3','4'});
set(gca,'FontSize',16,'LineWidth',2)
legend('boxoff')
ylabel('confidence','FontSize',20,'FontWeight','normal')
xlabel('distance [4=easiest]','FontSize',20,'FontWeight','normal')
title('pre-scan: model','FontWeight','normal');
box('off')

%% -----------------------------------------------------------------------
%% PREDICTION

%% loop through subjects
for i_sbj = 1:length(subjects)
             
    % load behavioural data
    
    % loop through scan runs
    for i_blk = 1:5     
        
    % load file
    file = [dataScanDir,fs,'s',num2str(subjects(i_sbj)),'_task_b',num2str(i_blk),'.mat'];
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
    
    % generate out-of-sample predictions
    % trial indices
    sindx = include;
    % variables
    c_con = (data.con(sindx).*10-4);
    c_coh = data.cohz(sindx);
    c_del = abs(data.deltaz(sindx))./360;
    c_rt1 = log(data.rt1(sindx)./1000);
    c_acc = data.acc(sindx)-.5;
    % predictor matrix
    X = [c_rt1;
         c_acc;
         c_coh;
         c_del;
         c_coh.*c_del;
         c_coh.*c_acc;
         c_del.*c_acc;
         c_coh.*c_del.*c_acc;
         ]';
    % generate prediction
    Yhat          = mnrval(fit.betas{i_sbj},X,'model','ordinal','link','logit');
    Yint          = sum(repmat([1:6],size(Yhat,1),1).*Yhat,2);
    % zero-pad
    zeropad= NaN(1,length(include)); j=0;
    for i= 1:length(zeropad); if include(i); j=j+1; zeropad(i)= Yint(j); end; end
    Yint= zeropad;
    % log prediction
    conmodel.Yint{i_sbj}= 7-Yint;
    
    % MODEL PREDICTION
    % wrong trials
    sindx = data.acc==0 & include & data.concat==1;
    % load relevant data
    c_con = (7-Yint(sindx))/10+.4;
    c_coh = data.cohz(sindx);
    c_del = abs(data.deltaz(sindx));
    u_coh = unique(c_coh);
    u_del = unique(c_del);    
    % loop through conditions
    j = 1; for l = 1:numel(u_del); scan.sim_con1w(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
    j = 2; for l = 1:numel(u_del); scan.sim_con2w(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
    % corret trials
    sindx = data.acc==1 & include & data.concat==1;
    % load relevant data
    c_con = (7-Yint(sindx))/10+.4;
    c_coh = data.cohz(sindx);
    c_del = abs(data.deltaz(sindx));
    u_coh = unique(c_coh);
    u_del = unique(c_del);    
    % loop through conditions
    j = 1; for l = 1:numel(u_del); scan.sim_con1c(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
    j = 2; for l = 1:numel(u_del); scan.sim_con2c(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
   
    % EMPIRICAL OBSERVATION
    % wrong trials
    sindx = data.acc==0 & include & data.concat==1;
    % load relevant data
    c_con = data.con(sindx);
    c_coh = data.cohz(sindx);
    c_del = abs(data.deltaz(sindx));
    u_coh = unique(c_coh);
    u_del = unique(c_del);    
    % loop through conditions
    j = 1; for l = 1:numel(u_del); scan.emp_con1w(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
    j = 2; for l = 1:numel(u_del); scan.emp_con2w(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
    % corret trials
    sindx = data.acc==1 & include & data.concat==1;
    % load relevant data
    c_con = data.con(sindx);
    c_coh = data.cohz(sindx);
    c_del = abs(data.deltaz(sindx));
    u_coh = unique(c_coh);
    u_del = unique(c_del);    
    % loop through conditions
    j = 1; for l = 1:numel(u_del); scan.emp_con1c(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
    j = 2; for l = 1:numel(u_del); scan.emp_con2c(i_sbj,l) = nanmean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
      
end

%% PLOT DATA
subplot(2,2,3)
fillsteplotblue(scan.emp_con1c,2,'-'); hold on
fillsteplotred(scan.emp_con2c,2,'-'); hold on
fillsteplotblue(scan.emp_con1w,2,'--'); hold on
fillsteplotnanred(scan.emp_con2w,2,'--'); hold on
ylim([.5 1]); xlim([0 3]);
set(gca,'YTick',.5:.1:1);
set(gca,'XTiCk',1:2);
set(gca,'XTickLabel',{'1','2'});
set(gca,'FontSize',16,'LineWidth',2)
legend('boxoff')
ylabel('confidence','FontSize',20,'FontWeight','normal')
xlabel('distance [2=easiest]','FontSize',20,'FontWeight','normal')
title('scan: empirical','FontWeight','normal');
box('off')
subplot(2,2,4)
fillsteplotblue(scan.sim_con1c,2,'-'); hold on
fillsteplotred(scan.sim_con2c,2,'-'); hold on
fillsteplotblue(scan.sim_con1w,2,'--'); hold on
fillsteplotnanred(scan.sim_con2w,2,'--'); hold on
ylim([.5 1]); xlim([0 3]);
set(gca,'YTick',.5:.1:1);
set(gca,'XTiCk',1:2);
set(gca,'XTickLabel',{'1','2'});
set(gca,'FontSize',16,'LineWidth',2)
legend('boxoff')
ylabel('confidence','FontSize',20,'FontWeight','normal')
xlabel('distance [4=easiest]','FontSize',20,'FontWeight','normal')
title('scan: model','FontWeight','normal');
box('off')

%% log prediction
save('conmodel_outOFsample_predictions.mat','conmodel');