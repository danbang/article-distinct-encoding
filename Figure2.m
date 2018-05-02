% Bang & Fleming (2018) Distinct encoding of decision confidence in human
% medial prefrontal cortex
%
% Visualises behavioural summary statistics for pre-scan and scan sessions and
% performs associated regression analyses
%
% Reproduces panels A-C in Figure 2
%
% Statistical tests can be reproduced by applying one-sample t-tests to
% regression betas: e.g. [H,P,CI,STATS]= ttest(prescan.beta.con,0)
%
% Posterior predictives from HDDM fit are overlaid as dots (see folder
% 'HDDM' for HDDM pipeline)
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
dataPrescanDir = [repoBase,fs,'Data',fs,'Behaviour',fs,'Prescan'];
dataScanDir = [repoBase,fs,'Data',fs,'Behaviour',fs,'Scan'];
dataHDDMDir = [repoBase,fs,'HDDM'];

% Add custom scripts
addpath('Custom');

%% -----------------------------------------------------------------------
%% PRESCAN

%% LOOP THROUGH SUBJECTS
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
    
    % Cardinality vector
    temp = [sort(linspace(1,45,45),'descend')./45 sort(linspace(1,44,44),'ascend')./45 1];
    card = repmat(temp,1,4);
    
    %% SUMMARY STATISTICS FOR VISUALISATION 
    
    % trial indices
    sindx = 1:length(data.acc);
    sindx = sindx(include);
    % re-assign relevant variables
    c_acc = data.acc(sindx);
    c_con = data.con(sindx);
    c_rt1 = data.rt1(sindx)/1000;
    c_coh = data.cohz(sindx);
    c_del = abs(data.deltaz(sindx));
    u_coh = unique(c_coh);
    u_del = unique(c_del);    
    % loop through conditions
    j = 1; for l = 1:numel(u_del); prescan.acc1(i_sbj,l) = mean(c_acc(c_coh ==u_coh(j) & c_del==u_del(l))); end
    j = 2; for l = 1:numel(u_del); prescan.acc2(i_sbj,l) = mean(c_acc(c_coh==u_coh(j) & c_del==u_del(l))); end   
    j = 1; for l = 1:numel(u_del); prescan.con1(i_sbj,l) = mean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end
    j = 2; for l = 1:numel(u_del); prescan.con2(i_sbj,l) = mean(c_con(c_coh==u_coh(j) & c_del==u_del(l))); end   
    j = 1; for l = 1:numel(u_del); prescan.rt11(i_sbj,l) = mean(c_rt1(c_coh==u_coh(j) & c_del==u_del(l))); end
    j = 2; for l = 1:numel(u_del); prescan.rt12(i_sbj,l) = mean(c_rt1(c_coh==u_coh(j) & c_del==u_del(l))); end   

    %% REGRESSION MODELS FOR STATISTCAL ANALYSIS
    
    %% Choice accuracy
    % trial indices
    sindx = find(include);
    % load relevant data
    c_acc = data.acc(sindx);
    c_cohcat = zscore(data.cohcat(sindx)-1.5);
    c_delcat = zscore(data.deltacat(sindx)-2.5);
    % predictor variables
    X = [c_cohcat;
        c_delcat;
        c_cohcat.*c_delcat;
        ]';
    % outcome variable
    Y = (c_acc)';
    % run regression
    [B,DEV,STATS] = glmfit(X,Y,'binomial','link','logit');
    % log output
    prescan.beta.acc(i_sbj,:) = B(2:end)';
    
    %% Confidence
    % trial indices
    sindx = find(include);
    % load relevant data
    c_con = data.con(sindx).*10-4;
    c_cohcat = zscore(data.cohcat(sindx)-1.5);
    c_delcat = zscore(data.deltacat(sindx)-2.5);
    % predictor variables
    X = [c_cohcat;
        c_delcat;
        c_cohcat.*c_delcat;
        ]';
    % outcome variable
    Y = 7-c_con;
    % run regression
    [B,~,STATS] = mnrfit(X,Y,'model','ordinal','link','probit');
    % log output
    prescan.beta.con(i_sbj,:) = B(end-size(X,2)+1:end)';
    
    %% Confidence: FULL
    % trial indices
    sindx = find(include);
    % load relevant data
    c_con = data.con(sindx).*10-4;
    c_acc = data.acc(sindx)-.5;
    c_cohcat = zscore(data.cohcat(sindx)-1.5);
    c_delcat = zscore(data.deltacat(sindx)-2.5);
    c_rt1 = zscore(log(data.rt1(sindx)./1000));
    c_mar = zscore(data.markini(sindx).*10-4);
    c_cho = data.choice(sindx)-.5;
    c_mu = data.meanz(sindx)./360;
    c_card = zscore(card(round(c_mu*360)));
    % predictor variables
    X = [c_cohcat;
        c_delcat;
        c_cohcat.*c_delcat;
        c_rt1;
        c_acc;
        c_mar;
        c_cho;
        c_card;
        ]';
    % outcome variable
    Y = 7-c_con;
    % run regression
    [B,~,STATS] = mnrfit(X,Y,'model','ordinal','link','probit');
    % log output
    prescan.beta.conFull(i_sbj,:) = B(end-size(X,2)+1:end)';
    
    %% Reaction time
    % trial indices
    sindx = find(include);
    % load relevant data
    c_rt1 = zscore(log(data.rt1(sindx)./1000));
    c_cohcat = zscore(data.cohcat(sindx)-1.5);
    c_delcat = zscore(data.deltacat(sindx)-2.5);
    % predictor variables
    X = [c_cohcat;
        c_delcat;
        c_cohcat.*c_delcat;
        ]';
    % outcome variable
    Y = (c_rt1)';
    % run regression
    [B,DEV,STATS] = glmfit(X,Y);
    % log output
    prescan.beta.rt1(i_sbj,:) = B(2:end)';
    
end

%% Load HDDM posterior predictives
file = [dataHDDMDir,fs,'hddm_posterior_predictives_prescan.mat'];
load(file);

%% -----------------------------------------------------------------------
%% FIGURES

%% Choice accuracy
figure('color',[1 1 1]);
fillsteplotblue(prescan.acc1,4); hold on
fillsteplotred(prescan.acc2,4); hold on
plot(mean(LC_acc),'ko','MarkerFaceColor',[51/255 153/255 255/255],'MarkerSize',16,'LineWidth',2); hold on;
plot(mean(HC_acc),'ko','MarkerFaceColor',[255/255 51/255 51/255],'MarkerSize',16,'LineWidth',2); hold on;
ylim([.5 1]); xlim([0 5]);
set(gca,'YTick',.5:.1:1);
set(gca,'XTiCk',1:4);
set(gca,'XTickLabel',{'1','2','3','4'});
set(gca,'FontSize',46,'LineWidth',4)
legend('boxoff')
ylabel('percent correct','FontSize',50,'FontWeight','normal')
xlabel('distance [4=easiest]','FontSize',50,'FontWeight','normal')
title('pre-scan','FontWeight','normal');
box('off')

%% Reaction time
figure('color',[1 1 1]);
fillsteplotblue(prescan.rt11,4,'-'); hold on
fillsteplotred(prescan.rt12,4,'-'); hold on
plot(mean(LC_rt1),'ko','MarkerFaceColor',[51/255 153/255 255/255],'MarkerSize',16,'LineWidth',2); hold on;
plot(mean(HC_rt1),'ko','MarkerFaceColor',[255/255 51/255 51/255],'MarkerSize',16,'LineWidth',2); hold on;
ylim([.9 2.1]); xlim([0 5]);
set(gca,'YTick',1:.25:2);
set(gca,'XTiCk',1:4);
set(gca,'XTickLabel',{'1','2','3','4'});
set(gca,'FontSize',46,'LineWidth',4)
legend('boxoff')
ylabel('seconds','FontSize',50,'FontWeight','normal')
xlabel('distance [4=easiest]','FontSize',50,'FontWeight','normal')
box('off')
title('pre-scan','FontWeight','normal');

%% Confidence
figure('color',[1 1 1]);
fillsteplotblue(prescan.con1,4,'-'); hold on
fillsteplotred(prescan.con2,4,'-'); hold on
ylim([.5 1]); xlim([0 5]);
set(gca,'YTick',.5:.1:1);
set(gca,'XTiCk',1:4);
set(gca,'XTickLabel',{'1','2','3','4'});
set(gca,'FontSize',46,'LineWidth',4)
legend('boxoff')
ylabel('confidence','FontSize',50,'FontWeight','normal')
xlabel('distance [4=easiest]','FontSize',50,'FontWeight','normal')
box('off')
title('pre-scan','FontWeight','normal');

%% -----------------------------------------------------------------------
%% SCAN

% Subjects
subjects= [1:9 11:13 15:20 22:35]; 

%% LOOP THROUGH SUBJECTS
for i_sbj= 1:length(subjects)
    
    %% HOUSE KEEEPING

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

    % Distance categories
    sindx= 1:length(data.acc);
    c_del= abs(data.deltaz(sindx));
    u_del= unique(c_del);
    for t= 1:length(sindx); data.deltacat(sindx(t))= sum(c_del(t)>=u_del); end
    
    % Cardinality vector
    temp = [sort(linspace(1,45,45),'descend')./45 sort(linspace(1,44,44),'ascend')./45 1];
    card = repmat(temp,1,4);
    
    %% SUMMARY STATISTICS FOR VISUALISATION 
    
    % trial indices
    sindx= 1:length(data.acc);
    sindx= sindx(include);
    % re-assign relevant variables
    c_acc= data.acc(sindx);
    c_con= data.con(sindx);
    c_concat= data.concat(sindx);
    c_rt1= data.rt1(sindx)/1000;
    c_coh = data.cohz(sindx);
    c_del = abs(data.deltaz(sindx));
    u_coh = unique(c_coh);
    u_del = unique(c_del);    
    % loop through conditions
    j= 1; for l= 1:numel(u_del); scan.acc1(i_sbj,l)= mean(c_acc(c_coh==u_coh(j) & c_del==u_del(l))); end
    j= 2; for l= 1:numel(u_del); scan.acc2(i_sbj,l)= mean(c_acc(c_coh==u_coh(j) & c_del==u_del(l))); end   
    j= 1; for l= 1:numel(u_del); scan.con1(i_sbj,l)= mean(c_con(c_coh==u_coh(j) & c_del==u_del(l) & c_concat)); end
    j= 2; for l= 1:numel(u_del); scan.con2(i_sbj,l)= mean(c_con(c_coh==u_coh(j) & c_del==u_del(l) & c_concat)); end   
    j= 1; for l= 1:numel(u_del); scan.rt11(i_sbj,l)= mean(c_rt1(c_coh==u_coh(j) & c_del==u_del(l))); end
    j= 2; for l= 1:numel(u_del); scan.rt12(i_sbj,l)= mean(c_rt1(c_coh==u_coh(j) & c_del==u_del(l))); end   

    %% REGRESSION MODELS FOR STATISTCAL ANALYSIS
    
    %% Choice accuracy
    % trial indices
    sindx= find(include);
    % load relevant data
    c_acc= data.acc(sindx);
    c_cohcat= zscore(data.cohcat(sindx)-1.5);
    c_delcat= zscore(data.deltacat(sindx)-2.5);
    % predictor variables
    X= [c_cohcat;
        c_delcat;
        c_cohcat.*c_delcat;
        ]';
    % outcome variable
    Y= (c_acc)';
    % run regression
    [B,DEV,STATS]= glmfit(X,Y,'binomial','link','logit');
    % log output
    scan.beta.acc(i_sbj,:)= B(2:end)';
    
    %% Confidence
    % trial indices
    sindx = find(include&data.concat);
    % load relevant data
    c_con    = data.con(sindx).*10-4;
    c_cohcat = zscore(data.cohcat(sindx)-1.5);
    c_delcat = zscore(data.deltacat(sindx)-1.5);
    % predictor variables
    X = [c_cohcat;
         c_delcat;
         c_cohcat.*c_delcat;
         ]';
    % outcome variable
    Y = 7-c_con;
    % run regression
    [B,~,STATS] = mnrfit(X,Y,'model','ordinal','link','probit');
    % log output
    scan.beta.con(i_sbj,:) = B(end-size(X,2)+1:end)';
    
    %% Confidence: FULL
    % trial indices
    sindx = find(include&data.concat);
    % load relevant data
    c_con    = data.con(sindx).*10-4;
    c_acc    = data.acc(sindx)-.5;
    c_cohcat = zscore(data.cohcat(sindx)-1.5);
    c_delcat = zscore(data.deltacat(sindx)-1.5);
    c_rt1    = zscore(log(data.rt1(sindx)./1000));
    c_mar    = zscore(data.markini(sindx).*10-4);
    c_cho    = data.choice(sindx)-.5;
    c_mu     = data.meanz(sindx)./360;
    c_card   = zscore(card(round(c_mu*360)));
    % predictor variables
    X = [c_cohcat;
         c_delcat;
         c_cohcat.*c_delcat;
         c_rt1;
         c_acc;
         c_mar;
         c_cho;
         c_card;
         ]';
    % outcome variable
    Y = 7-c_con;
    % run regression
    [B,~,STATS] = mnrfit(X,Y,'model','ordinal','link','probit');
    % log output
    scan.beta.conFull(i_sbj,:) = B(end-size(X,2)+1:end)';
    
    %% Choice accuracy
    % trial indices
    sindx = find(include);
    % load relevant data
    c_rt1    = zscore(log(data.rt1(sindx)./1000));
    c_cohcat = zscore(data.cohcat(sindx)-1.5);
    c_delcat = zscore(data.deltacat(sindx)-1.5);
    % predictor variables
    X = [c_cohcat;
         c_delcat;
         c_cohcat.*c_delcat;
         ]';
    % outcome variable
    Y = (c_rt1)';
    % run regression
    [B,DEV,STATS] = glmfit(X,Y);
    % log output
    scan.beta.rt1(i_sbj,:) = B(2:end)';
    
end

%% Load HDDM posterior predictives
file = [dataHDDMDir,fs,'hddm_posterior_predictives_scan.mat'];
load(file);

%% -----------------------------------------------------------------------
%% FIGURES

%% Choice accuracy
figure('color',[1 1 1]);
fillsteplotblue(scan.acc1,4); hold on
fillsteplotred(scan.acc2,4); hold on
plot(mean(LC_acc),'ko','MarkerFaceColor',[51/255 153/255 255/255],'MarkerSize',16,'LineWidth',2); hold on;
plot(mean(HC_acc),'ko','MarkerFaceColor',[255/255 51/255 51/255],'MarkerSize',16,'LineWidth',2); hold on;
ylim([.5 1]); xlim([0 3]);
set(gca,'YTick',.5:.1:1);
set(gca,'XTiCk',1:2);
set(gca,'XTickLabel',{'1','2'});
set(gca,'FontSize',46,'LineWidth',4)
legend('boxoff')
ylabel('percent correct','FontSize',50,'FontWeight','normal')
xlabel('distance [2=easiest]','FontSize',50,'FontWeight','normal')
title('scan','FontWeight','normal');
box('off')

%% Reaction time
figure('color',[1 1 1]);
fillsteplotblue(scan.rt11,4,'-'); hold on
fillsteplotred(scan.rt12,4,'-'); hold on
plot(mean(LC_rt1),'ko','MarkerFaceColor',[51/255 153/255 255/255],'MarkerSize',16,'LineWidth',2); hold on;
plot(mean(HC_rt1),'ko','MarkerFaceColor',[255/255 51/255 51/255],'MarkerSize',16,'LineWidth',2); hold on;
ylim([.9 2.1]); xlim([0 3]);
set(gca,'YTick',1:.25:2);
set(gca,'XTiCk',1:2);
set(gca,'XTickLabel',{'1','2'});
set(gca,'FontSize',46,'LineWidth',4)
legend('boxoff')
ylabel('seconds','FontSize',50,'FontWeight','normal')
xlabel('distance [2=easiest]','FontSize',50,'FontWeight','normal')
box('off')
title('scan','FontWeight','normal');

%% Confidence
figure('color',[1 1 1]);
fillsteplotblue(scan.con1,4,'-'); hold on
fillsteplotred(scan.con2,4,'-'); hold on
ylim([.5 1]); xlim([0 3]);
set(gca,'YTick',.5:.1:1);
set(gca,'XTiCk',1:2);
set(gca,'XTickLabel',{'1','2'});
set(gca,'FontSize',46,'LineWidth',4)
legend('boxoff')
ylabel('confidence','FontSize',50,'FontWeight','normal')
xlabel('distance [2=easiest]','FontSize',50,'FontWeight','normal')
box('off')
title('scan','FontWeight','normal');