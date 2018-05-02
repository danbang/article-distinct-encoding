% Bang & Fleming (2018) Distinct encoding of decision confidence in human
% medial prefrontal cortex
%
% Export data as .csv file for HDDM
%
% PRE-SCAN CONDITION NAMING CONVENTION
% #1 Condition: coherence 2, distance 1
% #2 Condition: coherence 2, distance 2
% #3 Condition: coherence 2, distance 3
% #4 Condition: coherence 2, distance 4
% #5 Condition: coherence 1, distance 1
% #6 Condition: coherence 1, distance 2
% #7 Condition: coherence 1, distance 3
% #8 Condition: coherence 1, distance 4
% 
% SCAN CONDITION NAMING CONVENTION
% #1 Condition: coherence 2, distance 1
% #2 Condition: coherence 2, distance 2
% #5 Condition: coherence 1, distance 1
% #6 Condition: coherence 1, distance 2
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

% Add custom scripts
addpath('Custom');

%% -----------------------------------------------------------------------
%% PRE-SCAN

% Initialise variables using HDDM naming convention
subj_idx= [];
response= [];
rt= [];
condition= [];

%% LOOP THROUGH SUBJECTS
for i_sbj = 1:length(subjects)
    
    %% HOUSE KEEEPING
    
    % Load file
    file = [dataPrescanDir,fs,'s',num2str(subjects(i_sbj)),'_task.mat'];
    load(file);
    
    % add stimulus time to RT1 and transform to seconds
    data.rt1= (data.rt1./1000)+1;
    
    % Include trials based on deviation from grand mean
    rt1 = log(data.rt1./1000);
    centre = mean(rt1);
    stdval = std(rt1)*2.5;
    include = (rt1>(centre-stdval))&(rt1<(centre+stdval));

    % delta categories
    sindx = 1:length(data.acc);
    c_del = abs(data.deltaz(sindx));
    u_del = unique(c_del);
    for t = 1:length(sindx); data.bodcat(sindx(t)) = sum(c_del(t)>=u_del); end
    
    % turn into conditions
    clear tmp
    for i = 1:length(data.cohcat);
        if data.cohcat(i)==2 && data.bodcat(i)==1
            tmp(i)=1;
        elseif data.cohcat(i)==2 && data.bodcat(i)==2
            tmp(i)=2;
        elseif data.cohcat(i)==2 && data.bodcat(i)==3
            tmp(i)=3;
        elseif data.cohcat(i)==2 && data.bodcat(i)==4
            tmp(i)=4;
        elseif data.cohcat(i)==1 && data.bodcat(i)==1
            tmp(i)=5;
        elseif data.cohcat(i)==1 && data.bodcat(i)==2
            tmp(i)=6;
        elseif data.cohcat(i)==1 && data.bodcat(i)==3
            tmp(i)=7;
        elseif data.cohcat(i)==1 && data.bodcat(i)==4
            tmp(i)=8;
        end
    end

    
    % concatenate
    subj_idx= [subj_idx; ones(sum(include),1).*i_sbj];
    response= [response; data.acc(include)'];
    rt= [rt; data.rt1(include)'];
    condition= [condition; tmp(include)'];
        
end

%% Create data matrix and labels
labels ={'subj_idx','response','rt','condition'};
values = [subj_idx response rt condition];

%% Export data matrix and labels
csvwrite_with_headers('hddm_input_prescan.csv',values,labels);

%% -----------------------------------------------------------------------
%% PRE-SCAN

% Initialise variables using HDDM naming convention
subj_idx= [];
response= [];
rt= [];
condition= [];

%% LOOP THROUGH SUBJECTS
for i_sbj = 1:length(subjects)
    
    %% HOUSE KEEEPING
    clear tmp;
    
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
    
    % re-assign
    data = tmp;
    
    % add stimulus time to RT1 and transform to seconds
    data.rt1= (data.rt1./1000)+1;
    
    % Include trials based on deviation from grand mean
    rt1 = (data.rt1./1000);
    centre = mean(rt1);
    stdval = std(rt1)*2.5;
    include = (rt1>(centre-stdval))&(rt1<(centre+stdval));

    % delta categories
    sindx = 1:length(data.acc);
    c_del = abs(data.deltaz(sindx));
    u_del = unique(c_del);
    for t = 1:length(sindx); data.bodcat(sindx(t)) = sum(c_del(t)>=u_del); end
    
    % turn into conditions
    clear tmp
    for i = 1:length(data.cohcat);
        if data.cohcat(i)==2 && data.bodcat(i)==1
            tmp(i)=1;
        elseif data.cohcat(i)==2 && data.bodcat(i)==2
            tmp(i)=2;
        elseif data.cohcat(i)==2 && data.bodcat(i)==3
            tmp(i)=3;
        elseif data.cohcat(i)==2 && data.bodcat(i)==4
            tmp(i)=4;
        elseif data.cohcat(i)==1 && data.bodcat(i)==1
            tmp(i)=5;
        elseif data.cohcat(i)==1 && data.bodcat(i)==2
            tmp(i)=6;
        elseif data.cohcat(i)==1 && data.bodcat(i)==3
            tmp(i)=7;
        elseif data.cohcat(i)==1 && data.bodcat(i)==4
            tmp(i)=8;
        end
    end

    
    % concatenate
    subj_idx= [subj_idx; ones(sum(include),1).*i_sbj];
    response= [response; data.acc(include)'];
    rt= [rt; data.rt1(include)'];
    condition= [condition; tmp(include)'];
        
end

%% Create data matrix and labels
labels ={'subj_idx','response','rt','condition'};
values = [subj_idx response rt condition];

%% Export data matrix and labels
csvwrite_with_headers('hddm_input_scan.csv',values,labels);