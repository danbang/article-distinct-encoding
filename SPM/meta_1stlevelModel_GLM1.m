function meta_1stlevelModel_GLM1(subjects)
% Specifies contrasts for 1st-level analysis
% Requires that model estimation be completed first by 
% meta_getonsetsModel......
%
% Dan Bang danbang.db@gmail.com 2018

% =========================== Set up contrasts ============================
% =========================================================================

%% Current directory;
cwd = pwd;

%% Model name
this_model = 'GLM1';

%% Scan parameters
n_block     = 5;
TR          = 3.36;
nslices     = 48;
hpcutoff    = 128; % hp filter
n_motion    = 12;
n_spike     = 14; 
n_noise     = n_motion+n_spike;  % covariates of no interest

%% Directories
fs               = filesep;
dir_spm          = 'path_spm';
dir_brain        = 'path_preprocessed_scan_data';
dir_epi          = 'Functional';
dir_run          = 'Run';
dir_spike        = 'Spike';
dir_behaviour    = 'path_behavioural_data';
dir_stats        = ['Functional_1st_level',fs,this_model];
dir_info         = 'path_subject_information';

%% Load subject details (see fil_subject_details)
load([dir_info,fs,'fil_subject_spike.mat']);

%% Add paths
addpath(genpath(dir_spm));
addpath([pwd,fs,'routines']);

%% Launch SPM
% spm fmri

%% loop through subjects
for i_s = 1 : length(name_subj)
    % Contrast names
    T.contrasts = {'high > low coh','high > low dis','high < low coh','high < low dis','interaction (+)','interaction (-)'};
   	% Initialise contrast index
    j= 1;
    % You need as many rows here as entries in the cell array above
    % high > low coherence
    contrasts_v = [-1 -1 +1 +1  0  0  0];  T.contrastVectors(j,:) = adaptive_block_contrast_spike(contrasts_v,subj{subjects(i_s)}.spike(1:n_block),n_block,n_motion,n_spike); j= j+1;
    % high > low absolute distance
    contrasts_v = [-1 +1 -1 +1  0  0  0];  T.contrastVectors(j,:) = adaptive_block_contrast_spike(contrasts_v,subj{subjects(i_s)}.spike(1:n_block),n_block,n_motion,n_spike); j= j+1;
    % high < low coherence
    contrasts_v = [+1 +1 -1 -1  0  0  0];  T.contrastVectors(j,:) = adaptive_block_contrast_spike(contrasts_v,subj{subjects(i_s)}.spike(1:n_block),n_block,n_motion,n_spike); j= j+1;
    % high < low absolute distance
    contrasts_v = [+1 -1 +1 -1  0  0  0];  T.contrastVectors(j,:) = adaptive_block_contrast_spike(contrasts_v,subj{subjects(i_s)}.spike(1:n_block),n_block,n_motion,n_spike); j= j+1;
    % interaction (+)
    contrasts_v = [+1 -1 -1 +1  0  0  0];  T.contrastVectors(j,:) = adaptive_block_contrast_spike(contrasts_v,subj{subjects(i_s)}.spike(1:n_block),n_block,n_motion,n_spike); j= j+1;
    % interaction (-)
    contrasts_v = [-1 +1 +1 -1  0  0  0];  T.contrastVectors(j,:) = adaptive_block_contrast_spike(contrasts_v,subj{subjects(i_s)}.spike(1:n_block),n_block,n_motion,n_spike); j= j+1;
    % Job
    jobs{1}.stats{1}.con.spmmat    = {[dir_brain,fs,name_subj{i_s},fs,dir_stats,fs,'SPM.mat']};
    j = 1; % job index
    %% loop through contrasts
    for cont_nr = 1:length(T.contrasts)
        % Specify tcontrasts
        contrasttype = T.contrasts{cont_nr};
        contr_input = T.contrastVectors(cont_nr,:);
        % setup job structure for contrasts
        jobs{1}.stats{1}.con.consess{j}.tcon.name           = contrasttype;
        jobs{1}.stats{1}.con.consess{j}.tcon.convec         = contr_input;
        jobs{1}.stats{1}.con.consess{j}.tcon.sessrep        = 'none';
        j=j+1; % update job index
    end  
    % if 1 then all existing contrasts are deleted
    jobs{1}.stats{1}.con.delete                         = 1;   
    % output
    outputDir = [dir_brain,fs,'s',num2str(subjects(i_sbj)),fs,dir_stats];
    cd (outputDir);
    % save and run job
    save contrasts.mat jobs
    disp(['RUNNING contrast specification for subject number  ', num2str(subjects(i_sbj))]);
    spm_jobman('run','contrasts.mat');
    disp(['Contrasts created for subject ', num2str(subjects(i_sbj))]);
    clear jobs T
    cd(cwd);
end   % end of subject loop

end
