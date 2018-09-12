% Bang & Fleming (2018) Distinct encoding of decision confidence in human
% medial prefrontal cortex
%
% Pseudo-script for ROI specification using leave-one-out cross-validation 
% procedure in order to avoid circularity in ROI analyses
%
% We do this in a fairly 'dumb' way: (a) run second-level analysis N times
% where N is number of subjects, each time leaving out the subject for whom
% we want to create an ROI mask --> (b) load a subject's second-level t-map 
% into xjview, select the cluster(s) of interest and save image as a mask 
% --> (c) use subject mask for ROI analyses.
%
% It should be possible to develop a more 'automated' procedure using coordinates 
% from second-level t-map for all subjects and the MarsBaR toolbox
%
% Requires that model estimation be completed first by 
% meta_getonsetsModel_...... and meta_1stlevelModel_......
% 
% Dan Bang danbang.db@gmail.com 2018

%% -----------------------------------------------------------------------
%% PREPARATION

% fresh memory
clear;

% Subjects
subjects = 1:10;

% =========================== Set up contrasts ============================
% =========================================================================

%% Current directory;
cwd = pwd;

%% Model name
this_model= 'GLM1';
this_ROI= 'pgACC';

%% subjects
k=0; for i_s=subjects; k=k+1; subjects{k}=['s',num2str(i_s)]; end

%% Directories
fs               = filesep;
dir_spm          = 'path_spm';
dir_brain        = 'path_preprocessed_scan_data';
dir_epi          = 'Functional';
dir_run          = 'Run';
dir_spike        = 'Spike';
dir_behaviour    = 'path_behavioural_data';
dir_stats        = ['Functional_1st_level',fs,this_model];
dir_sbjROI       = ['ROI',fs,this_ROI];
dir_info         = 'path_subject_information';

%% Load subject details (see fil_subject_details)
load([dir_info,fs,'fil_subject_details.mat']);

%% Add paths
addpath(genpath(dir_spm));

%% Current contrasts as per first-level analysis
contrastNames = {'CxD'};
conImages = {'con_0001.nii,1'};

%% Loop through subjects
for s = 1:length(subjects)
 
%% Loop through contrasts
for j = 1:length(contrastNames)
    
    cov = [];
    contrastFolder = [dir_brain_root,fs,subjects{s},fs,dir_sbjROI]; 
    
    % if the intended resultsFolder directory does not exits, make it and go into it
    if exist(contrastFolder,'dir') ~= 7
        mkdir(contrastFolder);
        cd(contrastFolder);
    else
        % if it does exist, go into it and delete its contents.
        % change this for more than one 2nd level test
        cd(contrastFolder);
        delete('*.*');
    end
    
    % setup job structure for 2nd level t-test of individual contrasts
    %---------------------------------------------------------------------
    AllSubjectsExcept1 = find(subjects~=subjects(s)); % all subjects except one
    k = 0;
    for n = AllSubjectsExcept1 % select con image from all subjects except one
        k = k + 1;
        conDir = [dir_brain_root,fs,subjects{n},fs,dir_stats];
        conImg = conImages{j};
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{k,1}    = [conDir,fs,conImg] ;
    end
        
    % Check if there is a covariate for this contrast, and compute
    if ~isempty(cov)
        for l = 1:length(cov)
            matlabbatch{1}.spm.stats.factorial_design.cov(l) = struct('c',cov{l}','cname',covName{l},'iCFI',1,'iCC',1);
        end
    else
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c',{},'cname',{},'iCFI',{},'iCC',{});
    end
    
    % Setup job
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = [];
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {[dir_spm,fs,'tpm',fs,'mask_ICV.nii']};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = [];
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = [];
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    matlabbatch{1}.spm.stats.factorial_design.dir = {contrastFolder};
    
    % run 2nd level test
    % ---------------------------------------------------------------------
    % save and run job
    save second_level.mat matlabbatch
    disp(['RUNNING second level test']);
    spm_jobman('run','second_level.mat');
    clear matlabbatch
    
    % Ensure implicit masking is switched off
    load SPM
    SPM.xM.I = 0;
    save SPM SPM
    
    jobs{1}.stats{1}.fmri_est.spmmat = {[contrastFolder,fs,'SPM.mat']};
    save estimate.mat jobs
    disp(['RUNNING model estimation'])
    spm_jobman('run','estimate.mat');
    clear jobs
    
    % Specify one-sample tcontrasts
    % ---------------------------------------------------------------------
    contrasttype = contrastNames{j};
    contr_input = [1];              % one-sample RFX significance
    
    % setup job structure for 2nd level contrast using newly created
    % SPM.mat
    jobs{1}.stats{1}.con.spmmat    = {[contrastFolder,fs,'SPM.mat']};
    jobs{1}.stats{1}.con.consess{1}.tcon.name           = contrasttype;
    jobs{1}.stats{1}.con.consess{1}.tcon.convec         = contr_input;
    jobs{1}.stats{1}.con.consess{1}.tcon.sessrep        = 'none';
    jobs{1}.stats{1}.con.delete                         = 1;
    
    % save and run job
    save contrasts.mat jobs
    disp(['RUNNING 2nd level test for 1st level contrast  ' contrastNames{j}]);
    spm_jobman('run','contrasts.mat');
    clear jobs
    
end

cd(cwd);

end