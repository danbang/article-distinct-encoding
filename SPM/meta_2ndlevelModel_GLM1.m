function meta_2ndlevelModel_GLM1(subjects)
% Specifies contrasts for 2nd-level analysis
% Requires that model estimation be completed first by 
% meta_getonsetsModel_...... and meta_1stlevelModel_......
%
% Dan Bang danbang.db@gmail.com 2018

% =========================== Set up contrasts ============================
% =========================================================================

%% Current directory;
cwd = pwd;

%% Model name
this_model = 'GLM1';

%% Directories
fs               = filesep;
dir_spm          = 'path_spm';
dir_brain        = 'path_preprocessed_scan_data';
dir_epi          = 'Functional';
dir_run          = 'Run';
dir_spike        = 'Spike';
dir_behaviour    = 'path_behavioural_data';
dir_stats        = ['Functional_1st_level',fs,this_model];
dir_stats2_base  = 'path_group_level_maps';
dir_stats2_name  = ['Functional_2nd_level',fs,this_model];
dir_info         = 'path_subject_information';

%% Load subject details (see fil_subject_details)
load([dir_info,fs,'fil_subject_details.mat']);

%% Add paths
addpath(genpath(dir_spm));

%% Current contrasts
contrastNames = {'high_beats_low_coh','high_beats_low_dis','low_beats_high_coh','low_beats_high_dis','interaction_pos','interaction_neg'};
conImages = {'con_0001.nii,1','con_0002.nii,1','con_0003.nii,1','con_0004.nii,1','con_0005.nii,1','con_0006.nii,1'};	% these correspond to the relevant contrast images in your first-level script

%% Loop through contrasts
for j = 1:length(contrastNames)
    
    cov = [];
    contrastFolder = [dir_stats2_base,fs,dir_stats2_name,fs,contrastNames{j}]; 
    
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
    for n = 1 : length(subjects)  % select con image from every subject
        conDir = [dir_brain,fs,'s',num2str(subjects(i_sbj)),fs,dir_stats];
        conImg = conImages{j};
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{n,1}    = [conDir,fs,conImg] ;
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
