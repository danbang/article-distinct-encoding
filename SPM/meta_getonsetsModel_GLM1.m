function meta_getonsetsModel_GLM1(subjects)
% Script constructs and estimates 1st-level design matrix
% Adds motion covariates (realignment + derivatives)
% Flexibly includes physiology covariates (Spike)
%
% Dan Bang danbang.db@gmail.com 2018

%% Current directory;
cwd = pwd;

%% Model name
this_model = 'GLM1';

%% Options for processing
compute_conditions  = 1;
construct_design    = 1;
estimate_design     = 1;
model_physiology    = 1;

%% Scan parameters
n_block     = 5;
TR          = 3.36;
nslices     = 48;
hpcutoff    = 128;  % hp filter
timeNorm    = 1000; % term for normalise times to seconds

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

%% Load physiology details
load([dir_info,fs,'fil_subject_spike.mat']);

%% Add paths
addpath(genpath(dir_spm));

%% loop through subjects
for n = 1:length(subjects)
           
    %% Create output directory
    outputDir=[dir_brain,fs,'s',num2str(subjects(i_sbj)),fs,dir_stats];
    % if it does not exits, make it and go into it
    if ~exist(outputDir,'dir'); mkdir(outputDir);
    % if it does exist, go into it and delete existing contrast files
    else cd(outputDir); delete('SPM.mat','*.img','*.hdr'); end
    
    %=========================================================================
    %% Get onsets in scans from behavioural datafiles and define images
    %======================================================================
    
    cd(cwd);   
    
    %% Functional data
    for k = 1:n_block;
        if compute_conditions
            %% Display
            disp(['Computing event onsets for subject ',num2str(subjects(i_sbj)),' -- session ',num2str(k)]);
            %% Behavioural data
            datafile = [dir_behaviour,'s',num2str(subjects(i_sbj)),'_task_b',num2str(k),'.mat'];
            load(datafile); hat_onsets = onsets;    
            %% Compute onsets and durations in SECONDS
            % Time constant for correction
            constant = slices.dummy.stamp+(TR/nslices)*timeNorm; % first slice after last slice of last dummy volume
            % Onsets
            stmOnsets  = []; for i_t = 1:length(data.trial); stmOnsets(i_t)  = (hat_onsets{i_t}.stm-constant)./timeNorm; end
            decOnsets  = []; for i_t = 1:length(data.trial); decOnsets(i_t)  = (hat_onsets{i_t}.dec-constant)./timeNorm; end
            decOffsets = []; for i_t = 1:length(data.trial); decOffsets(i_t) = (hat_onsets{i_t}.fbk1-constant)./timeNorm; end
            conOnsets  = []; for i_t = 1:length(data.trial); conOnsets(i_t)  = (hat_onsets{i_t}.con-constant)./timeNorm; end
            conOffsets = []; for i_t = 1:length(data.trial); conOffsets(i_t) = (hat_onsets{i_t}.fbk2-constant)./timeNorm; end
            % Confidence presses
            conidx = find(data.concat==1);
            keyPress      = []; for i_t = conidx; keyPress      = [keyPress numel(routes{i_t}.conf)-1]; end
            conidx = find(data.concat==1&data.acc==1);
            keyPress_Corr = []; for i_t = conidx; keyPress_Corr = [keyPress_Corr numel(routes{i_t}.conf)-1]; end
            conidx = find(data.concat==1&data.acc==0);
            keyPress_Wrng = []; for i_t = conidx; keyPress_Wrng = [keyPress_Wrng numel(routes{i_t}.conf)-1]; end
            % Durations
            dec_duration = decOffsets-stmOnsets;
            con_duration = conOffsets-conOnsets;
            % RT outlier
            rt1      = log(data.rt1./timeNorm);
            centre   = mean(rt1);
            std3     = std(rt1)*2.5;
            outlier1 = (rt1<(centre-std3))|(rt1>(centre+std3));
            data.concat = logical(data.concat);
            %% Define conditions
            % Define coherence by delta
            conditions = [5]; 
            for i_t = 2:length(data.trial); 
               %% OF INTEREST
               % low coherence, low distance
               if     data.cohcat(i_t) == 1 && data.bodcat(i_t) == 1 && outlier1(i_t)==0; conditions(i_t) = 1; 
               % low coherence, high distance
               elseif data.cohcat(i_t) == 1 && data.bodcat(i_t) == 2 && outlier1(i_t)==0; conditions(i_t) = 2; 
               % high coherence, low distance
               elseif data.cohcat(i_t) == 2 && data.bodcat(i_t) == 1 && outlier1(i_t)==0; conditions(i_t) = 3; 
               % high coherence, high distance
               elseif data.cohcat(i_t) == 2 && data.bodcat(i_t) == 2 && outlier1(i_t)==0; conditions(i_t) = 4; 
               %% OF NO INTEREST
               elseif outlier1(i_t)==1; conditions(i_t) = 5;
               end
            end
            %% Condition onsets
            cldlOnsets = []; cldlOnsets=stmOnsets(conditions==1); cldlDuration=dec_duration(conditions==1);
            cldhOnsets = []; cldhOnsets=stmOnsets(conditions==2); cldhDuration=dec_duration(conditions==2);
            chdlOnsets = []; chdlOnsets=stmOnsets(conditions==3); chdlDuration=dec_duration(conditions==3);
            chdhOnsets = []; chdhOnsets=stmOnsets(conditions==4); chdhDuration=dec_duration(conditions==4);
            wrngOnsets = []; wrngOnsets=stmOnsets(conditions==5); wrngDuration=dec_duration(conditions==5);
            %% Set up design
            % initialise
            names = []; onsets = []; durations = [];
            % create conditions
            i_c = 1;   names{i_c} = 'cldl';   onsets{i_c} = cldlOnsets;     durations{i_c} = cldlDuration;
            i_c = 2;   names{i_c} = 'cldh';   onsets{i_c} = cldhOnsets;     durations{i_c} = cldhDuration;      
            i_c = 3;   names{i_c} = 'chdl';   onsets{i_c} = chdlOnsets;     durations{i_c} = chdlDuration;
            i_c = 4;   names{i_c} = 'chdh';   onsets{i_c} = chdhOnsets;     durations{i_c} = chdhDuration;
            i_c = 5;   names{i_c} = 'noise';  onsets{i_c} = wrngOnsets;     durations{i_c} = wrngDuration;
            i_c = 6;   names{i_c} = 'rate';   onsets{i_c} = conOnsets(data.concat==1);      durations{i_c} = con_duration(data.concat==1);
            % parametric modulation
            i_c = 6;   pmod(i_c).name{1} = 'presses';  pmod(i_c).param{1} = zscore(keyPress);   pmod(i_c).poly{1} = 1;
            cd(outputDir);
            conditionFile = sprintf('conditions%d.mat',k);
            save(conditionFile, 'names', 'onsets', 'durations', 'pmod');    
        end
        
        %==========================================================================
        %% Construct design matrix
        %==========================================================================
        
        % Load files we have just created
        epiDir = [dir_brain,fs,'s',num2str(subjects(i_sbj)),fs,dir_epi,fs,dir_run,num2str(k)];
        conditionFile = sprintf('conditions%d.mat',k);       
        % Get text file with movement regressors and concatenate with
        % first derivative
        mCorrFile = spm_select('List',epiDir,'^rp_af.*\.txt$');
        M = textread([epiDir,fs,mCorrFile]);
        R = [M [zeros(1,6); diff(M)]];
        % Add physiology data to nuisance matrix
        if model_physiology
            % But only if data was recorded / is usable
            if subj{subjects(n)}.spike(k)
                if subjects(n)<9
                spikeFile = [dir_brain,fs,'s',num2str(subjects(i_sbj)),fs,dir_spike,fs,'s',num2str(subjects(i_sbj)),'_ex_b',num2str(k),'_R_session1'];
                else
                spikeFile = [dir_brain,fs,'s',num2str(subjects(i_sbj)),fs,dir_spike,fs,'s',num2str(subjects(i_sbj)),'_b',num2str(k),'_R_session1'];    
                end
                temp = load(spikeFile,'R');
                R = [R temp.R(1:size(R,1),:)];
            end
        end
        % Write matrix
        cd(outputDir);
        multiFile = sprintf('multireg%d.mat',k);
        save(multiFile, 'R');
        % Assign .mat file with onsets/names/pmods in to path
        conditionPath = [outputDir,fs,conditionFile];
        multiregPath = [outputDir,fs,multiFile];     
        % Get epi files for this session
        epiDir = [dir_brain,fs,'s',num2str(subjects(i_sbj)),fs,dir_epi,fs,dir_run,num2str(k)];
        % select scans and concatenate
        f      = spm_select('List',epiDir,'^swuaf.*\.nii$');     % Select smoothed normalised images
        files  = cellstr([repmat([epiDir fs],size(f,1),1) f]);
        % prepare job
        jobs{1}.stats{1}.fmri_spec.sess(k).scans = files; 
        jobs{1}.stats{1}.fmri_spec.sess(k).multi = {conditionPath};
        jobs{1}.stats{1}.fmri_spec.sess(k).multi_reg = {multiregPath};
        % high pass filter
        jobs{1}.stats{1}.fmri_spec.sess(k).hpf = hpcutoff;
        jobs{1}.stats{1}.fmri_spec.sess(k).cond = struct([]);
        jobs{1}.stats{1}.fmri_spec.sess(k).regress = struct([]);
        % clear temporary variables for next run
        f = []; files = [];   
    end
    
    %==========================================================================
    %======================================================================
    jobs{1}.stats{1}.fmri_spec.dir = {outputDir};
    % timing variables
    jobs{1}.stats{1}.fmri_spec.timing.units     = 'secs';
    jobs{1}.stats{1}.fmri_spec.timing.RT        = TR;
    jobs{1}.stats{1}.fmri_spec.timing.fmri_t    = nslices;
    jobs{1}.stats{1}.fmri_spec.timing.fmri_t0   = nslices/2;   
    % basis functions
    jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs = [0 0];
    % model interactions (Volterra) OPTIONS: 1|2 = order of convolution
    jobs{1}.stats{1}.fmri_spec.volt             = 1;
    % global normalisation
    jobs{1}.stats{1}.fmri_spec.global           = 'None';
    % explicit masking
    jobs{1}.stats{1}.fmri_spec.mask             = {[dir_spm,fs,'tpm/mask_ICV.nii']};
    % serial correlations
    jobs{1}.stats{1}.fmri_spec.cvi              = 'AR(1)';
    % no factorial design
    jobs{1}.stats{1}.fmri_spec.fact             = struct('name', {}, 'levels', {});
    
    
    %==========================================================================
    %% run model specification
    %==========================================================================
    
    if construct_design
        % save and run job
        cd(outputDir);
        save specify.mat jobs
        disp(['RUNNING model specification for subject ', num2str(subjects(i_sbj))]);
        spm_jobman('run','specify.mat');
        clear jobs
    end
    
    % Ensure implicit masking is switched off
    load SPM
    SPM.xM.TH = repmat(-Inf, length(SPM.xM.TH), 1);
    SPM.xM.I = 0;
    save SPM SPM
    
    %% Estimate
    % setup job structure for model estimation and estimate
    % ---------------------------------------------------------------------
    if estimate_design
        jobs{1}.stats{1}.fmri_est.spmmat = {[outputDir,fs,'SPM.mat']};
        save estimate.mat jobs
        disp(['RUNNING model estimation for subject ',num2str(subjects(i_sbj))])
        spm_jobman('run','estimate.mat');
        clear jobs
    end
end
cd(cwd);

end