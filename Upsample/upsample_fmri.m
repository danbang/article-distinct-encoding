% Bang & Fleming (2018) Distinct encoding of decision confidence in human
% medial prefrontal cortex
%
% Pseudo-script for upsamling fMRI activity for single-trial analysis
%
% We here assume that the same mask is used for all subjects 
% (e.g. anatomically defined ROI)
%
% Dan Bang danbang.db@gmail.com 2018

%% -----------------------------------------------------------------------
%% PREPARATION

% fresh memory
clear;

% Subjects
subjects = 1:10;

%% Current directory;
cwd = pwd;

%% Scan parameters
n_runs      = 5;    % number of functional runs
TR          = 3.36;
n_slices    = 48;
n_motion    = 12;   % motion covariates

%% Upsample information
uinfo.samp_reso     = .144;                  % sample resolution in miliseconds
uinfo.samp_rate     = uinfo.samp_reso./TR;   % sample rate
uinfo.window_length = 12;                    % time window in seconds from onset

%% Mask information
minfo.name = {'ventral striatum'};
minfo.file = {'vstriatum.img'};

% Paths [change 'repoBase' according to local setup]
fs = filesep;
dataBehaviour = 'POINT-TO-BEHAVIOURAL-DATA';
dataScan = 'POINT-TO-FMRI-DATA';
dirEPI = 'FUNCTIONAL-DATA-NAME';
dirMasks = 'POINT-TO-ROI-MASKS';

%% -----------------------------------------------------------------------
%% UPSAMPLE

%% Loop through ROIs
for i_roi = 1:length(minfo.file);
    
    %% Starting ROI
    fprintf(['===== Upsampling activity for ',minfo.name{i_roi},' ===== \n']);
        
    %% Clear variable for storage
    clear time_series;    
    
    %% loop through subjects
    for i_sbj = 1:length(subjects);
    
        %% Initialise variable for storage
        time_series = [];
        
        %% Load masks
        mask_path = [dirMasks,fs,minfo.file{i_roi}];
        
        %% Loop through runs
        for i_run = 1:n_runs;
            
            %% Load data with time logs
            data = [dirData,fs,'s',subjects(i_sbj),'_b',num2str(i_run),'.mat'];
            trialOnsets = REPLACE; % specify trial onsets (e.g. 1 second before stimulus onset)
            
            %% Scans
            % clear variables
            clear ts M R
            % data directory
            epiDir = [dirBrain,fs,'s',num2str(i_sbj),fs,dirEPI,fs,sess_prfx,num2str(i_run)];
            % select scans and concatenate
            tmp = spm_select('List', epiDir, '^swuaf.*\.nii$'); % Select smoothed normalised images
            files = cellstr([repmat([epiDir fs],size(tmp,1),1) tmp]); 
            % get scans
            V = spm_vol(files);
            % get and check mask (if different dimensions, use mask_reslice.m)
            Vmask   = spm_vol(deblank(mask_path));
            mask    = spm_read_vols(Vmask);
            if sum(V{1}.dim == Vmask.dim) ~= 3; error('mask and images have different dimensions'); end
            %% get time series
            for i = 1:length(files)
                img = spm_read_vols(V{i});
                dat = img(mask > 0);
                ts(i) = nanmean(dat(:));
            end
            %% remove confounds
            % scan specs for cosine basis set
            params.nscan    = length(files);
            params.nsess    = 1;
            params.tr       = TR;          
            % motion covariates + derivatives
            mCorrFile = spm_select('List',epiDir,'rp_af.*\.txt$');
            M         = textread([epiDir fs mCorrFile]);
            R         = [M [zeros(1,6); diff(M)]];
            % execute removal
            ts = ts(:);  % column vector must be entered into remove_confounds
            ts = remove_confounds(ts, R, params);
            ts = ts(:)'; % row vector for the future
            % re-sample and average timeseries at specified sample rate
            t      = 0:length(ts)-1;
            t_ups  = 0:spec.samp_rate:length(ts)-1;
            ts_ups = spline(t,ts,t_ups);    % interpolates time series
            % create trial-by-time series matrix for each scan run
            for i_trial = 1:length(trialOnsets);
                win_min  = ceil(trialOnsets(i_trial)./spec.samp_reso);       % window start in miliseconds
                win_max  = ceil(win_min+spec.window_length./spec.samp_reso); % window end in miliseconds
                if (size(ts_ups,2)-win_max)>0
                c_ts_ups = ts_ups(win_min:win_max);
                else 
                c_ts_ups = NaN(1,(win_max-win_min)+1);
                end
                time_series_run(i_trial,:) = c_ts_ups; % update
            end
            % concatenate run time series
            time_series = [time_series; time_series_run];
    
        end

        %% log subject data
        roi{i_roi}.time_series{i_sbj} = time_series;
          
    end
   
end

%% save output
save('upsampled_frmi','roi','uinfo','minfo');