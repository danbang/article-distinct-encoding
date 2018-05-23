function [output] = remove_confounds(input, R, params)
% Remove confounds from timeseries with spm_dctmtx

%% Remove confounds

n = fix((2*(params.nscan*params.tr))./128);
dct = spm_dctmtx(params.nscan*params.nsess,n);
% Motion covariates and session mean
dct = [dct R];
output=input-dct*pinv(dct)*input;