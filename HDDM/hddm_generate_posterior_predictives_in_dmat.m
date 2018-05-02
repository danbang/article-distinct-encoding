% Bang & Fleming (2018) Distinct encoding of decision confidence in human
% medial prefrontal cortex
%
% Plot aggregate posterior predictives after HDDM estimation of DDM
% parameters
%
% Requires the DMAT toolbox to run DDM simulations:
% https://ppw.kuleuven.be/okp/software/dmat/
%
% Dan Bang danbang.db@gmail.com 2018

%% -----------------------------------------------------------------------
%% PREPARATION

% fresh memory
clear; close all;

% Paths [change 'repoBase' according to local setup]
fs = filesep;
repoBase = [getDropbox(1),fs,'Ego',fs,'Matlab',fs,'ucl',fs,'sensory_vs_decision',fs,'Repository'];
dmatDir = [repoBase,fs,'Toolbox',fs,'dmat'];

% Add DMAT toolbox
addpath(dmatDir);

% DMAT specification
N = 1e4; % number of simulations per condition

%% -----------------------------------------------------------------------
%% PRE-SCAN

% load HDDM parameter estimates
load('hddm_parameters_group_mean_prescan.mat');

% assign parameters
% multiply by .1 for comptability with HDDM:
% http://ski.clps.brown.edu/hddm_docs/howto.html#compare-parameters-to-other-papers
v_drift= hddm.drift.*.1;
v_nondt= hddm.nondt; % subtract stimulus duration
v_bound= hddm.bound.*.1;
sv = hddm.sv.*.1;
st = hddm.st; 
sz = 0; % no variability in DDM bias
clear sim;

% loop through conditions
for i = 1:length(v_drift)
    
    % current parameters
    a = v_bound(i);
    t = v_nondt(i);
    v = v_drift(i);
    
    % simulate trial and save RT and accuracy 
    par = [a t sv a/2 sz st v];
    [rt1 acc] = simuldiff(par, N);
    sim.rt1(:,i) = rt1;
    sim.acc(:,i) = acc;
    
end

% organise simulation output 
% (L: low; H: high; C: coherence)
LC_rt1= sim.rt1(:,5:8);
HC_rt1= sim.rt1(:,1:4);
LC_acc= sim.acc(:,5:8);
HC_acc= sim.acc(:,1:4);

% save simulation output
save('hddm_posterior_predictives_prescan','LC_acc','HC_acc','LC_rt1','HC_rt1');

%% -----------------------------------------------------------------------
%% SCAN

% load HDDM parameter estimates
load('hddm_parameters_group_mean_scan.mat');

% assign parameters
% multiply by .1 for comptability with HDDM:
% http://ski.clps.brown.edu/hddm_docs/howto.html#compare-parameters-to-other-papers
v_drift= hddm.drift.*.1;
v_nondt= hddm.nondt; % subtract stimulus duration
v_bound= hddm.bound.*.1;
sv = hddm.sv.*.1;
st = hddm.st; 
sz = 0; % no variability in DDM bias
clear sim;

for i = 1:length(v_drift)
    
    % current parameters
    a = v_bound(i);
    t = v_nondt(i);
    v = v_drift(i);
    
    % simulate trial and save RT and accuracy 
    par = [a t sv a/2 sz st v];
    [rt1 acc] = simuldiff(par, N);
    sim.rt1(:,i) = rt1;
    sim.acc(:,i) = acc;
    
end

% organise simulation output 
% (L: low; H: high; C: coherence)
LC_rt1= sim.rt1(:,3:4);
HC_rt1= sim.rt1(:,1:2);
LC_acc= sim.acc(:,3:4);
HC_acc= sim.acc(:,1:2);

% save simulation output
save('hddm_posterior_predictives_scan','LC_acc','HC_acc','LC_rt1','HC_rt1');