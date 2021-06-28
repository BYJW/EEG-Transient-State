clear all
clc
close all

%%

% The script corresponding to the state inference and transfer analysis used in manuscript 'Spontaneous transient brain states in EEG source space of disorders of consciousness'
% an OSL and HMM toolbox are needed, which can be accessed from https://ohba-analysis.github.io/osl-docs/ and https://github.com/OHBA-analysis/HMM-MAR/wiki
% This script follows pipeline from Vidaurre et. al. Nat Commun, 2018: Vidaurre, D., Hunt, L. T., Quinn, A. J., Hunt, B. a. E., Brookes, M. J., Nobre, A. C. & Woolrich, M. W. 2018. Spontaneous cortical activity transiently organises into frequency specific phase-coupling networks. Nat Commun, 9, 2987.
% The original pipeline on MEG analysis could be found from https://github.com/OHBA-analysis/HMM-MAR/tree/master/examples
% 

%  edit by Yang Bai 2021-06-22


%% addin the toolboxs
scriptdir = '...'; % <- edit this line
addpath( fullfile( scriptdir , 'scripts' ) );
outdir = '...'; % <- edit this line
% Initialise OSL
addpath( fullfile(scriptdir,'toolboxes','osl','osl-core') );
osl_startup
% Add netlab and fmt
addpath( fullfile(osldir,'ohba-external','netlab3.3','netlab') );
addpath( fullfile(osldir,'ohba-external','fmt') );
% Add HMM-MAR to path
addpath(genpath( fullfile(scriptdir,'toolboxes','HMM-MAR-master') ));
% Add distibutionplot to path
addpath(genpath( fullfile(scriptdir,'toolboxes','distributionPlot') ));

%% reconstruct the data from all subjects and extract information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  above are calculated with single patient
%  next are calculated with all patient, need link all the data together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = [];  % HMM ready dataset
T_all=[];   % Length of each data
R_all=[];   % location of the single subject data in the concatenated dataset
T=[];
R=[];
for numsub=1:totalnum  %%%% concatenate all the single data into a matrix to do state inference
%%%%%%%%%%%%%%%%% load the parcellated sourcedata
load(fullfile(outdir,[' ... _paceldata_lcmv']));     
nsamples(numsub) = size(paceldata,2);
fs=250;
Ltrials=triallength*fs;
ntrials = nasamples(numsub)/Ltrials;

Ts = [];               % Length of continuous good segments
Rs = [];               % Indices of single run within data
for num=1:ntrials
pow=paceldata(:,(num-1)*Ltrials+1:num*Ltrials);
dat = ROInets.remove_source_leakage(pow, 'symmetric');
t = size(dat,2);
offset = sum(Ts);
Rs = cat(1,Rs,[offset+1 offset+size(dat,2)]);
Ts = cat(1,Ts,t);
data = cat(2,data,dat);
end
T=cat(1,T,Ts);
R=cat(1,R,Rs);
T_all(numsub)=nsamples(numsub);
R_all=cat(1,R_all,[sum(nsamples(1:numsub-1)),sum(nsamples)]);
end
%% Resolve dipole sign ambiguity  
x = data';
options_signflip = [];
options_signflip.maxlag = 7; 
options_signflip.verbose = 0;
flips = findflip(x,T_all,options_signflip);
data = flipdata(x,T_all,flips);
%% HMM inference
numstate=8;
embeddedlag = 7; 
K = numstate; Hz = 250; ndim = 42;

options = struct();
options.order = 0;
options.zeromean = 1;
options.covtype = 'full';
options.embeddedlags = -embeddedlag:embeddedlag;
options.pca = ndim*2;
options.K = K;
options.Fs = Hz;
options.verbose = 1;
options.onpower = 0; 
options.standardise = 1;
options.standardise_pc = options.standardise;
options.inittype = 'HMM-MAR';
options.cyc = 100;
options.initcyc = 10;
options.initrep = 5;
options.useParallel=0;  %%%%% use the parallel, can open it in real analysis
% % % stochastic options, refer to Guideline
options.BIGNinitbatch = 10;
options.BIGNbatch = 10;
options.BIGtol = 1e-7;
options.BIGcyc = 500;
options.BIGundertol_tostop = 5;
options.BIGdelay = 5;
options.BIGforgetrate = 0.7;
options.BIGbase_weights = 0.9;

[hmm, Gamma, ~, vpath] = hmmmar(data,T',options);

%% the transfer state analysis, using some scripts modified from HMMMAR toolbox
h=load('hmm');
h.hmm.P=[]; 
h.hmm.Pi=[];
h.hmm.train= rmfield( h.hmm.train, 'A' );
[Gamma, vpath, hmm]=statetransfer(data,T',options,h);

%% spectrum features
%%%%% padding the Gamma equal to the real data length,
pad_options = struct();
pad_options.embeddedlags = -7:7;
Gamma = padGamma(Gamma, T, pad_options);

N = length(T_all);
Hz=250;
options_mt = struct('Fs',Hz); 
options_mt.fpass = [1 40];  % band of frequency you're interested in
options_mt.tapers = [4 7]; % taper specification - leave it with default values
options_mt.p = 0; %0.01; % interval of confidence  
options_mt.win = 2 * Hz; % multitaper window
options_mt.to_do = [1 0]; % turn off pdc
options_mt.order = 0;
options_mt.embeddedlags = -7:7;

% average
fitmt = hmmspectramt(data,size(data,1),Gamma,options_mt);
% per subject
fitmt_subj = cell(N,1);
for ind=1:N
    subj_data = data(R_all(ind,1):R_all(ind,2),:);
    fitmt_subj{ind} = hmmspectramt(subj_data,R_all(ind,2)-R_all(ind,1),Gamma(R_all(ind,1):R_all(ind,2),:),options_mt);
    fitmt_subj{ind}.state = rmfield(fitmt_subj{ind}.state,'ipsd');
    fitmt_subj{ind}.state = rmfield(fitmt_subj{ind}.state,'pcoh');
    fitmt_subj{ind}.state = rmfield(fitmt_subj{ind}.state,'phase');
    disp(['Subject ' num2str(ind)])
end


%% Temporal features

scan_T = T_all; % Indexing individual scan sessions
subj_T = scan_T; % Indexing individal subjects
% Compute temporal stats
% Fractional Occupancy is the proportion of time spent in each state
FO = getFractionalOccupancy( Gamma, subj_T, 2);
% Interval Time is the time between subsequent visits to a state
IT = getStateIntervalTimes( Gamma, subj_T, []);
ITmerged = cellfun(@mean,IT);
% Life Times (or Dwell Times) is the duration of visits to a state
LT = getStateLifeTimes( Gamma, subj_T, []);
LTmerged = cellfun(@mean,LT); 

