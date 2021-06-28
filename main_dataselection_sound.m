clear all
close all


%%

% The script corresponding to the data selection used in manuscript 'Spontaneous transient brain states in EEG source space of disorders of consciousness'
% The details could be found in the supplementary information and Fig. S1
% variable 'data' is defined as a structure including resting-state data from 62 patients (column)
% an EEGLAB toolbox is needed, which can be accessed in https://sccn.ucsd.edu/eeglab/index.php
% the script used an algorithm named 'SOUND', details in Mutanen, T. P., Metsomaa, J., Liljander, S. & Ilmoniemi, R. J. 2018. Automatic and robust noise suppression in EEG and MEG: The SOUND algorithm. Neuroimage, 166, 135-151.
% signal quality indexes of patients was fitted by two guassian distribution (gauss_model_2fit) to determine the threshold discarding outliers


%  edit by Yang Bai 2021-06-22

%%

TOOLBOXPATH = ['...'];

datapath='...';  % <- edit this line
load(fullfile(datapath,'Groupdata'));

addpath(fullfile(TOOLBOXPATH, 'eeglab'))
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % launch eeglab


%%  running the sound algorithm for signal quaility eatimate

rf=0.1;
lambda=[];
i_ref=0;
max_iter=100;
LFM_sphere=[];
for numsub=1:size(data,2)
cleandata=data{1,numsub};
EEG=pop_importdata('dataformat','array','nbchan',62,'data','cleandata','setname','healthy','srate',250,'pnts',2000,'xmin',0);
EEG=pop_chanedit(EEG, 'lookup',fullfile(datapath,'BP62.locs'),'settype',{'' 'EEG'},'load',{fullfile(indir,'BP62.locs') 'filetype' 'autodetect'});

if isempty(LFM_sphere)
% Build the spherical lead field, using the theoretical electrode-locations
[LFM_sphere] = construct_spherical_lead_field(EEG);
end

[~, sigmas_S(numsub,:), lambda, x, i_ref, k] = SOUND(EEG.data, LFM_sphere,rf, lambda, i_ref, max_iter);

end


%% using gaussian model for data fitting
data=sigmas_S;
plotflg=1;  %%% plot the histogram and models
[ res ] = gauss_model_2fit( data , plotflag )
lowqualitysub=find(data>res.orig_th); %%% find out the subjects with lower signal quality

