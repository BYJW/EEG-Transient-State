clear all
clc
close all

TOOLBOXPATH = ['...'];
addpath((fullfile(TOOLBOXPATH,filesep, 'fieldtrip')));
ft_defaults % initialize fieldtrip

outdir = '...'; % <- edit this line
mripath='...';
patientname='...';


%%  --------------------------------  process the raw MRI images  ------------------------------------
mri = ft_read_mri(fullfile(mripath,patientname,'firstimage'),'dataformat','dicom');
%%%% choose a coordinate, use acpc
cfg          = [];
cfg.method   = 'interactive';
cfg.coordsys = 'acpc';   
mri          = ft_volumerealign(cfg, mri);   %%%%%%%   marker the landmark, nose left/right ear and positive Z point

%%% volume reslice
cfg            = [];
cfg.resolution = 1
cfg.dim        = [256 256 256];
mri            = ft_volumereslice(cfg, mri);

%%%% the volumesegment sometimes generate problem of overlapping or non-closed tissues, need to check and give correct parameters
cfg           = [];
cfg.output    = {'brain','skull','scalp'};  
% cfg.scalpthreshold =0.05;  %%% if did not extract right scalp
% cfg.scalpsmooth=10;
segmentedmri  = ft_volumesegment(cfg, mri);

%%%%%%  visualize the segment
seg_i = ft_datatype_segmentation(segmentedmri,'segmentationstyle','indexed');
cfg = [];
cfg.funparameter = 'seg';
cfg.funcolormap  = gray(4); % distinct color per tissue
cfg.location     = 'center';
cfg.atlas        = seg_i;
ft_sourceplot(cfg, seg_i);
%%%%%%  prepare mesch
cfg=[];
cfg.spmversion = 'spm12';
cfg.tissue={'brain','skull','scalp'};
cfg.numvertices = [3000 2000 1000];
bnd=ft_prepare_mesh(cfg,segmentedmri);
figure
hold on
ft_plot_mesh(bnd(1), 'edgecolor', 'none', 'facecolor', 'r')
ft_plot_mesh(bnd(2), 'edgecolor', 'none', 'facecolor', 'g')
ft_plot_mesh(bnd(3), 'edgecolor', 'none', 'facecolor', 'b')
alpha 0.3
view(132, 14)

%%%%%  prepare headmodel
cfg        = [];
cfg.spmversion = 'spm12';
cfg.method = 'dipoli';
headmodel  = ft_prepare_headmodel(cfg, bnd);  %%% dipoli

figure;
hold on
ft_plot_mesh(headmodel.bnd(1));
ft_plot_headmodel(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');


%%  --------- using the template_mri provided by fieldtrip to construct template grid -----------
cfg            = [];
cfg.spmversion = 'spm12';
cfg.resolution = 1
cfg.dim        = [256 256 256];
template_mri            = ft_volumereslice(cfg, template_mri);

cfg           = [];
cfg.output    = {'brain','skull','scalp'};  
template_segmentedmri  = ft_volumesegment(cfg, template_mri);

cfg=[];
cfg.spmversion = 'spm12';
cfg.tissue={'brain','skull','scalp'};
cfg.numvertices = [3000 2000 1000];
template_bnd=ft_prepare_mesh(cfg,template_segmentedmri);

cfg        = [];
cfg.spmversion = 'spm12';
cfg.method = 'dipoli';
template_headmodel  = ft_prepare_headmodel(cfg, template_bnd);

cfg = [];
cfg.xgrid  = -20:1:20;
cfg.ygrid  = -20:1:20;
cfg.zgrid  = -20:1:20;
cfg.unit   = 'cm';
cfg.tight  = 'yes';
cfg.inwardshift = -1.5;
cfg.headmodel   = template_headmodel;
template_grid   = ft_prepare_sourcemodel(cfg);
template_grid = ft_convert_units(template_grid,'mm');

%% ------------  using the template grid construct individual grid -----------------
cfg = [];
cfg.warpmni = 'yes';
cfg.nonlinear = 'yes';
cfg.unit = 'mm';
cfg.template = template_grid;
cfg.mri = mri;
cfg.spmversion = 'spm12';   % default is 'spm8'
cfg.spmmethod = 'new'      % default is 'old'
grid = ft_prepare_sourcemodel(cfg);

%%  --------------------------------  align electrodes and sensors -----------------------
%%%%%%  load the original elec
addpath(fullfile(TOOLBOXPATH, 'eeglab'))
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % launch eeglab
dataname='...'; %%% processed resting data
EEG = pop_loadset('filename',[dataname,'.set'],'filepath','...');
DOCdata = eeglab2fieldtrip(EEG,'preprocessing','none');
elec=DOCdata.elec;
%%%  align the electrode with brain ---- a manual way
elec = ft_convert_units(elec,'mm'); % should be the same unit as MRI
cfg = [];
cfg.method    = 'interactive';   
cfg.elec      = elec;
cfg.headshape = headmodel.bnd(1);  %%%%% bnd(1):scalp  bnd(2):skull  bnd(3):cortex
elec = ft_electroderealign(cfg);

cfg = [];
cfg.method = 'project'; % onto scalp surface
cfg.headshape = headmodel.bnd(1); % scalp surface
elec = ft_electroderealign(cfg, elec);

figure;
hold on
ft_plot_mesh(headmodel.bnd(3).pos,'vertexcolor','cortex');
% ft_plot_headmodel(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');
ft_plot_sens(elec,'style', '.k');

%% ------------------------------------- forward model -----------------------------------

headmodel = ft_convert_units(headmodel, 'mm');
grid = ft_convert_units(grid, 'mm');
elec = ft_convert_units(elec,'mm'); % should be the same unit as MRI

cfg = [];
cfg.sourcemodel = grid;   
cfg.headmodel= headmodel;
cfg.elec = elec;
cfg.singleshell.batchsize = 5000; % speeds up the computation
leadfield = ft_prepare_leadfield(cfg);


%%  --------------------  beamforming source analysis --------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
DOCdata.elec=elec;
%%%%
cfg=[];
cfg.reref='yes';
cfg.refchannel='all';
DOCdata = ft_preprocessing(cfg, DOCdata)

cfg = [];
cfg.trials = 1;
hdata = ft_preprocessing(cfg, DOCdata);

numtrail=size(DOCdata.trial);
data_conca=hdata.trial{1};
time_conca=hdata.time{1};
for i=2:numtrail
data_conca=cat(2,data_conca,DOCdata.trial{i});
time_conca=cat(2,time_conca,(DOCdata.time{i}+(i-1)*triallength*ones(1,length(hdata.time{1}))));
end
hdata.trial{1}=data_conca;
hdata.time{1}=time_conca;   
hdata.elec=elec;
%%%%% inverse soluation lcmv
cfg=[];
cfg.method='lcmv';
cfg.grid=leadfield;
cfg.headmodel=headmodel;
cfg.keepfilter = 'yes';
cfg.lcmv.keepfilter='yes';
cfg.elec = hdata.elec;
source=ft_sourceanalysis(cfg, hdata);
%%%% extract time curves of sources
cfg = [];
cfg.projectmom         = 'yes';
source_lcmv         = ft_sourcedescriptives(cfg,source);

Npos = size(source_lcmv.pos,1);
Ntime = length(source_lcmv.time);
source_lcmv.avg.momint = nan(Npos,Ntime);
idx = find(source_lcmv.inside);
for i = 1:length(idx)
    index = idx(i);
    source_lcmv.avg.momint(index,:) = source_lcmv.avg.mom{index};
end

save(fullfile(outdir,[patientname '_sourcepow_lcmv']),'sourcepow_lcmv','-v7.3');





