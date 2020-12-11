clear all
clc
close all


outdir = '...'; % <- edit this line


TOOLBOXPATH = ['...'];
addpath((fullfile(TOOLBOXPATH,filesep, 'fieldtrip')));
ft_defaults % initialize fieldtrip

%%
scriptdir = '...'; % <- edit this line
addpath( fullfile(scriptdir,'toolboxes','osl','osl-core') );
osl_startup

osldir='.../toolboxes/osl/';
p = parcellation(fullfile(osldir,'parcellations','fmri_d100_parcellation_with_3PCC_ips_reduced_2mm_ss5mm_ds8mm_adj.nii.gz')); 
anatomical.pos=p.template_coordinates;
fs=250;

%%
load(fullfile(outdir,[patientname '_sourcepow_lcmv']));

pow=NaN(prod(source_lcmv.dim),size(source_lcmv.avg.momint,2));
pow(find(source_lcmv.inside==1),:)=source_lcmv.avg.momint;

source_lcmv.avg.pow=pow; 
source_lcmv.time=[0:diff(source_lcmv.time([1,2])):size(pow,2)/fs];
source_lcmv.time(end)=[];
functional=source_lcmv;

cfg = [];
cfg.parameter = 'pow';
source2 = ft_sourceinterpolate(cfg, functional, anatomical);
sourcedata=source2.pow;

spatialBasis=p.parcelflag;
voxelData=sourcedata;
nParcels = ROInets.cols(spatialBasis);
goodSamples=[1:size(voxelData,2)];
top5pcInd = abs(spatialBasis) >= repmat(prctile(abs(spatialBasis), 95),[ROInets.rows(spatialBasis), 1]);
for iParcel = nParcels:-1:1,
    mapSign(iParcel) = sign(mean(spatialBasis(top5pcInd(:,iParcel), iParcel)));
end%for
scaledSpatialMaps = ROInets.scale_cols(spatialBasis, mapSign ./max(max(abs(spatialBasis), [], 1), eps));

temporalSTD = max(std(voxelData, [], 2), eps);

voxelWeightings = zeros(size(spatialBasis)); % Preallocate the nvoxels x nparcels weight matrix      

% find time-course for each spatial basis map
for iParcel = nParcels:-1:1, % allocate memory on the fly
    thisMap     = scaledSpatialMaps(:, iParcel);
    parcelMask  = logical(thisMap);

    weightedTS  = voxelData(find(parcelMask),:,:); %#ok Can't use logical indexing
    weightedTS  = reshape(weightedTS,[size(weightedTS,1),size(weightedTS,2)*size(weightedTS,3)]); % reshape to handle trials
    weightedTS  = weightedTS(:,goodSamples);
    weightedTS  = ROInets.scale_rows(weightedTS, thisMap(parcelMask));

    [U, S, V]   = ROInets.fast_svds(weightedTS, 1);
    clear weightedTS

    PCAscores   = S * V';
    maskThresh  = 0.5; % 0.5 is a decent arbitrary threshold chosen by Steve Smith and MJ after playing with various maps.
    thisMask    = thisMap(parcelMask) > maskThresh;   

    if any(thisMask), % the mask is non-zero
        relativeWeighting = abs(U(thisMask)) ./sum(abs(U(thisMask)));
        TSsign  = sign(mean(U(thisMask)));
        TSscale = dot(relativeWeighting, temporalSTD(thisMask));       
        nodeTS  = TSsign .*                               ...
                  (TSscale / max(std(PCAscores), eps)) .* ...      
                  PCAscores;

        % for Mark: this is the linear operator which is applied to
        % the voxel data to get nodeTS.
        voxelWeightings(parcelMask,iParcel) = TSsign .* (TSscale / max(std(PCAscores), eps)).* (U' .* thisMap(parcelMask)');

    else
        warning([mfilename ':EmptySpatialComponentMask'],          ...
                ['%s: When calculating ROI time-courses, ',        ...
                 'an empty spatial component mask was found for ', ...
                 'component %d. \n',                               ...
                 'The ROI will have a flat zero time-course. \n',  ...
                 'Check this does not cause further problems ',    ...
                 'with the analysis. \n'],                         ...
                 mfilename, iParcel);

        nodeTS = zeros(1, ROInets.cols(weightedTS));
        voxelWeightings(~thisMask, iParcel) = zeros(length(thisMask), 1);
    end%if
    nodeData(iParcel, :) = nodeTS;
end
paceldata=nodeData;
source2=[];sourcedata=[];nodeData=[];iParcel=[];

save(fullfile(outdir,[patientname '_paceldata_lcmv']),'paceldata','-v7.3');





