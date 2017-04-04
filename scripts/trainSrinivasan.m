clear all
% parameters region
numCols = 500;
cliprange = [135 134+numCols];
% define three different regions inside each B-Scan, for which seperate appearance models are trained
bscanregions = [1 200 300]';

% we define 17 regions across the volume, to reduce the amount of ground truth required
% for each region, we have the same appearance and shape model
% the central scan of the volume has its own region, most B-Scans are pooled in groups of four
regions = [0:4:28 30 31 33:4:57; 3:4:27 29 30 32 36:4:60];
regionsSelect = [2:16]; % the volumes provided by Srinivasan et al. cover slightly smaller 
BScanPositions = (mean(regions(:,regionsSelect)+1)-31)*0.1191; % geometry parameter, defining distance to central scan around the fovea in mm (can be used to fit volumes with different scan geometry)
sizeY = [6.1 6.0 6.1 6.1 5.7 5.8 5.6 5.9 5.9 6.0 6.3 6.2 5.8 6.2 6.0];

folderData = [getenv('OCT_DOCUMENTS') '/Datasets/OCT/Internet/Srinivasan/']; % set this to the folder where you stored the dataset of Srinivasan et al.
folderLabels = [getenv('OCT_CODE_DIR'),'/datafiles/modelFiles/labelsSrinivasan/'];
filterLabels = 1;

% All healthy volumes
for i = 1:15
	files(i).name = ['NORMAL' num2str(i)];
end
%[files labelIDs] = getLabelIDs(files,folderLabels,regions,regionsSelect);
regions = regions(:,regionsSelect);
labelIDs = zeros(length(files),size(regions,2));
for i = 1:length(files)
    [a filename] = fileparts(files(i).name);
    labelFiles = dir(sprintf('%s/%s_*_coordinates.mat',folderLabels,filename));
    filesBScans = dir([folderData filename '/TIFFs/8bitTIFFs/*.tif']);
    numBScans = length(filesBScans);
	dist =  sizeY(i)/(numBScans-1)*((1:numBScans) - 1 - (numBScans-1)/2);
	[~,region] = min(abs(bsxfun(@minus,BScanPositions,dist'))');

    for j = 1:length(labelFiles)
        [a filename] = fileparts(labelFiles(j).name);
        [a b c d e] = regexp(filename,'.*_(\d){1,2}_coordinates');
        id = str2num(e{1}{1})+1;
		% map BScan ID to 61 BScans coordinates
		regionMap = region(id); 
        labelIDs(i,regionMap) = id;
    end
end


%% ********** PARAMETERS *********
parameters = struct();
collector.options = struct('clip',1,'clipRange',cliprange,'labelIDs',labelIDs,'calcOnGPU',0,'printTimings',1,'BScanPositions',BScanPositions,'loadRoutineData','Srinivasan','loadRoutineLabels','LabelsFromLabelingTool');
collector.options = setCollectorDefaults(collector.options,parameters,files,folderData,folderLabels);
collectorTrn.name = @collectTrnData; collectorTrn.options = collector.options;

%% ******** TRAINING *********
trainFunc.name = @trainModels;
trainFunc.options.appearanceModel{1} = struct('centerPatches',1,'priorVolumesPaper',1,'BScanRegions',bscanregions);
trainFunc.options.appearanceModel{1}.preprocessing.patchLevel = {{@projToEigenspace,20}};
trainFunc.options.shapeModel = struct('numModes',12);

model = trainFunc.name(files,collectorTrn,parameters,trainFunc.options);
save([getenv('OCT_CODE_DIR') '/datafiles/model3D_Srinivasan.mat'],'model');
