% set folders that hold the data
folderData = '/path/to/data/';
folderLabels = '/path/to/ground/truth/';
files = dir([folderData '*.mat']); % fetch all files

% set up collector; we use the default parameters 
% use 'help setCollectorDefaults' for an overview about model parameters
parameters = struct();
% options that define how to collect ground truth and drawing patches from the OCT-scans 
collector.options.labelIDs = repmat(1,length(files),1);
collector.options.loadRoutineData = 'routineData';
collector.options.loadRoutineLabels = 'routineLabels'; 
% set the default values for the remaining parameters 
collector.options = setCollectorDefaults(collector.options,parameters,files,folderData,folderLabels);
% set the collector functions for training and testing
collectorTrn.name = @collectTrnData; collectorTrn.options = collector.options;
collectorTest.name = @collectTestData; collectorTest.options = collector.options;

% set the appearance model; we use the default options for both appearance models and shape prior
trainFunc.options.appearanceModel = struct('appearanceModel',@trainGlasso);
trainFunc.options.shapeModel = struct();
% train apparance and shape models using all files except the last 5 ones
model = trainModels(files(1:end-5),collectorTrn,parameters,trainFunc.options);

% setup predictor for appearance terms
testFunc.options.predAppearance = @predGlasso;
testFunc.options.plotting = 1;
testFunc.options.folderPlots = '/folder/to/store/plots/';
% make predictions for the last 5 files in the set
prediction = predVariational(files(end-4:end),collectorTest,parameters,model,testFunc.options);

% evaluate the results
results = signedUnsigned(prediction,struct());
