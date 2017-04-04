fileFolder = [getenv('OCT_CODE_DIR') '/datafiles/'];
exampleFile = dir([fileFolder '/exampleScans/circularScan.mat']);

% load the trained model provided with the package
load([fileFolder '/modelFiles/model2D']);

% set collector defaults
collector.options = model.params;
collector.options.folder_data = [fileFolder '/exampleScans/'];
collector.options.folder_labels = [fileFolder '/exampleScans/'];
collector.options.labelIDs = 1;

% segment the example scan
collectorTest.name = @collectTestData; collectorTest.options = collector.options;
testFunc.options.appearanceModel = model.options.appearanceModel
prediction = predVariational(exampleFile,collectorTest,struct(),model,testFunc.options);

% evaluate the prediction, calculate unsigned and signed error measures
results = signedUnsigned(prediction,struct());
str = [sprintf('The unsigned error for each layer in mum is:\n') sprintf('%d: %.3f\n',[(1:9); results.unsigned_per_layer'*3.87]) sprintf('Avg: %.3f',results.unsigned*3.87)];
disp(str);

% load the example scan for plotting
data = load([collector.options.folder_data 'circularScan.mat'],'B0'); B0 = sqrt(sqrt(data.B0));
% plot the scan and its segmentation, the ground truth and the scan itself
figure; imagesc(B0); colormap gray; hold on; plot(collector.options.columnsPred,prediction.prediction{1}');
