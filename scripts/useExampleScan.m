% set scriptFolder to the folder 'script' inside the package
scriptFolder = '/path/to/script/Folder/';
exampleFile = dir([scriptFolder 'circularScan.mat']);

% set collector defaults
collector.options = setCollectorDefaults(struct(),struct(),exampleFile,scriptFolder,scriptFolder);

% load the trained model provided with the package
load([scriptFolder 'circularScanModel']);
% calculate transition matrices
model.shapeModel = preCalcTransitionMatrices(collector,model.shapeModel);

% segment the example scan
collectorTest.name = @collectTestData; collectorTest.options = collector.options;
prediction = predVariational(exampleFile,collectorTest,struct(),model,struct());

% evaluate the prediction, calculate unsigned and signed error measures
results = signedUnsigned(prediction,struct());
str = [sprintf('The unsigned error for each layer in mum is:\n') sprintf('%d: %.3f\n',[(1:9); results.unsigned_per_layer'*3.87]) sprintf('Avg: %.3f',results.unsigned*3.87)];
disp(str);

% load the example scan for plotting
data = load([scriptFolder 'circularScan.mat'],'B0'); B0 = sqrt(sqrt(data.B0));
% plot the scan and its segmentation, the ground truth and the scan itself
plotBScan(B0,prediction.prediction{1},prediction.columnsPred,[]);
plotBScan(B0,prediction.trueLabels{1},[1:768],[]);
figure; imagesc(B0); colormap gray;
