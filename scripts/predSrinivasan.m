% Use the existing model file, trained on our own data to predict ground truth for the healthy volumes on the Srinivasan et al. dataset
% Download it from here: http://people.duke.edu/~sf59/Srinivasan_BOE_2014_dataset.htm 
% Since we are not allowed to publish our own datasets, this is a good way to explain the most important functionality of our model: training and prediction
% We will use the labels obtained here in another script (trainSrinivasan) to train a new model for 3-D volumesnew 

% load model file to segment volumetric data
load([getenv('OCT_CODE_DIR') '/datafiles/modelFiles/model3D']);
% load margins containing confidence measures expectations of the objective function
margins = load([getenv('OCT_CODE_DIR') '/datafiles/modelFiles/healthyMargins3D']);

collector.options = model.params; % use most standard parameters stored in the model file
collector.options.clip = 1; % cut the left and right boundary of each B-Scan to obtain the central part with width 500pxo
collector.options.folder_data = [getenv('OCT_DOCUMENTS') '/Datasets/OCT/Internet/Srinivasan/']; % set this folder to the directory where you unzipped the dataset
collector.options.loadRoutineData = 'Srinivasan'; % we specified a load-routine just for that dataset
collector.options.loadLabels = 0; % since we have no ground thruth; This can be used to print the segmentation error during the segmentation
collector.options.numRegionsPerVolume = 1; % since we have no ground thruth; This can be used to print the segmentation error during the segmentation

% set the folder that contains the data from Srinivasan et al. 
folder = collector.options.folder_data;
subFolders = {'NORMAL','AMD','DME'};

% volume diameter in y-direction in mm (information provided in the paper)
sizeY = [6.1 6.0 6.1 6.1 5.7 5.8 5.6 5.9 5.9 6.0 6.3 6.2 5.8 6.2 6.0; 5.9 5.9 4.6 4.5 4.5 4.5 4.5 4.4 4.4 4.5 4.5 4.5 4.5 4.3 4.5; 5.8 6.4 5.9 6.0 6.0 5.9 7.3 7.3 7.3 7.2 7.4 7.6 7.3 7.2 7.2];
testFunc.name = @predVariational; testFunc.options = struct('calcFuncVal',1,'plotting',0,'printTimings',0);
testFunc.options.appearanceModel = model.options.appearanceModel;

% split "standard" HDE Spectralis scans with 61 BScans into 17 regions that use the same shape- and appearance models
regions = [0:4:28 30 31 33:4:57; 3:4:27 29 30 32 36:4:60];


% only predict the healthy scans 
% (set j = 2 or j = 3 for AMD/DME scans; the model trained on healthy data will not segment them correctly)
% we are currently working on a model extension to handle pathological scans as well
% 
for j = 1:1
	for file = 1:15
		fprintf('******** FILE %d ******** \n',file);
		% files in folder 
		files = dir([folder subFolders{j} num2str(file) '/TIFFs/8bitTIFFs/*.tif']);
		numBScans = length(files);

		testFile.name = [subFolders{j} num2str(file)];
		% first determine the region for each B-Scan and then randomly pick one for each region 
		% divide the total length in y-direction by the number of B-Scans in the volume to obtain the distance between two B-Scans
		dist =  sizeY(j,file)/(numBScans-1)*((1:numBScans) - 1 - (numBScans-1)/2);
		[~,region] = min(abs(bsxfun(@minus,collector.options.BScanPositions,dist'))');

		for i = 1:length(collector.options.BScanPositions)
			% randomly select scan for the current region
			regionSelect = find(region==i);
			if length(regionSelect) > 0
				regionSelect = regionSelect(randperm(length(regionSelect),1));
				fprintf('********* SCAN %d *********\n',regionSelect);
							
				collector.options.labelIDs = regionSelect;
				collector.options.scanID = i;
				collectorTest.name = @collectTestData; collectorTest.options = collector.options;
				% perform the prediction
				testFunc.options.calcFuncVal = 0;
				predictionTmp = testFunc.name(testFile,collectorTest,struct(),model,testFunc.options);
				% store the prediction
%				predictionTmp = rmfield(predictionTmp,'q_c_singleton');
%				prediction(j,file,i) = predictionTmp;
		
				% plot the prediction
				if 0
					collector.options.labelID = regionSelect;
					B0 = loadData(testFile.name,collector.options);
					h = imagesc(B0); ax = gca; ax.ColorOrderIndex = 1;
					colormap gray; hold on; plot(collector.options.columnsPred,predictionTmp.prediction{1}');
					drawnow
					%h = imagesc(B0); colormap gray; hold on; plot(collector.options.columnsPred,predictionTmp.q_b{1}');

					% display confidence estimates (set testFunc.options.calcFuncVal = 1)
					%makeQualityOverlay(predictionTmp.prediction{1},squeeze(predictionTmp.funcVal.q_c_data),margins.mu(4,:),margins.sigma(4,:),collector.options.columnsPred);
				end
				% write label files		
				% add boundaries to bring it into the format [496 768]
				idxInterpolate = 1:500; idxInterpolate = idxInterpolate(~ismember(idxInterpolate,collector.options.columnsPred));

				clear interpolation
				interpolation(:,collector.options.columnsPred) = predictionTmp.prediction{1}; 
				for bound = 1:9
					interpolation(bound,idxInterpolate) = interp1(collector.options.columnsPred,predictionTmp.prediction{1}(bound,:),idxInterpolate);
				end
				interpolation = [zeros(9,134) interpolation zeros(9,134)];
				% map the label into the same range 1 to 61
				save([getenv('OCT_CODE_DIR') '/datafiles/modelFiles/labelsSrinivasan/' testFile.name '_' num2str(regionSelect-1) '_coordinates.mat'],'interpolation');
			end
		end
	end
end
