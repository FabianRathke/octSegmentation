function output = predAppearance(files,collector,models,options)
% predAppearance - fetches patches from test files, and predicts likelihood terms p(y|c) or p(c|y) using a user selected appearance model
%
% Syntax:
%   output = predAppearance(files,collector,models,options) 
%
% Inputs:
%   files     - [struct] files to fetch ground truth from
%   collector - [struct] variables that control the collection of training and test data; see setCollectorsDefault for a description of the variables
%     .options.columnsPred 
%     .options.labelIDs
%     .options.printTimings
%     .options.numRegionsPerBScan
%     .options.numRegionsPerVolume
%     .options.dataType
%     .options.preprocessing.patchLevel
%     .options.BScanRegions
%   models    - [struct] appearance models for all files and BScans; have to be compatible with the function called in 
%   options   - [struct] post-processing related to the inference approach
%     .normalizeDataTerm - [boolean] converts generative to discriminative terms (i.e. normalizes the sum of probabilities to 1 for each pixel)
%     .logOutput         - [boolean] outputs log-probabilities
%     .predAppearance    - [function handle] must be compatible to the appearance models
%
% Outputs:
%   output - [struct]
%     .trueLabels - [cell] pixel-wise labels for all B-Scans
%     .prediction - [cell] pixel-wise probabilities for all classes and B-Scans
%
% See also: predGlasso, predVariational

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 24-Nov-2016

% defaults are defined in local function, see below
options = setDefaults(options);

numClasses = collector.options.numLayers*2-1;

output.prediction = cell(length(files),collector.options.numRegionsPerVolume);
output.trueLabels = cell(length(files),collector.options.numRegionsPerVolume);

tic;

for i = 1:length(files)
	for regionVolume = 1:collector.options.numRegionsPerVolume
		% initialize results with the suitable data type; for GPU or CPU computations
		results = zeros(numClasses,length(collector.options.columnsPred)*collector.options.Y);
		if length(models) > 1
			resultsFirstModel = zeros(numClasses,length(collector.options.columnsPred)*collector.options.Y);
		end
		mixturePrior = ones(1,length(models))/length(models);
		% mixture of appearance models
		for model = 1:length(models)
			optionsApp = setAppearanceDefaults(options.appearanceModel{model},struct(),files,collector);
			fn = fieldnames(optionsApp);
			for j = 1:length(fn)
				eval(sprintf('collector.options.%s = optionsApp.%s;',fn{j},fn{j}));
			end

			% fetch patches
			collector.options.labelID = collector.options.labelIDs(i,regionVolume);
			a = tic; patches = collector.name(files(i),collector.options); timeFetchingPatches = toc(a); 
			if collector.options.printTimings
				fprintf('[Fetched patches]: %.3fs\n',timeFetchingPatches);
			end

			if collector.options.saveAppearanceTerms && collector.options.loadLabels
				output.trueLabels{i,regionVolume} = patches.classID;
			end

			% iterate over regions within the current B-scan
			for regionBScan = 1:collector.options.numRegionsPerBScan
				% select the subset of patches beloning to the current region
				idxSet = find(patches.idx(:,2) >= collector.options.BScanRegions(regionBScan,1) & patches.idx(:,2)<=collector.options.BScanRegions(regionBScan,2));
				if length(idxSet) > 0
					if (length(idxSet)*size(patches.data,2)*4/1024^3 > 1); nSplits = 10; else; nSplits = 1; end;
					for j = 1:nSplits
						idxSub = (1:length(idxSet)/nSplits) + (j-1)*length(idxSet)/nSplits;
						% fetch the current set of patches
						patchesCurrRegion = double(patches.data(idxSet(idxSub),:));
						
						% pre-processing on patch-level (for preprocessing on scan-level see loadData.m)
						for l = 1:length(collector.options.preprocessing.patchLevel)
							if  (size(optionsApp.preprocessing.patchLevel{l},1) == 1)
								patchesCurrRegion = optionsApp.preprocessing.patchLevel{l}{1}(patchesCurrRegion,optionsApp.preprocessing.patchLevel{l},models{model}(regionVolume,regionBScan,:));
							elseif (size(optionsApp.preprocessing.patchLevel{l},1) == collector.options.numRegionsPerBScan)
								 patchesCurrRegion = optionsApp.preprocessing.patchLevel{l}{l}{1}(patchesCurrRegion,optionsApp.preprocessing.patchLevel{l}{l},models{model}(regionVolume,regionBScan,:));
							else
								error(sprintf('Wrong number of options for preprosessing: %s\n',optionsApp.preprocessing.patchLevel{l}{1}));
							end
						end
				
						% predict appearance terms
						tmp = exp(options.predAppearance(patchesCurrRegion,models{model}(regionVolume,regionBScan,:)));
						results(:,idxSet(idxSub)) = results(:,idxSet(idxSub)) + mixturePrior(model)*tmp;
						if length(models) > 1 && model == 1
							resultsFirstModel(:,idxSet(idxSub)) = tmp;
						end
					end
				end
			end
		end
		output.predictionNonNormalized = double((results));
		if length(models) > 1 
			output.predictionFirstModel = double(resultsFirstModel./repmat(sum(resultsFirstModel),numClasses,1));
		end
		% post-processing (disciminative/generative, logarithmic/exponential output)
		if options.normalizeDataTerm
%			if ~options.logOutput
				output.prediction{i,regionVolume} = double(results./repmat(sum(results),numClasses,1));
%			else
%				results = exp(results)./repmat(sum(exp(results)),numClasses,1);
%				results(results==0) = realmin;
%				output.prediction{i,regionVolume} = double(log(results));
%			end
%		elseif ~options.normalizeDataTerm && ~options.logOutput
		else
			output.prediction{i,regionVolume} = double(exp(results));
%		else
%			% for numerical issues
%			output.prediction{i,regionVolume} = double(results-mean(max(results)));
		end
	end
end

if collector.options.printTimings
	fprintf('[Predicted appearance]: %.3fs\n',toc-timeFetchingPatches);
end

end

% set defaults
function options = setDefaults(options)
if (~isfield(options,'normalizeDataTerm'))
    options.normalizeDataTerm = true;
end

if (~isfield(options,'logOutput'))
    options.logOutput = false;
end

if (~isfield(options,'predAppearance'))
	options.predAppearance = @predGlasso;
end

end

