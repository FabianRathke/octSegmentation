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
%   options   - [struct]
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
% Last Revision: 10-Dec-2013

% defaults are defined in local function, see below
options = setAppearanceDefaults(options);
numClasses = size(models,3);

output.prediction = cell(length(files),collector.options.numRegionsPerVolume);
output.trueLabels = cell(length(files),collector.options.numRegionsPerVolume);

tic;
for i = 1:length(files)
	for regionVolume = 1:collector.options.numRegionsPerVolume
		% fetch patches
		collector.options.labelID = collector.options.labelIDs(i,regionVolume);
		a = tic; patches = collector.name(files(i),collector.options); timeFetchingPatches = toc(a); 
		if collector.options.printTimings
			fprintf('[Fetched patches]: %.3fs\n',timeFetchingPatches);
		end

		if collector.options.saveAppearanceTerms && collector.options.loadLabels
			output.trueLabels{i,regionVolume} = patches.classID;
		end
		% initialize results with the suitable data type; for GPU or CPU computations
		eval(sprintf('results = zeros(numClasses,size(patches.data,1),%s);',collector.options.dataType));

		% iterate over regions within the current B-scan
		for regionBScan = 1:collector.options.numRegionsPerBScan
			% select the subset of patches beloning to the current region
			idxSet = patches.idx(:,2) >= collector.options.BScanRegions(regionBScan,1) & patches.idx(:,2)<=collector.options.BScanRegions(regionBScan,2);
			if sum(idxSet) > 0
				% fetch the current set of patches
				patchesCurrRegion = double(patches.data(idxSet,:));
				
				% pre-processing on patch-level (for preprocessing on scan-level see loadData.m)
				for i = 1:length(collector.options.preprocessing.patchLevel)
					patchesCurrRegion = collector.options.preprocessing.patchLevel{i}{1}(patchesCurrRegion,collector.options.preprocessing.patchLevel{i},models(regionVolume,regionBScan,:));
				end
		
				% predict appearance terms
				results(:,idxSet) = options.predAppearance(patchesCurrRegion,models(regionVolume,regionBScan,:));
			end
		end

		% post-processing (disciminative/generative, logarithmic/exponential output)
		if options.normalizeDataTerm
			if ~options.logOutput
				output.prediction{i,regionVolume} = double(exp(results)./repmat(sum(exp(results)),numClasses,1));
			else
				results = exp(results)./repmat(sum(exp(results)),numClasses,1);
				results(results==0) = realmin;
				output.prediction{i,regionVolume} = double(log(results));
			end
		elseif ~options.normalizeDataTerm && ~options.logOutput
			output.prediction{i,regionVolume} = double(exp(results));
		else
			% for numerical issues
			output.prediction{i,regionVolume} = double(results-mean(max(results)));
		end
	end
end

if collector.options.printTimings
	fprintf('[Predicted appearance]: %.3fs\n',toc-timeFetchingPatches);
end

end

% set defaults
function options = setAppearanceDefaults(options)
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

