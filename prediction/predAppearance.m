function output = predAppearance(files,collector,params,models,options)
% predAppearance - fetches patches from test files, and predicts likelihood terms p(y|c) or p(c|y) using a user selected appearance model
%
% Syntax:
%   output = predAppearance(files,collector,params,models,options) 
%
% Inputs:
%   files     - [struct] files to fetch ground truth from
%   collector - [struct] variables that control the collection of training and test data; see setCollectorsDefault for a description of the variables
%     .storeLabels
%     .columnsPred 
%     .labelIDs
%     .printTimings
%     .numRegionsPerBScan
%     .dataType
%     .projToEigenspace
%   params    - [struct]
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
% TODO: Seperate pre-processing in seperate functions; controlled by the user

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 10-Dec-2013

options = setDefaultAppearance(options);
numClasses = length(models);

output.prediction = cell(length(files),collector.options.numRegionsPerVolume);
output.trueLabels = cell(length(files),collector.options.numRegionsPerVolume);

tic;
for i = 1:length(files)
	for regionVolume = 1:collector.options.numRegionsPerVolume
		% fetch patches
		collector.options.labelID = collector.options.labelIDs(i,regionVolume);
		a = tic; patches = collector.name(files(i),collector.options); timeFetchingPatches = toc(a); 
		if collector.options.printTimings
			fprintf('Fetching patches in %.2f\ns',timeFetchingPatches);
		end

		if collector.options.saveAppearanceTerms
			output.trueLabels{i,regionVolume} = patches.classID;
		end
		% initialize results with the suitable data type; for GPU or CPU computations
		eval(sprintf('results = zeros(numClasses,size(patches.data,1),%s);',collector.options.dataType));

		% iterate over regions within the current B-scan
		for regionBScan = 1:collector.options.numRegionsPerBScan
			% select the subset of patches beloning to the current region
			idxSet = find(collector.options.columnsPred>=collector.options.BScanRegions(regionBScan,1) & collector.options.columnsPred<=collector.options.BScanRegions(regionBScan,2)); 
			idxStart = 1 + (idxSet(1)-1)*collector.options.Y;
			idxEnd = idxSet(end)*collector.options.Y;
			numPatches = idxEnd-idxStart+1;
			% fetch current patches
			patchesCurrRegion = double(patches.data(idxStart:idxEnd,:));
	
			% pre-processing
			if collector.options.projToEigenspace
				patchesCurrRegion = (patchesCurrRegion - models(regionVolume,regionBScan,1).mu(1,1))*models(regionVolume,regionBScan,1).W;
			end

			% predict appearance terms
			results(:,idxStart:idxEnd) = options.predAppearance(patchesCurrRegion,models(regionVolume,regionBScan,:));
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
	fprintf('Predicted appearance in %.3fs\n',toc-timeFetchingPatches);
end

end

% set defaults
function options = setDefaultAppearance(options)
if (~isfield(options,'normalizeDataTerm'))
    options.normalizeDataTerm = true;
end

if (~isfield(options,'logOutput'))
    options.logOutput = false;
end

end

