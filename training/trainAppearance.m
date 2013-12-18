function models = trainAppearance(files,collector,params,options)
% trainAppearance - grabs training data using the collector defined in collector.name; calls the appearance model defined in options.appearanceModel
% 
% Syntax:
%   models = trainAppearance(files,collector,params,options)
%
% Inputs:
%   files     - [struct] output of the matlab function 'dir'
%   collector - [struct] holds the options for the collector grabbing the training patches, description given in setCollectorDefaults
%     .labelIDs
%     .numRegionsPerVolume
%     .numRegionsPerBScan
%     .projToEigenspace
%     .numModesAppearance
%   params    - [struct] holds params used for example in cross-validation 
%   options   - [struct] holds options used within the appearance function
%     .appearanceModel - [string] name of the function to call
%
% Outputs:
%   models - [struct] contains the appearance models, 3-dim struct
%
% See also: trainGlasso, trainShape

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 16-Dec-2013

tstart = tic;
models = struct();
% iterate over the regions of a volume; for 2-D B-Scans numRegionsPerVolume is 1
for regionVolume = 1:collector.options.numRegionsPerVolume
	% selects the IDs for the current region with the volume for all files for the collector
	collector.options.labelIDsCurrent = collector.options.labelIDs(:,regionVolume);
	% fetch training patches from all files
	trainData = collector.name(files,collector.options);

	% train seperate models for each Region
	for regionBScan = 1:collector.options.numRegionsPerBScan
		% pre-processing of the patches
		for i = 1:length(collector.options.preprocessing.patchLevel)
			[trainData(regionBScan).data appearanceModel] = collector.options.preprocessing.patchLevel{i}{1}(trainData(regionBScan),collector.options.preprocessing.patchLevel{i});
		end
		models = appendToModel(models,appearanceModel,regionVolume,regionBScan);

		% find class ids (i.e. classes representing boundary and layer patches)
		ids = unique(trainData(regionBScan).classID);

		% call the appearance model set be the user			
		appearanceModel = options.appearanceModel(trainData(regionBScan),ids,params,options);
	    models = appendToModel(models,appearanceModel,regionVolume,regionBScan);
	end
end
models(1).numClasses = length(ids);
printMessage(sprintf('... trained appearance models in %.2f s ...\n',toc(tstart)),1,collector.options.verbose);
end

function models = appendToModel(models,appearanceModel,regionVolume,regionBScan)

% grab fields and store them in the global models struct
names = fieldnames(appearanceModel);
for i = 1:length(names)
	for j = 1:length(appearanceModel)
		eval(sprintf('models(regionVolume,regionBScan,j).%s = appearanceModel(j).%s;',names{i},names{i}));
	end
end

end
