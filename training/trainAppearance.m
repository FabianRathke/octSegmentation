function models = trainAppearance(files,collector,params,options)
% trainAppearance - grabs training data using the collector defined in collector.name; calls the appearance model defined in options.appearanceModel
% 
% Syntax:
%   models = trainAppearance(files,collector,params,options)
%
% Inputs:
%   files     - [struct] output of the matlab function 'dir'
%   collector - [struct] holds the options for the collector grabbing the training patches, description given in setCollectorDefaults
%     .options.labelIDs
%     .options.numRegionsPerVolume
%     .options.numRegionsPerBScan
%     .options.preprocessing
%   params    - [struct] holds params used for example in cross-validation 
%   options   - [struct] holds options used within the appearance function
%     .appearanceModel   - [string] name of the function to call
%     .priorVolumesPaper - [boolean] assigns a specific non-uniform prior distribution to the appearance models; see documentation.pdf and the source code for details
%
% Outputs:
%   models - [struct] contains the appearance models, 3-dim struct
%
% See also: setAppearanceDefaults, trainGlasso, trainShape

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 05-May-2014

options = setAppearanceDefaults(options,params,files);

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
		% pre-processing of the patches [calls a arbitrary set of models defined in collector.options.preprocessing.patchLevel]
		for i = 1:length(collector.options.preprocessing.patchLevel)
			[trainData(regionBScan).data appearanceModel] = collector.options.preprocessing.patchLevel{i}{1}(trainData(regionBScan),collector.options.preprocessing.patchLevel{i});
			models = appendToModel(models,appearanceModel,regionVolume,regionBScan);
		end
	
		% find class ids (i.e. classes representing boundary and layer patches)
		ids = unique(trainData(regionBScan).classID);

		% call the appearance model set by the user	
		appearanceModel = options.appearanceModel(trainData(regionBScan),ids,params,options);
	    models = appendToModel(models,appearanceModel,regionVolume,regionBScan);
		
		if options.priorVolumesPaper
			models = setPriors(models,ids,collector,regionVolume,regionBScan);
		end
	end
end
models(1).numClasses = length(ids);
printMessage(sprintf('... trained appearance models in %.2f s ...\n',toc(tstart)),1,collector.options.verbose);

end % end function


function models = appendToModel(models,appearanceModel,regionVolume,regionBScan)
% grab fields and store them in the global models struct
names = fieldnames(appearanceModel);
for i = 1:length(names)
	for j = 1:length(appearanceModel)
		eval(sprintf('models(regionVolume,regionBScan,j).%s = appearanceModel(j).%s;',names{i},names{i}));
	end
end

end

% restores specific non-uniform priors for the appearance models as incorrectly used in our publication for the 3-D dataset
function models = setPriors(models,ids,collector,regionVolume,regionBScan)
for j = 1:length(ids)
	if j <= collector.options.numLayers
		models(regionVolume,regionBScan,j).prior = 1/2*(sum(log(eig(models(regionVolume,regionBScan,1).P))) - sum(log(eig(models(regionVolume,regionBScan,j).P))));
	else
		models(regionVolume,regionBScan,j).prior = 1/2*(sum(log(eig(models(regionVolume,regionBScan,2).P))) - sum(log(eig(models(regionVolume,regionBScan,j).P))));
	end
end

end
