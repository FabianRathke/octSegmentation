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
% Last Revision: 06-Dec-2013

tstart = tic;
% iterate over the regions of a volume; for 2-D B-Scans numRegionsPerVolume is 1
for regionVolume = 1:collector.options.numRegionsPerVolume
	% selects the IDs for the current region with the volume for all files for the collector
	collector.options.labelIDsCurrent = collector.options.labelIDs(:,regionVolume);
	% fetch training patches from all files
	trainData = collector.name(files,collector.options);

	% train seperate models for each Region
	for regionBScan = 1:collector.options.numRegionsPerBScan
		% find class ids (i.e. classes representing boundary and layer patches)
		ids = unique(trainData(regionBScan).classID);
		% optional projection onto a low-dimensional subspace
		if collector.options.projToEigenspace
			[V D] = eig(cov(trainData(regionBScan).data));
			W = V(:,end-collector.options.numModesAppearance+1:end);
			mu = mean(trainData(regionBScan).data,1);
			trainData(regionBScan).data = (trainData(regionBScan).data-repmat(mu,size(trainData(regionBScan).data,1),1))*W;
		end
		
		% for each class train an appearance model
		for k = 1:length(ids)
			M = squeeze(trainData(regionBScan).data(trainData(regionBScan).classID==ids(k),:));
			% call the appearance model set be the user			
			appearanceModel = options.appearanceModel(M,params,options);
			% grab fields and store them in the global models struct
			names = fieldnames(appearanceModel);
			for i = 1:length(names)
				eval(sprintf('models(regionVolume,regionBScan,k).%s = appearanceModel.%s;',names{i},names{i}));
			end
			models(regionVolume,regionBScan,k).ID = ids(i);
		end
		
		if collector.options.projToEigenspace
			models(regionVolume,regionBScan,1).W = W;
			models(regionVolume,regionBScan,1).mu = mu;
		end
	end
end

printMessage(sprintf('... trained appearance models in %.2f s ...\n',toc(tstart)),1,collector.options.verbose);
