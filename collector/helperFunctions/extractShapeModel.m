function [mu Sigma h] = extractShapeModel(files,collector)
% extractShapeModel - fetches the ground truth for all entries in files; calculates the mean shape and its covariance matrix
%
% Syntax:
%   [mu Sigma h] = extractShapeModel(files,collector)
%
% Inputs:
%   files     - [struct] files to fetch ground truth from
%   collector - [struct] variables that control the collection of training and test data; see setCollectorsDefault for their description
%     .numRegionsPerVolume
%     .labelIDs
%     .EdgesTrain
%     .columnsShape
%
% Outputs:
%   mu     - [array] mean of h
%   Sigma  - [matrix] covariance matrix of h
%   h      - [matrix] containts all shapes from files
%
% Calls: loadLabels
% See also: trainShape, setCollectorDefaults

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 06-Dec-2013

numColumnsShape = cellfun('length',collector.options.columnsShape);
numColPrev = 0;

% initialize data matrix
h = zeros(length(collector.options.EdgesTrain)*sum(numColumnsShape),length(files));

for i = 1:length(files)
	for regionVolume = 1:collector.options.numRegionsPerVolume
		idx = collector.options.columnsShape{regionVolume};
	    [a filename] = fileparts(files(i).name);
		% sets the ID for the current file (important for volumes)
	    collector.options.labelID = collector.options.labelIDs(i,regionVolume);

		interpolation = loadLabels(filename,collector.options);
		interpolation = interpolation(collector.options.EdgesTrain,idx);
	
		idxSave = (1:length(idx)*size(interpolation,1)) + sum(numColumnsShape(1:regionVolume-1))*size(interpolation,1);
		h(idxSave,i) = reshape(interpolation',1,size(interpolation,1)*length(idx));
	end
end

% calculate mean and covariance
mu = mean(h'); Sigma = cov(h');
