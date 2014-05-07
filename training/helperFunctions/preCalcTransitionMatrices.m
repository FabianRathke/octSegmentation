function model = preCalcTransitionMatrices(collector,model)
% preCalcTransitionMatrices - given a PPCA shape model, precalcs transition matrices used during the initialization of the inference approach 
%
% Syntax:
%   model = preCalcTransitionMatrices(collector,model)
% 
% Inputs:
%   collector - [struct] options struct, see setCollectorDefaults for details
%      .options.numRegionsPerVolume
%      .options.columnsShape
%      .options.EdgesTrain
%      .options.Y
%   model     - [struct] model struct, holds the PPCA shape model
%
% Outputs:
%   model - [struct] the shape model struct with a new field concatenated
%     .pTransV - [cell-array] for each 
%
% See also: trainShape, getCondTransMatrix, setCollectorDefaults
% Calls: getCondTransMatrix

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 02-May-2013

% some definitions
preCalcMat = tic;
numBounds = length(collector.options.EdgesTrain);
model.pTransV = cell(1,collector.options.numRegionsPerVolume);
numColumnsShape = cellfun('length',collector.options.columnsShape);
% for each pair of neighboring boundaries within one colume, calculate transition matrices from the learned shape prior
for region = 1:collector.options.numRegionsPerVolume
	numColumns = length(collector.options.columnsShape{region});
	pTransTmp = cell(numColumns,numBounds);
	for i = 1:numColumns
		for j = 2:numBounds
			idx_a = i + (j-2)*numColumns + numBounds*sum(numColumnsShape(1:region-1));
			idx_b = idx_a + numColumns;
			P = inv(model.WML([idx_a idx_b],:)*model.WML([idx_a idx_b],:)' + eye(2)*model.sigmaML);
			pTransTmp{i,j} = sparse(getCondTransMatrix([model.mu(idx_a) model.mu(idx_b)]',P,collector.options.Y));
		end
	end
	model.pTransV{region} = pTransTmp;
end

%printMessage(sprintf('... calculated transition matrices %.2f s ... \n',toc(preCalcMat)),1,collector.options.verbose);
