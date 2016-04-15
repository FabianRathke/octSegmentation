function model = preCalcTransitionMatrices(collector,model,eps,idxPredict)
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
% Last Revision: 04-Feb-2016

% some definitions
preCalcMat = tic;
numBounds = length(collector.options.EdgesTrain);
if collector.options.full3D
	model.pTransV = cell(1,collector.options.numRegionsPerVolume);
end
numColumnsShape = cellfun('length',collector.options.columnsShape);

% for each pair of neighboring boundaries within one colume, calculate transition matrices from the learned shape prior
% for a fully 3D model, we have one shape prior over all B-Scans
if collector.options.full3D
	for region = 1:collector.options.numRegionsPerVolume
		numColumns = length(collector.options.columnsShape{region});
		pTransTmp = cell(numColumns,numBounds);
		for i = 1:numColumns
			for j = 2:numBounds
				idx_a = i + (j-2)*numColumns + numBounds*sum(numColumnsShape(1:region-1));
				idx_b = idx_a + numColumns;
				P = inv(model.WML([idx_a idx_b],:)*model.WML([idx_a idx_b],:)' + eye(2)*model.sigmaML);
%				pTransTmp{i,j} = sparse(getCondTransMatrix([model.mu(idx_a) model.mu(idx_b)]',P,collector.options.Y));
				[iS jS sS numElements] = getCondTransMatrixC([model.mu(idx_a) model.mu(idx_b)]',P,int32(collector.options.Y),eps);
				pTransTmp{i,j} = sparse(iS(1:numElements),jS(1:numElements),sS(1:numElements),collector.options.Y,collector.options.Y);
			end
		end
		model.pTransV{region} = pTransTmp;
	end
% for 3D with 2D shape prior, we calculate a seperate model for each 2-D B-Scan
else
    for region = 1:collector.options.numRegionsPerVolume
        numColumns = length(collector.options.columnsShape{region});
        pTransTmp = cell(numColumns,numBounds);
        for i = 1:numColumns
			if nargin == 4
				idxPredictLocal = find(idxPredict(:,i)~=0);
				numBounds = length(idxPredictLocal);
			end
            for j = 2:numBounds
				if nargin == 4
					idx_a = i + (idxPredictLocal(j-1)-1)*numColumns;
					idx_b = i + (idxPredictLocal(j)-1)*numColumns;
				else
					idx_a = i + (j-2)*numColumns;
                	idx_b = idx_a + numColumns;
				end
                P = inv(model(region).WML([idx_a idx_b],:)*model(region).WML([idx_a idx_b],:)' + eye(2)*model(region).sigmaML);
                %pTransTmp{i,j} = sparse(getCondTransMatrix([model(region).mu(idx_a) model(region).mu(idx_b)]',P,collector.options.Y));
				[iS jS sS numElements] = getCondTransMatrixC([model(region).mu(idx_a) model(region).mu(idx_b)]',P,int32(collector.options.Y),eps);
                pTransTmp{i,j} = sparse(iS(1:numElements),jS(1:numElements),sS(1:numElements),collector.options.Y,collector.options.Y);
            end
        end
       	model(region).pTransV{1} = pTransTmp;
    end
end

printMessage(sprintf('... calculated transition matrices %.2f s ... \n',toc(preCalcMat)),1,collector.options.verbose);
