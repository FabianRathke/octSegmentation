function results = predGlassoGPU(patches,models,options)
% predGlassoGPU - GPU implementation of predGlasso; type 'help predGlasso' for more information
%
% Syntax:
%   results = predGlassoGPU(patches,models,options)

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 10-Dec-2013

numClasses = size(models,3);
numPatches = size(patches,1);
results = zeros(numClasses,numPatches,GPUsingle);

% init variables
dataGPU = GPUsingle(patches);
meanGPU = zeros(size(models(1).class_mean),GPUsingle);
PGPU = zeros(size(models(1).P),GPUsingle);

dataGPUTmp = zeros([numPatches,size(meanGPU,2)],GPUsingle);
				
for j = 1:numClasses
	meanGPU = GPUsingle(models(1,1,j).class_mean);
	eigenvaluesGPU = GPUsingle(models(1,1,j).eigenvalues);
	assign(0,dataGPUTmp, meanGPU, {ones(1,numPatches)},[1,1,size(meanGPU,2)]);
	dataGPUTmp = dataGPU - dataGPUTmp;
	PGPU = GPUsingle(full(models(1,1,j).P));
	results(j,:) = -(size(patches,2)/2)*log(2*pi) + 1/2*sum(log(models(1,1,j).eigenvalues)) - sum(dataGPUTmp*PGPU.*dataGPUTmp,2)/2;
end
