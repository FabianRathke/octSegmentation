function results = predGlasso(patches,models,options)
% predGlasso - evaluates Gaussian distributions for patches 
%
% Syntax:
%   results = predGlasso(patches,models,options)
%
% Inputs:
%   patches - [matrix] each row holds a patch
%   models  - [struct] Gaussian distributions for each appearance class
%   options - [struct] no options needed for predGlasso
%
% Outputs:
%   results - [matrix] the probabilities for all classes and patches
%
% See also: predAppearance, predGlassoGPU

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 10-Dec-2013

numClasses = size(models,3);
numPatches = size(patches,1);
results = zeros(numClasses,numPatches);

for j = 1:numClasses
	patches_tmp = double(patches) - models(1,1,j).class_mean(ones(1,numPatches),:);
    results(j,:) = -(size(patches_tmp,2)/2)*log(2*pi) + 1/2*sum(log(eig(models(1,1,j).P))) - sum(patches_tmp*(squeeze(models(1,1,j).P)).*patches_tmp,2)/2;
end

end
