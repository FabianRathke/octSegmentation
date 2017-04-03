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
% Last Revision: 28-Mar-2017

numClasses = size(models,3);
numPatches = size(patches,1);
results = zeros(numClasses,numPatches);

for j = 1:numClasses
	patches_tmp = bsxfun(@minus,patches,models(1,1,j).class_mean);
	% if P is a diagonal matrix
	if norm(diag(diag(models(1,1,j).P))-full(models(1,1,j).P)) == 0
		results(j,:) = -(size(patches_tmp,2)/2)*log(2*pi) + 1/2*sum(log(eig(models(1,1,j).P))) - sum(bsxfun(@times,patches_tmp.*patches_tmp,full(diag(models(1,1,j).P))'),2)*0.5;
	else
%		patches_tmp = double(patches_tmp);
		results(j,:) = -(size(patches_tmp,2)/2)*log(2*pi) + 1/2*sum(log(eig(models(1,1,j).P))) - sum(patches_tmp*(full(models(1,1,j).P)).*patches_tmp,2)*0.5;
	end
	if isfield(models(1,1,1),'prior')
		results(j,:) = results(j,:) + models(1,1,j).prior;
	end
end

end
