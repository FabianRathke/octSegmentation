function model = trainGlasso(data,params,options)
% trainGlasso - learns a Gaussian distribution as a appearanceModel; gets patches for one class
%
% Syntax:
%   model = trainGlasso(data,params,options)
%
% Inputs:
%   data    - [matrix](numPatches,dimData) patch information from the training set for one class
%   params  - [struct] 
%     .glasso - [float] controls the sparsity of the precision matrix, bigger values correspond to sparser matrices
%   options - [struct] options struct; can be used to control the behaviour of this function
%
% Outputs:
%   model - [struct] 
%
% See also: trainAppearance

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 06-Dec-2013

mu = mean(data);
S = cov(data); % empirical covariance matrix

% calls a C implementation of the glasso algorithm
[model.S model.P] = cglasso(S,params.glasso,1,[],[]);

% determine the sparsity of P
if (length(mu)^2-sum(sum(model.P==0)) < length(mu)^2/10)
	model.P = sparse(model.P);
end

eigenvalues = eig(model.P);
% in case of negative eigenvalues, make P positive definite
if min(eigenvalues)<0
	fprintf('Eigenvalues smaller than zero for model %d! (%.4f)\n',i,min(eigenvalues));
	model.P = model.P + diag((ones(1,patchSize)*abs(min(eigenvalues))+0.001));
end

model.class_mean = mu;
