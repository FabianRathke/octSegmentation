function model = trainGlasso(trainData,ids,params,options)
% trainGlasso - learns a Gaussian distributions as appearanceModels for each appearance class;
%               estimates are regularized by learning sparse precision matrices; details can be found
%               in the paper "Sparse inverse covariance estimation with the graphical lasso" by 
%               Friedman et al, 2008.
%
% Syntax:
%   model = trainGlasso(data,params,options)
%
% Inputs:
%   trainData - [matrix](numPatches,dimData) patch information from the training set for one class
%   ids - [array] a vector of classIDs indicating texture and bondary classes
%   params - [struct] 
%      .glasso - [float] controls the sparsity of the precision matrix, bigger values correspond to sparser matrices
%   options - [struct] options struct
%      .glasso - [float] same as above
%
% Outputs:
%   model - [struct] structure with one entry for each boundary/appearance; holds parameters of normal distribtions
%      .S          - [matrix] covarince matrix
%      .P          - [matrix] precision matrix (sparse)
%      .class_mean - [array] mean patch
%
% See also: trainAppearance, collectTrnData

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 04-May-2014

options = setDefaultOptions(params,options);

for k = 1:length(ids)
	M = squeeze(trainData.data(trainData.classID==ids(k),:));
	
	mu = mean(M);
	S = cov(M); % empirical covariance matrix

	% calls a C implementation of the glasso algorithm
	[model(k).S model(k).P] = cglasso(double(S),options.glasso,1,[],[]);

	% determine the sparsity of P
	if (length(mu)^2-sum(sum(model(k).P==0)) < length(mu)^2/10)
		model(k).P = sparse(model(k).P);
	end

	eigenvalues = eig(model(k).P);
	% in case of negative eigenvalues, make P positive definite
	if min(eigenvalues)<0
		fprintf('Eigenvalues smaller than zero for model(k) %d! (%.4f)\n',i,min(eigenvalues));
		model(k).P = model(k).P + diag((ones(1,patchSize)*abs(min(eigenvalues))+0.001));
	end

	model(k).class_mean = mu;
end

end

function options = setDefaultOptions(params,options)
	options = checkFields(options,params,0.01,'glasso');
end
