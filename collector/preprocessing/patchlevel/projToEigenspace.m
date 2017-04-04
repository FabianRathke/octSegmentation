function [dataReturn modelsToAppend] = projToEigenspace(data,options,models)
% projToEigenspace - projects patches onto the first n modes derived via PCA; returns modes W and mu inside modelsToAppend
%
% Syntax:
%   [dataReturn modelsToAppend] = projToEigenspace(data,options,models)
%
% Inputs:
%   data - [struct] field .data contains the patches (the output of fetchPatches)
%   options - [cell-array] first entry holds handle to the function; second entry the number of modes; optional entries are W and mu
%   models - [struct] appearance models structure; should only be passed during prediction 
%
% Outputs:
%   dataReturn - [matrix] patches projected onto the n modes
%   modelsToAppend - [struct] fields to be appended to the models struct in trainAppearance
%      .W  - [matrix] first n eigenvectors
%      .mu - [array] mean patch
%
% See also: trainAppearance, fetchPatches

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 24-Nov-2016

% training
if nargin < 3
    if length(options)==2
		[V D] = eig(cov(data.data));
		W = V(:,end:-1:end-options{2}+1);
		mu = mean(data.data,1);
	else
		W = options{3};
		mu = options{4};
	end
	data.data = bsxfun(@minus,data.data,mu)*W;
	
	dataReturn = data.data;
	modelsToAppend.W = W;
	modelsToAppend.mu = mu;
% prediction
else
 	dataReturn = bsxfun(@minus,data,models(1).mu)*models(1).W;
	modelsToAppend = [];
end


