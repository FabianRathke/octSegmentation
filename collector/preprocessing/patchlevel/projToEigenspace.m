function [dataReturn modelsToAppend] = projToEigenspace(data,options,models)
% projToEigenspace - projects patches onto the first n modes derived via PCA
%
% Syntax:
%   [dataReturn modelsToAppend] = projToEigenspace(data,options,models)
%
% Inputs:
%   data - [struct] field .data contains the patches (the output of fetchPatches)
%   options - [cell-array] first entry holds handle to the function; second entry the number of modes
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
% Last Revision: 04-May-2014

if nargin < 3
	[V D] = eig(cov(data.data));
	W = V(:,end-options{2}+1:end);
	mu = mean(data.data,1);
	data.data = (data.data-repmat(mu,size(data.data,1),1))*W;
	
	dataReturn = data.data;
	modelsToAppend.W = W;
	modelsToAppend.mu = mu;
else
 	dataReturn = (data - models(1).mu(ones(1,size(data,1)),:))*models(1).W;
	modelsToAppend = [];
end


