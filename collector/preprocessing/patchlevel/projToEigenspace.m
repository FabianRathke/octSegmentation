function [dataReturn modelsToAppend] = projToEigenspace(data,options,models)

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


