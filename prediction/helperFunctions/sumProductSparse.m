function [gamma epsilon] = sumProductSparse(pStart,pTrans,pObs,calcPairwise)
% sumProductSparse - performs sum-product for chain graphs
%
% Syntax:
%	[gamma epsilon] = sumProductSparse(pStart,pTrans,pObs,calcPairwise)
%
% Inputs:
%   pStart       - [array] marginal probabilities for p(x_1)
%   pTrans       - [cell] cell array where each entry is a (preferably sparse) transition matrix
%   pObs         - [matrix] data terms
%   calcPairWise - [boolean] triggers the calculation of pairwise marginals 
%
% Outputs:
%   gamma   - [matrix] uniwise marginals
%   epsilon - [matrix] pairwise marginals

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 12-Dec-2013

if nargin < 4
	calcPairwise = false;
end

warning off;
[K N] = size(pObs);

alpha = zeros(size(pObs));
beta = zeros(size(pObs));
gamma = zeros(size(pObs));
if calcPairwise
	epsilon = zeros(K^2,N-1);
end
c = size(1,N);

% do the forward message passing
alpha(:,1) = pStart'.*pObs(:,1);
c(1) = sum(alpha(:,1));
alpha(:,1) = alpha(:,1)/c(1);
for i = 2:N
	alpha(:,i) = pObs(:,i)'.*(alpha(:,i-1)'*pTrans{i});
	if i ==2
		log(alpha(:,i))';
	end
	% scale alpha
	c(i) = sum(alpha(:,i));
	alpha(:,i) = alpha(:,i)/c(i);
end

% do the backward message passing
beta(:,N) = 1;
gamma(:,N) = alpha(:,N).*beta(:,N);

for i = N-1:-1:1
	beta(:,i) = beta(:,i+1)'/c(i+1).*pObs(:,i+1)'*pTrans{i+1}';
	% calculate marginals p(z_n|X)
	gamma(:,i) = alpha(:,i).*beta(:,i);
	% calculate joint marginals p(z_n,z_n-1|X)
	if calcPairwise
		epsilon(:,i+1) = reshape(c(i+1).*alpha(:,i*ones(1,K)).*pObs(:,(i+1)*ones(K,1))'.*pTrans{i+1}.*beta(:,(i+1)*ones(1,K))',K^2,1);
		epsilon(:,i+1) = epsilon(:,i+1)/sum(epsilon(:,i+1));
	end
end

marginals = gamma;
if calcPairwise
	epsilon(epsilon<10^-25) = 0;
	epsilon = sparse(epsilon);
end

warning on;
