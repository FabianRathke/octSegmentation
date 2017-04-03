% calculate pairwise terms
%calcPairwiseTerms;

% terms of q_b
%K = (sigmaML^-1*eye(size(WML,1)) - sigmaML^-1*WML*M*WML');
%funcVal.Shape = 0.5*sum(sum(K'.*(q_b.mu*q_b.mu' -2*q_b.mu*models.shapeModel.mu')));

clear funcVal
if ~exist('predictionFuncVal')
	predictionFuncVal = prediction;
	%fprintf('DEBUG: No predictionFuncVal\n')
end

[funcVal.q_c_singleton funcVal.q_c_data] = calcFuncValsC(q_c.singleton,predictionFuncVal);
funcVal.q_c_singleton = reshape(funcVal.q_c_singleton,numBounds,numColumnsPred);
funcVal.q_c_data = reshape(funcVal.q_c_data,[1,numColumnsPred,numBounds]);

%% negative entropy of q_c (singleterm entropy)
%for volRegion = 1:numVolRegions
%	for j = 1:numColumnsPred
%		for k = 1:numBounds
%			tmp = q_c.singleton(:,k,j,volRegion);
%			funcVal.q_c_singleton(k,j,volRegion) = -sum(tmp(tmp~=0).*log(tmp(tmp~=0)));
%		end
%	end
%end
%
%% negative pairwise entropy (== mutual information)
%for volRegion = 1:numVolRegions
%	for j = 1:numColumnsPred
%		for k = 1:numBounds-1
%%			[I J] = find(q_c.pairwise{volRegion,j,k}~=0); I = unique(I); J = unique(J);
%%			tmp = q_c.pairwise{volRegion,j,k}(I,J).*log(q_c.pairwise{volRegion,j,k}(I,J)./(squeeze(q_c.singleton(I,k,j,volRegion))*squeeze(q_c.singleton(J,k+1,j,volRegion))'));
%%			funcVal.q_c_pairwise(volRegion,j,k) = full(sum(tmp(~isnan(tmp)&~isinf(tmp))));
%
%%			tmp = q_c.pairwise{volRegion,j,k}(I,J).*log(omegaTerms{volRegion,k+1,j}(I,J));
%%			funcVal.q_c_shape(volRegion,j,k+1) = - full(sum(tmp(~isnan(tmp)&~isinf(tmp)))); 
%		end
%	end
%end
%
%% Data & Shape Terms
%for volRegion = 1:numVolRegions
%	% unary terms
%	for j = 1:numColumnsPred
%		I = find(q_c.singleton(:,1,j,volRegion)~=0);
%		tmp = squeeze(q_c.singleton(I,1,j,volRegion))'.*squeeze(log(predictionFuncVal(I,1,j,volRegion)))';
%		funcVal.q_c_data(volRegion,j,1) = -full(sum(tmp(~isnan(tmp)&~isinf(tmp))));
%%		tmp = squeeze(q_c.singleton(volRegion,j,1,I))'.*log([omegaTerms{volRegion,1,j}(I)]);
%%		funcVal.q_c_shape(volRegion,j,1) = -full(sum(tmp(~isnan(tmp)&~isinf(tmp))));
%	end
%
%	for k = 2:numBounds
%		% pairwise terms
%		for j = 1:numColumnsPred
%			I = find(q_c.singleton(:,k,j,volRegion)~=0);
%			tmp = squeeze(q_c.singleton(I,k,j,volRegion))'.*squeeze(log(predictionFuncVal(I,k,j,volRegion)))';
%			funcVal.q_c_data(volRegion,j,k) = -full(sum(tmp(~isnan(tmp)&~isinf(tmp))));
%%			tmp = q_c.pairwise{volRegion,j,k-1}.*log(omegaTerms{volRegion,k,j});
%%			funcVal.q_c_shape(volRegion,j,k) = - full(sum(tmp(~isnan(tmp)&~isinf(tmp)))); 
%		end
%	end
%end
%funcVal.entropy = sum(sum(sum(funcVal.q_c_singleton))) + sum(sum(sum(funcVal.q_c_pairwise)));
%q_c_data = sum(sum(sum(funcVal.q_c_data)));
%q_c_shape = sum(sum(sum(funcVal.q_c_shape)));

%fprintf('Objective term values calculated and stored in output.funcVal\n');
