% calculate pairwise terms
calcPairwiseTerms;

% terms of q_b
K = (sigmaML^-1*eye(size(WML,1)) - sigmaML^-1*WML*M*WML');
funcVal.Shape = 0.5*sum(sum(K'.*(q_b.mu*q_b.mu' -2*q_b.mu*models.shapeModel.mu')));

% negative entropy of q_c (singleterm entropy)
for volRegion = 1:numVolRegions
	for j = 1:numColumnsPred
		for k = 1:numBounds
			tmp = q_c.singleton(volRegion,j,k,:);
			funcVal.q_c_singleton(volRegion,j,k) = sum(tmp(tmp~=0).*log(tmp(tmp~=0)));
		end
	end
end

% negative pairwise entropy (== mutual information)
for volRegion = 1:numVolRegions
	for j = 1:numColumnsPred
		for k = 1:numBounds-1
			tmp = q_c.pairwise{volRegion,j,k}.*log(q_c.pairwise{volRegion,j,k}./(squeeze(q_c.singleton(volRegion,j,k*ones(1,numRows),:))'.*squeeze(q_c.singleton(volRegion,j,k+1*ones(1,numRows),:))));
			funcVal.q_c_pairwise(volRegion,j,k) = full(sum(tmp(~isnan(tmp)&~isinf(tmp))));
		end
	end
end

% Data
for volRegion = 1:numVolRegions
	% unary terms
	for j = 1:numColumnsPred
		tmp = squeeze(q_c.singleton(volRegion,j,1,:))'.*log([omegaTerms{volRegion,1,j}]);
		funcVal.q_c_shape(volRegion,j,1) = -full(sum(tmp(~isnan(tmp)&~isinf(tmp))));
		tmp = squeeze(q_c.singleton(volRegion,j,1,:))'.*squeeze(prediction(:,1,j,volRegion))';
		funcVal.q_c_data(volRegion,j,1) = -full(sum(tmp(~isnan(tmp)&~isinf(tmp))));
	end

	for k = 2:numBounds
		% pairwise terms
		for j = 1:numColumnsPred
			tmp = squeeze(q_c.singleton(volRegion,j,k,:))'.*squeeze(prediction(:,k,j,volRegion))';
			funcVal.q_c_data(volRegion,j,k) = - full(sum(tmp(~isnan(tmp)&~isinf(tmp))));
			tmp = q_c.pairwise{volRegion,j,k-1}.*log(omegaTerms{volRegion,k,j});
			funcVal.q_c_shape(volRegion,j,k) = - full(sum(tmp(~isnan(tmp)&~isinf(tmp)))); 
		end
	end
end

entropy = sum(sum(sum(funcVal.q_c_singleton))) + sum(sum(sum(funcVal.q_c_pairwise)));
q_c_data = sum(sum(sum(funcVal.q_c_data)));
q_c_shape = sum(sum(sum(funcVal.q_c_shape)));

fprintf('Objective term values calculated and stored in output.funcVal\n');
