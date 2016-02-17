% calc $A^j_k = (K_{jj}^-1*K_{j,\setminus{j}})_{k,*}$
if collector.options.printTimings
	ticCalcShape = tic;
end

A_k = eval(sprintf('zeros(numBounds*numColumnsShapeTotal,%s);',collector.options.dataType));
K_jj_inverse = cell(numVolRegions,max(numColumnsShape));	
sigma_tilde_squared = zeros(1,numBounds*numColumnsShapeTotal);
idx_diag = 1:numBounds+1:numBounds^2;
idx_all = (1:numBounds*numColumnsShapeTotal);
	
% naive implementation for each image column
for volRegion = 1:numVolRegions
	for j = 1:numColumnsShape(volRegion)
		idx_j = (j:numColumnsShape(volRegion):numColumnsShape(volRegion)*numBounds) + sum(numColumnsShape(1:volRegion-1))*numBounds;
		idx_not_j = idx_all; idx_not_j(idx_j) = [];
		
		% use the CPU version of WML and M explicitly	
		K_jj_inverse{volRegion,j} = inv(eye(numBounds)*sigmaML^-1 - sigmaML^-1*WML(idx_j,:)*M*WML(idx_j,:)');
		sigma_tilde_squared(idx_j) = K_jj_inverse{volRegion,j}(idx_diag)/options.alpha;
	
		% calculating A_k^j: note that we do not need the sigma^-2I part, since we pull a part of K that does not include the diagonal
		if collector.options.calcOnGPU
			A_k(idx_j,idx_not_j) = GPUsingle(K_jj_inverse{volRegion,j})*-sigmaML^-1*WML(idx_j,:)*M*WML(idx_not_j,:)';
		else
			A_k(idx_j,idx_not_j) = K_jj_inverse{volRegion,j}*-sigmaML^-1*WML(idx_j,:)*prodWMT(:,idx_not_j);
		end
	end
end
clear idx_*
% cast to different data type if required
sigma_tilde_squared = eval(sprintf('%s(sigma_tilde_squared);',collector.options.dataTypeCast));

%factor = (1./sigma_tilde_squared');
%P_mu = (A_k.*factor(:,ones(1,numBounds*numColumnsShapeTotal)))'*A_k;

if collector.options.printTimings
	if collector.options.calcOnGPU
		GPUsync;
	end
	fprintf('[Initialized Shape Terms]: %.3fs \n',toc(ticCalcShape));
end
