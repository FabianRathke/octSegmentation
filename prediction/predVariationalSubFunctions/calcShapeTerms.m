% calc $A^j_k = (K_{jj}^-1*K_{j,\setminus{j}}p)_{k,*}$
if collector.options.printTimings
	ticCalcShape = tic;
end

%A_k = eval(sprintf('zeros(numBounds*numColumnsShapeTotal,%s);',collector.options.dataType));
sigma_tilde_squared = zeros(1,numBounds*numColumnsShapeTotal);
idx_diag = 1:numBounds+1:numBounds^2;

idxI = zeros(1,numColumnsShapeTotal*numBounds^2);
idxJ = zeros(1,numColumnsShapeTotal*numBounds^2);
s1 = zeros(1,numColumnsShapeTotal*numBounds^2); s2 = zeros(1,numColumnsShapeTotal*numBounds^2); s3 = zeros(1,numColumnsShapeTotal*numBounds^2);
idx_all = 1:numColumnsShapeTotal*numBounds;

% naive implementation for each image column
for volRegion = 1:numVolRegions
	for j = 1:numColumnsShape(volRegion)
		idx_j = (j:numColumnsShape(volRegion):numColumnsShape(volRegion)*numBounds) + sum(numColumnsShape(1:volRegion-1))*numBounds;
		
		% use the CPU version of WML and M explicitly
		tmp = eye(numBounds)*sigmaML^-1 - sigmaML^-1*WML(idx_j,:)*M*WML(idx_j,:)';
		K_jj_inverse = inv(tmp);
		sigma_tilde_squared(idx_j) = K_jj_inverse(idx_diag)/options.alpha;

		countJ = sum(numColumnsShape(1:volRegion-1)) + j; 
		idxRange = (countJ-1)*numBounds^2 + 1:countJ*numBounds^2;
		[idxI(idxRange) idxJ(idxRange)] = meshgrid(idx_j); 
		s1(idxRange) = K_jj_inverse(:);
		s3(idxRange) = tmp(:);	

		% calculating A_k^j: note that we do not need the sigma^-2I part, since we pull a part of K that does not include the diagonal
		s2(idxRange) = reshape(K_jj_inverse*-sigmaML^-1*WML(idx_j,:)*prodWMT(:,idx_j),1,[]);
		% DEBUG: CALC A_K
%		idx_not_j = idx_all; idx_not_j(idx_j) = [];
%		A_k(idx_j,idx_not_j) = K_jj_inverse*-sigmaML^-1*WML(idx_j,:)*prodWMT(:,idx_not_j);
	end
end
% cast to different data type if required
sigma_tilde_squared = eval(sprintf('%s(sigma_tilde_squared);',collector.options.dataTypeCast));

K_jj_inverse_block = sparse(idxI,idxJ,s1);
K_jj_block = sparse(idxI,idxJ,s3);
A_k_nonzero = sparse(idxI,idxJ,s2);
A_k_partial = -sigmaML^-1*K_jj_inverse_block*WML*M;

% DEBUG CALC \tilde{P}
%factor = (1./sigma_tilde_squared');
%P_mu = (A_k.*factor(:,ones(1,numBounds*numColumnsShapeTotal)))'*A_k;

if collector.options.printTimings
	if collector.options.calcOnGPU
		GPUsync;
	end
	fprintf('[Initialized Shape Terms]: %.3fs \n',toc(ticCalcShape));
end
