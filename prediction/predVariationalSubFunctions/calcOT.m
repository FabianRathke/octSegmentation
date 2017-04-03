if collector.options.printTimings
	ticCalcOT = tic;
end

%mu_diff = q_b.mu'-models.shapeModel.mu';
%product = sigmaML^-1*double(mu_diff') - sigmaML^-1*WML*M*(WML'*mu_diff');
%productB = product(idxProd) - (sigmaML^-1*mu_diff(idxProd) - sigmaML^-1*sum(WML(idxProd,:)'.*(M*squeeze(sum(reshape((WML(idxPrecB,:)'.*repmat(reshape(mu_diff(idxMu)',1,[])',1,size(M,1))')',[numBounds,numColumnsShapeTotal*numBounds,size(M,1)]),1))'),1))';
%% simple reshape does not work anymore since the number of columns per region can vary
%tmp = models.shapeModel.mu(idxProd) - sum(Sigma_a_a_.*productB(idxProd2),2);
%mu_a_b = [];
%for volRegion = 1:numVolRegions
%	mu_a_b = [mu_a_b; reshape(reshape(tmp((1:numColumnsShape(volRegion)*numBounds) + sum(numColumnsShape(1:volRegion-1))*numBounds),numBounds,[])',1,[])'];
%end
%mu_a_b = eval(sprintf('%s(mu_a_b)',collector.options.dataTypeCast));

% x was calculated in optQB
% mu_a_b is the mean of the second gaussian when calculating the expectation of p(b|c) with respect to q(b);
% covers terms p(b_{k,j}|b_{\setminus j})
productB = sigmaML^-1*x - sigmaML^-1*WML(idx_b,:)*(prodWMT(:,idx_b)*x) - K_jj_block(idx_b,idx_b)*x;
mu_a_b = zeros(numColumnsShape*numBounds,1);
mu_a_b(idx_b) = models.shapeModel.mu(idx_b) - K_jj_inverse_block(idx_b,idx_b)*productB;

%[condQB cmin cmax] = calcCondQB(mu_a_b,sigma_tilde_squared,numRows,numBounds,numColumnsShape,options.thresholdAccuracy);
[condQB cmin cmax] = calcCondQB(mu_a_b,sigma_tilde_squared,numRows,numBounds,numColumnsShape,numColumnsPred,int32(columnsPredShapeVec(1,:)-1),columnsPredShapeFactorVec(1,:),int32(columnsPredShapeVec(2,:)-1),columnsPredShapeFactorVec(2,:),options.thresholdAccuracy);

%tmp = X - mu_a_b(:,ones(1,numRows));
% evaluate 1-d gaussian for all columns
%condQB = repmat(1./sqrt(2*pi*sigma_tilde_squared),numRows,1).*exp((-0.5)*tmp.*tmp./sigma_tilde_squared(ones(1,numRows),:)')';
%condQB = (1./sqrt(2*pi*sigma_tilde_squared(ones(1,numRows),:)').*exp((-0.5)*((X - mu_a_b(:,ones(1,numRows)))).^2./sigma_tilde_squared(ones(1,numRows),:)'))';
% implicit cast to single precision
%condQB(condQB < options.thresholdAccuracy) = 0;

if collector.options.printTimings
   	if collector.options.calcOnGPU
   		GPUsync;
	end
	fprintf('[calcOT]: %.3f s\n',toc(ticCalcOT));
end
