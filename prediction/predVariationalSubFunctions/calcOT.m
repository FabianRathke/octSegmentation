if collector.options.printTimings
	ticCalcOT = tic;
end

mu_diff = q_b.mu'-models.shapeModel.mu';
product = sigmaML^-1*double(mu_diff')- sigmaML^-1*models.shapeModel.WML*models.shapeModel.M*(models.shapeModel.WML'*mu_diff');
productB = product(idxProd) - (sigmaML^-1*mu_diff(idxProd) - sigmaML^-1*sum(models.shapeModel.WML(idxProd,:)'.*(models.shapeModel.M*squeeze(sum(reshape((models.shapeModel.WML(idxPrecB,:)'.*repmat(reshape(mu_diff(idxMu)',1,[])',1,size(M,1))')',[numBounds,numColumnsShapeTotal*numBounds,size(M,1)]),1))'),1))';
% simple reshape does not work anymore since the number of columns per region can vary
tmp = models.shapeModel.mu(idxProd) - sum(Sigma_a_a_.*productB(idxProd2),2);
mu_a_b = [];
for volRegion = 1:numVolRegions
	mu_a_b = [mu_a_b; reshape(reshape(tmp((1:numColumnsShape(volRegion)*numBounds) + sum(numColumnsShape(1:volRegion-1))*numBounds),numBounds,[])',1,[])'];
end
mu_a_b = eval(sprintf('%s(mu_a_b)',collector.options.dataTypeCast));

condQB = (1./sqrt(2*pi*sigma_tilde_squared(ones(1,numRows),:)').*exp((-0.5)*((X - mu_a_b(:,ones(1,numRows)))).^2./sigma_tilde_squared(ones(1,numRows),:)'))';
% implicit cast to single precision
condQB(condQB < options.thresholdAccuracy) = 0;

if collector.options.printTimings
	GPUsync;
	fprintf('[calcOT]: %.3f s\n',toc(ticCalcOT));
end
