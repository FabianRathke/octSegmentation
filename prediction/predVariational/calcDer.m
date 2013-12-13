% calc terms described as $\tilde{P}_{j,k}$ in the paper
% P for Sigma is the same as P for mu

if collector.options.printTimings
	calcDerTic = tic;
end

p_mu = eval(sprintf('zeros(1,numColumnsShapeTotal*numBounds,%s);',collector.options.dataType));

% calc terms for k = 1
for volRegion = 1:numVolRegions
	idx = 1:numRows;
	idx = idx(ones(1,numColumnsShape(volRegion)),:);
	for k = 1:numBounds
		columns = (1:numColumnsShape(volRegion)) + (k-1)*numColumnsShape(volRegion) + sum(numColumnsShape(1:volRegion-1))*numBounds;
		if collector.options.calcOnGPU
			factor = GPUsingle(sum(squeeze(q_c.singleton(volRegion,columnsShapePred{volRegion},k,:)).*(idx - models.shapeModel.mu(columns,ones(1,numRows))),2))./sigma_tilde_squared(columns)';
		else
			factor = sum(squeeze(q_c.singleton(volRegion,columnsShapePred{volRegion},k,:)).*(idx - models.shapeModel.mu(columns,ones(1,numRows))),2)./sigma_tilde_squared(columns)';
		end
		p_mu = p_mu + factor'*A_k(columns,:);
	end
end

if collector.options.printTimings
	GPUsync;
	fprintf('[calcDer]: %.3fs\n',toc(calcDerTic));
end
