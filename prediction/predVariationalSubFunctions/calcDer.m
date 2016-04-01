% calc terms described as $\tilde{p}_{j,k}$ in the paper
% P for Sigma is the same as P for mu

if collector.options.printTimings
	calcDerTic = tic;
end

if volRegion == 1
	factor = ((reshape(squeeze(boundaries(1,:,:))',[],1)-models.shapeModel.mu)'./sigma_tilde_squared);
	p_mu = factor*A_k_partial*WML'- factor*A_k_nonzero;
else
	p_mu = eval(sprintf('zeros(1,numColumnsShapeTotal*numBounds,%s);',collector.options.dataType));
	for volRegion = 1:numVolRegions
		columns = (1:(numColumnsShape(volRegion)*numBounds)) + sum(numColumnsShape(1:volRegion-1))*numBounds;
		p_mu = p_mu + ((reshape(squeeze(boundaries(volRegion,:,:))',[],1)-models.shapeModel.mu(columns))'./sigma_tilde_squared)*A_k(columns,:);
	end
end

if collector.options.printTimings
    if collector.options.calcOnGPU
        GPUsync;
    end
	fprintf('[calcDer]: %.3fs\n',toc(calcDerTic));
end
