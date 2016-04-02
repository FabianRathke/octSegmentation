if collector.options.printTimings
	optQBTic = tic;
end

% calc terms described as $\tilde{p}_{j,k}$ in the paper
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

% conjugate gradient using the low-rank decomposition of \Sigma
x = zeros(length(p_mu),1);
r = -p_mu';
p=r;
rsold=r'*r;

%prodA = WML*M; % moved definition to initSegmentation.m
sigmaMLSquared = sigmaML^-2;


for i=1:10000000
    y = (A_k_partial*(WML'*p) - A_k_nonzero*p)./sigma_tilde_squared';
    Ap = WML*(A_k_partial'*y) - A_k_nonzero'*y + sigmaML^-1*p - sigmaML^-1*prodWM*(WML'*p);

    alpha=rsold/(p'*Ap);
    x=x+alpha*p;
    r=r-alpha*Ap;
    rsnew=r'*r;
    if sqrt(rsnew)<1e-10
          break;
    else
    %   fprintf('%.10f\n',rsnew);
    end
    p=r+rsnew/rsold*p;
    rsold=rsnew;
end
q_b.mu = x + models.shapeModel.mu;

if collector.options.printTimings
   if collector.options.calcOnGPU
        GPUsync;
    end
	fprintf('[optQB]: %.3f s\n',toc(optQBTic));
end

if options.plotting
	for volRegion = 1:numVolRegions
		idx = (1:numBounds*numColumnsShape(volRegion)) + sum(numColumnsShape(1:volRegion-1))*numBounds;
		toPlot = single(reshape(q_b.mu(idx),numColumnsShape(volRegion),numBounds));
		fileSaveName = sprintf('%s/qb_%d/%s_%d.eps',folderName,iter,filename,collector.options.labelIDs(volRegion));
		eval(sprintf('plotBScan(B%d,toPlot,collector.options.columnsShape{volRegion},fileSaveName)',collector.options.labelIDs(volRegion)))
	end
end
