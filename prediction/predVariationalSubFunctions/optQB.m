if collector.options.printTimings
	optQBTic = tic;
end

% conjugate gradient: Ax = b (K+\tilde{P}(\bar{mu}-\mu) =  \tilde{p}
x = conjgrad(P_mu + sigmaML^-1*eval(sprintf('eye(size(WML,1),%s)',collector.options.dataType)) - sigmaML^-1*WML*M*WML',-p_mu',eval(sprintf('zeros(1,length(p_mu),%s)',collector.options.dataType))');
q_b.mu = single(x) + models.shapeModel.mu;

% alternative conjugate gradient for very large shape prior matrices, where the explicit representation is to expensive in terms of memory
% inactive has to be checked and updated
if 0
	x = zeros(length(p_mu),1);
	r=single(-p_mu');
	p=r;
	rsold=r'*r;
	prodA = sigma_ml^-1*W_ml*M;

	matProd = sigma_ml^-1*M*W_ml';
	idx_all = (1:numBounds*numColumnsShapeTotal);

	for i=1:10000000
		Ap = zeros(length(p_mu),1);
		for volRegion = 1:numVolRegions
			for j = 1:numColumnsShape(volRegion)
				idx_j = (j:numColumnsShape(volRegion):numColumnsShape(volRegion)*numBounds) + sum(numColumnsShape(1:volRegion-1))*numBounds;
				idx_not_j = idx_all; idx_not_j(idx_j) = [];
				
				tmp = (1./sigma_tilde_squared(idx_j)').*(K_jj_inverse{volRegion,j}*(W_ml(idx_j,:)*(matProd(:,idx_not_j)*p(idx_not_j))));
				Ap(idx_not_j,1) = Ap(idx_not_j,1) + sum(matProd(:,idx_not_j)'*W_ml(idx_j,:)'*K_jj_inverse{volRegion,j}.*tmp(:,ones(1,length(idx_not_j)))',2);
			end
		end
		% add the contribution of the precision matrix of p(b)
		Ap = Ap + sigma_ml^-1*p - prodA*(W_ml'*p);
	 
		alpha=rsold/(p'*Ap);
		x=x+alpha*p;
		r=r-alpha*Ap;
		rsnew=r'*r;
		if sqrt(rsnew)<1e-10
			  break;
		else
		%	fprintf('%.10f\n',rsnew);
		end
		p=r+rsnew/rsold*p;
		rsold=rsnew;
	end
	q_b.mu = x + models.shapeModel.mu;
end

if collector.options.printTimings
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
