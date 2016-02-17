if collector.options.printTimings
	optQBTic = tic;
end

% conjugate gradient: Ax = b ==> (K+\tilde{P})(\bar{mu}-\mu) =  \tilde{p}
%x = conjgrad(P_mu + sigmaML^-1*eval(sprintf('eye(size(WML,1),%s)',collector.options.dataType)) - sigmaML^-1*WML*M*WML',-p_mu',eval(sprintf('zeros(1,length(p_mu),%s)',collector.options.dataType))');
%q_b.mu = single(x) + models.shapeModel.mu;

% alternative conjugate gradient for very large shape prior matrices, where the explicit representation is to expensive in terms of memory --> use low-rank representation
x = zeros(length(p_mu),1);
r = -p_mu';
p=r;
rsold=r'*r;

%prodA = WML*M; % moved definition to initSegmentation.m
sigmaMLSquared = sigmaML^-2;

for i=1:10000000
	Ap = zeros(length(p_mu),1);
	prodC = prodWMT*p;
	for volRegion = 1:numVolRegions
		bj = zeros(size(M,1),numColumnsShape(volRegion));
		for j = 1:numColumnsShape(volRegion)
			idx_j = (j:numColumnsShape(volRegion):numColumnsShape(volRegion)*numBounds) + sum(numColumnsShape(1:volRegion-1))*numBounds;
			WMLj = K_jj_inverse{volRegion,j}*WML(idx_j,:);
			bj(:,j) = sigmaMLSquared*WMLj'*diag(1./sigma_tilde_squared(idx_j))*WMLj*(prodC-prodWMT(:,idx_j)*p(idx_j));
			Ap(idx_j) = -prodWM(idx_j,:)*bj(:,j);
		end
	end
	Ap = Ap + prodWM*sum(bj,2) + sigmaML^-1*p - sigmaML^-1*prodWM*(WML'*p);
	%Ap = sum(repmat((A_k*p)./sigma_tilde_squared',1,size(A_k,2)).*A_k)';%
	% add the contribution of the precision matrix of p(b)
%	Ap = Ap + sigmaML^-1*p - prodA*(WML'*p);
 
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
