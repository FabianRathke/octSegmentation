if collector.options.printTimings
	optQBTic = tic;
end

% DEBUG: direct implementation
%x = (P_mu(idx_b,idx_b) + sigmaML^-1*eye(length(idx_b)) - sigmaML^-1*WML(idx_b,:)*M*WML(idx_b,:)')\p_mu(idx_b);

% at the very first iteration, only trust A-Scans with a good initialization
if 0 && iter == 1 && options.calcFuncVal && exist('margins')
	columnWiseFuncVal = sum(squeeze(funcVal.q_c_data(volRegion,columnsShapePred{volRegion},:))');
	boundFilter = [sum(margins.mu(4,:)) + sum(margins.sigma(4,:)*2)];
	% check for a systematic offset for the funcVals (for example if there is a very dark scan, the trained appearance models are not representative for the funcVal)
	if min(columnWiseFuncVal) > boundFilter
		columnWiseFuncVal = columnWiseFuncVal - min(columnWiseFuncVal);
	end

	idxFail = columnWiseFuncVal - [sum(margins.mu(4,:)) + sum(margins.sigma(4,:)*2)]>0;
	tmp= 1:numColumnsShape;
	idx_b_tmp = reshape([repmat(tmp(~idxFail),numBounds,1) + repmat([0:8]'*numColumnsShape,1,sum(~idxFail))]',1,[]);
   	idx_a_tmp = 1:numColumnsShape*numBounds;
	idx_a_tmp(idx_b_tmp) = [];
	tmp = reshape(squeeze(boundaries(1,:,columnsShapePred{volRegion}))',1,[])';

	fprintf('%d positions OK, %d positions failed\n',length(idx_b_tmp),length(idx_a_tmp));
	if length(idx_a_tmp) > 0	
		muAB = models.shapeModel.mu(idx_a_tmp) + inv(sigmaML^-1*eye(length(idx_a_tmp)) - sigmaML^-1*WML(idx_a_tmp,:)*M*WML(idx_a_tmp,:)')*sigmaML^-1*WML(idx_a_tmp,:)*M*(WML(idx_b_tmp,:)'*(tmp(idx_b_tmp)-models.shapeModel.mu(idx_b_tmp))); 
		tmp(idx_a_tmp) = muAB;
	%	factor = (tmp-models.shapeModel.mu)./sigma_tilde_squared';
	%	p_mu(idx_b) = WML(idx_b,:)*(A_k_partial(idx_b_tmp,:)'*factor(idx_b_tmp)) - A_k_nonzero(idx_b_tmp,idx_b)'*factor(idx_b_tmp);
		boundaries(1,:,columnsShapePred{volRegion}) = reshape(tmp,[],numBounds)';
		clear idx_tmp_b;
	end
end


if 0 &&(iter == 1 && exist('textureQuality'))
	faultyColumns = unique(columnsPredShape{1}(1,find(textureQuality(collector.options.columnsPred) < 5)));
	idxDelete = reshape((repmat(faultyColumns,numBounds,1) + repmat((0:(numBounds-1))'*numColumnsShape,1,length(faultyColumns)))',1,[]);
	% temporarily change idx_b
	idx_b_save = idx_b;
	idx_b = setdiff(idx_b,idxDelete);
	plot(find(textureQuality(collector.options.columnsPred) < 5),reshape(reshape(squeeze(boundaries(1,:,find(textureQuality(collector.options.columnsPred) < 5)))',1,[])',[],9),'o');
	prodWMSave = prodWM;
	MSave = M;
	M = inv(WML(idx_b,:)'*WML(idx_b,:) + sigmaML*eye(size(WML,2)));
	prodWM = WML*M;
	calcShapeTerms;
	%	% DEBUG plot
	if 0
		eval(sprintf('BTmp = loadData(files(file).name,collector.options);',collector.options.labelID));
		figure; imagesc(BTmp(:,collector.options.columnsPred)); colormap gray; hold on;
		plot(reshape(reshape(squeeze(boundaries(1,:,:))',1,[])',[],9),'*');
	end
end

% calc terms described as $\tilde{p}_{j,k}$ in the paper
if volRegion == 1
	factor = (reshape(squeeze(boundaries(1,:,columnsShapePred{1}))',[],1)-models.shapeModel.mu)./sigma_tilde_squared';
	p_mu = zeros(numColumnsShape*numBounds,1);
	p_mu(idx_b) = WML(idx_b,:)*(A_k_partial(idx_b,:)'*factor(idx_b)) - A_k_nonzero(idx_b,idx_b)'*factor(idx_b);
%    p_mu = factor'*A_k_partial*WML'- factor*A_k_nonzero;
else
	error('Reimplement optQB.m for 3_D');
end

% conjugate gradient using the low-rank decomposition of \Sigma
x = zeros(length(idx_b),1);
r = -p_mu(idx_b);
p=r;
rsold=r'*r;
sigmaMLSquared = sigmaML^-2;
sigma_tilde_squared_inv = 1./sigma_tilde_squared;

for i=1:100
%  	y = (A_k_partial(idx_b,:)*(WML(idx_b,:)'*p) - A_k_nonzero(idx_b,idx_b)*p).*sigma_tilde_squared_inv(idx_b)';
%    Ap = WML(idx_b,:)*(A_k_partial(idx_b,:)'*y) - A_k_nonzero(idx_b,idx_b)*y + sigmaML^-1*p - sigmaML^-1*prodWM(idx_b,:)*(WML(idx_b,:)'*p);
  	y = (A_k_partial*(WML'*p) - A_k_nonzero*p).*sigma_tilde_squared_inv';
    Ap = WML*(A_k_partial'*y) - A_k_nonzero*y + sigmaML^-1*p - sigmaML^-1*prodWM*(WML'*p);

    alpha=rsold/(p'*Ap);
    x=x+alpha*p;
    r=r-alpha*Ap;
    rsnew=r'*r;
    if sqrt(rsnew)<1e-10
          break;
    end
    p=r+rsnew/rsold*p;
    rsold=rsnew;
end
q_b.mu = zeros(numColumnsShape*numBounds,1);
q_b.mu(idx_b) = x + models.shapeModel.mu(idx_b);

% use conditional mean estimate to fill missing indices
if 0 && (iter == 1) && exist('textureQuality')
	prodWM = prodWMSave;
	M = MSave;
	calcShapeTerms;
    fprintf('DEBUG: Texture Quality %d positions OK, %d positions failed\n',length(idx_b),length(idxDelete));
    if length(idxDelete) > 0
        muAB = models.shapeModel.mu(idxDelete) + inv(sigmaML^-1*eye(length(idxDelete)) - sigmaML^-1*WML(idxDelete,:)*M*WML(idxDelete,:)')*sigmaML^-1*WML(idxDelete,:)*M*(WML(idx_b,:)'*(q_b.mu(idx_b)-models.shapeModel.mu(idx_b)));
        q_b.mu(idxDelete) = muAB;
	end
	idx_b = idx_b_save;
end

if 0
    eval(sprintf('BTmp = loadData(files(file).name,collector.options);',collector.options.labelID));
    figure; imagesc(BTmp(:,models.shapeModel.columnsShape)); colormap gray; hold on;
    plot(reshape(q_b.mu(idx_b),[],9));
end

vec = q_b.mu(idx_b)-models.shapeModel.mu(idx_b);
% Mahalonbis distance
%distance = sqrt(sigmaML^-1*vec'*vec - sigmaML^-1*(vec'*WML)*M*(WML'*vec));
z = M*(WML(idx_b,:)'*vec);
%fprintf('Distance of E[q_b] to shape prior: %.2f\n',norm(z));

if collector.options.printTimings
   if collector.options.calcOnGPU
        GPUsync;
    end
	fprintf('[optQB]: %.3f s\n',toc(optQBTic));
end

if options.plotting
	for volRegion = 1:numVolRegions
		idx = (1:numBounds*numColumnsShape(volRegion)) + sum(numColumnsShape(1:volRegion-1))*numBounds;
		toPlot = single(reshape(q_b.mu(idx),numColumnsShape(volRegion),numBounds))';
		fileSaveName = sprintf('%s/qb_%d/%s_%d.eps',folderName,iter,filename,collector.options.labelIDs(volRegion));
		eval(sprintf('plotBScan(B%d,toPlot,collector.options.columnsShape{volRegion},fileSaveName,idxPredictFull)',collector.options.labelIDs(volRegion)))
	end
end
