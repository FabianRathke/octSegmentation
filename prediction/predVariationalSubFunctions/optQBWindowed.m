if collector.options.printTimings
	optQBTic = tic;
end

% DEBUG: direct implementation
p_mu = zeros(numColumnsShape*numBounds,1);
p_mu = p_mu + A_k'*((reshape(squeeze(boundaries(:,:,columnsShapePred{1}))',[],1)-models.shapeModel.mu)./sigma_tilde_squared');

b = -p_mu;
if iter == 1
	x = zeros(length(p_mu),1);
end
delta = 3;

n = numColumnsShape*numBounds;
% Newton
%tic; NewtonLeastSquares; toc;
Ab = zeros(size(b));
for i = 1:numWindows
	Ab(idx_b_Windowed{i}) = -2*(AWindowed{i}*b(idx_b_Windowed{i}));
end
BFGSLLeastSquares;
q_b.mu = zeros(numColumnsShape*numBounds,1);
q_b.mu(idx_b) = x + models.shapeModel.mu(idx_b);

% Mahalonbis distance
vec = q_b.mu(idx_b)-models.shapeModel.mu(idx_b);
%distance = sqrt(sigmaML^-1*vec'*vec - sigmaML^-1*(vec'*WML)*M*(WML'*vec));
for i = 1:numWindows
	z{i} = MWindowed{i}*(WMLWindowed{i}(idx_b_Windowed{i},:)'*vec(idx_b_Windowed{i}));
end

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
