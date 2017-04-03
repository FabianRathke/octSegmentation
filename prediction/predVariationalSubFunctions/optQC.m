tic;

% calculate parameters of a Gaussian density, that is the product of two Gaussian densities (3.2.11 and 3.2.12 in the thesis)
% naming convention due to matrix cookbook (371) 
%Sigma_c = zeros(1,numColumnsShape*(numBounds-1));
%Sigma_c(idxB-numColumnsShape) = 1./(factorsPrecAVec(idxB-numColumnsShape) + 1./sigma_tilde_squared(idxB));

%mu_c = zeros(numRows,numColumnsShape*(numBounds-1));
%mu_c(:,idxB-numColumnsShape) = repmat(Sigma_c(idxB-numColumnsShape).*(1./sigma_tilde_squared(idxB).*mu_a_b(idxB)'),numRows,1) + repmat(Sigma_c(idxB-numColumnsShape).*factorsPrecAVec(idxB-numColumnsShape),numRows,1).*mu_a_b2(idxB-numColumnsShape,:)';

%c_c = (1./sqrt(2*pi*(factorsPrecAVec(ones(1,numRows),:).^-1) + sigma_tilde_squared(ones(1,numRows),idxB))'.*exp(-0.5*(mu_a_b2 - mu_a_b(idxB,ones(1,numRows))).^2.*(factorsPrecAVec(ones(1,numRows),:)'.^-1 + sigma_tilde_squared(ones(1,numRows),idxB)').^-1))';
%tmp = 1./(1./factorsPrecAVec(idxB-numColumnsShape) + sigma_tilde_squared(idxB));
%tmp2 = (mu_a_b2(idxB-numColumnsShape,:) - mu_a_b(idxB,ones(1,numRows)))';
%c_c = zeros(numRows,numColumnsShape*(numBounds-1));
%c_c(:,idxB-numColumnsShape) = repmat(sqrt(tmp)/sqrt(2*pi),numRows,1).*exp(-0.5*tmp2.*tmp2.*repmat(tmp,numRows,1));
%c_c(c_c < options.thresholdAccuracy) = 0;

if ~isfield(options,'doNotPredict')
	% call the c-function in case we need not to calculate pairwise marginals
	timeCFunc = 0;
	% old slower version of the C-code
	%a = tic; q_cc = optQCC(double(condQB),prediction,double(Sigma_c.^-1),double(mu_c),double(c_c),double(mu_a_b2'),numColumnsPred,numColumnsShape,columnsPredShapeVec(1,:),columnsPredShapeFactorVec(1,:),columnsPredShapeVec(2,:),columnsPredShapeFactorVec(2,:),double(factorsPrecAVec),hashTable); 
%	a = tic; [q_cc boundaries] = optQCCFaster(double(condQB),prediction,double(mu_a_b2),numColumnsPred,numColumnsShape,int32(columnsPredShapeVec(1,:)),columnsPredShapeFactorVec(1,:),int32(columnsPredShapeVec(2,:)),columnsPredShapeFactorVec(2,:),double(factorsPrecAVec),hashTable,cmin,cmax);
%	a = tic; [q_cc boundaries] = optQCCAvgVals(double(condQB),prediction,double(mu_a_b2),double(factorsPrecAvg),hashTable,numColumnsPred,cmin,cmax); 
	a = tic; [q_cc boundaries] = optQCC(condQB,prediction,mu_a_b2,factorsPrecAvg,hashTable,numColumnsPred,cmin,cmax); 
	boundaries = permute(boundaries,[3 1 2]);
	q_c.singleton = reshape(q_cc,[numRows,numBounds,numColumnsPred,numVolRegions]);
	timeCFunc = toc(a);
else
	% old matlab implementation; implicitly uses OmegaMatrices (see calcFuncVal for code to calculate these matrices)
	alpha = zeros(numRows,numBounds);
	beta = zeros(numRows,numBounds);
	c = size(1,numBounds);
	for volRegion = 1:numVolRegions
%		mu_a_b = mu_a_b2((1:numColumnsShape(volRegion)*(numBounds-1)) + sum(numColumnsShape(1:volRegion-1))*(numBounds-1),:);
		mu_a_b = mu_a_b2;
	
		for j = 1:numColumnsPred
			colFacA = columnsPredShapeFactor{volRegion}(1,j); colFacB = columnsPredShapeFactor{volRegion}(2,j);
			colA = columnsPredShape{volRegion}(1,j); colB = columnsPredShape{volRegion}(2,j);
			alpha(:) = 0; beta(:) = 0; c(:) = 0;
			idxPredict = find(~options.doNotPredict(:,j)); numBoundsLocal = length(idxPredict);

			if length(idxPredict > 0)
				
				pObs = squeeze(prediction(:,idxPredict,j,volRegion));
				numPrevCols = (0:numBounds-1)*numColumnsShape(volRegion) + sum(numColumnsShape(1:volRegion-1))*numBounds;
				% idx not considering the first boundary
				numPrevColsWithout = (0:numBounds-2)*numColumnsShape(volRegion) + sum(numColumnsShape(1:volRegion-1))*(numBounds-1);
				numPrevCols = numPrevCols(idxPredict); numPrevColsWithout = numPrevColsWithout(idxPredict(2:end)-1);
		
				% do the forward message passing
				alpha(:,1) = (colFacA*condQB(:,colA+numPrevCols(1)) + (colFacB*condQB(:,colB+numPrevCols(1)))).*pObs(:,1);
				c(1) = sum(alpha(:,1));
				alpha(:,1) = alpha(:,1)/c(1);
				for i = 2:numBoundsLocal
					idxNonZeroA = find((condQB(:,colA+numPrevCols(i)) + condQB(:,colB+numPrevCols(i)))~=0);
					idxNonZeroB = find(alpha(:,i-1)~=0);
		
					% directly calculates the marginals conditioned on all other boundaries
					alpha(idxNonZeroA,i) = pObs(idxNonZeroA,i).*((alpha(idxNonZeroB,i-1)*colFacA.*c_c(idxNonZeroB,colA + numPrevColsWithout(i-1)))'*exp(-0.5*Sigma_c(colA + numPrevColsWithout(i-1))^-1*((idxNonZeroA(:,ones(1,length(idxNonZeroB)))' - mu_c(idxNonZeroB,ones(1,length(idxNonZeroA))*(colA + numPrevColsWithout(i-1)))).^2)) + ...
					(alpha(idxNonZeroB,i-1)*colFacB.*c_c(idxNonZeroB,colB+numPrevColsWithout(i-1)))'*exp(-0.5*Sigma_c(colB+numPrevColsWithout(i-1))^-1*((idxNonZeroA(:,ones(1,length(idxNonZeroB)))' - mu_c(idxNonZeroB,ones(1,length(idxNonZeroA))*(colB + numPrevColsWithout(i-1)))).^2)))';
		
					% scale alpha
					c(i) = sum(alpha(:,i));
					alpha(:,i) = alpha(:,i)/c(i);
				end
		
				% do the backward message passing
				beta(:,numBoundsLocal) = ones(1,numRows);
		
				for i = numBoundsLocal-1:-1:1
					idxNonZeroA = (condQB(:,colA+numPrevCols(i+1)) + condQB(:,colB+numPrevCols(i+1)))~=0;
					idxNonZeroB = beta(:,i+1)~=0;
					idxFinal = find(logical(idxNonZeroA.*idxNonZeroB));
					idxB_ = find(alpha(:,i)~=0);
		
					beta(idxB_,i) = colFacA/c(i+1)*(beta(idxFinal,i+1).*pObs(idxFinal,i+1).*condQB(idxFinal,(colA+numPrevCols(i+1))))'*exp(-0.5*factorsPrecA{volRegion}((idxPredict(i)-1)*numColumnsShape(volRegion)+colA)*(idxFinal(:,ones(1,length(idxB_)))' - mu_a_b(ones(1,length(idxFinal))*((idxPredict(i)-1)*numColumnsShape(volRegion)+colA),idxB_)').^2)' ...
					+ colFacB/c(i+1)*(beta(idxFinal,i+1).*pObs(idxFinal,i+1).*condQB(idxFinal,(colB+numPrevCols(i+1))))'*exp(-0.5*factorsPrecA{volRegion}((idxPredict(i)-1)*numColumnsShape(volRegion)+colB)*(idxFinal(:,ones(1,length(idxB_)))' - mu_a_b(ones(1,length(idxFinal))*((idxPredict(i)-1)*numColumnsShape(volRegion)+colB),idxB_)').^2)';
				end

				tmp = sum(alpha(:,1:numBoundsLocal).*beta(:,1:numBoundsLocal),1);
				q_c.singleton(:,idxPredict,j) = alpha(:,1:numBoundsLocal).*beta(:,1:numBoundsLocal)./tmp(ones(1,numRows),:);
				% calc pairwise terms in order to calculate the value of the energy function
				if options.calcFuncVal && 0
					for i = 2:numBounds
						q_c.pairwise{volRegion,j,i-1} = sparse(c(i)^-1.*alpha(:,(i-1)*ones(1,numRows))'.*pObs(i*ones(numRows,1),:).*omegaTerms{volRegion,i,j}.*beta(i*ones(1,numRows),:));
						q_c.pairwise{volRegion,j,i-1} = q_c.pairwise{volRegion,j,i-1}/sum(sum(q_c.pairwise{volRegion,j,i-1}));
					end
				end
			end
		end
	end
end
if collector.options.printTimings
	if ~isfield(options,'doNotPredict')
		fprintf('[optQC]: %.3f s (C-Func %.3f)\n',toc,timeCFunc);
	else
		fprintf('[optQC]: %.3f s\n',toc);
	end
end

if options.plotting
   for volRegion = 1:numVolRegions
		toPlot = squeeze(sum(q_c.singleton(:,:,:,volRegion).*repmat((1:numRows)',[1,numBounds,numColumnsPred])));
		fileSaveName = sprintf('%s/qc_%d/%s_%d.eps',folderName,iter,filename,collector.options.labelIDs(volRegion));
		eval(sprintf('plotBScan(B%d,toPlot(:,columnsShapePred{volRegion}),collector.options.columnsShape{volRegion},fileSaveName,idxPredictFull)',collector.options.labelIDs(volRegion)));
	end
end

