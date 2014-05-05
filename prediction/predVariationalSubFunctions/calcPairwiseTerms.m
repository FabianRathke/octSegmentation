% we do not need pairwise terms during the inference part
% but we require them for evaluating the different terms of the objective function

% in principle we could adapt the c-implementation of sum-product to directly calculate the pairwise terms for faster calculation
% but due to time constraints, we just stick with the matlab implementation (at the moment)

% first we calculate omega terms
omegaTerms = cell(numVolRegions,numBounds,numColumnsPred);
condQB2 = sparse(full(condQB'));
% interpolate between two shape columns
if size(columnsPredShape{1},1)==2
	for volRegion = 1:numVolRegions
		numPrevCols = sum(numColumnsShape(1:volRegion-1))*numBounds;
		for j = 1:numColumnsPred
			omegaTerms{volRegion,1,j} = columnsPredShapeFactor{volRegion}(1,j)*condQB2(columnsPredShape{volRegion}(1,j)+numPrevCols,:) + (columnsPredShapeFactor{volRegion}(2,j)*condQB2(columnsPredShape{volRegion}(2,j)+numPrevCols,:));
		end

		for k = 2:numBounds
			numPrevCols = (k-1)*numColumnsShape(volRegion) + sum(numColumnsShape(1:volRegion-1))*numBounds;
			for j = 1:numColumnsPred
				omegaTerms{volRegion,k,j} = columnsPredShapeFactor{volRegion}(1,j)*condQB2((columnsPredShape{volRegion}(1,j)+numPrevCols)*ones(1,numRows),:).*models.shapeModel.pTransV{volRegion}{columnsPredShape{volRegion}(1,j),k} + columnsPredShapeFactor{volRegion}(2,j)*condQB2((columnsPredShape{volRegion}(2,j)+numPrevCols)*ones(1,numRows),:).*models.shapeModel.pTransV{volRegion}{columnsPredShape{volRegion}(2,j),k};
			end
		end
	end
% choose the closest shape column
else
	for volRegion = 1:numVolRegions
		for j = 1:numColumnsPred
			idx = columnsPredShape(j) + (volRegion-1)*numColumnsShape*numBounds;
			omegaTerms{volRegion,1,j} = condQB2(idx,:);
		end

		for k = 2:numBounds
			for j = 1:numColumnsPred
				idx = (k-1)*numColumnsShape + columnsPredShape(j) + (volRegion-1)*numColumnsShape*numBounds;
				omegaTerms{volRegion,k,j} = condQB2(idx*ones(1,numRows),:).*models.shapeModel.pTransV{columnsPredShape(j)}{k}(end-numRows+1:end,end-numRows+1:end);
			end
		end
	end
end


% matlab implementation of sum-product; explicitly requires omega terms (c-implementation does not but calculates them on the fly)
% can also be used to check the correctness of alternative c-implementations
for volRegion = 1:numVolRegions
	mu_a_b = mu_a_b2((1:numColumnsShape(volRegion)*(numBounds-1)) + sum(numColumnsShape(1:volRegion-1))*(numBounds-1),:);

	for j = 1:numColumnsPred
		pObs = squeeze(prediction(:,:,j,volRegion))';

		alpha = zeros(numBounds,numRows);
		beta = zeros(numBounds,numRows);
		c = size(1,numBounds);
		numPrevCols = (0:numBounds-1)*numColumnsShape(volRegion) + sum(numColumnsShape(1:volRegion-1))*numBounds;
		% idx not considering the first boundary
		numPrevColsWithout = (0:numBounds-2)*numColumnsShape(volRegion) + sum(numColumnsShape(1:volRegion-1))*(numBounds-1);

		% do the forward message passing
		alpha(1,:) = (columnsPredShapeFactor{volRegion}(1,j)*condQB(:,columnsPredShape{volRegion}(1,j)+numPrevCols(1)) + (columnsPredShapeFactor{volRegion}(2,j)*condQB(:,columnsPredShape{volRegion}(2,j)+numPrevCols(1)))).*pObs(1,:)';
		c(1) = sum(alpha(1,:));
		alpha(1,:) = alpha(1,:)/c(1);
		for i = 2:numBounds
			idxNonZeroA = find((condQB(:,columnsPredShape{volRegion}(1,j)+numPrevCols(i)) + condQB(:,columnsPredShape{volRegion}(2,j)+numPrevCols(i)))~=0);
			idxNonZeroB = find(alpha(i-1,:)~=0);

			% directly calculates the marginals conditioned on all other boundaries
			alpha(i,idxNonZeroA) = pObs(i,idxNonZeroA).*((alpha(i-1,idxNonZeroB)*columnsPredShapeFactor{volRegion}(1,j).*c_c(idxNonZeroB,columnsPredShape{volRegion}(1,j)+numPrevColsWithout(i-1))')*exp(-0.5*Sigma_c(columnsPredShape{volRegion}(1,j)+numPrevColsWithout(i-1))^-1*((idxNonZeroA(:,ones(1,length(idxNonZeroB)))' - mu_c(idxNonZeroB,ones(1,length(idxNonZeroA))*(columnsPredShape{volRegion}(1,j)+numPrevColsWithout(i-1)))).^2)) + ...
			(alpha(i-1,idxNonZeroB)*columnsPredShapeFactor{volRegion}(2,j).*c_c(idxNonZeroB,columnsPredShape{volRegion}(2,j)+numPrevColsWithout(i-1))')*exp(-0.5*Sigma_c(columnsPredShape{volRegion}(2,j)+numPrevColsWithout(i-1))^-1*((idxNonZeroA(:,ones(1,length(idxNonZeroB)))' - mu_c(idxNonZeroB,ones(1,length(idxNonZeroA))*(columnsPredShape{volRegion}(2,j)+numPrevColsWithout(i-1)))).^2)));

			% scale alpha
			c(i) = sum(alpha(i,:));
			alpha(i,:) = alpha(i,:)/c(i);
		end

		% do the backward message passing
		beta(end,:) = ones(1,numRows);

		for i = numBounds-1:-1:1
			idxNonZeroA = (condQB(:,columnsPredShape{volRegion}(1,j)+numPrevCols(i+1)) + condQB(:,columnsPredShape{volRegion}(2,j)+numPrevCols(i+1)))~=0;
			idxNonZeroB = beta(i+1,:)~=0;
			idxFinal = find(logical(idxNonZeroA.*idxNonZeroB'));
			idxB_ = find(alpha(i,:)~=0);

			beta(i,idxB_) = columnsPredShapeFactor{volRegion}(1,j)/c(i+1)*(beta(i+1,idxFinal).*pObs(i+1,idxFinal).*condQB(idxFinal,(columnsPredShape{volRegion}(1,j)+numPrevCols(i+1)))')*exp(-0.5*factorsPrecA{volRegion}((i-1)*numColumnsShape(volRegion)+columnsPredShape{volRegion}(1,j))*(idxFinal(:,ones(1,length(idxB_)))' - mu_a_b(ones(1,length(idxFinal))*((i-1)*numColumnsShape(volRegion)+columnsPredShape{volRegion}(1,j)),idxB_)').^2)' ...
			+ columnsPredShapeFactor{volRegion}(2,j)/c(i+1)*(beta(i+1,idxFinal).*pObs(i+1,idxFinal).*condQB(idxFinal,(columnsPredShape{volRegion}(2,j)+numPrevCols(i+1)))')*exp(-0.5*factorsPrecA{volRegion}((i-1)*numColumnsShape(volRegion)+columnsPredShape{volRegion}(2,j))*(idxFinal(:,ones(1,length(idxB_)))' - mu_a_b(ones(1,length(idxFinal))*((i-1)*numColumnsShape(volRegion)+columnsPredShape{volRegion}(2,j)),idxB_)').^2)';
		end

		tmp = sum(alpha.*beta,2);
		q_c.singleton(volRegion,j,:,:) = alpha.*beta./tmp(:,ones(1,numRows));
		% calc pairwise terms in order to calculate the value of the energy function
		for i = 2:numBounds
			q_c.pairwise{volRegion,j,i-1} = sparse(c(i)^-1.*alpha((i-1)*ones(1,numRows),:)'.*pObs(i*ones(numRows,1),:).*omegaTerms{volRegion,i,j}.*beta(i*ones(1,numRows),:));
			q_c.pairwise{volRegion,j,i-1} = q_c.pairwise{volRegion,j,i-1}/sum(sum(q_c.pairwise{volRegion,j,i-1}));
		end
	end
end

