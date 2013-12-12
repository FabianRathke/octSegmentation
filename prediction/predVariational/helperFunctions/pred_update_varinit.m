idxPrec = zeros(numColumnsShapeTotal*numBounds,numBounds);
idxMu = zeros(numColumnsShapeTotal*numBounds,numBounds);
idxProd = zeros(numColumnsShapeTotal*numBounds,1);
idxProd2 = zeros(numColumnsShapeTotal*numBounds,numBounds);
Sigma_a_a_ = zeros(numColumnsShapeTotal*numBounds,numBounds);
idxPrecA = []; idxPrecB = [];
for volRegion = 1:numVolRegions
	for i = 1:numColumnsShape(volRegion)
		%idx_a = numColumnsShape*[0:numBounds-1]+i;
		numPreviousColumns = sum(numColumnsShape(1:volRegion-1))*numBounds;
		% where to look in the precision matrix
		idx_a = (i:numColumnsShape(volRegion):numColumnsShape(volRegion)*numBounds) + numPreviousColumns;
		% where to save, just stock it into one vector one after another
		a = ((i-1)*numBounds+1 + numPreviousColumns):(i*numBounds + numPreviousColumns);
		[A B] = meshgrid(idx_a);
		idxPrecB = [idxPrecB; reshape(B,[],1)];

		idxMu(a,:) = idx_a(ones(1,numBounds),:);
		idxProd(a) = idx_a;
		idxProd2(a,:) = a(ones(1,numBounds),:);
		Sigma_a_a_(a,:) = inv(eye(numBounds)*sigmaML^-1 - sigmaML^-1*models.shapeModel.WML(idx_a,:)*models.shapeModel.M*models.shapeModel.WML(idx_a,:)');
	end
end


