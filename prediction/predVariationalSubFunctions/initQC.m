%function [q_c output] = initQC(file,collector,models,options,q_c)
if ~isfield(options,'appearance')
	% fetch appearance Terms
	predictionA = predAppearance(files(file),collector,models.appearanceModel,options);
else
	predictionA.prediction = options.appearance;
end

prediction = zeros(numRows,numClasses,numColumnsPred,numVolRegions);

for volRegion = 1:numVolRegions
	% load ground truth and the scan itself
	collector.options.labelID = collector.options.labelIDs(file,volRegion);
	% load labels if set by the user
	if collector.options.loadLabels
		output.trueLabels{file,volRegion} = loadLabels(files(file).name,collector.options);
	end

	if options.plotting
		eval(sprintf('B%d = loadData(files(file).name,collector.options);',collector.options.labelID));
	end
	prediction(:,:,:,volRegion) = permute(reshape(predictionA.prediction{volRegion},[numClasses,collector.options.Y,numColumnsPred]),[2 1 3]);
end
% we only need the subset of boundary classes
prediction = prediction(:,length(collector.options.LayersPred) + (1:length(collector.options.EdgesPred)),:,:);
% renormalize
%prediction = prediction./repmat(sum(prediction,2),[1 numBounds 1 1]);
% reduce the precision
prediction(prediction < options.thresholdAccuracy) = 0;

% outputs pixelwise appearance terms too
if collector.options.saveAppearanceTerms
	output.appearanceTerms.prediction{file} = single(predictionA.prediction{1});
	if collector.options.loadLabels
		output.appearanceTerms.trueLabels{file} = single(predictionA.trueLabels{1});
	end
end
clear predictionA;

% use old matlab code when conditioning on partial segmentation
initqc = tic;

if isfield(options,'segmentation')
	if ~isfield(models.shapeModel,'pTransV')
		models.shapeModel = preCalcTransitionMatrices(collector,models.shapeModel,10^-20);
	end
%	for volRegion = 1:numVolRegions
	columnsToSegment = find(sum(options.idxRecalc)>0);

	for j = 1:length(columnsToSegment)
		i = columnsToSegment(j);
		idxA = sum(numColumnsShape(1:volRegion-1))*numBounds + i;
		pObs = squeeze(prediction(:,:,columnsShapePred{volRegion}(i),volRegion));
		fixBounds = find(~options.idxRecalc(:,columnsToSegment(j)));
		for k = 1:length(fixBounds)
			pObs(:,fixBounds(k)) = zeros(numRows,1);
			segm = options.segmentation(fixBounds(k),columnsToSegment(j));
			pObs([floor(segm) ceil(segm)],fixBounds(k)) = [ceil(segm)-segm segm-floor(segm)];
		end

		variance = sum(models.shapeModel.WML(idxA,:).^2) + models.shapeModel.sigmaML;
		% calculate probabilities for first boundary
		pStart = 1/sqrt(2*pi*variance)*exp(-0.5*(1/variance)*((1:numRows) - models.shapeModel.mu(idxA)).^2);
		q_c.singleton(volRegion,columnsShapePred{volRegion}(i),:,:) = sumProductSparse(pStart,models.shapeModel.pTransV{volRegion}(i,:),pObs)';
	end
else
	if 0
		if ~isfield(models.shapeModel,'pTransV')
			models.shapeModel = preCalcTransitionMatrices(collector,models.shapeModel,10^-20);
		end
		for volRegion = 1:numVolRegions
			% initialization only has to be made for columns relvant for updating the q(b) distribution
			for i = 1:numColumnsShape(volRegion)
				idxA = sum(numColumnsShape(1:volRegion-1))*numBounds + i;
				pObs = squeeze(prediction(:,:,columnsShapePred{volRegion}(i),volRegion));
	
				variance = sum(models.shapeModel.WML(idxA,:).^2) + models.shapeModel.sigmaML;
				% calculate probabilities for first boundary
				pStart = 1/sqrt(2*pi*variance)*exp(-0.5*(1/variance)*((1:numRows) - models.shapeModel.mu(idxA)).^2);
				q_c.singleton(volRegion,columnsShapePred{volRegion}(i),:,:) = sumProductSparse(pStart,models.shapeModel.pTransV{volRegion}(i,:),pObs)';
			end
		end
	end
	% C version; segments all columns for one BScan; needs as input the indices of the first boundary inside the shape model
	for volRegion = 1:numVolRegions
		% the -1 is C-indexing
		idxA = sum(numColumnsShape(1:volRegion-1))*numBounds + (1:numColumnsShape(volRegion)) - 1;
		q_c.singleton(volRegion,columnsShapePred{volRegion},:,:) = permute(reshape(sumProductSparseC(prediction(:,:,columnsShapePred{volRegion},volRegion),models.shapeModel(volRegion).mu,models.shapeModel(volRegion).WML,models.shapeModel(volRegion).sigmaML,int32(idxA),hashTable),[numRows,numBounds,numColumnsShape(volRegion)]),[3 2 1]);
		% column-wise C-code (used only for testing)
		%	offset = numColumnsShape(volRegion)*[0:numBounds-2];
		% initialization only has to be made for columns relvant for updating the q(b) distribution
		%	for i = 1:numColumnsShape(volRegion)
		%		% the -1 is C-indexing
		%		idxA = sum(numColumnsShape(1:volRegion-1))*numBounds + i - 1;
		%		idxA = [idxA+offset; idxA+offset+numColumnsShape(volRegion)];
		%		pObs = squeeze(prediction(:,:,columnsShapePred{volRegion}(i),volRegion));
		%		
		%		q_c.singleton2(volRegion,columnsShapePred{volRegion}(i),:,:) = reshape(sumProductSparseBackup(pObs,models.shapeModel(volRegion).mu,models.shapeModel(volRegion).WML,models.shapeModel(volRegion).sigmaML,int32(idxA),hashTable),[numRows numBounds])';
		%	end
	end
end
q_c.singleton(q_c.singleton < options.threshold_q_c) = 0;

if collector.options.printTimings
	fprintf('[Initialized q_c]: %.3fs\n',toc(initqc));
end



