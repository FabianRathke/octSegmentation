% fetch appearance Terms
predictionA = predAppearance(files(file),collector,models.appearanceModel,options);
prediction = zeros(numRows,numClasses,numColumnsPred,numVolRegions);

for volRegion = 1:numVolRegions
    % load ground truth and the scan itself
	collector.options.labelID = collector.options.labelIDs(file,volRegion);
	% load labels if set by the user
	if collector.options.loadLabels
		output.trueLabels{file,volRegion} = loadLabels(files(file).name,collector.options);
	else % otherwise provide some dummy labels
		output.trueLabels{file,volRegion} = zeros(collector.options.numLayers-1,collector.options.X);
	end

	if options.plotting
		eval(sprintf('B%d = loadData(files(file).name,collector.options);',collector.options.labelID));
	end
    prediction(:,:,:,volRegion) = permute(reshape(predictionA.prediction{volRegion},[numClasses,collector.options.Y,numColumnsPred]),[2 1 3]);
end
% we only need the subset of boundary classes
prediction = prediction(:,length(collector.options.LayersPred) + (1:length(collector.options.EdgesPred)),:,:);
% reduce the precision
prediction(prediction < options.thresholdAccuracy) = 0;

% outputs pixelwise appearance terms too
if collector.options.saveAppearanceTerms
    output.appearanceTerms.prediction{file} = single(predictionA.prediction{1});
    output.appearanceTerms.trueLabels{file} = single(predictionA.trueLabels{1});
end

clear predictionA;
boundaries = zeros(numVolRegions,numBounds,numColumnsPred);

initqc = tic;
%if 0
%	if ~isfield(models.shapeModel,'pTransV')
%		models.shapeModel = preCalcTransitionMatrices(collector,models.shapeModel);
%	end
%	for volRegion = 1:numVolRegions
%		% initialization only has to be made for columns relvant for updating the q(b) distribution
%		for i = 1:numColumnsShape(volRegion)
%			idxA = sum(numColumnsShape(1:volRegion-1))*numBounds + i;
%			pObs = squeeze(prediction(:,:,columnsShapePred{volRegion}(i),volRegion));
%
%			variance = sum(models.shapeModel.WML(idxA,:).^2) + models.shapeModel.sigmaML;
%			% calculate probabilities for first boundary
%			pStart = 1/sqrt(2*pi*variance)*exp(-0.5*(1/variance)*((1:numRows) - models.shapeModel.mu(idxA)).^2);
%			q_c.singleton(volRegion,columnsShapePred{volRegion}(i),:,:) = sumProductSparse(pStart,models.shapeModel.pTransV{volRegion}(i,:),pObs)';
%		end
%	end
%end

% c version
for volRegion = 1:numVolRegions
	% the -1 is C-indexing
    idxA = sum(numColumnsShape(1:volRegion-1))*numBounds + (1:numColumnsShape(volRegion)) - 1;
    q_c.singleton(volRegion,columnsShapePred{volRegion},:,:) = permute(reshape(sumProductSparseC(prediction(:,:,columnsShapePred{volRegion},volRegion),models.shapeModel(volRegion).mu,models.shapeModel(volRegion).WML,models.shapeModel(volRegion).sigmaML,int32(idxA),hashTable),[numRows,numBounds,numColumnsShape(volRegion)]),[3 2 1]);
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
q_c.singleton(q_c.singleton < options.threshold_q_c) = 0;

if collector.options.printTimings
	fprintf('[Initialized q_c]: %.3fs\n',toc(initqc));
end

