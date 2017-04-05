checkPredictionGlobal(files(file).name,collector.options.labelIDs(file,volRegion));

global isStoredInGlobal
boundsPred = [ones(numColumnsShape,1) ones(numColumnsShape,1)*numRows]';
boundsPredAll = [ones(numColumnsPred,1) ones(numColumnsPred,1)*numRows]';
% Calculate appearance terms
% Use globally stored pre-calculated appearance terms
if isStoredInGlobal.prediction
	% check if all columns that are to be predicted are contained in the globally saved prediction
	if sum(~ismember(collector.options.columnsPred,predictionGlobal.columns{collector.options.labelIDs(file,volRegion)})) == 0
		fprintf('Reusing appearance terms stored in global variable\n');
		% calculate indices
		prediction = predictionGlobal.data{collector.options.labelIDs(file,volRegion)}(:,:,ismember(predictionGlobal.columns{collector.options.labelIDs(file,volRegion)},collector.options.columnsPred));
		if isfield(predictionGlobal,'dataFuncVal') & options.calcFuncVal
			predictionFuncVal = predictionGlobal.dataFuncVal{collector.options.labelIDs(file,volRegion)}(:,:,ismember(predictionGlobal.columns{collector.options.labelIDs(file,volRegion)},collector.options.columnsPred));
		end
		if isfield(predictionGlobal,'textureQuality')
			textureQuality = predictionGlobal.textureQuality{collector.options.labelIDs(file,volRegion)};
		end
		if isfield(predictionGlobal,'boundsPred') && ~isempty(predictionGlobal.boundsPred{collector.options.labelIDs(file,volRegion)})
			boundsPredAll = round(predictionGlobal.boundsPred{collector.options.labelIDs(file,volRegion)}(:,collector.options.columnsPred) + [-50*ones(1,numColumnsPred); +50*ones(1,numColumnsPred)]);
			boundsPredAll(1,:) = max(1,boundsPredAll(1,:));
			boundsPredAll(2,:) = min(numRows,boundsPredAll(2,:));
			boundsPred = boundsPredAll(:,ismember(collector.options.columnsPred,collector.options.columnsShape{1}));
		end
	else
		fprintf('Reusing of appearance terms failed, calculating them anew! (Information is not stored for all image columns).\n');
		isStoredInGlobal.prediction = 0;
	end
% No global variable: Calculate them
else
	fprintf('No cached appearance terms found, calculating them\n');
	if ~isfield(options,'appearance') && ~isStoredInGlobal.prediction
		% use prediction on a subset of columns to restrict the area that has to be taken into consideration
		if length(collector.options.columnsPred) > 100 && isfield(collector.options,'margins') && 0
			volRegion = 1;
			if numVolRegions > 1
				error('Not implemented yet for 3D models');
			end
			columnsPred = collector.options.columnsPred; % save for later use
			subVec = round(linspace(1+10,numColumnsPred-10,15)); % the subset of columns used for the initial prediction
			collector.options.columnsPred = collector.options.columnsPred(subVec);

			% fetch appearance Terms
			predictionSparse = predAppearance(files(file),collector,models.appearanceModel,options);
			prediction = permute(reshape(predictionSparse.prediction{volRegion},[numClasses,collector.options.Y,length(collector.options.columnsPred)]),[2 1 3]);
			prediction = prediction(:,length(collector.options.LayersPred) + (1:length(collector.options.EdgesPred)),:);

			% predict q_c
			idxA = columnsPredShape{volRegion}(1,subVec)-1;
	%		idxA = (1:length(subVec))-1; idxB = repmat(columnsPredShape{volRegion}(1,subVec),numBounds,1)' + repmat((0:(numBounds-1))*numColumnsShape,length(subVec),1);
			boundsPred = [zeros(1,length(idxA)); ones(1,length(idxA))*(numRows-1)];
			q_c_init = permute(reshape(sumProductSparseC(prediction(:,:,:),models.shapeModel(volRegion).mu,models.shapeModel(volRegion).WML,models.shapeModel(volRegion).sigmaML,int32(idxA),hashTable,int32(boundsPred)),[numRows,numBounds,length(subVec)]),[3 2 1]);

			% obtain estimated boundaries
			z = size(q_c_init);
			boundsInit = squeeze(sum(permute(q_c_init,[2 3 1]).*repmat(1:numRows,[numBounds,1,length(subVec)]),2))';

			% obtain quality of estimate
			for k = 1:numBounds
				for j = 1:length(subVec)
					I = find(q_c_init(j,k,:)>=10^-30);
					tmp = squeeze(q_c_init(j,k,I))'.*squeeze(log(prediction(I,k,j)))';
					q_c_data(j,k) = -full(sum(tmp(~isnan(tmp)&~isinf(tmp))));
				end
			end

			% which entries do we trust?
			idxTrust = q_c_data - repmat(collector.options.margins.mu(4,:),length(subVec),1)-repmat(collector.options.margins.sigma(4,:),length(subVec),1)*1 < 0;

			% if several b-scan columns map to the same shape prior column we select a unique one here
			[C IA] = unique(columnsShapePred{1}(subVec));
			idxTrust = idxTrust(IA,:); boundsInit = boundsInit(IA,:);
			idxGiven = repmat(C',1,numBounds)+repmat((0:(numBounds-1))*numColumnsShape,length(C),1);
			idxGiven = idxGiven(idxTrust);
			% get indices from the first and last boundary
			idxCond = setdiff([[1:numColumnsShape] [1:numColumnsShape]+numColumnsShape*(numBounds-1)],[C(idxTrust(:,1)),C(idxTrust(:,end))+numColumnsShape*(numBounds-1)]); 
			% obtain posterior mean
			Sigma_b_b = sigmaML*eye(length(idxGiven)) + WML(idxGiven,:)*WML(idxGiven,:)'; 
			Sigma_a_b = WML(idxCond,:)*WML(idxGiven,:)'; 
			mu_a_b = models.shapeModel.mu(idxCond) + Sigma_a_b*inv(Sigma_b_b)*(boundsInit(idxTrust)-models.shapeModel.mu(idxGiven));
			idxGiven = [C(idxTrust(:,1)),C(idxTrust(:,end))+numColumnsShape];	
			idxCond(idxCond>numColumnsShape)= idxCond(idxCond>numColumnsShape) - numColumnsShape*(numBounds-2);

			boundsFull(idxCond) = mu_a_b;
			boundsFull(setdiff(1:numColumnsShape*2,idxCond)) = [boundsInit(idxTrust(:,1),1); boundsInit(idxTrust(:,end),end)];
			boundsPred = round(reshape(boundsFull,numColumnsShape,2) + [-25*ones(numColumnsShape,1) 25*ones(numColumnsShape,1)]);
			collector.options.columnsPred = columnsPred;

			% create idxSet that contains positions of patches for the reduced image region
			idxSet = zeros(2,sum(boundsPred(:,2)-boundsPred(:,1)+1));
			tmp = [arrayfun(@(x,y,z) ones(1,y-x+1)*z,boundsPred(:,1),boundsPred(:,2),collector.options.columnsPred','UniformOutput',false)];
			idxSet(2,:) = [tmp{:}];
			tmp = [arrayfun(@(x,y) colon(x,y),boundsPred(:,1),boundsPred(:,2),'UniformOutput',false)];
			idxSet(1,:) = [tmp{:}];
			collector.options.idxSet = idxSet;
			
			% fetch appearance terms
			prediction = predAppearance(files(file),collector,models.appearanceModel,options);
			predictionA.prediction{1} = zeros(size(prediction.prediction{1},1),numColumnsPred*collector.options.Y);
			% change the index over columns from the real B-Scan to columnsPred
			idxChange = zeros(1,collector.options.X);
			idxChange(collector.options.columnsPred) = 1:numColumnsPred;
			idxSet(2,:) = idxChange(idxSet(2,:));
			IND = sub2ind([collector.options.Y numColumnsPred],idxSet(1,:),idxSet(2,:));
			predictionA.prediction{1}(:,IND) = prediction.prediction{1};
		else
			% predict all rows in each column
			predictionA = predAppearance(files(file),collector,models.appearanceModel,options);
	
			isfieldName = 0;
%			for j = 1:length(models.appearanceModel)
				if isfield(models.appearanceModel{1},'prior')
					isfieldName = 1;
				end
%			end
					
			if options.calcFuncVal && isfieldName
				modelTmp = models.appearanceModel(1);
				modelTmp{1} = rmfield(modelTmp{1},'prior');
				predictionFuncValA = predAppearance(files(file),collector,modelTmp,options);
			end
		end
	elseif isfield(options,'appearance')
		predictionA.prediction = options.appearance;
	end

	prediction = zeros(numRows,numClasses,numColumnsPred,numVolRegions);
	for volRegion = 1:numVolRegions
		prediction(:,:,:,volRegion) = permute(reshape(predictionA.prediction{volRegion},[numClasses,collector.options.Y,numColumnsPred]),[2 1 3]);
	end
	% we only need the subset of boundary classes
	prediction = prediction(:,length(collector.options.LayersPred) + (1:length(collector.options.EdgesPred)),:,:);
	% reduce the precision
	%prediction(prediction < options.thresholdAccuracy) = 0;

	if isfield(predictionA,'predictionFirstModel')
		predictionFuncValA.prediction = predictionA.predictionFirstModel;
	end

	if exist('predictionFuncValA')
		predictionFuncVal = zeros(numRows,numClasses,numColumnsPred,numVolRegions);
		for volRegion = 1:numVolRegions
			if iscell(predictionFuncValA.prediction)
				predictionFuncVal(:,:,:,volRegion) = permute(reshape(predictionFuncValA.prediction{volRegion},[numClasses,collector.options.Y,numColumnsPred]),[2 1 3]);
			else
				predictionFuncVal(:,:,:,volRegion) = permute(reshape(predictionFuncValA.prediction,[numClasses,collector.options.Y,numColumnsPred]),[2 1 3]);
			end
		end
		% we only need the subset of boundary classes
		predictionFuncVal = predictionFuncVal(:,length(collector.options.LayersPred) + (1:length(collector.options.EdgesPred)),:,:);
		% reduce the precision
		predictionFuncVal(predictionFuncVal < options.thresholdAccuracy) = 0;
	end
	
	% returns pixelwise appearance terms to the user as part of the output struct
	if collector.options.saveAppearanceTerms
		output.appearanceTerms.prediction{file} = single(predictionA.prediction{1});
		if collector.options.loadLabels
			output.appearanceTerms.trueLabels{file} = single(predictionA.trueLabels{1});
		end
	end
	clear predictionA;
end

% calc conditional means for transition matrices: the mean of boundary k conditioned on the preceding boundary k-1 in the same image column
factorMuAB2 = [factorsPrecAVec(idxB-numColumnsShape).^-1.*factorsPrecBVec(idxB-numColumnsShape)]';
mu_a_b2 = calcMuAB2(models.shapeModel.mu,factorMuAB2,numRows,numBounds,numColumnsShape,numColumnsPred,int32(columnsPredShapeVec(1,:)-1),columnsPredShapeFactorVec(1,:),int32(columnsPredShapeVec(2,:)-1),columnsPredShapeFactorVec(2,:),int32(boundsPredAll-1));
%clear factorsPrecBVec;
%if isfield(options,'doNotPredictShape')
%    A = 1:numRows;
%    mu_a_b2 = zeros(numColumnsShape*(numBounds-1),numRows);
%    mu_a_b2(idxB-numColumnsShape,:) = models.shapeModel.mu(idxB,ones(1,numRows)) - factorMuAB2(:,ones(1,numRows)).*(A(ones(1,length(idxA)),:)-models.shapeModel.mu(idxA,ones(1,numRows)));
%end
factorsPrecAvg = reshape((repmat(columnsPredShapeFactorVec(1,:),numBounds-1,1).*factorsPrecAVec(repmat(columnsPredShapeVec(1,:),numBounds-1,1) + repmat((0:numColumnsShape:numColumnsShape*(numBounds-2))',1,numColumnsPred)) + repmat(columnsPredShapeFactorVec(2,:),numBounds-1,1).*factorsPrecAVec(repmat(columnsPredShapeVec(2,:),numBounds-1,1) + repmat((0:numColumnsShape:numColumnsShape*(numBounds-2))',1,numColumnsPred)))',[],1);

% set the appearance models for layers in the white areas to almost zero
if strcmp(collector.options.loadRoutineData,'Srinivasan')
	collector.options.labelID = collector.options.labelIDs(file,volRegion);
	eval(sprintf('BTmp = loadData(files(file).name,collector.options);',collector.options.labelID));
	
	% find connected components of pixels with label 1
	D = bwconncomp(BTmp(:,collector.options.columnsPred)==1);
	idxSet = find(arrayfun(@(x) length(x{1}),D.PixelIdxList) > 1000);
	for jConn = 1:length(idxSet)
		for kk = 1:numBounds
			predictionTmp = squeeze(prediction(:,kk,:));
			predictionTmp(D.PixelIdxList{idxSet(jConn)}) = 10^-200;
			prediction(:,kk,:) = predictionTmp;
		end
	end
end

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
end

%if ~isfield(predictionGlobal,'data')
%	predictionGlobal.data = cell(1,collector.options.labelIDs(file,volRegion));
%end

% saves the appearance terms in a global variable --> can be reused in subsequent runs and requires less user interaction 
if collector.options.makePredGlobal && ~isStoredInGlobal.prediction %isempty(predictionGlobal.data{collector.options.labelIDs(file,volRegion)})
	fprintf('Store appearance terms in global variable\n');
	predictionGlobal.data{collector.options.labelIDs(file,volRegion)} = prediction;
	predictionGlobal.columns{collector.options.labelIDs(file,volRegion)} = collector.options.columnsPred;
	% save the filename in order to prevent using global data for the wrong file
	predictionGlobal.filename = hash(files(file).name,'MD5');
	if exist('predictionFuncVal')
		predictionGlobal.dataFuncVal{collector.options.labelIDs(file,volRegion)} = predictionFuncVal;
	end
end


initqc = tic;

if isfield(options,'doNotPredict')
	% old Matlab code
	models.shapeModel = preCalcTransitionMatrices(collector,models.shapeModel,10^-20,~options.doNotPredict);
	for volRegion = 1:numVolRegions
		% initialization only has to be made for columns relvant for updating the q(b) distribution
		for i = 1:numColumnsShape(volRegion)
			%idxPredict = find(~options.doNotPredict(:,i));
			idxPredict= 1:9;
			if length(idxPredict) > 0
				%idxA = sum(numColumnsShape(1:volRegion-1))*numBounds + i;
				pObs = squeeze(prediction(:,idxPredict,columnsShapePred{volRegion}(i),volRegion));
				idxA = i + (idxPredict(1)-1)*numColumnsShape;
				variance = sum(models.shapeModel.WML(idxA,:).^2) + models.shapeModel.sigmaML;
				% calculate probabilities for first boundary
				pStart = 1/sqrt(2*pi*variance)*exp(-0.5*(1/variance)*((1:numRows) - models.shapeModel.mu(idxA)).^2);
				
				clear pTransTmp;
				for layer = 1:8
					idxA = i + (layer-1)*numColumnsShape;
					idxB = i + layer*numColumnsShape;
					P = inv(WML([idxA idxB],:)*WML([idxA idxB],:)' + eye(2)*sigmaML);
					[iS jS sS numElements] = getCondTransMatrixC(models.shapeModel.mu([idxA idxB]),P,int32(numRows),10^-20,1);
					pTransTmp{layer+1} = sparse(iS(1:numElements),jS(1:numElements),sS(1:numElements),numRows,numRows);
				end
				q_c.singleton(:,idxPredict,columnsShapePred{volRegion}(i),volRegion) = sumProductSparse(pStart,pTransTmp,pObs);
			end
		end
	end
else
	if numVolRegions > 1
		warning('Variable boundsPred has to be updated to multiple regions');
	end
	% C version; segments all columns for one BScan; needs as input the indices of the first boundary inside the shape model
    if isfield(options,'q_c_singleton')
		q_c.singleton = options.q_c_singleton;
	else
		for volRegion = 1:numVolRegions
			% the -1 is C-indexing
			idxA = sum(numColumnsShape(1:volRegion-1))*numBounds + (1:numColumnsShape(volRegion)) - 1;
			q_c.singleton(:,:,columnsShapePred{volRegion},volRegion) = sumProductSparseC(prediction(:,:,columnsShapePred{volRegion},volRegion),models.shapeModel(volRegion).mu,WML,sigmaML,int32(idxA),hashTable,int32(boundsPred)-1);
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
end
q_c.singleton(q_c.singleton < options.threshold_q_c) = 0;

if collector.options.printTimings
	fprintf('[Initialized q_c]: %.3fs\n',toc(initqc));
end
