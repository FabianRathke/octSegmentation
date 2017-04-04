ColorOrder = [0 0 1.0000;	
	 0 0.5000 0;
	1.0000         0         0;
	0    0.7500    0.7500;
	0.7500         0    0.7500;
	0.7500    0.7500         0;
	0.8000    0.8000    0.8000];

options = setSegmentationDefaults(options,params);

% if available, use params specified by the trained model
fieldNames = fieldnames(models.params);
for i = 1:length(fieldNames)
	setfield(collector.options,fieldNames{i},getfield(models.params,fieldNames{i}));
%	eval(sprintf('collector.options.%s = models.params.%s;',fieldNames{i},fieldNames{i}));
end

% name used in plots containing current parameter combination
names = fieldnames(params);
filename = '';
for i = 1:length(names)
    filename = [filename sprintf('_%s=%.6f',names{i},getfield(params,names{i}))];
end

% some standard variables
%numVolRegions = length(collector.options.labelIDs);
numVolRegions = 1;
% quick hack --> change for final version: check which model components are required and modify the model to fit these requirements;
% better --> just save an array of model indices and leave the model itself untouched
modelsTmp.appearanceModel = cell(size(models.appearanceModel));
if length(models.shapeModel) > 1
	for i = 1:numVolRegions
		% directly determine the scanID (i.e. the slice of the 3-D shape prior that is selected)
		if ~isfield(collector.options,'scanID')
			% calculate distance of i'th Bscan from the center scan in the volume (which is at zero px)
	%		dist = double((collector.options.labelIDs(i) - (floor(double(collector.options.numBScansPerVolume)/2)+1))*collector.options.distBScans/collector.options.clipFactor);
			% DEBUG? removed collector.options.clipFactor --> has nothing do to with z-direction
			dist = double((collector.options.labelIDs(i) - (floor(double(collector.options.numBScansPerVolume)/2)+1))*collector.options.distBScans);
			% find closest B-scan in the shape model
			[~,scanID] = min(abs(models.params.BScanPositions-dist));
			fprintf('DEBUG: Distance to central scan is %.2fmm, picked %d''th component of shape model\n',dist,scanID);
		else
			scanID = collector.options.scanID(i);
		end
		for model = 1:length(models.appearanceModel)
			modelTmp.appearanceModel{model}(i,1:size(models.appearanceModel{model},2),:) = models.appearanceModel{model}(scanID,1:size(models.appearanceModel{model},2),:);
		end
		modelTmp.shapeModel(i) = models.shapeModel(scanID);
		output.modelSelect{1}(i) = scanID;
	end
	% set modified model as the new model
	models = modelTmp;
else
	output.modelSelect{1}(1) = 1;
end
numBounds = length(collector.options.EdgesPred);

% adding simple shape structure to model pathologies
if isfield(options,'addDrusenMode')
	if options.addDrusenMode
		models.shapeModel.WML = [models.shapeModel.WML options.modesToAdd];
	end
end

if length(models.shapeModel) > 1
	error('Reimplement full 3D shape prior');
end

% ************** DETERMINE CONNECTION BETWEEN QB AND QC ******************
collector.options.columnsShape = cell(1,numVolRegions);
% cut shape model to subregion defined by columnsPred
for i = 1:numVolRegions
	minMaxColumns = [min(collector.options.columnsPred) max(collector.options.columnsPred)];
	subsetColumns = find(models.shapeModel(i).columnsShape>=minMaxColumns(1) & models.shapeModel(i).columnsShape<=minMaxColumns(2));
	numColumnsShape = length(models.shapeModel.columnsShape);
	subsetColumnsFull = reshape((repmat(subsetColumns,numBounds,1)+repmat((0:numBounds-1)'.*numColumnsShape,1,length(subsetColumns)))',1,[]);

	models.shapeModel(i).mu = models.shapeModel(i).mu(subsetColumnsFull);
	models.shapeModel(i).WML = models.shapeModel(i).WML(subsetColumnsFull,:);
	models.shapeModel(i).columnsShape = models.shapeModel(i).columnsShape(subsetColumns);
%	numColumns = length(models.shapeModel.columnsShape);

	collector.options.columnsShape{i} = models.shapeModel(i).columnsShape;
end	
	
numRows = collector.options.Y;
numClasses = length(collector.options.EdgesPred) + length(collector.options.LayersPred);
numColumnsPred = length(collector.options.columnsPred);
for i = 1:numVolRegions
	numColumnsShape(i) = length(collector.options.columnsShape{i});
	% the mapping prediction-columns to each shape-columns
	for column = 1:numColumnsShape(i)
		[a b] = min(abs(collector.options.columnsShape{i}(column) - collector.options.columnsPred));
		columnsShapePred{i}(column) = b;
	end
	% the mapping shape-columns to prediction-columns
	for column=1:numColumnsPred
		[a b] = sort(abs(collector.options.columnsPred(column) - collector.options.columnsShape{i}));
		columnsPredShape{i}(:,column) = b(1:2);
		columnsPredShapeFactor{i}(:,column) = a(2:-1:1)./(a(1)+a(2));
	end
end

columnsPredShapeVec = [columnsPredShape{:}];
columnsPredShapeFactorVec = [columnsPredShapeFactor{:}];
numColumnsShapeTotal = sum(numColumnsShape);
% ******************** END *****************************************

output.prediction = cell(length(files),numVolRegions);
output.trueLabels = cell(length(files),numVolRegions);
options.positions = 1:numRows;

% indices for leaving out boundaries, used when accessing mu and W from the shapePrior
if isfield(options,'doNotPredict')
	% first, move from predColumns to shapeColumns
	options.doNotPredictShape = options.doNotPredict(:,columnsShapePred{1}(1,:));
    idx_a = find(reshape(options.doNotPredictShape',1,[]));
    idx_b = find(reshape(~options.doNotPredictShape',1,[]));
	idxPredictFull = ~options.doNotPredictShape;
else
	idx_a = [];
    idx_b = 1:numColumnsShape*numBounds;
	idxPredictFull = ones(numBounds,numColumnsShape);
end

% save for later use, save as short form for better code readability
WML = models.shapeModel.WML; sigmaML = models.shapeModel.sigmaML;
if isfield(options,'variance')
	sigmaML = sigmaML*options.variance;
end
if isfield(options,'windowSize')
	WML = [WML [options.windowModes{:}]];
end
M = inv(WML(idx_b,:)'*WML(idx_b,:) + sigmaML*eye(size(WML,2)));
prodWM = WML*M;
prodWMT = prodWM';


% calculate factors for mu_a_b where a and b are neighbouring boundaries in the same column : 2x2 conditional gaussian distributions
factorsPrecA = cell(1,numVolRegions);
factorsPrecB = cell(1,numVolRegions);
for volRegion = 1:numVolRegions
	if isfield(options,'doNotPredictShape')
		% for each image row detect preceding boundaries
		idxA = []; idxB = [];
		for j = 2:numBounds
			for i = 1:numColumnsShape
				if options.doNotPredictShape(j,i) == 0 && sum(~options.doNotPredictShape(1:j-1,i)) > 0
					idxB = [idxB (j-1)*numColumnsShape + i];
					idxA = [idxA (find(options.doNotPredictShape(1:j-1,i)==0,1,'last')-1)*numColumnsShape + i];
				end
			end
		end
	else
		idxA = (1:numColumnsShape(volRegion)*(numBounds-1)) + sum(numColumnsShape(1:volRegion-1))*numBounds;
		idxB = idxA + numColumnsShape(volRegion);
	end
	

	% calculate all entries of symmetric covariance matrix
	d = sum(WML(idxB,:).*WML(idxB,:),2) + sigmaML;
	b = sum(WML(idxB,:).*WML(idxA,:),2);
 	a = sum(WML(idxA,:).*WML(idxA,:),2) + sigmaML;
	
	% manual inversion of a 2x2 matrix
	tmp = 1./(a.*d - b.*b);
	factorsPrecA{volRegion} = zeros(1,numColumnsShape*(numBounds-1)); factorsPrecA{volRegion}(idxB-numColumnsShape) = (tmp.*a)';
	factorsPredB{volRegion} = zeros(1,numColumnsShape*(numBounds-1)); factorsPrecB{volRegion}(idxB-numColumnsShape) = (-tmp.*b)';
end

factorsPrecAVec = [factorsPrecA{:}];
factorsPrecBVec = [factorsPrecB{:}];

% for the c-sum-product
if ~isempty(predictionGlobal)
	if isfield(predictionGlobal,'hashTable')
		hashTable = predictionGlobal.hashTable;
	else
		hashTable = sort(exp(-1000:0.001:0),'descend');
	end
else
	hashTable = sort(exp(-1000:0.001:0),'descend');
end

