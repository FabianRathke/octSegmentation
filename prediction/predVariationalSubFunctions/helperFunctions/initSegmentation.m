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
end

% name used in plots containing current parameter combination
names = fieldnames(params);
filename = '';
for i = 1:length(names)
    filename = [filename sprintf('_%s=%.6f',names{i},getfield(params,names{i}))];
end

% some standard variables
%numVolRegions = collector.options.numRegionsPerVolume;
numVolRegions = length(collector.options.labelIDs);
% quick hack --> change for final version: check which model components are required and modify the model to fit these requirements;
% better --> just save an array of model indices and leave the model itself untouched
if length(models.shapeModel) > 1
	for i = 1:numVolRegions
		% calculate distance of i'th Bscan from the center scan in the volume (is at zero px)
		dist = double((collector.options.labelIDs(i) - (floor(collector.options.numBScansPerVolume/2)+1))*collector.options.distBScans);
		% find closest B-scan in the shape model
		[~,scanID] = min(abs(models.params.BScanPositions-dist));
		fprintf('DEBUG: Distance to central scan is %.fpx, picked %d''th component of shape model\n',dist,scanID);
		modelTmp.appearanceModel(i,:,:) = models.appearanceModel(scanID,:,:);
		modelTmp.shapeModel(i) = models.shapeModel(scanID);
	end
	% set modified model as the new model
	models = modelTmp;
end

if length(models.shapeModel) > 1
	error('Reimplement full 3D shape prior');
end

% ************** DETERMINE CONNECTION BETWEEN QB AND QC ******************
numBounds = length(collector.options.EdgesPred);
collector.options.columnsShape = cell(1,numVolRegions);
for i = 1:numVolRegions
	% cut shape model to subregion defined by columnsPred
	minMaxColumns = [min(collector.options.columnsPred) max(collector.options.columnsPred)];
	subsetColumns = find(models.shapeModel(i).columnsShape>=minMaxColumns(1) & models.shapeModel(i).columnsShape<=minMaxColumns(2));
	%  
	numColumns = length(models.shapeModel.columnsShape);
	subsetColumnsFull = reshape((repmat(subsetColumns,numBounds,1)+repmat((0:numBounds-1)'.*numColumns,1,length(subsetColumns)))',1,[]);

	models.shapeModel(i).mu = models.shapeModel(i).mu(subsetColumnsFull);
	models.shapeModel(i).WML = models.shapeModel(i).WML(subsetColumnsFull,:);
	models.shapeModel(i).columnsShape = models.shapeModel(i).columnsShape(subsetColumns);

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

% save for later use, save as short form for better code readability
WML = models.shapeModel.WML; sigmaML = models.shapeModel.sigmaML;
M = inv(WML'*WML + sigmaML*eye(size(WML,2)));
models.shapeModel.M = M;
WML = eval(sprintf('%s(WML)',collector.options.dataTypeCast));
M = eval(sprintf('%s(M)',collector.options.dataTypeCast));
prodWM = WML*M;
prodWMT = prodWM';

% calculate factors for mu_a_b where a and b are neighbouring boundaries in the same column
factorsPrecA = cell(1,numVolRegions);
factorsPrecB = cell(1,numVolRegions);
for volRegion = 1:numVolRegions
	idxB = (1:numColumnsShape(volRegion)*(numBounds-1)) + sum(numColumnsShape(1:volRegion-1))*numBounds;
	idxA = idxB + numColumnsShape(volRegion);

	d = sum(WML(idxA,:).*WML(idxA,:),2) + sigmaML;
	b = sum(WML(idxA,:).*WML(idxB,:),2);
 	a = sum(WML(idxB,:).*WML(idxB,:),2) + sigmaML;

 
	tmp = 1./(a.*d - b.*b);
	factorsPrecA{volRegion} = (tmp.*a)';
	factorsPrecB{volRegion} = (-tmp.*b)';

%	for j = 1:length(idxB)
%		P = inv(models.shapeModel.WML([idxB(j) idxA(j)],:)*models.shapeModel.WML([idxB(j) idxA(j)],:)' + sigmaML*eye(2));
%		factorsPrecA{volRegion}(j) = P(2,2);
%		factorsPrecB{volRegion}(j) = P(2,1);
%	end
end

factorsPrecAVec = [factorsPrecA{:}];
factorsPrecBVec = [factorsPrecB{:}];

% for the c-sum-product
hashTable = sort(exp(-1000:0.001:0),'descend');


% initialize 
%pred_update_varinit;

% auslagerung optQC
% calculate mu_a_b for all regions at once
idxNotB = []; idxNotA = [];
for volRegion = 1:numVolRegions
	idxNotB = [idxNotB (1:numColumnsShape(volRegion))+sum(numColumnsShape(1:volRegion-1))*numBounds];
	idxNotA = [idxNotA (1:numColumnsShape(volRegion))+sum(numColumnsShape(1:volRegion-1))*numBounds+(numBounds-1)*numColumnsShape(volRegion)];
end
idxA = 1:numColumnsShapeTotal*numBounds; idxB = idxA;
idxB(idxNotB) = []; idxA(idxNotA) = [];

A = 1:numRows;
mu_a_b2 = (models.shapeModel.mu(idxB,ones(1,numRows)) - factorsPrecAVec(ones(1,numRows),:)'.^-1.*factorsPrecBVec(ones(1,numRows),:)'.*(A(ones(1,length(idxA)),:)-models.shapeModel.mu(idxA,ones(1,numRows))));
%mu_a_b2 =  eval(sprintf('%s(mu_a_b2)',collector.options.dataTypeCast));
%factorsPrecAVec = eval(sprintf('%s(factorsPrecAVec)',collector.options.dataTypeCast));

% auslagerung calcOT
X = 1:numRows;
X = X(ones(1,numColumnsShapeTotal*numBounds),:);
if collector.options.calcOnGPU
	X = colon(1,1,numRows,GPUsingle);
	X = X(ones(1,numColumnsShapeTotal*numBounds),:);
end
