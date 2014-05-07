ColorOrder = [0 0 1.0000;	
	 0 0.5000 0;
	1.0000         0         0;
	0    0.7500    0.7500;
	0.7500         0    0.7500;
	0.7500    0.7500         0;
	0.8000    0.8000    0.8000];

options = setSegmentationDefaults(options,params);

% name used in plots containing current parameter combination
names = fieldnames(params);
filename = '';
for i = 1:length(names)
    filename = [filename sprintf('_%s=%.6f',names{i},getfield(params,names{i}))];
end

% some standard variables
numVolRegions = collector.options.numRegionsPerVolume;
numBounds = length(collector.options.EdgesPred);
numRows = collector.options.Y;
numClasses = length(collector.options.EdgesPred) + length(collector.options.LayersPred);
numColumnsPred = length(collector.options.columnsPred);
for i = 1:length(collector.options.columnsShape)
	numColumnsShape(i) = length(collector.options.columnsShape{i});
	% find the closest prediction-column to each shape-column
	for column = 1:numColumnsShape(i)
		[a b] = min(abs(collector.options.columnsShape{i}(column) - collector.options.columnsPred));
		columnsShapePred{i}(column) = b;
	end
	% find the to closest Shape columns for each Pred column for interpolation
	for column=1:numColumnsPred
		[a b] = sort(abs(collector.options.columnsPred(column) - collector.options.columnsShape{i}));
		columnsPredShape{i}(:,column) = b(1:2);
		columnsPredShapeFactor{i}(:,column) = a(2:-1:1)./(a(1)+a(2));
	end
end

columnsPredShapeVec = [columnsPredShape{:}];
columnsPredShapeFactorVec = [columnsPredShapeFactor{:}];
numColumnsShapeTotal = sum(numColumnsShape);

output.prediction = cell(length(files),numVolRegions);
output.trueLabels = cell(length(files),numVolRegions);
options.positions = 1:numRows;

% save for later use, save as short form for better code reading
WML = models.shapeModel.WML; sigmaML = models.shapeModel.sigmaML;
M = inv(WML'*WML + sigmaML*eye(size(WML,2)));
models.shapeModel.M = M;
WML = eval(sprintf('%s(WML)',collector.options.dataTypeCast));
M = eval(sprintf('%s(M)',collector.options.dataTypeCast));

% calculate factors for mu_a_b where a and b are neighbouring boundaries in the same column
factorsPrecA = cell(1,numVolRegions);
factorsPrecB = cell(1,numVolRegions);
for volRegion = 1:numVolRegions
	idxB = (1:numColumnsShape(volRegion)*(numBounds-1)) + sum(numColumnsShape(1:volRegion-1))*numBounds;
	idxA = idxB + numColumnsShape(volRegion);
	
	for j = 1:length(idxB)
		P = inv(models.shapeModel.WML([idxB(j) idxA(j)],:)*models.shapeModel.WML([idxB(j) idxA(j)],:)' + sigmaML*eye(2));
		factorsPrecA{volRegion}(j) = P(2,2);
		factorsPrecB{volRegion}(j) = P(2,1);
	end
end

factorsPrecAVec = [factorsPrecA{:}];
factorsPrecBVec = [factorsPrecB{:}];

% for the c-sum-product
hashTable = sort(exp(-10000:0.01:0),'descend');


% initialize 
pred_update_varinit;

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
mu_a_b2 =  eval(sprintf('%s(mu_a_b2)',collector.options.dataTypeCast));
factorsPrecAVec = eval(sprintf('%s(factorsPrecAVec)',collector.options.dataTypeCast));

% auslagerung calcOT
X = 1:numRows;
X = X(ones(1,numColumnsShapeTotal*numBounds),:);
if collector.options.calcOnGPU
	X = colon(1,1,numRows,GPUsingle);
	X = X(ones(1,numColumnsShapeTotal*numBounds),:);
end
