function results = signedUnsigned(prediction,options,collector)

options = setEvalDefaults(options);

numFiles = size(prediction.prediction,1);
numRegions = size(prediction.prediction,2);
numBounds = size(prediction.prediction{1},1);
if options.interpolation
	numColumns = size(prediction.trueLabels{1},2);
	columnsInterp = setdiff(1:numColumns,prediction.columnsPred);
else
	numColumns = size(prediction.prediction{1},2);
end
div = numBounds*numColumns;

for i = 1:numFiles
	for j = 1:numRegions
		if options.interpolation
			for k = 1:numBounds
				predInter(k,columnsInterp) = interp1(prediction.columnsPred,prediction.prediction{i,j}(k,:),columnsInterp,'spline','extrap');
				predInter(k,prediction.columnsPred) = prediction.prediction{i,j}(k,:);
				predInterInit(k,columnsInterp) = interp1(prediction.columnsPred,prediction.prediction_init{i,j}(k,:),columnsInterp,'spline','extrap');
				predInterInit(k,prediction.columnsPred) = prediction.prediction_init{i,j}(k,:);
			end
			results.quadratic(i,j) = sum(sum(abs(predInter-prediction.trueLabels{i,j}).^2,2))/div;
			results.unsigned(i,j) = sum(sum(abs(predInter-prediction.trueLabels{i,j}),2))/div;
			results.unsigned_per_layer(:,i,j) = sum(abs(predInter-prediction.trueLabels{i,j}),2)/numColumns;
			results.unsigned_per_column(i,j,:,:) = abs(predInter-prediction.trueLabels{i,j});
			results.unsigned_init(i,j) = sum(sum(abs(predInterInit-prediction.trueLabels{i,j}),2))/div;
			results.unsigned_init_per_layer(:,i,j) = sum(abs(predInterInit-prediction.trueLabels{i,j}),2)/numColumns;
			results.signed(i,j) = sum(sum(predInter-prediction.trueLabels{i,j},2))/div;
			results.signed_per_layer(:,i,j) = sum(predInter-prediction.trueLabels{i,j},2)/numColumns;
		else
			results.quadratic(i,j) = sum(sum(abs(prediction.prediction{i,j}-prediction.trueLabels{i,j}(:,prediction.columnsPred)).^2,2))/div;
			results.unsigned(i,j) = sum(sum(abs(prediction.prediction{i,j}-prediction.trueLabels{i,j}(:,prediction.columnsPred)),2))/div;
			results.unsigned_per_layer(:,i,j) = sum(abs(prediction.prediction{i,j}-prediction.trueLabels{i,j}(:,prediction.columnsPred)),2)/numColumns;
			results.unsigned_per_column(i,j,:,:) = abs(prediction.prediction{i,j}-prediction.trueLabels{i,j}(:,prediction.columnsPred));
			results.unsigned_init(i,j) = sum(sum(abs(prediction.prediction_init{i,j}-prediction.trueLabels{i,j}(:,prediction.columnsPred)),2))/div;
			results.unsigned_init_per_layer(:,i,j) = sum(abs(prediction.prediction_init{i,j}-prediction.trueLabels{i,j}(:,prediction.columnsPred)),2)/numColumns;
			results.signed(i,j) = sum(sum(prediction.prediction{i,j}-prediction.trueLabels{i,j}(:,prediction.columnsPred),2))/div;
			results.signed_per_layer(:,i,j) = sum(prediction.prediction{i,j}-prediction.trueLabels{i,j}(:,prediction.columnsPred),2)/numColumns;
		end
	end
end

results.performance = -mean(mean(results.unsigned));
