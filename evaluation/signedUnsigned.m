function results = signedUnsigned(prediction,options,collector)
% signedUnsigned - given ground truth and predictions both contained in the struct prediction, returns the signed and unsigned error for the whole scan, for each boundary respectively or for each image column seperately.
%  
% Syntax:
%   options = signedUnsigned(prediction,options,collector)
%
% Inputs:
%   prediction - [struct] contains fields prediction and trueLabels which are cell arrays with one entry for each file
%   options    - [struct] options struct
%     .interpolation - [boolean] shall unpredicted columns be interpolated and included in the error measure?
%   collector  - [struct] not needed for this function
%
% Outputs:
%  results - [struct] contains fields for unsigned and signed error
%
% See also: predEvalDefaults

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 30-Apr-2014

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
		if isfield(prediction,'trueLabels')
			% interpolated unpredicted columns and calculate error measure for the whole scan
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
			% without interpolation only obtain error measure of columns for which predictions where obtained
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
end

results.performance = -mean(mean(results.unsigned));
