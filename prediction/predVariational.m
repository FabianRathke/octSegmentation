function output = predVariational(files,collector,params,models,options)
% predVariational - performs variational inference given appearance models and a shape prior model; collects test data as specified in files and collector
%
% Syntax:
%   output = predVariational(files,collector,params,models,options) 
%
% Inputs:
%   files     - [struct] files for which predictions will be made
%   collector - [struct] variables that control the collection of training and test data; see setCollectorsDefault for a description of the variables
%      .options.labelIDs
%      .options.numRegionsPerVolume
%   params    - [struct]
%   models    - [struct] appearance models for all files and BScans; have to be compatible with the function called in trainModels.m 
%   options   - [struct]
%      .plotting       - [boolean] plots are made of 
%      .folderPlots    - [string]
%      .calcFuncVal    - [boolean]
%      .detailedOutput - [boolean] 
%      .iterations     - [int]
%      .threshold      - [float]
%      .threshold_q_c  - [float] 
%
% Outputs:
%   output - [struct]
%     .trueLabels - [cell] pixel-wise labels for all B-Scans
%     .prediction - [cell] predicted segmentations for all B-Scans
%     .prediction_init - [cell] segmentation after initialization for all B-Scan
%     .columnsPred - [array] columns for which predictions are produced
%
% See also: predAppearance, setSegmentationDefaults

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 05-May-2014

initSegmentation;

for file = 1:length(files)
	fprintf('... make predictions for ** %s **\n',files(file).name);
	
	if options.plotting
		% check folder to store plots for current file
		folderName = checkForFolder(sprintf('%s/%s',options.folderPlots,files(file).name));
	end		

	boundaries = zeros(numVolRegions,numBounds,numColumnsPred);
	q_c.singleton = zeros(numVolRegions,numColumnsPred,numBounds,numRows);

	% the pairwise terms of q_c are only needed, when calculating the mutual information
	if options.calcFuncVal
		q_c.pairwise = cell(numVolRegions,numColumnsPred,numBounds-1);
		for region = 1:numVolRegions
			for j = 1:numColumnsPred
				for k = 1:numBounds-1
					q_c.pairwise{region,j,k} = sparse(numRows,numRows);
				end
			end
		end
	end

	initMarginalsWithHMM;

	error_init = 0;
	for volRegion = 1:collector.options.numRegionsPerVolume
		q_c_plot = squeeze(sum(permute(squeeze(q_c.singleton(volRegion,:,:,:)),[2 3 1]).*repmat(1:numRows,[numBounds,1,numColumnsPred]),2));

		old_boundaries(volRegion,:,:) = q_c_plot;
		output.prediction_init{file,volRegion} = q_c_plot;
		if collector.options.loadLabels
			error_init = error_init + sum(sum(abs(output.prediction_init{file,volRegion}(:,columnsShapePred{volRegion})-output.trueLabels{file,volRegion}(:,collector.options.columnsShape{volRegion}))))/(numColumnsShape(volRegion)*numBounds);
			fprintf('Initial unsigned error: %.3fpx\n',error_init/numVolRegions);
		else 
			error_init = 0;
		end

		if options.plotting
			fileSaveName = sprintf('%s/init/qc_0%s_%d.eps',folderName,filename,collector.options.labelIDs(file,volRegion));
			eval(sprintf('plotBScan(B%d,q_c_plot(:,columnsShapePred{volRegion}),collector.options.columnsShape{volRegion},fileSaveName)',collector.options.labelIDs(file,volRegion)))
		end
	end
	
	% calculates the constant contributions of the shape prior to the omega-matrices
	calcShapeTerms;

	for iter = 1:options.iterations
		fprintf('Iteration %d\n',iter);
		% save last iterations prediction for convergence criteria
		if iter > 1
			old_boundaries = boundaries;
		end

		% P and p terms
		calcDer;
		% optimize q_b
		optQB;
		% intializes omega matrices 
		calcOT;
		% optimize q_c
		optQC;

		change_per_iteration(iter) = 0;
		for volRegion = 1:numVolRegions
			boundaries(volRegion,:,:) = squeeze(sum(permute(squeeze(q_c.singleton(volRegion,:,:,:)),[2 3 1]).*repmat(1:numRows,[numBounds,1,numColumnsPred]),2));
			change_per_iteration(iter) =  change_per_iteration(iter) + mean(mean(mean(abs(boundaries(volRegion,:,columnsShapePred{volRegion})-old_boundaries(volRegion,:,columnsShapePred{volRegion})))));
		end
		fprintf('Mean change: %.3fpx\n',change_per_iteration(iter));

		for volRegion = 1:numVolRegions
			if options.detailedOutput
				idx = (1:numBounds*numColumnsShape(volRegion)) + numBounds*sum(numColumnsShape(1:volRegion-1));
				q_b_error(iter,volRegion,:) = sum(abs(single(reshape(q_b.mu(idx),numColumnsShape(volRegion),numBounds)) - output.trueLabels{file,volRegion}(:,collector.options.columnsShape{volRegion})'))/numColumnsShape(volRegion);
				q_c_error(iter,volRegion,:) = sum(abs(squeeze(sum(permute(squeeze(q_c.singleton(volRegion,:,:,:)),[2 3 1]).*repmat(1:numRows,[numBounds,1,numColumnsPred]),2)) -output.trueLabels{file,volRegion}(:,collector.options.columnsPred)),2)/numColumnsPred;
			end
			unsigned_error(iter,volRegion) = sum(sum(abs(squeeze(boundaries(volRegion,:,:)) - output.trueLabels{file,volRegion}(:,collector.options.columnsPred))))/(numColumnsPred*numBounds); 	

			if isnan(unsigned_error(iter,:))
				fprintf('NaN detected, aborting prediction\n');
				boundaries(iter,:,:) = old_boundaries;
				break;
			end	
		end
		fprintf('Unsigned error: %.3fpx\n',mean(unsigned_error(iter,:)));

		% after convergence, save the final result and quit the optimization
		if (iter == options.iterations || (iter > 1 && change_per_iteration(iter) < options.threshold))
			output.columnsPred = collector.options.columnsPred;
			for volRegion = 1:numVolRegions
				output.prediction{file,volRegion} = squeeze(boundaries(volRegion,:,:));
			end
			if options.detailedOutput
				output.change_per_iteration{file} = change_per_iteration; change_per_iteration = [];
				output.unsigned_error{file} = unsigned_error; unsigned_error = [];		
				output.q_c_error{file} = q_c_error; q_c_error = [];
				output.q_b_error{file} = q_b_error; q_b_error = [];
				output.q_b{file} = q_b.mu;
			end
			% calc the terms of the objective function
			if options.calcFuncVal
				calcFuncVals;
				output.funcVal(file) = funcVal;
			end

			% clean after mex-functions
			clear mex

			break
		end
	end
end
