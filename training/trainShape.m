function models = trainShape(files,collector,params,options)
% trainShape - collects ground truth and estimates a Gaussian distribution regularized by PPCA; alternatively loads a shape stored on disk
% 
% Syntax:
%   models = trainShape(files,collector,params,options)
%
% Inputs:
%   files     - [struct] output of the matlab function 'dir'
%   collector - [struct] holds the options for the collector grabbing the training patches, description given in setCollectorDefaults
%     .options.columnsShape
%     .options.numRegionsPerVolume
%     .options.Y
%     .options.verbose
%     .options.EdgesTrain
%   params    - [struct] holds params used for example in cross-validation 
%   options   - [struct] holds options used within the appearance function
%     .numModes      - [int] number of eigenvectors to be used
%     .loadShape     - [boolean] loads a stored shape model
%     .loadShapeName - [string] the filename of the shape model to load
%     .saveShape     - [boolean] saves the calculated shape model as mat-file
%     .saveShapeName - [string] the filename for the model to store
%     .shapeFolder   - [string] folder where shape models are stored
%
% Outputs:
%   models - [struct] contains the shape Model
%     .mu      - [array] the mean shape vector
%     .WML     - [matrix] contains the first n reweighted eigenvectors of the emperical covariance matrix
%     .sigmaML - [float] the variance contained in the remaining eigenvalues not contained in WML
%     .pTransV - [cell] contains a sparse transition matrix for each image column and boundary pair for the initialization step
%
% Calls: extractShapeModel
% See also: trainShape

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 04-Feb-2016

% generates a default filename out of the set of the trainings set filenames
h = hash([files.name],'MD5');
options = setShapeDefaults(options,params,files);

% loads the previously calculated shape prior
if options.loadShape
	tic;
	if isfield(options,'loadShapeName') h = options.loadShapeName; end

	if exist(sprintf('%s/%s.mat',options.shapeFolder,h),'file')
		load(sprintf('%s/%s.mat',options.shapeFolder,h));
	else
		error(sprintf('File %s does not exist, specify a valid filename (via options.loadShapeName).',h));
	end

	% compare number of B-Scans in the shape model with the number set in the options
	if size(models.columnsShape,1) ~= collector.options.numRegionsPerVolume
		error(sprintf('Model containts %d regions, but %d are set via options',size(models.columnsShape,1),collector.options.numRegionsPerVolume));
	end

	% in case the transition matrices are not stored in the model to save storage
	if ~isfield(models,'pTransV')
		models = preCalcTransitionMatrices(collector,models);
	end

	printMessage(sprintf('... loaded saved shape prior model in %.2f s ... \n',toc),1,collector.options.verbose);
else
	printMessage(sprintf('... train shape prior model (may take several minutes) ... \n'),1,collector.options.verbose); tic;
	% fetch the ground truth for files
	[mu,Sigma,data] = extractShapeModel(files,collector);
	
	for j = 1:length(mu)
		% perform PPCA (see the book 'Pattern Recognition' by Bishop for details)
		[V D] = eig(Sigma{j}); % eigendecomposition of the covariance matrix 
		D = diag(D);
		% estimate for sigma^2 is the average of the disregarded eigenvalues
		sigmaML = 1/(length(D)-options.numModes)*sum(D(1:end-options.numModes));
		% select modes accoding to the options.numModes highest eigenvalues
		WML = V(:,end:-1:end-options.numModes+1)*sqrt(diag(D(end:-1:end-options.numModes+1)) - eye(options.numModes)*sigmaML);
		% add parameters to the model struct
		models(j).mu = mu{j}'; models(j).WML = WML; models(j).sigmaML = sigmaML;
	    if collector.options.full3D
			models.columnsShape = collector.options.columnsShape;
		else
			models(j).columnsShape = collector.options.columnsShape{j};
		end

		if options.returnShapeData
			models(j).data = data;
		end
	end
	clear D V Sigma

	% pre-calculate the shape prior for the columnwise graphical models to speed up prediction
%	models = preCalcTransitionMatrices(collector,models); 

	% save the trained shape model as matfile
	if options.saveShape
		h = options.saveShapeName;
		% save all model components
		if ~options.saveCompressedModel
			save(sprintf('%s/%s.mat',options.shapeFolder,h),'models','-v7.3');
		% reduce the required storage by not saving the transition matrices but calculate them upon loading the model
		else
			for j = 1:length(models)
				pTransV{j} = models(j).pTransV;
			end
			models = rmfield(models,'pTransV');
			save(sprintf('%s/%s.mat',options.shapeFolder,h),'models','-v7.3');
			for j = 1:length(models)
				models(j).pTransV = pTransV{j};
			end
		end
		printMessage(fprintf('Stored shape model for later usage in %s\n',sprintf('%s/%s.mat',options.shapeFolder,h)),1,collector.options.verbose);
	end

	printMessage(sprintf('... trained shape prior model in %.2f s ... \n',toc),1,collector.options.verbose);
end

end
