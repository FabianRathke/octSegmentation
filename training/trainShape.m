function models = trainShape(files,collector,params,options)
% trainShape - collects ground truth and estimates a Gaussian distribution regularized by PPCA; alternatively loads a shape stored on disk
% 
% Syntax:
%   models = trainShape(files,collector,params,options)
%
% Inputs:
%   files     - [struct] output of the matlab function 'dir'
%   collector - [struct] holds the options for the collector grabbing the training patches, description given in setCollectorDefaults
%     .columnsShape
%     .
%   params    - [struct] holds params used for example in cross-validation 
%   options   - [struct] holds options used within the appearance function
%     .numModes         - [int] number of eigenvectors to be used
%     .loadShape        - [boolean]
%     .loadShapeName    - [string]
%     .saveShape        - [boolean]
%     .saveShapeName    - [string]
%     .shapePriorFolder - [string] folder where shape models are stored
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
% Last Revision: 10-Dec-2013

% generates a default filename out of the set of the trainings set filenames
h = hash([files.name],'MD5');
options = setDefaultShape(options,params,files);

if options.loadShape
	tic;
	if isfield(options,'loadShapeName') h = options.loadShapeName; end

	if exist(sprintf('%s/%s.mat',options.shapePriorFolder,h),'file')
		load(sprintf('%s/%s.mat',options.shapePriorFolder,h));
	else
		error(sprintf('File %s does not exist, specify a valid filename (via options.loadShapeName).',h));
	end

	% compare number of B-Scans in the shape model with the number set in the options
	if size(models.columnsShape,1) ~= collector.options.numRegionsPerVolume
		error(sprintf('Model containts %d regions, but %d are set via options',size(models.columnsShape,1),collector.options.numRegionsPerVolume));
	end

	printMessage(sprintf('... loaded saved shape prior model in %.2f s ... \n',toc),1,collector.options.verbose);
else
	% fetch the ground truth for files
	[tmp mu Sigma] = extractShapeModel(files,collector);
	printMessage(sprintf('... train shape prior model (can take several minutes) ... \n'),1,collector.options.verbose);
	tic;
	% ******* PPCA *********, see Tipping et al, 1999
	[V D] = eig(Sigma); % eigendecomposition of the covariance matrix 
	D = diag(D);
	% estimate for sigma^2 is the average of the disregarded eigenvalues
	sigmaML= 1/(length(D)-options.numModes)*sum(D(1:end-options.numModes));
	% select modes accoding to the options.numModes highest eigenvalues
	WML = V(:,end:-1:end-options.numModes+1)*sqrt(diag(D(end:-1:end-options.numModes+1)) - eye(options.numModes)*sigmaML);
	% add parameters to the model struct
	models.mu = mu'; models.WML = WML; models.sigmaML = sigmaML; models.columnsShape = collector.options.columnsShape;
	clear D V Sigma
	% ********** end PPCA **********

	numBounds = length(collector.options.EdgesTrain); 
	models.pTransV = cell(1,collector.options.numRegionsPerVolume);
	numColumnsShape = cellfun('length',collector.options.columnsShape);
	% pre-calculate the shape prior for the columnwise graphical models to speed up prediction
	for region = 1:collector.options.numRegionsPerVolume
		numColumns = length(collector.options.columnsShape{region});
		pTransTmp = cell(numColumns,numBounds);
		for i = 1:numColumns
			for j = 2:numBounds
				idx_a = i + (j-2)*numColumns + numBounds*sum(numColumnsShape(1:region-1));
				idx_b = idx_a + numColumns;
				P = inv(models.WML([idx_a idx_b],:)*models.WML([idx_a idx_b],:)' + eye(2)*models.sigmaML);
				pTransTmp{i,j} = sparse(getCondTransMatrix([models.mu(idx_a) models.mu(idx_b)]',P,collector.options.Y,options));
			end
		end
		models.pTransV{region} = pTransTmp;
	end
	if options.saveShape
		if isfield(options,'saveShapeName') h = options.saveShapeName; end
		save(sprintf('%s/%s.mat',options.shapePriorFolder,h),'models','-v7.3')
	end

	printMessage(sprintf('... trained shape prior model in %.2f s ... \n',toc),1,collector.options.verbose);
end
