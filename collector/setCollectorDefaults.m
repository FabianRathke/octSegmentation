function options = setCollectorDefaults(options,params,files)
% setCollectorDefaults - sets defauls values for variables used in the collectors; variables from params, used in a CV, can overwrite the standard values e.g. for the patch size
%  
% Details for each variable and its default value can be found within the source code
%
% Syntax:	
%   options = setCollectorDefaults(options,params,files)
%     .testFullImage
%     .width
%     .height
%
% Inputs:
%   options - [struct] options struct
%   params  - [struct] holds params used for example in cross-validation 
%   files   - [struct] output of Matlabs dir function
%
% Outputs:
%   options - [struct] options struct with default values set
%
% See also: collectTestData, collectTrnData, loadLabels, loadData, fetchPatches

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 06-Dec-2013

% patch width and height
options = checkFields(options,params,15,'width');
options = checkFields(options,params,15,'height');

% edges and layers used for training and prediction
if ~isfield(options,'EdgesTrain') options.EdgesTrain = 1:9; end
if ~isfield(options,'LayersTrain') options.LayersTrain = 1:10; end

if ~isfield(options,'EdgesPred') options.EdgesPred = 1:9; end
if ~isfield(options,'LayersPred') options.LayersPred = 1:10; end

% shall patches be projected onto a low-dimensional manifold determined by PCA?
if ~isfield(options,'projToEigenspace') options.projToEigenspace = true; end
% dimension of the low-dimensional manifold
if ~isfield(options,'numModesAppearance') options.numModesAppearance = 20; end

% fetches patches for the complete test scan; if false grabs only a subset 
if ~isfield(options,'testFullImage') options.testFullImage = true; end

% width and height of a B-Scan
if ~isfield(options,'X')
	if isfield(options,'clipRange')
		options.X = options.clipRange(2) - options.clipRange(1) + 1;
	else
		options.X = 768;
	end
end
if ~isfield(options,'Y') options.Y = 496; end

options.numFiles = length(files);
% the number of regions per file/volume, 1 ==  2-D Scan, > 1 == 3-D Volume
if ~isfield(options,'numRegionsPerVolume') options.numRegionsPerVolume = 1; end
% indicates how single scans of a volume are to be identified (are then used in functions loadLabels and loadData)
if ~isfield(options,'labelIDs') options.labelIDs = zeros(options.numFiles,options.numRegionsPerVolume); end

% can be used to divide each B-Scan into regions, which use seperate appearance models
if ~isfield(options,'BScanRegions') options.BScanRegions = [1 options.X]; end
if ~isfield(options,'numRegionsPerBScan') options.numRegionsPerBScan = size(options.BScanRegions,1); end


% add noise to patches (in order to simulate noisier B-Scans)
if ~isfield(options,'addNoise') options.addNoise = 0; end

% default behavior is to perform all calculations on the CPU
if ~isfield(options,'calcOnGPU') options.calcOnGPU = 0; end
% the datatype for those variables that are moved onto the GPU has to be variable
if options.calcOnGPU
	options.dataType = 'GPUSingle';
else
	options.dataType = '''double''';
end

% most probably deprecated (05-12-2013)
%if ~isfield(options,'columns')
%	options.columns = round(linspace(1,options.X,options.X));
%end

% number of patches to draw for training; per class and file
if ~isfield(options,'numPatches') options.numPatches = 30; end
% substract the patch mean from each patch; make appearance terms less vulnerable to variations of intensity between and within OCT scans
if ~isfield(options,'centerPatches') options.centerPatches = 1; end

% enables to clip the width of B-Scans; i.e. useful for volumnes to remove the part containing the nerve head
if ~isfield(options,'clip') options.clip = 0; else
	if ~isfield(options,'clipRange') options.clipRange = [1 options.X]; end
end

% which columns are part of the shape prior p(b) and are used for q_b (allows for sparse representations; intermediate columns are interpolated);
% can be set for each region within the volume separately
if ~isfield(options,'columnsShape') 
	options.columnsShape = cell(1,options.numRegionsPerVolume);
	for i = 1:options.numRegionsPerVolume
		options.columnsShape{i} = round(linspace(1,options.X,options.X/2));
	end
end
% which columns are to be predicted, i.e. are used in q_c
if ~isfield(options,'columnsPred') options.columnsPred = round(linspace(1,options.X,options.X/2)); end


% patches are drawn from random positions within their layer; other option is 'middle', patches are drawn from the center of their respective layer
if ~isfield(options,'patchPosition') options.patchPosition = 'random'; end

% if one segments volumes, this has to be set to 1
if ~isfield(options,'ThreeD') options.ThreeD = 0; end

% exist labels for the file, if true the collector collects ground trouth for later evaluation
if ~isfield(options,'storeLabels') options.storeLabels = 1; end

% the amount of information output
if ~isfield(options,'verbose') options.verbose = 0; end

% print Timings during prediction
if ~isfield(options,'printTimings') options.printTimings = 0; end
options.numLayers = length(options.LayersTrain);

if ~isfield(options,'saveAppearanceTerms') options.saveAppearanceTerms = 0; end
