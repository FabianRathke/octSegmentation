function options = setCollectorDefaults(options,params,files,folderData,folderLabels)
% setCollectorDefaults - sets defauls values for variables used for collecting training and/or testdata and some global variables used for trainin and prediction
%  
% Syntax:	
%   options = setCollectorDefaults(options,params,files,folderData,folderLabels)
%
% Inputs:
%   options      - [struct] options struct
%       .clip                - [boolean] determines whether to clip parts of the scan before segmentation; useful for example to remove the nerve head. Default: [false]
%       .clipRange           - [array](2) left and right boundary for clipping. Default: [1 scan-width]
%		.clipFactor			 - [int] the interval in y-direction (i.e. 2 means that every second column is grabbed)
%       .loadRoutineData     - [string] user-defined routine to load a B-Scan; see loadData.m for details. Default: ['spectralisMat'] 
%       .loadRoutineLabels   - [string] user-defined routine to load ground truth; see loadLabels.m for details. Default: ['LabelsFromLabelingTool']
%       .labelIDs            - [array](numFiles,numRegionsPerVolume) holds ids to identify the scans used in a 3-D volume 
%       .numRegionsPerVolume - [int] for 3-D volumes, the number of B-Scans; for 2-D set = 1
%       .calcOnGPU           - [boolean] move parts of the calculation onto the GPU. Default: [false]
%       .columnsShape        - [cell-array](1,numRegionsPerVolume) holds an array for each region; indicates which columns of the B-Scan are used for the shape-prior (intermediate columns are interpolated). Default: columnsShape{i} = [1:2:scan-width]
%       .columnsPred         - [array](2) determines for which columns predictions are made. Will be applied to each scan in the volume. Default: [1:2:scan-width]
%       .BScanPositions      - [double] position in SLO scan in px distance from central region (to position shape prior and select correct subregion)
%       .saveAppearanceTerms - [boolean] return appearance model predictions for each pixel. Default: [0]
%       .printTimings        - [boolean] print cpu/gpu timings of the different modules. Default: [false]
%       .verbose             - [int] the amount of printed information while runnning the programm (0 (nothing) - 2 (maximal)). Default: [1]
%   params       - [struct] holds params used for example in cross-validation 
%   files        - [struct] training set, output of Matlabs dir function; if empty some defaults will not be set
%   folderData   - [string] path with mat-files containing the OCT scans 
%   folderLabels - [string] path with mat-files containing ground truth
%
% Outputs:
%   options - [struct] options struct with default values set
%
% See also: collectTestData, collectTrnData, loadLabels, loadData, fetchPatches

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Major Revision: 15-Nov-2016

options.folder_data = folderData;
options.folder_labels = folderLabels;

if ~isfield(options,'loadLabels') options.loadLabels = 1; end
if ~isfield(options,'mirrorBScan') options.mirrorBScan = ''; end

% the amount of information output
if ~isfield(options,'verbose') options.verbose = 1; end

% define loading routines for labels and data sets
if ~isfield(options,'loadRoutineData') options.loadRoutineData = 'spectralisMat'; end
if ~isfield(options,'loadRoutineLabels') options.loadRoutineLabels = 'LabelsFromLabelingTool'; end

% the number of regions per file/volume, 1 ==  2-D Scan, > 1 == 3-D Volume
if ~isfield(options,'numRegionsPerVolume') 
	if ~isfield(options,'labelIDs')
		options.numRegionsPerVolume = 1;
	else
		options.numRegionsPerVolume = size(options.labelIDs,2);
	end
end
% ids for single scans of a volume  (used in functions loadLabels and loadData)
if ~isfield(options,'labelIDs') options.labelIDs = ones(length(files),options.numRegionsPerVolume); end

% enables to clip the width of B-Scans; i.e. useful for volumnes to remove the part containing the nerve head
if ~isfield(options,'clip') 
	options.clip = 0; 
else
	if options.clip
		if ~isfield(options,'clipRange') 
			error('Please specify the clip range in options.clipRange');
		end
	end
end
if ~isfield(options,'clipFactor')
	options.clipFactor = 1;
end

if length(files) > 0
	% pull a sample scan and set its dimensions
	options.labelID = options.labelIDs(1);
	options.preprocessing = struct();
	B0 = loadData(files(1).name,options);
	if isfield(options,'clipRange')
		options.X = options.clipRange(2) - options.clipRange(1) + 1;
	else
		options.X = size(B0,2);
	end
	options.Y = size(B0,1);

	% pull sample ground truth
	segmentation = loadLabels(files(1).name,options);
	options = rmfield(options,'labelID');

	numBoundaries = size(segmentation,1); numLayers = numBoundaries + 1;
	% edges and layers used for training and prediction
	options.EdgesTrain = 1:numBoundaries;
	options.LayersTrain = 1:numLayers;
	options.numLayers = length(options.LayersTrain);
	options.EdgesPred = 1:numBoundaries;
	options.LayersPred = 1:numLayers;

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
end

% default behavior is to perform all calculations on the CPU
if ~isfield(options,'calcOnGPU') options.calcOnGPU = 0; end
if options.calcOnGPU
	GPUstart;
end

% the datatype for those variables that are moved onto the GPU has to be variable
if options.calcOnGPU
	options.dataType = 'GPUsingle';
	options.dataTypeCast = 'GPUsingle';
else
	options.dataType = '''double''';
	options.dataTypeCast = 'double';
end

% print Timings during prediction
if ~isfield(options,'printTimings') options.printTimings = 0; end

% whether to return the training data for the shape model, i.e. all training segmentations
if ~isfield(options,'returnShapeData') options.returnShapeData = 0; end

if ~isfield(options,'full3D') options.full3D = 0; end

% params of shape prior for volumes 
if options.numRegionsPerVolume > 1
	if ~isfield(options,'BScanPositions') options.BScanPositions = 0; end
	if ~isfield(options,'BscansSelect') options.BscansSelect = 0; end
end

% the results struct will also contain the pixel wise probabilities for each appearance model
if ~isfield(options,'saveAppearanceTerms') options.saveAppearanceTerms = 0; end

% assign the appearance terms to a global variable --> can be reused in subsequent runs
if ~isfield(options,'makePredGlobal') options.makePredGlobal = 0; end

