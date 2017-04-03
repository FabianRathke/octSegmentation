function options = setAppearanceDefaults(options,params,files,collector)
% setAppearanceDefaults - sets the default values for variables that are used to train the appearance models
%
% Syntax:
%   options = setAppearanceDefaults(options,params)
%
% Inputs:
%   options - [struct] options struct
%       .appearanceModel     - [function-handle] points to the method used for training the appearance models
%       .width             	 - [int] patch width in px (has to be odd number). Default: [15]
%       .height            	 - [int] patch height in px. Default: [15]
%       .preprocessing       - [struct] struct with fields 'patchLevel' and 'scanLevel'; for each field cell array of methods and their parameters; Default: [options.preprocessing.patchLevel = {{@projToEigenspace,20}}]
%       .BScanRegions        - [array]([1,2],numBScanRegions) can be used to divide each B-Scan into regions with seperate appearance models; left and right boundaries of each reagion. Default: [[1 numColumns]], that is one appearance model for the whole scan
%       .numRegionsPerBScan  - [int] number of regions in each B-Scan, automatically determined from BScanRegions
%       .numPatches          - [int] number of patches to draw from each file for each appearance class. Default: [30]
%       .patchPosition       - [string] draw patches for layer-classes randomly ('random') or from the middle ('middle') for each column. Default: ['middle']
%       .centerPatches       - [boolean] subtract mean of each patch. Default: [true]
%       .priorVolumesPaper   - [boolean] assigns a specific non-uniform prior distribution to the appearance models; see documentation.pdf and the source code for details
%       .Patches3D			 - [boolean] uses voxes instead of 2-D patches
%		.depth				 - [int] voxel depth in z-direction
%   params  - [struct] holds params used for example in cross-validation 
%
% Outputs:
%   options - [struct] options struct augmented by default values for unset fields
%
% See also: trainShape, setCollectorDefaults, setShapeDefaults

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 15-Nov-2016

if ~isfield(options,'appearanceModel')
	options.appearanceModel = @trainGlasso;
end

% patch height and width
options = checkFields(options,params,15,'width');
options = checkFields(options,params,15,'height');

if ~isfield(options,'preprocessing') options.preprocessing = struct(); end
% preprocessing on patch-Level (performed in trainAppearance and predictAppearance)
if ~isfield(options.preprocessing,'patchLevel')
    options.preprocessing.patchLevel = {{@projToEigenspace,20}};
end

% preprocessing on scan-Level (performed on loadData) 
if ~isfield(options.preprocessing,'scanLevel') 
    options.preprocessing.scanLevel = {}; 
end 

% can be used to divide each B-Scan into regions, which use seperate appearance models
if isfield(options,'BScanRegions')
	% of only one column is provided (specifying left limits of each region), add another column with right limits
	if size(options.BScanRegions,2) == 1
		options.BScanRegions = [options.BScanRegions [options.BScanRegions(2:end)-1; collector.options.X]];
	end
end
if ~isfield(options,'BScanRegions') options.BScanRegions = [1 collector.options.X]; end
options.numRegionsPerBScan = size(options.BScanRegions,1);

% number of patches to draw for training; per class and file
if ~isfield(options,'numPatches') options.numPatches = 50; end
% substract the patch mean from each patch; make appearance terms less vulnerable to variations of intensity between and within OCT scans
if ~isfield(options,'centerPatches') options.centerPatches = 1; end
% patches are drawn from the center of their layer
if ~isfield(options,'patchPosition') options.patchPosition = 'middle'; end

if ~isfield(options,'priorVolumesPaper')
	options.priorVolumesPaper = 0;
end

if ~isfield(options,'Patches3D')
	options.Patches3D = 0;
end

if ~isfield(options,'depth') & options.Patches3D
	options.depth = 5;
end

if ~options.Patches3D
	options.depth = 1;
end
