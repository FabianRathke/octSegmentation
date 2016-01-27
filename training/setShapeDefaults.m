function options = setShapeDefaults(options,params,files)
% setShapeDefaults - sets the default values for variables that are used to train the shape prior
%
% Syntax:
%   options = setShapeDefaults(options,params)
%
% Inputs:
%   options   - [struct] options struct
%     .numModes            - [int] number of modes used in PPCA (default: 20)
%     .shapeFolder         - [string] points to the folder where shape models are loaded from and saved to. (default: '$HOME/shapemodels')
%     .loadShape           - [boolean] if true, shape model is loaded from shapePriorFolder and not trained (default: false)
%     .saveShape           - [boolean] if true, trained shape model is saved in shapePriorFolder (default: false)
%     .loadShapeName       - [string] filename of the shape model loaded (default: hash from concatenated filenames)
%     .saveShapeName       - [string] filename of the shape model saved (default: hash from concatenated filenames)
%     .saveCompressedModel - [boolean] if true, only the less storage intensive parts of the model are stored; the other parts are recalculated upon loading (default: true)
%	  .returnShapeData     - [boolean] if true, adds the segmentation training data to the model struct (default: false)
%   params    - [struct] holds params used for example in cross-validation 
%     .numModes   - same as above
%
% Outputs:
%   options - [struct] options struct augmented by default values for unset field
%
% See also: trainShape, setCollectorDefaults, setAppearanceDefaults

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 30-Apr-2013

options = checkFields(options,params,20,'numModes');

if ~isfield(options,'shapeFolder')
	options.shapeFolder = [getuserdir,'/shapemodels'];
	% check if folder exists
	if ~exist(options.shapeFolder)
		mkdir(options.shapeFolder);
	end
end

if ~isfield(options,'loadShape')
	options.loadShape = 0;
end

if ~isfield(options,'saveShape')
	options.saveShape = 0;
end

if options.loadShape
	if ~isfield(options,'loadShapeName')
		options.loadShapeName = hash([files.name],'MD5');
	end
end

if options.saveShape
	if ~isfield(options,'saveShapeName')
		options.saveShapeName = hash([files.name],'MD5');
	end
end

if ~isfield(options,'saveCompressedModel')
	options.saveCompressedModel = true;
end

if ~isfield(options,'returnShapeData')
	options.returnShapeData = 0;
end
