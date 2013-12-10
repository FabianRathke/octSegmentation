function options = setDefaultShape(options,params,files)
% setDefaultShape - sets the default values for variables that are used to train the shape prior
%
% Syntax:
%   options = setDefaultShape(options,params)
%
% Inputs:
%   options   - [struct] options struct
%     .triu             - [int] number of diagonals from the transition matrices that are set to zero; parameter for the Matlab function triu
%								this enforces the ordering constraint when doing inference for q_c; boundaries will be at least optinos.triu pixels away from each other
%     .shapePriorFolder - [string] points to the folder where shape models are loaded from and saved to; default is $HOME/shapemodels [not tested for Windows]
%     .loadShape        - [boolean] if true, shape model is loaded from shapePriorFolder and not trained
%     .saveShape        - [boolean] if true, shape model is saved in shapePriorFolder
%     .filenameShape    - [string] filename of the shape loaded or saved
%   params    - [struct] holds params used for example in cross-validation 
%
% Outputs:
%   options - [struct] options struct augmented by the new field
%
% See also: trainShape, setCollectorDefault, setAppearanceDefault

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 06-Dec-2013

options = checkFields(options,params,20,'numModes');

if ~isfield(options,'triu')
	options.triu = 1;
end

if ~isfield(options,'shapePriorFolder')
	options.shapePriorFolder = [getenv('HOME'),'/shapemodels'];
	% check if folder exists
	if ~exist(options.shapePriorFolder)
		mkdir(options.shapePriorFolder);
	end
end

if ~isfield(options,'loadShape')
	options.loadShape = 0;
end

if ~isfield(options,'saveShape')
	options.saveShape = 0;
end

if options.loadShape || options.saveShape
	if ~isfield(options,'filenameShape')
		options.filenameShape = hash([files.name],'MD5');
	end
end

