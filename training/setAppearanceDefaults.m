function options = setAppearanceDefaults(options,params,files)
% setAppearanceDefaults - sets the default values for variables that are used to train the appearance models
%
% Syntax:
%   options = setAppearanceDefaults(options,params)
%
% Inputs:
%   options   - [struct] options struct
%     .appearanceModel - [function-handle] points to the method used for training the appearance models
%   params    - [struct] holds params used for example in cross-validation 
%
% Outputs:
%   options - [struct] options struct augmented by default values for unset fields
%
% See also: trainShape, setCollectorDefaults, setShapeDefaults

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 29-Apr-2014

if ~isfield(options,'appearanceModel')
	options.appearanceModel = @trainGlasso;
end
