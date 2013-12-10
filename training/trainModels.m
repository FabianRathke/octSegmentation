function models = trainModels(files,collector,params,options)
% trainModels - calls the functions that train appearance and shape model, wrapped into a single function to be used with the cross-validation script
%
% Syntax:
%   models = trainModels(files,collector,params,options)
% 
% Inputs:
%   files     - [struct] list of files used for trainin
%	collector - [struct] options for the collector
%   params    - [struct]
%   options   - [struct] options for the functions training appearance and shape models
%     .appearanceModel - [struct] for details --> trainAppearance 
%     .shapeModel      - [struct] for details --> trainShape
%
% Outputs:
%   models - [struct]
%    .appearanceModel - [struct] for details --> trainAppearance
%    .shapeModel      - [struct] for details --> trainShape
%
% See also: trainAppearance, trainShape
% Calls: trainAppearance, trainShape

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 06-Dec-2013

models.appearanceModel = trainAppearance(files,collector,params,options.appearanceModel);
models.shapeModel = trainShape(files,collector,params,options.shapeModel);

