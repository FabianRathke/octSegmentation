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
% Last Major Revision: 03-Feb-2016

if isfield(options,'appearanceModel')
	if isstruct(options.appearanceModel)
		appearanceModel{1} = options.appearanceModel;
		options = rmfield(options,'appearanceModel');
		options.appearanceModel = appearanceModel;
	end
	if length(options.appearanceModel) > 1
		fprintf('Learn mixture of appearance models\n');
	end
	for j = 1:length(options.appearanceModel)
		[models.appearanceModel{j} options.appearanceModel{j}] = trainAppearance(files,collector,params,options.appearanceModel{j});
	end
	models.options.appearanceModel = options.appearanceModel;
end
if isfield(options,'shapeModel')
	models.shapeModel = trainShape(files,collector,params,options.shapeModel);
end

models.params = collector.options;
