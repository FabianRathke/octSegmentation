function options = setSegmentationDefaults(options,params)
% setSegmentationDefaults - parameters that are used during the inference
%  
% Syntax:
%   options = setSegmentationDefaults(options,params)
%
% Inputs:
%   options - [struct] options struct
%     .iterations          - [int] max number of iterations of alternating optimizations between q_c and q_b (default: 30)
%     .plotting            - [boolean] trigger to plot segmentation results for each iteration (default: 0)
%     .threshold           - [float] determines the relative ...
%     .threshold_q_c       - [float] 
%     .thresholdAccuracy   - [float] 
%     .alpha               - [float] smaller values increase the area around the mean of q_b that is taken into consideration when inferring q_c (default: 0.1) 
%   params  - [struct] holds params used for example in cross-validation 
%
% Outputs:
%   options - [struct] options struct with default values set
%
% See also: predVariational

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 10-Dec-2013

if (~isfield(options,'iterations'))
    options.iterations = 30;
end

if ~isfield(options,'plotting')
    options.plotting = 0;
end

if ~isfield(options,'threshold_q_c')
	options.threshold_q_c = 10^-6;
end

options = checkFields(options,params,0.1,'alpha');

if ~isfield(options,'threshold')
	options.threshold = 0.01;
end

if ~isfield(options,'calcFuncVal')
	options.calcFuncVal = 0;
end

if ~isfield(options,'thresholdAccuracy')
	options.thresholdAccuracy = realmin('single');
end
