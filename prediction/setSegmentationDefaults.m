function options = setSegmentationDefaults(options,params)
% setSegmentationDefaults - parameters that are used during the variational inference part in predVariational.m
%  
% Syntax:
%   options = setSegmentationDefaults(options,params)
%
% Inputs:
%   options - [struct] options struct
%     .iterations      - [int] max number of iterations of alternating optimizations between q_c and q_b (default: 30)
%     .plotting        - [boolean] plot segmentation results for each iteration (default: false)
%     .folderPlots     - [string] the folder where the plots are stored. (default: '$HOME/plots')
%     .threshold       - [float] the level of relative change which determines convergence (default: 0.01)
%     .threshold_q_c   - [float] the treshold, below which probabilities in q_c are set to zero, used to speed up calculations (default: 10^-6) 
%     .alpha           - [float] smaller values increase the area around the mean of q_b that is taken into consideration when inferring q_c (default: 0.1)
%     .calcFuncValue   - [boolean] if true outputs terms, that can be used for pathology detection and to access the segmentation uncertainty (see the corresponding publication for detailts) (default: false)
%     .detailedOutput  - [boolean] triggers the ouput of more detailed results (the q_b mode as well as errors for q_b and q_c for each iteration) (default: false)
%     .onlyInitialize  - [boolean] just executes the initialization of the q_c distribution and returns afterwards
%   params  - [struct] holds params used for example in cross-validation 
%     .alpha           - [float] same as above
%
% Outputs:
%   options - [struct] options struct with default values set
%
% See also: predVariational

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 29-Apr-2014

if (~isfield(options,'iterations'))
    options.iterations = 100;
end

if ~isfield(options,'plotting')
    options.plotting = 0;
end

if options.plotting == 1
	if ~isfield(options,'folderPlots')
    	options.folderPlots = [getuserdir,'/plots'];
	end
	if ~exist(options.folderPlots)
    	mkdir(options.folderPlots);
	end
end

if ~isfield(options,'threshold_q_c')
	options.threshold_q_c = 10^-6;
end

options = checkFields(options,params,0.1,'alpha');

if ~isfield(options,'threshold')
	options.threshold = 0.05;
end

if ~isfield(options,'calcFuncVal')
	options.calcFuncVal = 0;
end

%if ~isfield(options,'thresholdAccuracy')
	options.thresholdAccuracy = realmin('single');
%end

if ~isfield(options,'detailedOutput')
	options.detailedOutput = 0;
end

if ~isfield(options,'printTimings')
	options.printTimings = 0;
end

if ~isfield(options,'onlyInitialize')
	options.onlyInitialize = 0;
end
