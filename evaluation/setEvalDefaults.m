function options = setEvalDefaults(options)
% setEvalDefaults - parameters that are needed to compare the prediction against ground truth
%  
% Syntax:
%   options = setEvalDefaults(options,params)
%
% Inputs:
%   options - [struct] options struct
%     .interpolation  - [boolean] in case we only obtained predictions for a subset of image columns, interpolate the rest? if false, error measures will only be calculated for column subset (default: true)
%
% Outputs:
%   options - [struct] options struct with default values set
%
% See also: signedUnsigned

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 30-Apr-2014  

if ~isfield(options,'interpolation')
    options.interpolation = 1;
end

end
