function options = checkFields(options,params,value,fieldname)
% checkFields - checks whether fieldname exists in options, if not checks if it exists in params and sets this value; elseif value is used
%
% Syntax:
%   options = checkFields(options,params,value,fieldname)
%
% Inputs:
%   options   - [struct] options struct
%   params    - [struct] holds params used for example in cross-validation 
%   value     - [numeric/string] value that options.fieldname is set to, in case params.fieldname does not exists
%   fieldname - [string] name of the field to be set on options
%
% Outputs:
%   options - [struct] options struct augmented by the new field
%
% See also: setCollectorDefault

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 05-Dec-2013

if ~isfield(options,fieldname)
    % if the option is set as a parameter used for example in cross-validation, set it here
	if isfield(params,fieldname)
        eval(sprintf('options.%s = params.%s;',fieldname,fieldname));
	% otherwise use the standard value
    else
		eval(sprintf('options.%s = value;',fieldname));
    end
end

