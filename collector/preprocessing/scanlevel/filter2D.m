function B0 = filter2D(B0,options)
% filter2D - Preprocessing Module - filters the image B0 with the convolution matrix defined in options{2}
% 
% Syntax:
%   B0 = filter2D(B0,options)
%
% Inputs:
%   B0		- [matrix] OCT gray scan
%   options - [cell array]
%
% Outputs:
%   B0		- [matrix] OCT scan after filtering
%
% See also: trainAppearance, loadData

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 25-Feb-2014

B0 = conv2(B0,options{2},'same');
