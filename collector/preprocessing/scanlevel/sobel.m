function B0 = filter2D(B0,options)
% sobel - Preprocessing Module - applies the sobel filter --> edge detection
% 
% Syntax:
%   B0 = sobel(B0,options)
%
% Inputs:
%   B0		- [matrix] OCT gray scan
%   options - [cell array] - not used
%
% Outputs:
%   B0		- [matrix] OCT scan after filtering
%
% See also: trainAppearance, loadData

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 20-Mar-2016

% example code from http://de.mathworks.com/help/coder/examples/edge-detection-on-images.html
k = [1 2 1; 0 0 0; -1 -2 -1];
H = conv2(B0,k, 'same');
V = conv2(B0,k','same');
B0 = sqrt(H.*H + V.*V);
