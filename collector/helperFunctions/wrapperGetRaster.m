function Raster = wrapperGetRaster(border,dim)
% loadData - calls getRaster.c, takes a matrix of boundary positions and converts this into a pixel-wise labeling for the B-scan
%
% Inputs:
%   border - [matrix] holds segmentation of B-scan
%   dim    - [array](2) dimension of the B-scan
%
% Outputs:
%   Raster - [matrix] a pixel-wise labeling of the B-scan 
%
% See also: getRaster.c, fetchPatches, collectTrnData, collectTestData

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 05-Dec-2013

[a idx] = sort(border(:,1));
border = round(border(idx,:) - 0.5);
borders = [zeros(1,dim(2)); border; ones(1,dim(2))*dim(1)];

Raster = getRaster(double(borders),dim);

