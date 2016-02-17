function B0_aug = getAugmentedImage(B0,p)
% getAugmentedImage  - augments image to allow selection of patches from the scan-boundary
%
% Syntax:
%   getAugmentedImage(B0,p)
%
% Inputs:
%   B0 	   - [matrix] the B-Scan
%   p      - [array](2) patch-size
%
% Outputs:
%   B0_aug - [matrix] B-Scan augmented at the boundaries 
%
% See also: fetchPatches, collectTrnData, collectTestData

% TODO: At the moment the B-scan is trivially extended on its left and right border; for circular scans a more sophisticaed approach would be possible (i.e. continue at the opposite sides);
% 	for volumes, where clipping is used, the augmented scan could obtain the boundary information from the cutted parts; the trivial way of augmenting the B-Scan increases the error on the boundaries measurably; 
%	but since this affects only a small area the total error increases by only <= 1%
%
% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 05-Dec-2013

[y x] = size(B0);
pw = (p-1)/2;
B_middle =  [B0(:,1:pw(2)) B0 B0(:,x-pw(2)+1:end)]; 

% B-scan augmented by neighboring borders
B0_aug = [B0(1:pw(1),1:pw(2)) B0(1:pw(1),:) B0(1:pw(1),x-pw(2)+1:end); B_middle; B0(y-pw(1)+1:end,1:pw(2)) B0(y-pw(1)+1:end,:) B0(y-pw(1)+1:end,x-pw(2)+1:end)]; 
