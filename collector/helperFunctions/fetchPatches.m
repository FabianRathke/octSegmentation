function patches = fetchPatches(filename,idxSet,options)
% fetchPatches - fetches patches of size options.width x options.height from the B-Scan in filename
% 
% Inputs:
%   filename - [string] filename of the matfile to load
%   options  - [struct]
%      .centerPatches - [boolean] substract the mean of each patch
%      .height        - [int] patch height in px (has to odd number)
%      .width         - [int] patch width in px
%
% Outputs:
%   patches - [struct]
%      .data - [matrix] contains the patches
%
% See also: collectTrnData, collectTestData
% Calls: loadLabels, loadData, wrapperGetRaster, getAugmentedImage

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 05-Dec-2013

%tic;
pw = ([options.height options.width]-1)/2;
patchSize = options.width*options.height;
B0 = loadData(filename,options);

[Y X] = size(B0);

% increase the scan to draw patches from the border
B0Aug = getAugmentedImage(B0,[options.height options.width]);

% fetch indices to draw patches at positions idxSet
A = repmat([1:options.height]',[options.width 1])';
A = A(ones(1,size(idxSet,2)),:);
a = idxSet(1,:)'-pw(1)-1;
a = a(:,ones(1,size(A,2)));
A = A + a;

B = reshape(repmat(1:options.width,[options.height 1]),1,[]);
B = B(ones(1,size(idxSet,2)),:);
b = idxSet(2,:)'-pw(2)-1;
b = b(:,ones(1,size(B,2)));
B = B + b;

indices = sub2ind(size(B0Aug),A,B);
%fprintf('Overhead when fetching patches in %.3f\n',toc);

patches.data = B0Aug(indices);

% center all patches
if options.centerPatches
	mean_tmp = mean(patches.data')';
	patches.data = patches.data - mean_tmp(:,ones(1,patchSize));
end

if options.saveAppearanceTerms && options.loadLabels
	interpolation = loadLabels(filename,options);
	Classes = wrapperGetRaster(interpolation(options.EdgesTrain,:),[options.Y size(interpolation,2)]);
%	patches.classID = Classes(sub2ind(size(B0),idxSet));
	patches.classID = Classes;
end
