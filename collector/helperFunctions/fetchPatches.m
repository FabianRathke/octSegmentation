function patches = fetchPatches(filename,idxSet,options)
% fetchPatches - fetches patches of size options.width x options.height from the B-Scan in filename
% 
% Inputs:
%   filename - [string] filename of the matfile to load
%   options  - [struct]
%      .centerPatches - [boolean] substract the mean of each patch
%      .height        - [int] patch height in px (has to odd number)
%      .width         - [int] patch width in px
%      .Patches3D     - [boolean] whether to fetch 3D-patches 
%      .depth         - [int] patch depth in px (for 3D only)
%
% Outputs:
%   patches - [struct]
%      .data - [matrix] contains the patches
%	   .positions - [matrix] containts the row/column indices
%
% See also: collectTrnData, collectTestData
% Calls: loadLabels, loadData, wrapperGetRaster, getAugmentedImage

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 28-Mar-2017

checkPredictionGlobal(filename);
global predictionGlobal;
global isStoredInGlobal;
%tic;
pw = ([options.height options.width]-1)/2;
patchSize = options.width*options.height*options.depth;
if options.Patches3D
	labelIDs = options.labelID + [-(options.depth-1)/2:(options.depth-1)/2];
	for j = 1:options.depth
		options.labelID = labelIDs(j);
		B0{j} = loadData(filename,options);
	end
else
	if isStoredInGlobal.data & isfield(predictionGlobal,'BScans') & ~isempty(predictionGlobal.BScans{options.labelID})
		B0{1} = predictionGlobal.BScans{options.labelID};
	else
		B0{1} = loadData(filename,options);
	end
end

[Y X Z] = size(B0{1});

% increase the scan to draw patches from the border
for j = 1:length(B0)
	for i = 1:Z
		B0Aug{j}(:,:,i) = getAugmentedImage(squeeze(B0{j}(:,:,i)),[options.height options.width]);
	end
end

% fetch indices to draw patches at positions idxSet
A = int32(repmat([1:options.height]',[options.width 1])'-pw(1)-1);
%A = A(ones(1,size(idxSet,2)),:);
%A = A + idxSet(ones(1,size(A,2)),:)';
A = bsxfun(@plus,A,idxSet(1,:)');

B = int32(reshape(repmat(1:options.width,[options.height 1]),1,[])-pw(2)-2);
%B = B(ones(1,size(idxSet,2)),:);
%B = B + idxSet(ones(1,size(B,2))*2,:)';
B = bsxfun(@plus,B,idxSet(2,:)');

% indices = sub2ind(size(B0Aug),A,B);
%indices = A + (B-1)*size(B0Aug{1},1);
indices = A + B*size(B0Aug{1},1);
%fprintf('Overhead when fetching patches in %.3f\n',toc);

if length(B0Aug) == 1 & Z==1
	patches.data = B0Aug{j}(indices);
else
	patches.data = [];
	for j=1:length(B0Aug)
		for i = 1:Z
			idxTmp = indices + (i-1)*X*Y;
			patches.data = [patches.data B0Aug{j}(idxTmp)];
		end
	end
end

% center all patches
if options.centerPatches
	patches.data = bsxfun(@minus,patches.data,mean(patches.data,2));
end

if options.saveAppearanceTerms && options.loadLabels
	interpolation = loadLabels(filename,options);
	Classes = wrapperGetRaster(interpolation(options.EdgesTrain,:),[options.Y size(interpolation,2)]);
%	patches.classID = Classes(sub2ind(size(B0),idxSet));
	patches.classID = Classes;
end
