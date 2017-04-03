function trnData = collectTrnData(files,options)
% trnData - draws patches from a set of training images for training the apperance model 
% 
% Syntax:
%   trnData = collectTrnData(files,options)
%
% Inputs:
%   files - [struct] list of files
%     .name	- [string] the filename
%   options - [struct] holds options for fetching patches
%     .BScanRegions   - [matrix](2xn) n subdivisions of the B-Scan; each row holds x and y coordinates; 
%     .LayersTrain    - [array] indices of layers to train
%     .EdgesTrain     - [array] indices of edges to train
%     .X              - [int] width of the BScan
%     .Y              - [int] height of the BScan
%     .centerPatches  - [boolean] points to the folder of filename
%     .height         - [int] patch height
%     .width          - [int] patch width
%     .numPatches     - [int] number of patches to be drawn per file and appearance class 
%
% Outputs:
%   trnData - [struct]
%     .data    - [matrix] the actual patches drawn for training
%     .classID - [array] classID for each patch; assigns IDs to different layers (1:numLayers) and boundaries (1:numBoundaries)
%     .type    - [array] 1 or 2 for each patch; encoding patches of layers (=1) and boundaries (=2)
%     .idx     - [idx] the position of each patch (1-d indexing)
%
% See also: collectTestData, setCollectorDefaults

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 28-Jan-2016

pw = ([options.height options.width]-1)/2; % number of pixels left and right of the center column in each patch

numClasses = length(options.LayersTrain) + length(options.EdgesTrain);
numPatches = options.numPatches*length(files)*numClasses;

for region = 1:size(options.BScanRegions,1)
%	trnData(region).data = zeros(numPatches,options.height*options.width*options.depth);
	trnData(region).data = [];
	trnData(region).classID = zeros(1,numPatches,'int8');
	trnData(region).type = zeros(1,numClasses,'int8');
	trnData(region).idx = zeros(numPatches,2,'int16');
end

helper = 1:options.Y;
numPatchesAdded = zeros(1,size(options.BScanRegions,1));

for i = 1:length(files)
	[a filename] = fileparts(files(i).name);
	options.labelID = options.labelIDsCurrent(i);
 	interpolation = loadLabels(filename,options);
	for region = 1:size(options.BScanRegions,1)
		% draw patches from within layers	
		idxRand = zeros(2,options.numPatches*numClasses); % x and y coordinates of patches to draw, coordinates indicate center of the patch
		setToZero = []; % in case we can not provide numPatches for that region we delete these samples to zero
		counter = 1;
		if length(options.LayersTrain) > 0
			% calculate Raster; i.e. the pixel-wise labeling of the B-Scan
			Classes = wrapperGetRaster(interpolation(options.EdgesTrain,:),[options.Y size(interpolation,2)]);
			
			for j = 1:length(options.LayersTrain)
				idxSave = (1:options.numPatches) + (counter-1)*options.numPatches;
				counter = counter + 1;

				% applicable columns
			    idxCols = options.BScanRegions(region,1):options.BScanRegions(region,2);
				% find indices corresponding to the selected patch position
				if strcmp(options.patchPosition,'random')
					idxTmp = find(Classes(:,idxCols) == options.LayersTrain(j)) + (options.BScanRegions(region,1)-1)*options.Y;
					randTmp = randperm(length(idxTmp));
					[A B] = ind2sub(size(Classes),idxTmp(randTmp(1:options.numPatches)));
					idxRand(1,idxSave) = A;
					idxRand(2,idxSave) = B;	
				elseif strcmp(options.patchPosition,'middle')
					% first filter the set of usable columns: find columns with at least one pixels thick layers;
					idxCols = idxCols(find(sum((Classes(:,idxCols) == options.LayersTrain(j)).*ones(size(Classes(:,idxCols))))>0));
					% the patch drawn from the image must lie inside the image
					idxCols = idxCols(round(sum((Classes(:,idxCols) == options.LayersTrain(j)).*helper(ones(1,length(idxCols)),:)')./sum(Classes(:,idxCols) == options.LayersTrain(j)))>round(pw(1)/2));
					random = randperm(length(idxCols));
				
					if length(random) < options.numPatches
						idxRand(2,idxSave(1:length(random))) = idxCols;
						idxRand(1,idxSave(1:length(random))) = round(sum((Classes(:,idxRand(2,idxSave(1:length(random)))) == options.LayersTrain(j)).*helper(ones(1,length(random)),:)')./sum(Classes(:,idxRand(2,idxSave(1:length(random)))) == options.LayersTrain(j)));
						setToZero = [setToZero idxSave(length(random)+1:options.numPatches)];
					else
						idxRand(2,idxSave) = idxCols(random(1:options.numPatches));		
						idxRand(1,idxSave) = round(sum((Classes(:,idxRand(2,idxSave)) == options.LayersTrain(j)).*helper(ones(1,options.numPatches),:)')./sum(Classes(:,idxRand(2,idxSave)) == options.LayersTrain(j)));
					end
				end
				trnData(region).type(j) = 1;
			end
		end

		% draw patches representing edges
		if length(options.EdgesTrain) > 0 
			for j = 1:length(options.EdgesTrain);
				idxSet = options.BScanRegions(region,1):options.BScanRegions(region,2);
				idxSave = (1:options.numPatches) + (counter-1)*options.numPatches;
				counter = counter + 1;
				idxSet = idxSet(find(interpolation(options.EdgesTrain(j),idxSet)>pw(1)+1));
				random = randperm(length(idxSet));
				idxRand(2,idxSave) = idxSet(random(1:options.numPatches));
				idxRand(1,idxSave) = round(interpolation(options.EdgesTrain(j),idxRand(2,idxSave)));
						
				trnData(region).type(j+length(options.LayersTrain)) = 2;
			end
		end

		%idxSave = (1:(options.numPatches*numClasses)) + (i-1)*(options.numPatches*numClasses);
		numPatchesToAdd = options.numPatches*numClasses-length(setToZero);
		idxSave = (1:numPatchesToAdd) + numPatchesAdded(region);
	   	numPatchesAdded(region) = numPatchesAdded(region) + numPatchesToAdd;

		idxRand(:,setToZero) = [];
		
		trnData(region).idx(idxSave,:) = idxRand';
		
		% transfer indices from B0 to B0Aug
		idxRand = idxRand + repmat([pw(1); pw(2)],1,numPatchesToAdd);
		patches = fetchPatches(filename,int32(idxRand),options);

		classIDs = reshape(repmat([options.LayersTrain (1:length(options.EdgesTrain)) + options.numLayers],options.numPatches,1),1,options.numPatches*numClasses);
		classIDs(setToZero) = [];
		trnData(region).classID(idxSave) = classIDs; 
%		trnData(region).data(idxSave,:) = patches.data;
		trnData(region).data = [trnData(region).data; patches.data];
	end
end

end %function
