function testData = collectTestData(file,options)
% trnData - draws patches from a set of training images for training the apperance model 
% 
% Inputs:
%   files - [struct] list of files
%       .name - [string] the filename
%   options -[struct]
%       .BScanRegions  - [matrix](2xn) n subdivisions of the B-Scan; each row holds x and y coordinates; 
%       .LayersPred    - [array] indices of layers to predict; should coincide with .LayersTrain
%       .EdgesPred     - [array] indices of edges to predict; should coincide with .EdgesTrain
%       .Y             - [int] height of the BScan
%       .height        - [int] patch height
%       .width         - [int] patch width
%       .testFullImage - [boolean] if true fetches all patches in each file, otherwise only a subset (not implemented)
%
% Outputs:
%   testData - [struct]
%     .data     - 
%     .idx 
%     .classID
%
% See also: collectTrnData
% Calls: fetchPatches

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 06-Dec-2013

pw = ([options.height options.width]-1)/2;
[a filename] = fileparts(file.name);

numBounds = length(options.EdgesPred); numLayers = length(options.LayersPred);

% take the full image as test set
if options.testFullImage
	[a b] = meshgrid(options.columnsPred+pw(2),1+pw(1):options.Y+pw(1));
	idxSet = [reshape(b,1,options.Y*length(options.columnsPred)); reshape(a,1,options.Y*length(options.columnsPred))];
	testData = fetchPatches(filename,idxSet,options);
	testData.idx = int16(idxSet' - pw(ones(1,size(idxSet,2)),:));

	% set class to boundary classes where transitions between layers occur (change of index)
	if options.saveAppearanceTerms
		tmp = testData.classID(2:end)-testData.classID(1:end-1);
		idx = find(tmp==1);
		bounds = repmat((numLayers+1):(numLayers+numBounds),1,length(options.columnsPred));
		testData.classID(idx+1) = testData.classID(idx)+length(options.LayersPred);
		testData.classID(idx) = testData.classID(idx+1);
	end
else % choose a pixel subset for testing
	error('Old version deleted, implement newly using the implementation from collectorTrnData');
end
