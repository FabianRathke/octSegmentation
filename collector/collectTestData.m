function testData = collectTestData(file,options)
% collectTestData - draws patches from a set of training images for training the apperance model 
%
% Syntax:
%   testData = collectTestData(file,options)
% 
% Inputs:
%   files - [struct] list of files
%       .name - [string] the filename
%   options - [struct] holds options for fetching patches
%       .BScanRegions - [matrix](2xn) n subdivisions of the B-Scan; each row holds x and y coordinates; 
%       .LayersPred   - [array] indices of layers to predict; should coincide with .LayersTrain
%       .EdgesPred    - [array] indices of edges to predict; should coincide with .EdgesTrain
%       .Y            - [int] height of the BScan
%       .height       - [int] patch height
%       .width        - [int] patch width
%
% Outputs:
%   testData - [struct] holds different matrices containing patches from the test scan
%     .data    - [matrix]
%     .idx     - [matrix] 
%     .classID - [array]
%
% See also: collectTrnData, setCollectorDefaults

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 02-May-2014

pw = ([options.height options.width]-1)/2;
%[a filename] = fileparts(file.name);

% determines which columns and which rows to fetch from the B-scan
if isfield(options,'idxSet')
	idxSet = options.idxSet;
	idxSet(1,:) = idxSet(1,:) + pw(1); 
	idxSet(2,:) = idxSet(2,:) + pw(2);
else
	[a b] = meshgrid(int32(options.columnsPred+pw(2)),int32(1+pw(1):options.Y+pw(1)));
	idxSet = [reshape(b,1,options.Y*length(options.columnsPred)); reshape(a,1,options.Y*length(options.columnsPred))];
end
testData = fetchPatches(file.name,idxSet,options);
testData.idx = int16(idxSet') - int16(pw(ones(1,size(idxSet,2)),:));

% set class to boundary classes where transitions between layers occur (change of index)
%numBounds = length(options.EdgesPred); numLayers = length(options.LayersPred);
%if options.saveAppearanceTerms
%	tmp = testData.classID(2:end)-testData.classID(1:end-1);
%	idx = find(tmp==1);
%	bounds = repmat((numLayers+1):(numLayers+numBounds),1,length(options.columnsPred));
%	testData.classID(idx+1) = testData.classID(idx)+length(options.LayersPred);
%	testData.classID(idx) = testData.classID(idx+1);
%end
