function B0 = loadData(filename,options)
% loadData - loads B-Scan as defined in filename, folder is given by options.folder_data;
% 
% The loading logic for a given dataset, controlled by options.loadRoutineData has to be implemented here individually for each dataset. The output should be a matrix holding a B-Scan. If the file holds several B-Scans options.labelID can be used to pick the correct one.
%
% Inputs:
%	filename - [string] filename of the matfile to load
%	options  - [struct] controller options
%     .folder_data     - [string] points to the folder of filename
%     .labelID 	       - [int] Spectralis B-Scans are labeled B0,B1,..., labelID is ID of the B-Scan to load
%     .clip		       - [boolean] indicates whether to clip a B-Scan at the left and right border
%     .clipRange	   - [array](2) defines the left and right border of the B-Scan after clipping
%     .loadRoutineData - [string] the user-defined routine that is used to load the data
%
% Outputs:
%	B0 - [matrix] the B-Scan 
%
% See also: loadLabels, fetchPatches, collectTrnData, collectTestData

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 26-Feb-2014

if strcmp(options.loadRoutineData,'spectralis')
	load([options.folder_data filename '.mat'],sprintf('B%d',options.labelID));
	eval(sprintf('B0 = B%d;',options.labelID));
	B0(B0>10) = 0;
	B0 = sqrt(sqrt(B0));
elseif strcmp(options.loadRoutineData,'AMDDataset')
	load([options.folder_data filename],'images');
	B0 = squeeze(images(:,:,options.labelID));
else
	error('Please specify a valid routine for fetching data in collector.options.loadRoutineData');
end

printMessage(sprintf('Loaded data for %s and region %d.\n',filename,options.labelID),2,options.verbose);

if options.clip
	B0 = B0(:,options.clipRange(1):options.clipRange(2));
end

% apply scan-wise preprocessing
for i = 1:length(options.preprocessing.scanLevel)
	B0 = options.preprocessing.scanLevel{i}{1}(B0,options.preprocessing.scanLevel{i});
end

end
