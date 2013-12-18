function B0 = loadData(filename,options)
% loadData - loads B-Scan as defined in filename, folder is given by options.folder_data;
% 
% The loading logic for a given dataset, controlled by options.dataset has to be implemented here. The output should be a matrix holding the B-Scan
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
% Last Revision: 13-Dec-2013

if strcmp(options.loadRoutineData,'spectralis')
	if isfield(options,'labelID')
		load([options.folder_data filename '.mat'],sprintf('B%d',options.labelID));
		eval(sprintf('B0 = B%d;',options.labelID));
		printMessage(sprintf('Loaded data for %s and region %d.\n',filename,options.labelID),2,options.verbose);
	else
	   load([options.folder_data filename '.mat'],'B0');
	end  
	B0(B0>10) = 0;
	B0 = sqrt(sqrt(B0));
elseif strcmp(options.loadRoutineData,'AMDDataset')
	load([options.folder_data filename],'images');
	B0 = squeeze(images(:,:,options.labelID));
end

if options.clip
	B0 = B0(:,options.clipRange(1):options.clipRange(2));
end

for i = 1:length(options.preprocessing.scanLevel)
	B0 = options.preprocessing.scanLevel{i}{1}(B0,options.preprocessing.scanLevel{i});
end

%if options.addNoise
%	fprintf('Add speckle noise with variance .2f to the image\n',options.addNoise);
%	B0 = imnoise(B0,'speckle',options.addNoise);
%	B0 = B0 + sqrt(12*options.addNoise)*B0.*(rand(size(B0))-.5);
%	B0 = max(0,min(B0,1));
%end

end
