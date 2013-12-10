function B0 = loadData(filename,options)
% loadData - loads spectralis mat files from filename, folder is given by options.folder_data
% 
% Inputs:
%	filename - [string] filename of the matfile to load
%	options  - [struct]
%		.folder_data - [string] points to the folder of filename
%		.labelID 	 - [int] Spectralis B-Scans are labeled B0,B1,..., labelID is ID of the B-Scan to load
%		.clip		 - [boolean] indicates whether to clip a B-Scan at the left and right border
%		.clipRange	 - [array](2) defines the left and right border of the B-Scan after clipping
%
% Outputs:
%	B0 - [matrix] the B-Scan 
%
% See also: loadLabels, fetchPatches, collectTrnData, collectTestData

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 05-Dec-2013

if isfield(options,'labelID')
   	load([options.folder_data filename '.mat'],sprintf('B%d',options.labelID));
   	eval(sprintf('B0 = B%d;',options.labelID));
	printMessage(sprintf('Loaded data for %s and region %d.\n',filename,options.labelID),2,options.verbose);
else
   load([options.folder_data filename '.mat'],'B0');
end  

% The spectralis scans are sometimes set to values > 1, correct this
B0(B0>1) = 0;

if options.clip
	B0 = B0(:,options.clipRange(1):options.clipRange(2));
end

%if options.addNoise
%	fprintf('Add speckle noise with variance .2f to the image\n',options.addNoise);
%	B0 = imnoise(B0,'speckle',options.addNoise);
%	B0 = B0 + sqrt(12*options.addNoise)*B0.*(rand(size(B0))-.5);
%	B0 = max(0,min(B0,1));
%end

end
