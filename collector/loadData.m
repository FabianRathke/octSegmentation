function B0 = loadData(filename,options)
% loadData - loads B-Scan as defined in filename, folder is given by options.folder_data;
% 
% The loading logic for a given dataset, controlled by options.loadRoutineData has to be implemented here individually for each dataset. The output should be a matrix holding a B-Scan. If the file holds several B-Scans options.labelID can be used to pick the correct one.
%
% Inputs:
%	filename - [string] filename of the matfile to load
%	options  - [struct] collector options
%     .folder_data     	- [string] points to the folder of filename
%     .labelID 	       	- [int] Spectralis B-Scans are labeled B0,B1,..., labelID is ID of the B-Scan to load
%     .clip		       	- [boolean] indicates whether to clip a B-Scan at the left and right border
%     .clipRange	   	- [array](2) defines the left and right border of the B-Scan after clipping
%     .clipFactor	   	- [int] the interval in y-direction (i.e. 2 means that every second column is grabbed)
%     .loadRoutineData 	- [string] the user-defined routine that is used to load the data
%	  .mirrowBScan		- [string] defines a trigger to mirrow the scan (left/right eye)
%
% Outputs:
%	B0 - [matrix] the B-Scan 
%
% See also: loadLabels, fetchPatches, collectTrnData, collectTestData

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 27-Feb-2016

[pathstr,name,ext] = fileparts(filename);

if strcmp(options.loadRoutineData,'spectralisMat')
	if isempty(ext)
		ext = '.mat';
	end
	load([options.folder_data name ext],sprintf('B%d',options.labelID-1),'ScanPosition');
	eval(sprintf('B0 = B%d;',options.labelID-1));
	B0(B0>10) = 0;
	B0 = sqrt(sqrt(B0));
	if length(options.mirrorBScan) > 0
		if findstr(ScanPosition,options.mirrorBScan)
			B0 = B0(:,end:-1:1);
		end
	end
elseif strcmp(options.loadRoutineData,'spectralisVol')
	if isempty(ext)
		ext = '.vol';
	end
	% which slice of the volume should be loaded
	optionsVolImport = struct('BScansSelect',options.labelID, 'verbose',0);
	[BScanData fileHeader] = HDEVolImporter(options.folder_data,[name ext],optionsVolImport);
	if (~iscell(BScanData))
		error(sprintf('Loading %s%s failed',options.folder_data,filename));
	end
	B0 = BScanData{1};
	B0(B0>10) = 0;
	B0(isnan(B0)) = 0;
	B0 = sqrt(sqrt(B0));
	% if left eye, rotate the scan
	if length(options.mirrorBScan) > 0
		if findstr(fileHeader.ScanPosition,options.mirrorBScan)
			B0 = B0(:,end:-1:1);
		end
	end
elseif strcmp(options.loadRoutineData,'AMDDataset')
	if isempty(ext)
		ext = '.mat';
	end
	load([options.folder_data name ext],'images');
	B0 = squeeze(images(:,:,options.labelID));
else
	error('Please specify a valid routine for fetching data in collector.options.loadRoutineData');
end

printMessage(sprintf('Loaded data for %s and region %d.\n',filename,options.labelID),2,options.verbose);

if options.clip
	B0 = B0(:,options.clipRange(1):options.clipFactor:options.clipRange(2));
end

% apply scan-wise preprocessing
for i = 1:length(options.preprocessing.scanLevel)
	B0 = options.preprocessing.scanLevel{i}{1}(B0,options.preprocessing.scanLevel{i});
end

end
