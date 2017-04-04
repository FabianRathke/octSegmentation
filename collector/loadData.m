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
%	  .mirrowBScan		- [string] defines a trigger to mirrow the scan --> the model is trained on left eyes (left/right eye)
%
% Outputs:
%	B0 - [matrix] the B-Scan 
%
% See also: loadLabels, fetchPatches, collectTrnData, collectTestData

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 20-Nov-2016

[pathstr,name,ext] = fileparts(filename);

if strcmp(options.loadRoutineData,'spectralisMat')
	if isempty(ext)
		ext = '.mat';
	end
	load([options.folder_data name ext],sprintf('B%d',options.labelID-1));
	vars = whos('-file',[options.folder_data name ext]);
	if ismember('ScanPosition', {vars.name})
		load([options.folder_data name ext],'ScanPosition');
	end
	eval(sprintf('B0 = B%d;',options.labelID-1));
	B0(B0>1) = 0.1;
    B0(isnan(B0)) = 0.1;
	B0 = sqrt(sqrt(B0));
	if length(options.mirrorBScan) > 0
		fprintf('Rotate scan\n');
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
	B0(B0>10) = 0.1;
	B0(isnan(B0)) = 0.1;
	B0 = sqrt(sqrt(B0));
	% if left eye, rotate the scan
	if length(options.mirrorBScan) > 0
		if findstr(fileHeader.ScanPosition,options.mirrorBScan)
			B0 = B0(:,end:-1:1);
		end
	end
elseif strcmp(options.loadRoutineData,'TianDataset')
	% left/right eye has to be set
	%rotate = [0 1 1 0 0 0 0 0 0 1];
	load([options.folder_data name '.mat'],'volumedataClipped')
	B0 = single(squeeze(volumedataClipped(:,:,options.labelID)))/255;

	%	s = regexp(filename,'(\d{1,2})','tokens');
	%filenumber = str2num(s{1}{1});
	%if rotate(filenumber)
	%	B0 = B0(:,end:-1:1);
	%end
elseif strcmp(options.loadRoutineData,'ChiuDME')
	% left/right eye has to be set
	rotate = [0 0 0 0 1 0 1 0 1 0];
%	offset = [0 0 0 -4 -2 0 1 -2 0 0]; % correct offset of central scan (determined by the labels, see 7.2 of Chiu et al. 2015)
	load([options.folder_data name '.mat'],'images')
	s = regexp(filename,'(\d{1,2})','tokens');
	filenumber = str2num(s{1}{1});
	B0 = squeeze(images(:,:,options.labelID))/255;
%	B0(B0==1) = rand(1,sum(B0(:)==1))*0.2;
	B0(B0==1) = 0.2;

	if rotate(filenumber)
		B0 = B0(:,end:-1:1);
	end
elseif strcmp(options.loadRoutineData,'ChiuAMD')
	% left/right eye has to be set
	rotate = [[0 0 0 1 1]; [0 0 0 0 0]; [0 0 0 0 0]; [0 0 0 0 0];];
	load([options.folder_data name '.mat'],'images')
	s = regexp(filename,'(\d{1,2})','tokens');
	group = str2num(s{1}{1}); volume = str2num(s{2}{1});
	tmp = squeeze(images(:,:,options.labelID))/255;
	[X Y] = meshgrid(1:512,1:1000);
	% convert to HDE Spectralis geometrics --> in order to fit our learned appearance and shape models
	[X1 Y1] = meshgrid((1:423)*1.21,(1:568)*1.76);
	B0 = [zeros(496-423,568); interp2(X,Y,tmp',X1,Y1)'];

	if rotate(group,volume)
		B0 = B0(:,end:-1:1);
	end
elseif strcmp(options.loadRoutineData,'TianDatasetP')
	name = 'jbio201500239-sup-0003-Data-S1';
	load([options.folder_data name '.mat'],'volumedata');
	s = regexp(filename,'(\d{1,2})','tokens');
	filenumber = str2num(s{1}{1});
	B0 = single(squeeze(volumedata(:,:,(filenumber-1)*5 + options.labelID)))/255;

	% left/right eye has to be set
	rotate = [1 0 1 0 0 0 0 0 0 0];
	if rotate(filenumber)
		B0 = B0(:,end:-1:1);
	end
elseif strcmp(options.loadRoutineData,'Srinivasan')
	% AMD
	rotateAMD = [ones(1,10) zeros(1,4) 1];
	% DME
   	rotateDME = [1 1 0 1 1 0 1 1 0 0 1 0 0 1 0];
	% NORMAL
    rotateNORMAL = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
	if options.labelID < 10
		B0 = single(imread([options.folder_data filename sprintf('/TIFFs/8bitTIFFs/%02.0f.tif',options.labelID)]))/255;
	else
		B0 = single(imread([options.folder_data filename sprintf('/TIFFs/8bitTIFFs/%03.0f.tif',options.labelID)]))/255;
	end
	numCols = size(B0,2);
  	if numCols == 1024
    	options.clipFactor = 2;
        options.clipRange = [round((numCols-1000)/2)*ones(1,2) + [1 1000]];
    else
    	options.clipFactor = 1;
        options.clipRange = [round((numCols-500)/2)*ones(1,2) + [1 500]];
    end

	s = regexp(filename,'(\d{1,2})','tokens');
	filenumber = str2num(s{1}{1});
	s = regexp(filename,'([A-Z]{1,})','tokens');

	if eval(sprintf('rotate%s(%d)',s{1}{1},filenumber))
		B0 = B0(:,end:-1:1);
	end
else
	error('Please specify a valid routine for fetching data in collector.options.loadRoutineData');
end

printMessage(sprintf('Loaded data for %s and region %d.\n',filename,options.labelID),2,options.verbose);

%addpath('/export/home/frathke/Documents/Code/Libs/Filters/Denoising/BM3D');
%tic; [~,B0] = BM3D(1,B0,50); toc;

if options.clip
	B0 = B0(:,options.clipRange(1):options.clipFactor:options.clipRange(2));
end

% apply scan-wise preprocessing
if isfield(options.preprocessing,'scanLevel')
    for i = 1:length(options.preprocessing.scanLevel)
        B0 = options.preprocessing.scanLevel{i}{1}(B0,options.preprocessing.scanLevel{i});
    end
end

end
