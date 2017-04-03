function [BScanData, fileHeader, BScanHeader, BScanSeg, SLO] = HDEVolImporter(folder,filename,options)

options = setDefaultOptions(options);

HEADERSIZE = 2048;

if options.verbose > 0
	fprintf('Start reading...\n');
	fprintf('Filename: %s\n',filename);
	fprintf('Folder: %s\n',folder);
end

[pathstr,name,ext] = fileparts(filename);
% in case there is no extension, set it to .vol
if length(ext) == 0
	filename = [filename '.vol'];
elseif ~strcmp(ext,'.vol')
	disp(sprintf('%s is not a vol-file.',filename));
	BScanHeader = -1; BScanSeg = []; BScanData = []; SLO = []; fileHeader = [];
	return
end

file = fopen([folder filesep filename],'rb');
if (file==-1)
	disp(sprintf('%s could not be read.',filename));
	BScanHeader = -1; BScanSeg = []; BScanData = []; SLO = []; fileHeader = [];
	return
end

fileHeaderStruct = 	{{'Version','c',0}, ... 		% Version identifier: HSF-OCT-xxx, xxx = version number of the file format, Current version: xxx = 102
					{'SizeX','i',12},  ... 			% Number of A-Scans in each B-Scan, i.e. the width of each B-Scan in pixel
					{'NumBScans','i',16}, ... 		% Number of B-Scans in OCT scan
					{'SizeZ','i',20}, ... 			% Number of samples in an A-Scan, i.e. the Height of each B-Scan in pixel
					{'ScaleX','d',24}, ... 			% Width of a B-Scan pixel in mm
					{'Distance','d',32}, ... 		% Distance between two adjacent B-Scans in mm
					{'ScaleZ','d',40}, ... 			% Height of a B-Scan pixel in mm
					{'SizeXSlo','i',48}, ...  		% Width of the SLO image in pixel
					{'SizeYSlo','i',52}, ... 		% Height of the SLO image in pixel
					{'ScaleXSlo','d',56}, ... 		% Width of a pixel in the SLO image in mm
					{'ScaleYSlo','d',64}, ...		% Height of a pixel in the SLO image in mm
					{'FieldSizeSlo','i',72}, ... 	% Horizontal field size of the SLO image in dgr
					{'ScanFocus','d',76}, ...		% Scan focus in dpt
					{'ScanPosition','c',84}, ... 	% Examined eye (zero terminated string). "OS" for left eye; "OD" for right eye
					{'ExamTime','i',88}, ... 		% Examination time. The structure holds an unsigned 64-bit date and time value and represents the number of 100-nanosecond units	since the beginning of January 1, 1601.
					{'ScanPattern','i',96}, ...		% Scan pattern type: 0 = Unknown pattern, 1 = Single line scan (one B-Scan only), 2 = Circular scan (one B-Scan only), 3 = Volume scan in ART mode, 4 = Fast volume scan, 5 = Radial scan (aka. star pattern)
					{'BScanHdrSize','i',100}, ...	% Size of the Header preceding each B-Scan in bytes
					{'ID','c',104}, ...				% Unique identifier of this OCT-scan (zero terminated string). This is identical to the number <SerID> that is part of the file name. Format: n[.m] n and m are numbers. The extension .m exists only for ScanPattern 1 and 2. Examples: 2390, 3433.2
					{'ReferenceID','c',120}, ...	% Unique identifier of the reference OCT-scan (zero terminated string). Format: see ID, This ID is only present if the OCT-scan is 	part of a progression otherwise this string is empty. For the reference scan of a progression ID and ReferenceID are identical.
					{'PID','i',136}, ...			% Internal patient ID used by HEYEX.
					{'PatientID','c',140}, ...		% User-defined patient ID (zero terminated string).
					{'Padding','c',161}, ...		% To align next member to 4-byte boundary.
					{'DOB','date',164}, ... 		% Patient's date of birth
					{'VID','i',172}, ...			% Internal visit ID used by HEYEX.
					{'VisitID','c',176}, ...		% User-defined visit ID (zero terminated string). This ID can be defined in the Comment-field of the Diagnosis-tab of the Examination Data dialog box. The VisitID must be defined in the first row of the comment field. It has to begin with an "#" and ends with any white-space character. It can contain up to 23 alpha-numeric characters (excluding the "#").
					{'VisitDate','date',200}, ...	% Date the visit took place. Identical to the date of an examination tab in HEYEX.
					{'GridType','i',208}, ...		% Type of grid used to derive thickness data. 0 No thickness data available, >0 Type of grid used to derive thickness  values. Seeter "Thickness Grid"	for more details on thickness data, Thickness data is only available for ScanPattern 3 and 4.
					{'GridOffset','i',212}, ...		% File offset of the thickness data in the file. If GridType is 0, GridOffset is 0.
					{'GridType1','i',216}, ...		% Type of a 2nd grid used to derive a 2nd set of thickness data.
					{'GridOffset1','i',220}, ...	% File offset of the 2 nd thickness data set in the file.
					{'ProgID','c',224}, ...			% Internal progression ID (zero terminated string). All scans of the same progression share this ID.
					{'Spare','c',258}}; 			% Spare bytes for future use. Initialized to 0.	   

% the number of data fields depending on the version
versionSize = {{'100',19}, {'101',26}, {'102',28}, {'103',31}};

% read verion
version = fread(file,12,'int8=>char');
if options.verbose > 0
	fprintf('Version: %s\n',version);
end
version = version(9:fileHeaderStruct{2}{3}-1)';

numDataFields = 0;
for i = 1:length(versionSize)
	if strcmp(version,versionSize{i}{1})
		numDataFields = versionSize{i}{2};
	end
end

if (numDataFields == 0)
	warning(sprintf('Current version %s is newer than last documented version %s, default to this last version',version,versionSize{end}{1}));
	numDataFields = versionSize{end}{2};
end

% point file to the beginning 
frewind(file);

% read file header
fileHeader = readHeader(file,fileHeaderStruct,numDataFields);

if options.verbose > 0
	if strcmp(fileHeader.ScanPosition(1:2),'OD')
		fprintf('Right eye examined\n');
	else
		fprintf('Left eye examined\n');
	end
end

% point file to end of header
fseek(file, HEADERSIZE, 'bof');

% read SLO image, dimensions encoded in SizeXSLO and SizeYSLO; each pixel gray value between 0 and 255
SLO = reshape(fread(file,fileHeader.SizeXSlo*fileHeader.SizeYSlo,'uint8'),fileHeader.SizeXSlo,fileHeader.SizeYSlo)';

% read B0 scans
BScanHeaderLength = fileHeader.BScanHdrSize; % size in bytes
BScanDim = [fileHeader.SizeX,fileHeader.SizeZ]; % dimension in pixels

BScansNum = fileHeader.NumBScans; % number of BScans

BScanHeaderStruct = {{'Version','c',0}, ... 		% Version identifier (zero terminated string). Version Format: "HSF-BS-xxx, xxx = version number of the B-Scan header format. Current version: xxx = 103
					{'BScanHdrSize','i',12}, ...	% Size of the B-Scan header in bytes. It is identical to the same value of the file header.
					{'StartX','d',16}, ...		% X-Coordinate of the B-Scan's start point in mm.
					{'StartY','d',24}, ...		% Y-Coordinate of the B-Scan's start point in mm.
					{'EndX','d',32}, ...			% X-Coordinate of the B-Scan's end point in mm. For circle scans, this is the X-Coordinate of the circle's center point.
					{'EndY','d',40}, ...			% Y-Coordinate of the B-Scan's end point in mm. For circle scans, this is the Y-Coordinate of the circle's center point.
					{'NumSeg','i',48}, ...		% Number of segmentation vectors
					{'OffSeg','i',52}, ...		% Offset of the array of segmentation vectors relative to the beginning of this B-Scan header.
					{'Quality','f',56}, ...		% Image quality measure. If this value does not exist, its value is set to INVALID.
					{'Shift','i',60}, ...			% Horizontal shift (in # of A-Scans) of the classification band against the segmentation lines (for circular scan only).
					{'IVTrafo','f',64}, ...		% Intra volume transformation matrix. The values are only available for volume and radial scans and if alignment is turned off, otherwise the values are initialized to 0.
					{'Spare','c',88}};			% Spare bytes for future use.

versionSize = {{'100',8},{'101',9},{'102',10},{'103',11}};

numDataFields = 0;
for i = 1:length(versionSize)
    if strcmp(version,versionSize{i}{1})
        numDataFields = versionSize{i}{2};
    end
end

if (numDataFields == 0)
    warning(sprintf('Current version %s is newer than last documented version %s, default to this last version (may produce errors)',version,versionSize{end}{1}));
    numDataFields = versionSize{end}{2};
end

BScanHeader = cell(BScansNum,1);
BScanSeg = cell(BScansNum,1);
BScanData = cell(BScansNum,1); 

currPosition = ftell(file);


for n = 1:BScansNum
	% read B-Scan header
	BScanHeader{n} = readHeader(file,BScanHeaderStruct,numDataFields);

	fseek(file, currPosition+BScanHeader{n}.OffSeg, 'bof');
	% read segmentation (number of segmentation lines * number of A-Scans)
	BScanSeg{n} = reshape(fread(file,BScanHeader{n}.NumSeg*fileHeader.SizeX,'float'),fileHeader.SizeX,BScanHeader{n}.NumSeg);
	BScanSeg{n}(BScanSeg{n} == realmax('single')) = NaN;

	% read B-Scan itself
	fseek(file, currPosition+BScanHeader{n}.BScanHdrSize, 'bof');
	BScanData{n} = reshape(fread(file,fileHeader.SizeX*fileHeader.SizeZ,'float'),fileHeader.SizeX,fileHeader.SizeZ);
	BScanData{n}(BScanData{n} == realmax('single')) = NaN;
	BScanData{n} = BScanData{n}';
	
	currPosition = ftell(file);
end

if isfield(options,'BScansSelect')
	BScanHeader = BScanHeader(options.BScansSelect);
	BScanSeg = BScanSeg(options.BScansSelect);
	BScanData = BScanData(options.BScansSelect);
	BScansNum = length(options.BScansSelect);
end
	
if options.plotBScans
	for n = 1:BScansNum
		figure; imagesc(sqrt(sqrt(BScanData{n}))); title(sprintf('B-Scan %d',n)); colormap gray;
	end	
end

% print aera covered by BScans
sizeX = norm([BScanHeader{end}.EndX BScanHeader{end}.EndY]-[BScanHeader{end}.StartX BScanHeader{end}.StartY]);
sizeY = norm([BScanHeader{end}.StartX BScanHeader{end}.StartY]-[BScanHeader{1}.StartX BScanHeader{1}.StartY]);
if options.verbose > 0
	fprintf('Number of BScans: %d, resolution: %d px x %d px\n',length(BScanHeader),fileHeader.SizeX,fileHeader.SizeZ);
	fprintf('Size SLO Scan: %.d x %.d; mm per Pixel: %.4f, %.4f\n',fileHeader.SizeXSlo,fileHeader.SizeYSlo,fileHeader.ScaleXSlo,fileHeader.ScaleYSlo);
	fprintf('Area covered (X x Y): %.2f mm x %.2f mm = %.2f mm^2\n',sizeX,sizeY,sizeX*sizeY);
	fprintf('Area covered in px (X x Y): %.2f px x %.2f px = %.2f px^2\n',sizeX/fileHeader.ScaleXSlo,sizeY/fileHeader.ScaleYSlo,sizeX*sizeY/(fileHeader.ScaleXSlo*fileHeader.ScaleYSlo));
	fprintf('Distance between B-Scans: %.4f mm, %.4f px\n',sizeY/(fileHeader.NumBScans-1),sizeY/fileHeader.ScaleYSlo/(fileHeader.NumBScans-1));
end
%fileHeader.distanceBScans = fileHeader.Distance/fileHeader.ScaleYSlo;
fileHeader.distanceBScans = fileHeader.Distance;

if options.plotSLO
	figure; imagesc(sqrt(sqrt(SLO))); t = title(['SLO ' filename]); colormap gray;
	set(t,'Interpreter','none'); 
	xlabel('X');
	ylabel('Y');
	hold on;

	% plot positions of BScans
	for n = 1:BScansNum
		hLine = line([BScanHeader{n}.StartX/fileHeader.ScaleXSlo BScanHeader{n}.EndX/fileHeader.ScaleXSlo],[BScanHeader{n}.StartY/fileHeader.ScaleYSlo BScanHeader{n}.EndY/fileHeader.ScaleYSlo]);
		if n == floor(BScansNum/2)+1
			set(hLine,'Color','red');
		end

		if mod(n,2)
			text(BScanHeader{n}.EndX/fileHeader.ScaleXSlo+10,BScanHeader{n}.EndY/fileHeader.ScaleYSlo,num2str(n),'Color','white');
		else
			text(BScanHeader{n}.EndX/fileHeader.ScaleXSlo+50,BScanHeader{n}.EndY/fileHeader.ScaleYSlo,num2str(n),'Color','white');
		end
	end
end

end

function options = setDefaultOptions(options)
	if ~isfield(options,'plotSLO')
		options.plotSLO = 0;
	end

	if ~isfield(options,'plotBScans')
		options.plotBScans = 0;
	end

	if ~isfield(options,'verbose')
		options.verbose = 1;
	end
end

function header = readHeader(file,headerStruct,numDataFields)
	for i = 1:numDataFields
		fieldLength = headerStruct{i+1}{3}-headerStruct{i}{3};
		if (strcmp(headerStruct{i}{2}(1),'c')) % string
			value = fread(file,fieldLength,'int8=>char')';
		elseif (strcmp(headerStruct{i}{2}(1),'i')) % integer
			numInts = fieldLength/4; value = [];
			for j = 1:numInts
				value(j) = fread(file,1,'int32');
			end
		elseif (strcmp(headerStruct{i}{2}(1),'f')) % float
			numFloats = fieldLength/4; value = [];
			for j = 1:numFloats
				value(j) = fread(file,1,'float');
			end
		elseif (strcmp(headerStruct{i}{2}(1),'d') || strcmp(headerStruct{i}{2}(1),'date')) % double or date
			numDoubles = fieldLength/8; value=[];
			for j = 1:numDoubles
				value(j) = fread(file,1,'double');
			end
		end

%		header(i,:) = {headerStruct{i}{1},value};
		header.(headerStruct{i}{1}) = value;
	end
end
