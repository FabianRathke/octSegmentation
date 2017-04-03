function [BScanData, fileHeader, BScanHeader, BScanSeg, SLO] = HDEMatImporter(folder,filename,options)

options = setDefaultOptions(options);

if options.verbose > 0
	fprintf('Start reading...\n');
	fprintf('Filename: %s\n',filename);
end

[pathstr,name,ext] = fileparts(filename);
% in case there is no extension, set it to .vol
if length(ext) == 0
	filename = [filename '.mat'];
elseif ~strcmp(ext,'.mat')
	disp(sprintf('%s is not a mat-file.',filename));
	BScanHeader = -1; BScanSeg = []; BScanData = []; SLO = []; fileHeader = [];
	return
end

if (~exist([folder filesep filename]))
	disp(sprintf('%s could not be read.',filename));
	BScanHeader = -1; BScanSeg = []; BScanData = []; SLO = []; fileHeader = [];
	return
else
	load([folder filesep filename]);
end
if exist('SLO')
	SLO = double(SLO);
else
	SLO = [];
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
					{'ProgID','c',224}};			% Internal progression ID (zero terminated string). All scans of the same progression share this ID.

% check which variables are available
fileHeader = struct();
for i = 1:length(fileHeaderStruct)
	string = fileHeaderStruct{i}{1};
	if exist(string)
		eval(sprintf('fileHeader.%s = %s;',string,string));
	end
end

fieldNames = fieldnames(fileHeader);
for i = 1:length(fieldNames)
	val = getfield(fileHeader,fieldNames{i});
	if isinteger(val)
		fileHeader = setfield(fileHeader,fieldNames{i},int32(val));
	end
end

% get the number of BScans
varnames = who;
BScansNum = sum(cell2mat(regexpi(varnames,'^B\d{1,}$')));
if BScansNum == 1
	fileHeader.NumBScans = BScansNum;
end
if ~isfield(fileHeader,'ScanPosition')
	fileHeader.ScanPosition = 'OS';
end

BScanSeg = cell(BScansNum,1);
BScanData = cell(BScansNum,1); 
BScanHeader = cell(BScansNum,1);

headerVars = {'StartX','StartY','EndX','EndY'};

for n = 1:BScansNum
	BScanSeg{n} = eval(sprintf('B%dseg;',n-1));
	BScanData{n} = eval(sprintf('B%d',n-1));
  	BScanSeg{n}(BScanSeg{n} == realmax('single')) = NaN;
	BScanData{n}(BScanData{n} == realmax('single')) = NaN;

	for i = 1:length(headerVars)
		if BScansNum > 1
			eval(sprintf('BScanHeader{%d}.%s = B%d%s;',n,headerVars{i},n-1,headerVars{i}));
		end
	end
end

[fileHeader.SizeZ fileHeader.SizeX] = size(BScanData{1});

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

fileHeader.distanceBScans = 0;
if BScansNum > 1 && options.verbose > 0
	sizeX = norm([BScanHeader{end}.EndX BScanHeader{end}.EndY]-[BScanHeader{end}.StartX BScanHeader{end}.StartY]);
	sizeY = norm([BScanHeader{end}.StartX BScanHeader{end}.StartY]-[BScanHeader{1}.StartX BScanHeader{1}.StartY]);
	if options.verbose > 0
		fprintf('Number of BScans: %d, resolution: %d px x %d px\n',length(BScanHeader),fileHeader.SizeX,fileHeader.SizeZ);
		fprintf('Size SLO Scan: %.d x %.d; mm per Pixel: %.4f, %.4f\n',fileHeader.SizeXSlo,fileHeader.SizeYSlo,fileHeader.ScaleXSlo,fileHeader.ScaleYSlo);
		fprintf('Area covered (X x Y): %.2f mm x %.2f mm = %.2f mm^2\n',sizeX,sizeY,sizeX*sizeY);
		fprintf('Area covered in px (X x Y): %.2f px x %.2f px = %.2f px^2\n',sizeX/fileHeader.ScaleXSlo,sizeY/fileHeader.ScaleYSlo,sizeX*sizeY/(fileHeader.ScaleXSlo*fileHeader.ScaleYSlo));
		fprintf('Distance between B-Scans: %.4f mm, %.4f px\n',sizeY/double(fileHeader.NumBScans-1),sizeY/fileHeader.ScaleYSlo/double(fileHeader.NumBScans-1));
	end
%	fileHeader.distanceBScans = sizeY/fileHeader.ScaleYSlo/double(fileHeader.NumBScans-1);
	fileHeader.distanceBScans = sizeY/double(fileHeader.NumBScans-1);
	fileHeader.sizeY = sizeY;
	fileHeader.sizeX = sizeX;
end

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
		options.verbose = 0;
	end
end
