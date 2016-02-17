function octGUI()
%margins = load('~/healthyMargins.mat');
close all;

% define global variables --> nested functions have access to these variables;
hLine = [];
currentBScan = 0; numBScans = 0;
BScanHeader = ''; BScanData = ''; BScanSeg = ''; SLO = ''; fileHeader = ''; predictions = ''; funcVal = '';
model = ''; collector = '';
fileNameScan = ''; pathNameScan = '';

f = figure('Position', [100 100 1200 700],'Tag','mainWindow');

tgroup = uitabgroup('Parent', f);
tab1 = uitab('Parent', tgroup, 'Title', 'Viewer');
tab2 = uitab('Parent', tgroup, 'Title', 'Model params');

% Position of BScan window
BScanBottom = 100;
BScanLeft = 325;
BScanWidth = 800;
BScanHeight = 500;

% Position of SLO window
SLOBottom = 100;
SLOLeft = 40;
SLOWidth = 250;
SLOHeight = 250;
% axis
axisSLO = axes('Parent',tab1,'Units','Pixels','Position',[SLOLeft,SLOBottom,SLOWidth,SLOHeight]); set(axisSLO,'YTickLabel',[],'XTickLabel',[]);
axisBScan = axes('Parent',tab1,'Units','Pixels','Position',[BScanLeft,BScanBottom,BScanWidth,BScanHeight]); set(axisBScan,'YTickLabel',[],'XTickLabel',[],'Tag','axisBScan');

% Control bar to the left
BBottom = SLOBottom + SLOHeight;

% open scan dialog
openFile = uicontrol('Parent',tab1,'Style','pushbutton','String','Open Scan','Position',[60,BBottom + 20,75,25],'Callback',{@openFile_Callback},'Tag','openFile');
selectFileType = uicontrol('Parent',tab1,'Style','popupmenu','String',{'.vol','.mat'},'Position',[145,BBottom + 20,50,25],'Tag','selectFileType');
loadModel = uicontrol('Parent',tab1,'Style','pushbutton','String','Load Model','Position',[60,BBottom+60,100,25],'Enable','off','Tag','loadModel','Callback',{@loadModel_Callback});
% button to segment the currently display BScan
segmentBScan = uicontrol('Parent',tab1,'Style','pushbutton','String','Segment','Position',[60,BBottom+100,80,25],'Enable','off','Tag','segmentBScan','Callback',{@segmentBScan_Callback});
qualityBScan = uicontrol('Parent',tab1,'Style','checkbox','Position',[150,BBottom+100,25,25],'Enable','off','Value',0','Tag','qualityBScan');
segmentAllBScans = uicontrol('Parent',tab1,'Style','pushbutton','String','Segment All','Position',[60,BBottom + 140,100,25],'Enable','off','Tag','segmentAllBScans','Callback',{@segmentBScan_Callback});

% buttons to circle through B-Scans in a volume
nextBScan = uicontrol('Parent',tab1,'Style','pushbutton','String','>','Position',[BScanLeft+BScanWidth/2+50,BScanBottom-75,25,25],'Enable','off','Tag','nextBScan','Callback',{@switchBScan_Callback});
prevBScan = uicontrol('Parent',tab1,'Style','pushbutton','String','<','Position',[BScanLeft+BScanWidth/2-50,BScanBottom-75,25,25],'Enable','off','Tag','prevBScan','Callback',{@switchBScan_Callback});
% text field showing current BScan number
selectBScan = uicontrol('Parent',tab1,'Style','popupmenu','String','','Position',[BScanLeft+BScanWidth/2-10,BScanBottom-75,50,25],'Tag','selectBScan','Callback',{@switchBScan_Callback},'Visible','off');

% text field that shows status reports
statusText = uicontrol('Parent',tab1,'Style','text','String','Waiting for action...','Position',[40,10,500,25],'Tag','statusText','HorizontalAlignment','Left');

% buttons to switch what to plot
showDataQuality = uicontrol('Parent',tab1,'Style','pushbutton','String','DataTerm','Position',[BScanLeft BScanBottom+BScanHeight+10 75 25],'Enable','off','Callback',{@showDataQuality_Callback}); 
showShapeQuality = uicontrol('Parent',tab1,'Style','pushbutton','String','ShapeTerm','Position',[BScanLeft+85 BScanBottom+BScanHeight+10 100 25],'Enable','off','Callback',{@showShapeQuality_Callback}); 
noQuality = uicontrol('Parent',tab1,'Style','pushbutton','String','Pred','Position',[BScanLeft+195 BScanBottom+BScanHeight+10 70 25],'Enable','off','Callback',{@noQuality_Callback}); 
noPred = uicontrol('Parent',tab1,'Style','pushbutton','String','NoPred','Position',[BScanLeft+275 BScanBottom+BScanHeight+10 70 25],'Enable','on','Callback',{@noPred_Callback}); 

function openFile_Callback(hObject,eventdata,handles)
	h = findobj('Tag','selectFileType');
	String = get(h,'String');
	Value = get(h,'Value');
	if length(pathNameScan) > 0
		[FileName,PathName,FilterIndex] = uigetfile(String{Value},'Load Scan',pathNameScan);
	else
		[FileName,PathName,FilterIndex] = uigetfile(String{Value},'Load Scan');
	end


	if FileName ~= 0
		 fileNameScan = FileName; pathNameScan = PathName;

		[BScanData, BScanHeader, BScanSeg, SLO, fileHeader] = HDEVolImporter(PathName,FileName,struct('plotSLO',0,'plotBScans',0));
		numBScans = fileHeader.NumBScans;
		updateStatus(sprintf('%s loaded.',FileName));
		set(findobj('Tag','loadModel'),'Enable','on');

		% plot SLO scan
		axes(axisSLO); reset(axisSLO);
		imagesc(sqrt(sqrt(SLO))); t = title('SLO'); colormap gray;
		set(t,'Interpreter','none');
		hold on;
		predictions = cell(1,numBScans); funcVal = cell(1,numBScans);
		% plot positions of BScans
		for n = 1:numBScans
			hLine(n) = line([BScanHeader{n}.StartX/fileHeader.ScaleXSlo BScanHeader{n}.EndX/fileHeader.ScaleXSlo],[BScanHeader{n}.StartY/fileHeader.ScaleYSlo BScanHeader{n}.EndY/fileHeader.ScaleYSlo]);
			if n == 1
				set(hLine(n),'Color','red');
			else
				set(hLine(n),'Color','blue');
			end

			if mod(n,2)
				text(BScanHeader{n}.EndX/fileHeader.ScaleXSlo+10,BScanHeader{n}.EndY/fileHeader.ScaleYSlo,num2str(n),'Color','white');
			else
				text(BScanHeader{n}.EndX/fileHeader.ScaleXSlo+50,BScanHeader{n}.EndY/fileHeader.ScaleYSlo,num2str(n),'Color','white');
			end
		end

		% plot first B-Scan
		axes(axisBScan); reset(axisBScan);
		imagesc(sqrt(sqrt(BScanData{1}))); title(sprintf('B-Scan %d',1)); colormap gray;
		set(findobj('Tag','nextBScan'),'Enable','on');
		StringVals = cell(1,numBScans);
		for i = 1:numBScans
			StringVals{i} = num2str(i);
		end
		set(findobj('Tag','selectBScan'),'String',StringVals,'Value',1,'Visible','on');
		currentBScan = 1;
	end
end

% switch to the next BScan
function switchBScan_Callback(hObject,eventdata)
	Caller = get(hObject,'Tag');
	set(hLine(currentBScan),'Color','blue');
	if strcmp(Caller,'nextBScan')
		currentBScan = currentBScan + 1;
	elseif strcmp(Caller,'prevBScan')
		currentBScan = currentBScan - 1;
	else
		currentBScan = get(hObject,'Value');
	end
	set(hLine(currentBScan),'Color','red');
	set(findobj('Tag','selectBScan'),'Value',currentBScan);

	% plot BScan
	axes(axisBScan); reset(axisBScan);
	imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap gray; hold on;
	if prod(size(predictions{currentBScan})) > 0
		plot(collector.options.columnsPred+collector.options.clipRange(1),predictions{currentBScan});
	end

	% modify status of next and previous buttons
	set(findobj('Tag','prevBScan'),'Enable','on'); set(findobj('Tag','nextBScan'),'Enable','on'); 	
	if currentBScan == numBScans
		set(findobj('Tag','nextBScan'),'Enable','off');
	elseif currentBScan == 1
		set(findobj('Tag','prevBScan'),'Enable','off');
	end

	if prod(size(predictions{currentBScan})) > 0
		set(noQuality,'Enable','on');
	else
		set(noQuality,'Enable','off');
	end

	if prod(size(funcVal{currentBScan})) > 0
		set(showDataQuality,'Enable','on'); set(showShapeQuality,'Enable','on');
	else
		set(showDataQuality,'Enable','off'); set(showShapeQuality,'Enable','off');
	end
end

% segment BScan
function segmentBScan_Callback(hObject,eventdata)
	files.name = fileNameScan;

	sizeDiff = (fileHeader.SizeX-model.params.X)
	if sizeDiff < 0
		error('Scan width is smaller than for the trained model, aborting segmentation');
	end

	if strcmp(get(hObject,'Tag'),'segmentAllBScans')
		toSegment = 1:numBScans;
	elseif strcmp(get(hObject,'Tag'),'segmentBScan')
		toSegment = currentBScan;
	end
	collector.options = struct('numBScansPerVolume',fileHeader.NumBScans,'distBScans',fileHeader.distanceBScans,'folder_labels','','loadLabels',0,'clipRange',[1 model.params.X]+sizeDiff/2,'folder_data',pathNameScan,'loadRoutineData','spectralisVol','numRegionsPerVolume',1);
	fieldNames = fieldnames(model.params);
	for i = 1:length(fieldNames)
		if ~isfield(collector.options,fieldNames{i})
			eval(sprintf('collector.options.%s = model.params.%s;',fieldNames{i},fieldNames{i}));
		end
	end	

	testFunc.name = @predVariational; testFunc.options = struct('calcFuncVal',get(qualityBScan,'Value'));

	for i = 1:length(toSegment)
		collector.options.labelIDs = toSegment(i);
		collectorTest.name = @collectTestData; collectorTest.options = collector.options;

		updateStatus(sprintf('Segmenting BScan %d.',toSegment(i))); 
		set(f,'Pointer','watch');
		drawnow

		prediction = testFunc.name(files,collectorTest,struct(),model,testFunc.options);
		predictions{toSegment(i)} = prediction.prediction{1};
		set(noQuality,'Enable','on');

		if get(qualityBScan,'Value')
			set(showDataQuality,'Enable','on'); set(showShapeQuality,'Enable','off');			
			funcVal{toSegment(i)} = prediction.funcVal;
		end
		updateStatus(sprintf('Fnished segmenting BScan %d.',toSegment(i))); drawnow
	
		% fit plot to clip
		if toSegment(i) == currentBScan
			axes(axisBScan); reset(axisBScan);
			imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap gray;
			hold on; 
			plot(collector.options.columnsPred+collector.options.clipRange(1),predictions{currentBScan});
		end
	end
	set(f,'Pointer','arrow');
end

function showDataQuality_Callback(hObject,eventdata)
	axes(axisBScan); reset(axisBScan);
  	imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap gray; hold on;
	plot(collector.options.columnsPred+collector.options.clipRange(1),predictions{currentBScan});

	makeQualityOverlay(predictions{currentBScan},squeeze(funcVal{currentBScan}.q_c_data),margins.mu(4,:),margins.sigma(4,:),collector.options.columnsPred+collector.options.clipRange(1));
end

function showShapeQuality_Callback(hObject,eventdata)
	axes(axisBScan); reset(axisBScan);
  	imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap gray; hold on;
	plot(collector.options.columnsPred+collector.options.clipRange(1),predictions{currentBScan});

	makeQualityOverlay(predictions{currentBScan},squeeze(funcVal{currentBScan}.q_c_shape),margins.mu(3,:),margins.sigma(3,:),collector.options.columnsPred+collector.options.clipRange(1));
end

function noQuality_Callback(hObject,eventdata)
	axes(axisBScan); reset(axisBScan);
  	imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap gray; hold on;
	plot(collector.options.columnsPred+collector.options.clipRange(1),predictions{currentBScan});
end

function noPred_Callback(hObject,eventdata)
	axes(axisBScan); reset(axisBScan);
  	imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap gray;
end

% load model
function loadModel_Callback(hObject,eventdata)
	[FileName,PathName,FilterIndex] = uigetfile('*.mat');

    if FileName ~= 0
		tmp = load([PathName FileName],'model');
		model = tmp.model;
		if length(fieldnames(model)) > 0
			updateStatus('Model loaded.');
			set(findobj('Tag','segmentBScan'),'Enable','on');
			set(findobj('Tag','segmentAllBScans'),'Enable','on');
			set(findobj('Tag','qualityBScan'),'Enable','on');
		else
			updateStatus('Loading of model failed.');
		end
	end
end

function updateStatus(string)
    h = findobj('Tag','statusText'); set(h,'String',string);
	drawnow
end

end
