function octGUI()
close all;

% define global variables --> nested functions have access to these variables;
hLine = [];
currentBScan = 0; numBScans = 0; currentAppLayer = 0; numModes = 0; currentMode = 0; modesPerView = 4; shapeModel = 1; z = [];
BScanHeader = ''; BScanData = ''; BScanSeg = ''; SLO = ''; fileHeader = ''; predictions = ''; funcVal = ''; fileExt = '';
model = ''; collector = '';
fileNameScan = ''; pathNameScan = ''; pathNameGT = '';
% initialize buttons above the BScan axis
showDataQuality = ''; showShapeQuality = ''; showPred = ''; showPredInit = ''; noPred = ''; showAppearance = ''; switchAppearance = ''; showPredNew = ''; showQB = ''; showGT = ''; buttons = '';
% initialize buttons above the SLO axis
openFile = ''; selectFileType = ''; loadModel = ''; segmentBScan = ''; resegmentBScan = ''; addMode = ''; fileNameText = ''; nextBScan = ''; prevBScan = ''; selectBScan = ''; statusText = ''; selectGTDir = '';


% *****************************************************
% ********************* LAYOUT ************************
% *****************************************************

f = figure('Position', [100 100 1250 750],'Tag','mainWindow');

tgroup = uitabgroup('Parent', f);
tab1 = uitab('Parent', tgroup, 'Title', 'Viewer');
tab2 = uitab('Parent', tgroup, 'Title', 'Scan params');
tab3 = uitab('Parent', tgroup, 'Title', 'Model params');
tab4 = uitab('Parent', tgroup, 'Title', 'Shape mode viewer');

% Position of BScan axes
BScanBottom = 100; BScanLeft = 350; BScanWidth = 850; BScanHeight = 500; leftPos = BScanLeft;
% position of mode axes
ModeBottom = 310; ModeLeft = 300; ModeWidth = 650; ModeHeight = 350; subModeHeight = 150; subModeWidth = 250;

% Position of SLO axes
SLOBottom = 100; SLOLeft = 40; SLOWidth = 250; SLOHeight = 250;
% axis
axisSLO = axes('Parent',tab1,'Units','Pixels','Position',[SLOLeft,SLOBottom,SLOWidth,SLOHeight]); set(axisSLO,'YTickLabel',[],'XTickLabel',[]);
axisBScan = axes('Parent',tab1,'Units','Pixels','Position',[BScanLeft,BScanBottom,BScanWidth,BScanHeight]); set(axisBScan,'YTickLabel',[],'XTickLabel',[],'Tag','axisBScan');
% Control bar to the left
BBottom = SLOBottom + SLOHeight;

axisModeComp = axes('Parent',tab4,'Units','Pixels','Position',[ModeLeft,ModeBottom,ModeWidth,ModeHeight]); set(axisModeComp,'YTickLabel',[],'XTickLabel',[],'Tag','axisModeComp')
for i = 1:modesPerView
	posLeft = 50 + (i-1)*(subModeWidth+50);
	posBottom = 75;
	axisMode(i) = axes('Parent',tab4,'Units','Pixels','Position',[posLeft,posBottom,subModeWidth,subModeHeight]); set(axisMode(i),'YTickLabel',[],'XTickLabel',[]);
	decreaseMode(i) = uicontrol('Parent',tab4,'Style','pushbutton','String','<','Position',[posLeft+subModeWidth/2-60,posBottom-55,25,25],'Enable','on','Callback',{@changeComposition_Callback},'Tag',sprintf('DecMode%d',i));
	selectModeComp(i) = uicontrol('Parent',tab4,'Style','edit','Units','Pixels','Position',[posLeft+subModeWidth/2-30,posBottom-55,60,25],'String',0,'Callback',{@updateComposition_Callback}); 
	increaseMode(i) = uicontrol('Parent',tab4,'Style','pushbutton','String','>','Position',[posLeft+subModeWidth/2+35,posBottom-55,25,25],'Enable','on','Callback',{@changeComposition_Callback},'Tag',sprintf('IncMode%d',i));
% text field showing current Mode number
end
  
nextMode = uicontrol('Parent',tab4,'Style','pushbutton','String','>','Position',[ModeLeft+ModeWidth/2+40,ModeBottom-55,25,25],'Enable','off','Tag','nextMode','Callback',{@switchMode_Callback});
prevMode = uicontrol('Parent',tab4,'Style','pushbutton','String','<','Position',[ModeLeft+ModeWidth/2-60,ModeBottom-55,25,25],'Enable','off','Tag','prevMode','Callback',{@switchMode_Callback});
% text field showing current Mode number
selectMode = uicontrol('Parent',tab4,'Style','popupmenu','String','','Position',[ModeLeft+ModeWidth/2-20,ModeBottom-55,50,25],'Tag','selectMode','Callback',{@switchMode_Callback},'Visible','off');
resetMode = uicontrol('Parent',tab4,'Style','pushbutton','String','Reset','Position',[50,500,100,25],'Enable','on','Tag','resetMode','Callback',{@resetMode_Callback});
projMode = uicontrol('Parent',tab4,'Style','pushbutton','String','Proj. Prediction','Position',[50,460,150,25],'Enable','on','Tag','resetMode','Callback',{@projMode_Callback});


openFile = uicontrol('Parent',tab1,'Style','pushbutton','String','Open Scan','Position',[60,BBottom + 20,75,25],'Callback',{@openFile_Callback},'Tag','openFile');
selectFileType = uicontrol('Parent',tab1,'Style','popupmenu','String',{'.vol','.mat'},'Position',[145,BBottom + 20,50,25],'Tag','selectFileType');
selectGTDir = uicontrol('Parent',tab1,'Style','pushbutton','String','Select GT Folder','Position',[60,BBottom + 60,140,25],'Callback',{@selectGTDir_Callback});
loadModel = uicontrol('Parent',tab1,'Style','pushbutton','String','Load Model','Position',[60,BBottom+100,100,25],'Enable','off','Tag','loadModel','Callback',{@loadModel_Callback});
% button to segment the currently display BScan
segmentBScan = uicontrol('Parent',tab1,'Style','pushbutton','String','Segment','Position',[60,BBottom+140,80,25],'Enable','off','Tag','segmentBScan','Callback',{@segmentBScan_Callback});

%resegmentBScan = uicontrol('Parent',tab1,'Style','pushbutton','String','Resegment','Position',[60,BBottom + 140,100,25],'Enable','off','Tag','resegmentBScan','Callback',{@resegmentBScan_Callback});
%addMode = uicontrol('Parent',tab1,'Style','pushbutton','String','Add Mode','Position',[60,BBottom + 180,100,25],'Enable','off','Callback',{@addMode_Callback});

fileNameText = uicontrol('Parent',tab1,'Style','text','String','No file loaded...','Position',[60 BBottom+270,250,25],'Tag','fileNameText','HorizontalAlignment','Left');
modelNameText = uicontrol('Parent',tab1,'Style','text','String','No model loaded...','Position',[60 BBottom+240,250,25],'Tag','modelNameText','HorizontalAlignment','Left');
unsignedErrorText = uicontrol('Parent',tab1,'Style','text','String','Error:','Position',[60 BBottom+210,250,25],'Tag','unsignedErrorText','HorizontalAlignment','Left');

%segmentAllBScans = uicontrol('Parent',tab1,'Style','pushbutton','String','Segment All','Position',[60,BBottom + 180,100,25],'Enable','off','Visible','off','Tag','segmentAllBScans','Callback',{@segmentBScan_Callback});

% buttons to circle through B-Scans in a volume
nextBScan = uicontrol('Parent',tab1,'Style','pushbutton','String','>','Position',[BScanLeft+BScanWidth/2+50,BScanBottom-75,25,25],'Enable','off','Tag','nextBScan','Callback',{@switchBScan_Callback});
prevBScan = uicontrol('Parent',tab1,'Style','pushbutton','String','<','Position',[BScanLeft+BScanWidth/2-50,BScanBottom-75,25,25],'Enable','off','Tag','prevBScan','Callback',{@switchBScan_Callback});
% text field showing current BScan number
selectBScan = uicontrol('Parent',tab1,'Style','popupmenu','String','','Position',[BScanLeft+BScanWidth/2-10,BScanBottom-75,50,25],'Tag','selectBScan','Callback',{@switchBScan_Callback},'Visible','off');

% text field that shows status reports
statusText = uicontrol('Parent',tab1,'Style','text','String','Waiting for action...','Position',[40,10,500,25],'Tag','statusText','HorizontalAlignment','Left');

% buttons to show various segmentation results for the current BScan
buttons = {{'showDataQuality','DataTerm',75},{'showShapeQuality','ShapeTerm',100},{'showPred','Pred',75},{'showPredNew','PredNew',20,'checkbox'},{'showPredInit','PredInit',100},{'showQB','QB',50},{'noPred','BScan',85},{'showGT','GT',50},{'showAppearance','AppTerms',85}};

for i = 1:length(buttons)
    if length(buttons{i}) == 3
        eval(sprintf('%s = uicontrol(''Parent'',tab1,''Style'',''pushbutton'',''String'',''%s'',''Position'',[%d BScanBottom+BScanHeight+10 %d 25],''Enable'',''off'',''Callback'',{@%s_Callback});',buttons{i}{1},buttons{i}{2},leftPos,buttons{i}{3},buttons{i}{1}));
    else
        eval(sprintf('%s = uicontrol(''Parent'',tab1,''Style'',''checkbox'',''Position'',[%d BScanBottom+BScanHeight+10 %d 25],''Enable'',''off'',''Value'',0);',buttons{i}{1},leftPos,buttons{i}{3}));
    end
    leftPos = leftPos + 10 + buttons{i}{3};
end
switchAppearance = uicontrol('Parent',tab1,'Style','popupmenu','String',cellstr(num2str((1:19)')),'Position',[leftPos,BScanBottom+BScanHeight+8,50,25],'Enable','off','Callback',{@switchAppearance_Callback});
handlesVec = cellfun(@(c) c(1), buttons); handlesVec{end+1} = 'switchAppearance';



% *****************************************************
% ****************** CALLBACKS ************************
% *****************************************************

function selectGTDir_Callback(hObject,eventdata)
	if length(pathNameScan) > 0
		folderGT = uigetdir(pathNameScan);
	else
		folderGT = uigetdir();
	end

	if folderGT ~=0
		pathNameGT = folderGT;
		loadGT();
	end
end

function loadGT()
	if length(pathNameGT) > 0 && length(fileNameScan) > 0
		[pathstr,name,ext] = fileparts(fileNameScan);
		files = dir([pathNameGT filesep name '*']);
		if length(files) > 0
			for i = 1:length(files)
				interpolation = load([pathNameGT filesep files(i).name],'interpolation');
				% select labelID from filename
				[mat tok] = regexp(files(i).name,'_(\d{1,})_coordinates','match','tokens');
				BScanSeg{str2num(tok{1}{1})+1} = interpolation.interpolation;
			end
		end
	end
end

function openFile_Callback(hObject,eventdata,handles)
	h = findobj('Tag','selectFileType');
	String = get(h,'String');
	Value = get(h,'Value');
	if length(pathNameScan) == 0
		pathNameScan = '/home/fabian/Documents/DatasetsWork/drusendaten';
	end

	if length(pathNameScan) > 0
		[FileName,PathName,FilterIndex] = uigetfile(String{Value},'Load Scan',pathNameScan);
	else
		[FileName,PathName,FilterIndex] = uigetfile(String{Value},'Load Scan');
	end

	if FileName ~= 0
		fileNameScan = FileName; pathNameScan = PathName;
		set(fileNameText,'String',FileName);
		[pathstr,name,ext] = fileparts(fileNameScan); fileExt = ext;
		if strcmp(ext,'.vol')
			[BScanData, fileHeader, BScanHeader, BScanSeg, SLO] = HDEVolImporter(PathName,FileName,struct('plotSLO',0,'plotBScans',0));
		elseif strcmp(ext,'.mat')
			[BScanData, fileHeader, BScanHeader, BScanSeg, SLO] = HDEMatImporter(PathName,FileName,struct('plotSLO',0,'plotBScans',0));
		end
		tBScan = printDataInTable(tab2,fileHeader);
		tBScan.Position(3) = tBScan.Extent(3); tBScan.Position(4) = tBScan.Extent(4);           
		numBScans = fileHeader.NumBScans;
		updateStatus(sprintf('%s loaded.',FileName));
		set(findobj('Tag','loadModel'),'Enable','on');

		loadGT();

		% mirrow BScans for Volumes
		if findstr(fileHeader.ScanPosition,'OS') & numBScans > 1
			for i = 1:numBScans
				BScanData{i} = BScanData{i}(:,end:-1:1);
				BScanSeg{i} = BScanSeg{i}(:,end:-1:1);
			end
		end

		% load variables for evaluating the error term
		if numBScans > 1
			margins = load([getenv('OCT_CODE_DIR') 'datafiles/healthyMargins3D.mat']);
		else
			margins = load([getenv('OCT_CODE_DIR') 'datafiles/healthyMargins2D.mat']);
		end



		% plot SLO scan
		axes(axisSLO); reset(axisSLO); imagesc(sqrt(sqrt(SLO))); colormap gray; hold on;
		predictions = cell(1,numBScans);
		% plot positions of BScans for volumes
		if numBScans > 1
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
		end

		setButtonProps(handlesVec,{'Enable','off'});	

		% plot first B-Scan
		axes(axisBScan); reset(axisBScan);
		imagesc(sqrt(sqrt(BScanData{1}))); colormap(axisBScan,'gray');
		set(findobj('Tag','nextBScan'),'Enable','off');	set(findobj('Tag','prevBScan'),'Enable','off');
		if numBScans > 1
			set(findobj('Tag','nextBScan'),'Enable','on');
		end
		StringVals = cell(1,numBScans);
		for i = 1:numBScans
			StringVals{i} = num2str(i);
		end
		set(findobj('Tag','selectBScan'),'String',StringVals,'Value',1,'Visible','on');
		set(noPred,'Enable','on'); set(showGT,'Enable','on');
		currentBScan = 1;
	end
end

% switch to the next BScan
function switchBScan_Callback(hObject,eventdata)
	Caller = get(hObject,'Tag');
	if numBScans > 1 
		set(hLine(currentBScan),'Color','blue');
	end
	if strcmp(Caller,'nextBScan')
		currentBScan = currentBScan + 1;
	elseif strcmp(Caller,'prevBScan')
		currentBScan = currentBScan - 1;
	else
		currentBScan = get(hObject,'Value');
	end
	if numBScans > 1
		set(hLine(currentBScan),'Color','red');
	end
	set(findobj('Tag','selectBScan'),'Value',currentBScan);

	% plot BScan
	axes(axisBScan); reset(axisBScan);
	imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap(axisBScan,'gray'); hold on;
	if isfield(predictions{currentBScan},'prediction') > 0
		plot(collector.options.columnsPred+collector.options.clipRange(1)-1,predictions{currentBScan}.prediction{1});
	end

	% modify status of next and previous buttons
	set(findobj('Tag','prevBScan'),'Enable','on'); set(findobj('Tag','nextBScan'),'Enable','on'); 	
	if currentBScan == numBScans
		set(findobj('Tag','nextBScan'),'Enable','off');
	end
	if currentBScan == 1
		set(findobj('Tag','prevBScan'),'Enable','off');
	end

	set(showPredNew,'Value',0); set(showPred,'Enable','off'); set(showPredNew,'Enable','off'); set(resegmentBScan,'Enable','off'); set(showQB,'Enable','off'); set(showPredInit,'Enable','off'); set(showAppearance,'Enable','off'); set(addMode,'Enable','off');
	if isfield(predictions{currentBScan},'prediction')
		set(showPred,'Enable','on'); set(resegmentBScan,'Enable','on'); set(showQB,'Enable','on'); set(showPredInit,'Enable','on');
		if length(predictions{currentBScan}.prediction) > 1
			set(showPredNew,'Enable','on'); set(addMode,'Enable','on');
		end
	end

	if isfield(predictions{currentBScan},'funcVal')
		set(showDataQuality,'Enable','on'); %set(showShapeQuality,'Enable','on');
	else
		set(showDataQuality,'Enable','off'); set(showShapeQuality,'Enable','off');
	end
end

function showQB_Callback(hObject,eventdata)
	displaySegmentation(predictions{currentBScan}.q_b,collector.options.columnsPred+collector.options.clipRange(1)-1);
end

function showGT_Callback(hObject,eventdata)
	displaySegmentation(BScanSeg{currentBScan},1:fileHeader.SizeX);
end

function addMode_Callback(hObject,eventdata)
	% add the difference between the new segmentation and the
	diff = predictions{currentBScan}.prediction{1}-predictions{currentBScan}.prediction{2};
	for i = 1:9
    	if sum(diff(i,:)~=0) > 5
        	PP = splinefit(1:250,diff(i,:),round(sum(diff(i,:)~=0)/2));
	        diffSmooth(i,:) = ppval(PP,1:250);
    	else
	        diffSmooth(i,:) = zeros(1,250);
	    end
	end
	dist = (currentBScan - (floor(collector.options.numBScansPerVolume/2)+1))*collector.options.distBScans;
    % find closest B-scan in the shape model
    [~,scanID] = min(abs(model.params.BScanPositions-dist));
	
	model.shapeModel(scanID).WML = [model.shapeModel(scanID).WML reshape(diffSmooth',[],1)*0.5];
	updateStatus(sprintf('Added mode to %d''th shape model.',scanID)); drawnow
	predictions{currentBScan}.prediction(2) = []; predictions{currentBScan}.funcVal(2) = []; set(addMode,'Enable','off'); set(showPredNew,'Value',0); set(showPredNew,'Enable','off');
	set(hObject,'Tag','segmentBScan');
	segmentBScan_Callback(hObject,eventdata);
end

% segment BScan
function segmentBScan_Callback(hObject,eventdata)
	if strcmp(get(hObject,'Tag'),'segmentAllBScans')
		toSegment = 1:numBScans;
	elseif strcmp(get(hObject,'Tag'),'segmentBScan')
		toSegment = currentBScan;
	end
	pred = startSegmentation(struct(),struct(),toSegment);
    
	for i = 1:length(pred)
	    predictions{toSegment(i)}.prediction{1} = pred{i}.prediction;
        predictions{toSegment(i)}.predictionInit = pred{i}.predictionInit;
        predictions{toSegment(i)}.appearance = pred{i}.appearance;
		predictions{toSegment(i)}.q_b = pred{i}.q_b;
        predictions{toSegment(i)}.funcVal{1} = pred{i}.funcVal;
		if size(BScanSeg{toSegment(i)},1) == size(predictions{toSegment(i)}.prediction{1},1)
			predictions{toSegment(i)}.unsignedError = mean(mean(abs(predictions{toSegment(i)}.prediction{1}-BScanSeg{toSegment(i)}(:,collector.options.columnsPred+collector.options.clipRange(1) - 1))));
		end
		% if the currently displayed BSCan was segmented, display the prediction result
		if toSegment(i) == currentBScan
			set(showPred,'Enable','on'); set(showPredInit,'Enable','on'); set(resegmentBScan,'Enable','on'); set(showQB,'Enable','on');
			set(showAppearance,'Enable','on'); currentAppLayer = 1; set(switchAppearance,'Enable','on');
			displaySegmentation(predictions{currentBScan}.prediction{1},collector.options.columnsPred+collector.options.clipRange(1)-1);
			set(showDataQuality,'Enable','on'); set(showShapeQuality,'Enable','off');
			% display unsigned error
			if isfield(predictions{currentBScan},'unsignedError')
				set(unsignedErrorText,'String',sprintf('Error: %.2f px | %.2f mum',predictions{currentBScan}.unsignedError,predictions{currentBScan}.unsignedError*3.87));
			end
		end
	end
end

function resegmentBScan_Callback(hObject,eventdata)
	toSegment = currentBScan;

	q_c_data = squeeze(predictions{toSegment}.funcVal{1}.q_c_data);
	idxHealthy =  (q_c_data < repmat(margins.mu(4,:),250,1) + repmat(margins.sigma(4,:)*5,250,1))';
	% recalculate all boundaries 6-9
	idxHealthy(6:9,sum(~idxHealthy(6:9,:))>0) = 0;
	optionsAdd.segmentation = predictions{toSegment}.prediction{1};
	optionsAdd.idxRecalc = ~idxHealthy;
	optionsAdd.onlyInitialize = 1;
	optionsAdd.appearance = predictions{toSegment}.appearance;

	pred = startSegmentation(struct(),optionsAdd,toSegment);

	predictions{toSegment}.prediction{2} = predictions{toSegment}.prediction{1};
	predictions{toSegment}.prediction{2}(~idxHealthy) = pred{1}.predictionInit(~idxHealthy);
	predictions{toSegment}.funcVal{2} = pred{1}.funcVal;
	set(showPredNew,'Enable','on','Value',1); set(addMode,'Enable','on');
	displaySegmentation(predictions{toSegment}.prediction{2},collector.options.columnsPred+collector.options.clipRange(1)-1);
end

function pred = startSegmentation(collectorAdd,optionsAdd,toSegment)
	files.name = fileNameScan;
	sizeDiff = (int32(fileHeader.SizeX)-model.params.X);
	if sizeDiff < 0
		updateStatus('Scan width is smaller than for the trained model, aborting segmentation');
		return
	end

	collector.options = struct('numBScansPerVolume',fileHeader.NumBScans,'folder_labels','','loadLabels',0,'clipRange',double([1 int32(model.params.X)]+sizeDiff/2),'folder_data',pathNameScan,'numRegionsPerVolume',1,'saveAppearanceTerms',1,'printTimings',1);
	if isfield(fileHeader,'distanceBScans')
		collector.options.distBScans = fileHeader.distanceBScans;
	end
	if numBScans > 1
		collector.options.mirrorBScan = 'OS';
	else
		collector.options.mirrorBScan = '';
	end
	collector.options.loadRoutineData = ['spectralis' upper(fileExt(2)) fileExt(3:4)];
	testFunc.name = @predVariational; testFunc.options = struct('calcFuncVal',1);

	% add params stored in the model file, in case they are not defined yet
	fieldNames = fieldnames(model.params);
	for i = 1:length(fieldNames)
		if ~isfield(collector.options,fieldNames{i})
			eval(sprintf('collector.options.%s = model.params.%s;',fieldNames{i},fieldNames{i}));
		end
	end
	% add all options that are added by the calling function
	fieldNames = fieldnames(collectorAdd);
	for i = 1:length(fieldNames)
		eval(sprintf('collector.options.%s = collectorAdd.%s;',fieldNames{i},fieldNames{i}));
	end
	fieldNames = fieldnames(optionsAdd);
	for i = 1:length(fieldNames)
		eval(sprintf('testFunc.options.%s = optionsAdd.%s;',fieldNames{i},fieldNames{i}));
	end
		
	for i = 1:length(toSegment)
		collector.options.labelIDs = toSegment(i);
		collectorTest.name = @collectTestData; collectorTest.options = collector.options;

		updateStatus(sprintf('Segmenting BScan %d.',toSegment(i))); 
		set(f,'Pointer','watch');
		drawnow

		prediction = testFunc.name(files,collectorTest,struct(),model,testFunc.options);
		pred{i}.prediction = prediction.prediction{1};
		pred{i}.predictionInit = prediction.prediction_init{1};
		pred{i}.appearance = prediction.appearanceTerms.prediction;
		pred{i}.q_b = prediction.q_b{1};
		set(showDataQuality,'Enable','on'); set(showShapeQuality,'Enable','off');			
		pred{i}.funcVal = prediction.funcVal;
		updateStatus(sprintf('Fnished segmenting BScan %d.',toSegment(i))); drawnow
	end
	set(f,'Pointer','arrow');
end

function showAppearance_Callback(hObject,eventdata)
	displayAppearance();
end

function switchAppearance_Callback(hObject,eventdata)
	currentAppLayer = get(hObject,'Value'); 
	displayAppearance();
end

function displayAppearance()
	axes(axisBScan); reset(axisBScan); cla(axisBScan);
	h1 = imagesc(sqrt(sqrt(BScanData{currentBScan}(:,collector.options.columnsPred+collector.options.clipRange(1)-1)))); colormap(axisBScan,'gray'); hold on;
	app = reshape(predictions{currentBScan}.appearance{1}(currentAppLayer,:),fileHeader.SizeZ,[]);
	% normalize the maximum of each column to one
	app = app./repmat(max(app),fileHeader.SizeZ,1);
	h2 = imagesc(app);
	set(h2,'AlphaData',app);
	colormap(axisBScan,[gray(64); jet(64)])
	set(h1,'CData',get(h1,'CData')/2)
	set(h2,'CData',(get(h2,'CData')+1)/2)
end

function h = displaySegmentation(prediction,columns)
	axes(axisBScan); reset(axisBScan);
	imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap(axisBScan,'gray'); hold on;
	h = plot(columns,prediction);
end

function showDataQuality_Callback(hObject,eventdata)
	h = displaySegmentation(predictions{currentBScan}.prediction{get(showPredNew,'Value')+1},collector.options.columnsPred+collector.options.clipRange(1)-1);
	set(h(:),'Color',[0.8 0.8 1]);
	makeQualityOverlay(predictions{currentBScan}.prediction{get(showPredNew,'Value')+1},squeeze(predictions{currentBScan}.funcVal{get(showPredNew,'Value')+1}.q_c_data),margins.mu(4,:),margins.sigma(4,:),collector.options.columnsPred+collector.options.clipRange(1)-1);
end

function showShapeQuality_Callback(hObject,eventdata)
	h = displaySegmentation(predictions{currentBScan}.prediction{get(showPredNew,'Value')+1},collector.options.columnsPred+collector.options.clipRange(1)-1);
	set(h(:),'Color',[0.8 0.8 1]);
	makeQualityOverlay(predictions{currentBScan}.prediction{get(showPredNew,'Value')+1},squeeze(predictions{currentBScan}.funcVal{get(showPredNew,'Value')+1}.q_c_shape),margins.mu(3,:),margins.sigma(3,:),collector.options.columnsPred+collector.options.clipRange(1)-1);
end

function showPredInit_Callback(hObject,eventdata)
	displaySegmentation(predictions{currentBScan}.predictionInit,collector.options.columnsPred+collector.options.clipRange(1)-1);
end
function showPred_Callback(hObject,eventdata)
	displaySegmentation(predictions{currentBScan}.prediction{get(showPredNew,'Value')+1},collector.options.columnsPred+collector.options.clipRange(1)-1);
end

function noPred_Callback(hObject,eventdata)
	axes(axisBScan); reset(axisBScan);
  	imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap(axisBScan,'gray');
end

% load model
function loadModel_Callback(hObject,eventdata)
	[FileName,PathName,FilterIndex] = uigetfile('*.mat','Load Model','/home/fabian/Documents/Arbeit/Code/MyCode/OCT/datafiles');

    if FileName ~= 0
		tmp = load([PathName FileName],'model');
		model = tmp.model;
		tBScan = printDataInTable(tab3,model.params);
		tBScan.Position(3) = tBScan.Extent(3); tBScan.Position(4) = tBScan.Extent(4);           
		if length(fieldnames(model)) > 0
			updateStatus('Model loaded.');
			set(findobj('Tag','segmentBScan'),'Enable','on');
			set(findobj('Tag','segmentAllBScans'),'Enable','on');

			% show first mode in the respective tab
			numModes = size(model.shapeModel(1).WML,2); z = zeros(numModes,1);
			currentMode = 1; 
			stringVals = cellstr(num2str((1:ceil(numModes/modesPerView))'));
			set(selectMode,'String',stringVals,'Value',1,'Visible','on');
			set(nextMode,'Enable','on');
			set(modelNameText,'String',FileName);
			switchMode(currentMode);
			updateComposition_Callback();
		else
			updateStatus('Loading of model failed.');
		end
	end
end

function switchMode_Callback(hObject,eventdata)
	Caller = get(hObject,'Tag');
 	if strcmp(Caller,'nextMode')
        currentMode = currentMode + 1;
    elseif strcmp(Caller,'prevMode')
        currentMode = currentMode - 1;
    else
        currentMode = get(hObject,'Value');
    end
	set(selectMode,'Value',currentMode);
	switchMode(currentMode);	
end

function switchMode(currentMode)
	startMode = (currentMode-1)*modesPerView;
	for i = 1:modesPerView
		axes(axisMode(i)); reset(axisMode(i)); set(selectModeComp(i),'String',0)
		if i+startMode <= numModes
			plot(reshape(model.shapeModel(shapeModel).WML(:,i+startMode),length(model.shapeModel(shapeModel).columnsShape),[]));
			set(selectModeComp(i),'String',sprintf('%.2f',z(i+startMode)));
		end
	end
	set(nextMode,'Enable','on'); set(prevMode,'Enable','on');
	if currentMode == ceil(numModes/modesPerView)
		set(nextMode,'Enable','off');
	end
	if currentMode == 1
		set(prevMode,'Enable','off');
	end

end

function changeComposition_Callback(hObject,eventdata)
	startMode = (currentMode-1)*modesPerView;
	id = regexp(get(hObject,'Tag'),'\d','match');
	id = str2num(id{1});
	if strfind(get(hObject,'Tag'),'DecMode')
		z(id+startMode) = z(id+startMode) - 1;
	else
		z(id+startMode) = z(id+startMode) + 1;
	end
	set(selectModeComp(id),'String',z(id+startMode));
	updateComposition_Callback();
end

function updateComposition_Callback(hObject,eventdata)
	startMode = (currentMode-1)*modesPerView;
	for i = 1:modesPerView
		z(i+startMode) = str2num(get(selectModeComp(i),'String'));
	end	
	% calculate the composition based on the current z
	comp = model.shapeModel(shapeModel).mu + model.shapeModel.WML*z;
	axes(axisModeComp); reset(axisModeComp);
	imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap(axisBScan,'gray'); hold on;
	plot(collector.options.clipRange(1) - 1 + collector.options.columnsShape{1},reshape(comp,length(model.shapeModel(shapeModel).columnsShape),[])); set(gca,'YDir','reverse');
%	plot(reshape(comp,length(model.shapeModel(shapeModel).columnsShape),[])); set(gca,'YDir','reverse');
end

function resetMode_Callback(hObject,eventdata)
	z = zeros(numModes,1);
	set(selectModeComp,'String',0);
	updateComposition_Callback;
end


function projMode_Callback(hObject,eventdata)
	if isfield(predictions{currentBScan},'prediction')
		M = model.shapeModel(shapeModel).WML'*model.shapeModel(shapeModel).WML + eye(numModes)*model.shapeModel(shapeModel).sigmaML;
		z = inv(M)*(model.shapeModel(shapeModel).WML'*(reshape(predictions{currentBScan}.prediction{1}',[],1)-model.shapeModel(shapeModel).mu));
		for i = 1:modesPerView
			set(selectModeComp(i),'String',sprintf('%.2f',z(i + (currentMode-1)*modesPerView)));
		end
		comp = model.shapeModel(shapeModel).mu + model.shapeModel.WML*z;
		axes(axisModeComp); reset(axisModeComp);
		imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap(axisBScan,'gray'); hold on;
		plot(collector.options.clipRange(1) - 1 + collector.options.columnsShape{1},reshape(comp,length(model.shapeModel(shapeModel).columnsShape),[])); set(gca,'YDir','reverse');
	end
end

function setButtonProps(buttons,props)
for i = 1:size(props,1)
    for j = 1:length(buttons)
        eval(sprintf('set(%s,''%s'',''%s'');',buttons{j},props{i,1},props{i,2}));
    end
end

end

function updateStatus(string)
    h = findobj('Tag','statusText'); set(h,'String',string);
	drawnow
end

end
