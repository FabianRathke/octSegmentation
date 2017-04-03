function octGUI()
if length(getenv('OCT_CODE_DIR')) == 0
	error('Set system variable "OCT_CODE_DIR" to the toolbox directory via setenv');
end
close all;

% define global variables --> nested functions have access to these variables;
hLine = [];
currentBScan = 0; numBScans = 0; currentAppLayer = 0; numModes = 0; currentMode = 0; modesPerView = 4; shapeModel = 1; z = []; columnsToPlot = []; numBounds = 0;
BScanHeader = ''; BScanData = ''; BScanSeg = ''; SLO = ''; fileHeader = ''; predictions = ''; funcVal = ''; fileExt = '';
model = ''; collector = ''; margins = '';
fileNameScan = ''; pathNameScan = ''; pathNameGT = '';
% initialize buttons above the BScan axis
showDataQuality = ''; showDataInit = ''; showEntropy = ''; showPred = ''; showPredInit = ''; noPred = ''; showAppearance = ''; switchAppearance = ''; showQB = ''; showGT = ''; buttons = '';
showMask = ''; addToMask = ''; removeFromMask = ''; resetMask = '';
% initialize buttons above the SLO axis
openFile = ''; selectFileType = ''; loadModel = ''; segmentBScan = ''; fileNameText = ''; nextBScan = ''; prevBScan = ''; selectBScan = ''; statusText = ''; selectGTDir = '';
% segmentation mask
mask = '';
plotCounter = 1;

% *****************************************************
% ********************* LAYOUT ************************
% *****************************************************

f = figure('Position', [100 100 1250 750],'Tag','mainWindow');

tgroup = uitabgroup('Parent', f);
tab1 = uitab('Parent', tgroup, 'Title', 'Viewer');
tab2 = uitab('Parent', tgroup, 'Title', 'Scan params');
tab3 = uitab('Parent', tgroup, 'Title', 'Model params');
tab4 = uitab('Parent', tgroup, 'Title', 'Shape mode viewer');
tab5 = uitab('Parent', tgroup, 'Title', 'BScan filter');

% Position of BScan axes
BScanBottom = 100; BScanLeft = 350; BScanWidth = 860; BScanHeight = 500; leftPos = BScanLeft;

% Position of SLO axes
SLOBottom = 100; SLOLeft = 40; SLOWidth = 250; SLOHeight = 250;
% axis
axisSLO = axes('Parent',tab1,'Units','Pixels','Position',[SLOLeft,SLOBottom,SLOWidth,SLOHeight]); set(axisSLO,'YTickLabel',[],'XTickLabel',[]);
axisBScan = axes('Parent',tab1,'Units','Pixels','Position',[BScanLeft,BScanBottom,BScanWidth,BScanHeight]); set(axisBScan,'YTickLabel',[],'XTickLabel',[],'Tag','axisBScan');
% Control bar to the left
BBottom = SLOBottom + SLOHeight;

% ##################### Shape Mode Viewer ####################
% initialize buttons
projMode = ''; prodQBMode = ''; selectMode = ''; resetMode = '';

% position of mode axes
ModeBottom = 300; ModeLeft = 300; ModeWidth = 800; ModeHeight = 400; subModeHeight = 150; subModeWidth = 250;

axisModeComp = axes('Parent',tab4,'Units','Pixels','Position',[ModeLeft,ModeBottom,ModeWidth,ModeHeight]); set(axisModeComp,'YTickLabel',[],'XTickLabel',[],'Tag','axisModeComp')
for i = 1:modesPerView
	posLeft = 50 + (i-1)*(subModeWidth+50);
	posBottom = 70;
	axisMode(i) = axes('Parent',tab4,'Units','Pixels','Position',[posLeft,posBottom,subModeWidth,subModeHeight]); set(axisMode(i),'YTickLabel',[],'XTickLabel',[]);
	decreaseMode(i) = uicontrol('Parent',tab4,'Style','pushbutton','String','<','Position',[posLeft+subModeWidth/2-60,posBottom-55,25,25],'Enable','on','Callback',{@changeComposition_Callback},'Tag',sprintf('DecMode%d',i));
	selectModeComp(i) = uicontrol('Parent',tab4,'Style','edit','Units','Pixels','Position',[posLeft+subModeWidth/2-30,posBottom-55,60,25],'String',0,'Callback',{@updateComposition_Callback}); 
	increaseMode(i) = uicontrol('Parent',tab4,'Style','pushbutton','String','>','Position',[posLeft+subModeWidth/2+35,posBottom-55,25,25],'Enable','on','Callback',{@changeComposition_Callback},'Tag',sprintf('IncMode%d',i));
% text field showing current Mode number
end
  
nextMode = uicontrol('Parent',tab4,'Style','pushbutton','String','>','Position',[ModeLeft+ModeWidth/2+40,ModeBottom-55,25,25],'Enable','off','Tag','nextMode','Callback',{@switchMode_Callback});
prevMode = uicontrol('Parent',tab4,'Style','pushbutton','String','<','Position',[ModeLeft+ModeWidth/2-60,ModeBottom-55,25,25],'Enable','off','Tag','prevMode','Callback',{@switchMode_Callback});
% text field showing current Mode number
selectMode = uicontrol('Parent',tab4,'Style','popupmenu','String','','Position',[ModeLeft+ModeWidth/2-20,ModeBottom-55,50,25],'Callback',{@switchMode_Callback},'Visible','off');
resetMode = uicontrol('Parent',tab4,'Style','pushbutton','String','Reset','Position',[50,500,100,25],'Enable','off','Callback',{@resetMode_Callback});
projMode = uicontrol('Parent',tab4,'Style','pushbutton','String','Proj. Pred','Position',[50,460,120,25],'Enable','off','Tag','prediction','Callback',{@projMode_Callback});
projInitMode = uicontrol('Parent',tab4,'Style','pushbutton','String','Proj. PredInit','Position',[50,420,120,25],'Enable','off','Tag','predictionInit','Callback',{@projMode_Callback});
projQBMode = uicontrol('Parent',tab4,'Style','pushbutton','String','Proj. QB','Position',[50,380,120,25],'Enable','off','Tag','q_b','Callback',{@projMode_Callback});
zGrad = uicontrol('Parent',tab4,'Style','pushbutton','String','Gradient z','Position',[50,340,100,25],'Enable','off','Callback',{@zGrad_Callback});
normZText = uicontrol('Parent',tab4,'Style','text','String','||z|| = 0','Position',[50,530,100,25],'HorizontalAlignment','Left');
% ##############################################################

% ##################### LEFT PANEL CONTROLS #####################
openFile = uicontrol('Parent',tab1,'Style','pushbutton','String','Open Scan','Position',[60,BBottom + 20,75,25],'Callback',{@openFile_Callback},'Tag','openFile');
selectFileType = uicontrol('Parent',tab1,'Style','popupmenu','String',{'.vol','.mat'},'Position',[145,BBottom + 20,50,25],'Tag','selectFileType');
selectGTDir = uicontrol('Parent',tab1,'Style','pushbutton','String','Select GT Folder','Position',[60,BBottom + 60,140,25],'Callback',{@selectGTDir_Callback});
loadModel = uicontrol('Parent',tab1,'Style','pushbutton','String','Load Model','Position',[60,BBottom+100,100,25],'Enable','off','Tag','loadModel','Callback',{@loadModel_Callback});
% button to segment the currently display BScan
segmentBScan = uicontrol('Parent',tab1,'Style','pushbutton','String','Segment','Position',[60,BBottom+140,80,25],'Enable','off','Tag','segmentBScan','Callback',{@segmentBScan_Callback});
% set the alpha value
setAlpha = uicontrol('Parent',tab1,'Style','edit','String','0.1','Position',[150,BBottom+140,50,25],'Enable','off','Tag','setAlpha');
% set factor for the variance
setVariance = uicontrol('Parent',tab1,'Style','edit','String','1','Position',[210,BBottom+140,50,25],'Enable','off','Tag','setVariance');

% plot
printFigure = uicontrol('Parent',tab1,'Style','pushbutton','String','printFig','Position',[60,BBottom+180,100,25],'Enable','on','Callback',{@printFigure_Callback});

fileNameText = uicontrol('Parent',tab1,'Style','text','String','No file loaded...','Position',[60 BBottom+270,250,25],'Tag','fileNameText','HorizontalAlignment','Left');
modelNameText = uicontrol('Parent',tab1,'Style','text','String','No model loaded...','Position',[60 BBottom+240,250,25],'Tag','modelNameText','HorizontalAlignment','Left');
unsignedErrorText = uicontrol('Parent',tab1,'Style','text','String','Error:','Position',[60 BBottom+210,250,25],'Tag','unsignedErrorText','HorizontalAlignment','Left');
% ##############################################################

%segmentAllBScans = uicontrol('Parent',tab1,'Style','pushbutton','String','Segment All','Position',[60,BBottom + 180,100,25],'Enable','off','Visible','off','Tag','segmentAllBScans','Callback',{@segmentBScan_Callback});

% buttons to circle through B-Scans in a volume
nextBScan = uicontrol('Parent',tab1,'Style','pushbutton','String','>','Position',[BScanLeft+BScanWidth/2+50,BScanBottom-75,25,25],'Enable','off','Tag','nextBScan','Callback',{@switchBScan_Callback});
prevBScan = uicontrol('Parent',tab1,'Style','pushbutton','String','<','Position',[BScanLeft+BScanWidth/2-50,BScanBottom-75,25,25],'Enable','off','Tag','prevBScan','Callback',{@switchBScan_Callback});
% text field showing current BScan number
selectBScan = uicontrol('Parent',tab1,'Style','popupmenu','String','','Position',[BScanLeft+BScanWidth/2-10,BScanBottom-75,50,25],'Tag','selectBScan','Callback',{@switchBScan_Callback},'Visible','off');

% text field that shows status reports
statusText = uicontrol('Parent',tab1,'Style','text','String','Waiting for action...','Position',[40,10,600,25],'Tag','statusText','HorizontalAlignment','Left');

% buttons to show various segmentation results for the current BScan
buttons = {{'showDataQuality','DataTerm',75},{'showDataInit','DataTermInit',100},{'showEntropy','EntropyTerm',100},{'showPred','Pred',75},{'showPredInit','PredInit',100},{'showQB','QB',50},{'noPred','BScan',85},{'showGT','GT',50},{'showAppearance','AppTerms',85}};

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

% ############### MASK CONTROL #######################
leftPos = BScanLeft;
showMask = uicontrol('Parent',tab1,'Style','pushbutton','String','ShowMask','Position',[leftPos,BScanBottom+BScanHeight+50,75,25],'Enable','off','Callback',{@showMask_Callback}); handlesVec{end+1} = 'showMask';
addToMask = uicontrol('Parent',tab1,'Style','pushbutton','String','+','Position',[leftPos+85,BScanBottom+BScanHeight+50,40,25],'Enable','off','Callback',{@addToMask_Callback}); handlesVec{end+1} = 'addToMask';
removeFromMask = uicontrol('Parent',tab1,'Style','pushbutton','String','-','Position',[leftPos+130,BScanBottom+BScanHeight+50,40,25],'Enable','off','Callback',{@removeFromMask_Callback}); handlesVec{end+1} = 'removeFromMask';
resetMask = uicontrol('Parent',tab1,'Style','pushbutton','String','ResetMask','Position',[leftPos+185,BScanBottom+BScanHeight+50,80,25],'Enable','off','Callback',{@resetMask_Callback}); handlesVec{end+1} = 'resetMask';
autoMask = uicontrol('Parent',tab1,'Style','pushbutton','String','AutoMask','Position',[leftPos+275,BScanBottom+BScanHeight+50,80,25],'Enable','off','Callback',{@autoMask_Callback}); handlesVec{end+1} = 'autoMask';
% set left and right boundary for segmentation
setLeftBound = uicontrol('Parent',tab1,'Style','edit','String','','Position',[leftPos+365,BScanBottom+BScanHeight+50,50,25],'Enable','off'); handlesVec{end+1} = 'setLeftBound';
setRightBound = uicontrol('Parent',tab1,'Style','edit','String','','Position',[leftPos+425,BScanBottom+BScanHeight+50,50,25],'Enable','off'); handlesVec{end+1} = 'setRightBound';

addModeControl = uicontrol('Parent',tab1,'Style','pushbutton','String','AddMode','Position',[leftPos+500,BScanBottom+BScanHeight+50,80,25],'Enable','off','Callback',{@addMode_Callback}); handlesVec{end+1} = 'addModeControl';
rmMode = uicontrol('Parent',tab1,'Style','pushbutton','String','DelMode','Position',[leftPos+590,BScanBottom+BScanHeight+50,80,25],'Enable','off','Callback',{@rmMode_Callback}); handlesVec{end+1} = 'rmMode';
updateModes = uicontrol('Parent',tab1,'Style','pushbutton','String','UpdateMode','Position',[leftPos+700,BScanBottom+BScanHeight+50,80,25],'Enable','off','Callback',{@updateMode_Callback}); handlesVec{end+1} = 'updateModes';

% ####################################################

handlesVec{end+1} = 'projMode'; handlesVec{end+1} = 'projQBMode'; handlesVec{end+1} = 'projInitMode'; handlesVec{end+1} = 'zGrad';


% *****************************************************
% ****************** CALLBACKS ************************
% *****************************************************

function updateMode_Callback(hObject,eventdata)
	model.shapeModel(shapeModel(currentBScan)).WML = model.shapeModel(shapeModel(currentBScan)).mu - reshape(predictions{currentBScan}.prediction',1,[])';
    numModes = size(model.shapeModel(shapeModel(currentBScan)).WML,2);
    stringVals = cellstr(num2str((1:ceil(numModes/modesPerView))'));
    set(selectMode,'String',stringVals,'Value',1,'Visible','on');

end

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
			fprintf('Loading external ground truth found in GT folder');

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
			[BScanData, fileHeader, BScanHeader, BScanSeg, SLO] = HDEVolImporter(PathName,FileName,struct('plotSLO',0,'plotBScans',0,'verbose',1));
		elseif strcmp(ext,'.mat')
			[BScanData, fileHeader, BScanHeader, BScanSeg, SLO] = HDEMatImporter(PathName,FileName,struct('plotSLO',0,'plotBScans',0,'verbose',1));
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
				BScanSeg{i} = BScanSeg{i}(end:-1:1,:);
			end
		end

		% load variables for evaluating the error term
		if numBScans > 1
			margins = load([getenv('OCT_CODE_DIR') 'datafiles/healthyMargins3D.mat']);
		else
			margins = load([getenv('OCT_CODE_DIR') 'datafiles/healthyMargins2D.mat']);
		end
	
		shapeModel = ones(1,numBScans);


		% plot SLO scan
		axes(axisSLO); reset(axisSLO); imagesc(sqrt(sqrt(SLO))); colormap gray; hold on;
		predictions = cell(1,numBScans);
		mask = cell(1,numBScans);
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
		plot(columnsToPlot,predictions{currentBScan}.prediction);
	end

	% modify status of next and previous buttons
	set(findobj('Tag','prevBScan'),'Enable','on'); set(findobj('Tag','nextBScan'),'Enable','on'); 	
	if currentBScan == numBScans
		set(findobj('Tag','nextBScan'),'Enable','off');
	end
	if currentBScan == 1
		set(findobj('Tag','prevBScan'),'Enable','off');
	end

	if length(model) > 0
		numModes = size(model.shapeModel(shapeModel(currentBScan)).WML,2);
	end

	setButtonProps(handlesVec,{'Enable','off'});
	set(noPred,'Enable','on'); set(showGT,'Enable','on');
	if isfield(predictions{currentBScan},'prediction')
		setButtonProps(handlesVec,{'Enable','on'});
		if length(model) > 0
			projMode_Callback(hObject,eventdata);
		end
	else
		if length(model) > 0
			resetMode_Callback();
		end
	end
end

function printFigure_Callback(hObject,eventdata)
	hfig = figure;
	hax_new = copyobj(axisBScan, hfig);
	set(gca, 'Units', 'normalized', 'Position', [.1 .1 .85 .8]);
	set(get(gca,'Parent'),'Position',[200 50 700 600]);
    set(gcf, 'PaperPositionMode', 'auto');
	set(gca,'YTickLabels','');
	set(gca,'XTickLabels','');
	symbols = ['a':'z' 'A':'Z' '0':'9'];
%	string = symbols(randi(length(symbols),1,10));

	print(hfig,['plots/figure_' num2str(plotCounter) '.eps'],'-depsc2');
	plotCounter = plotCounter + 1;
%	close(hfig);
end

function showQB_Callback(hObject,eventdata)
	displaySegmentation(predictions{currentBScan}.q_b,columnsToPlot);
end

function showGT_Callback(hObject,eventdata)
	maskSave = mask{currentBScan};
	mask{currentBScan} = zeros(size(mask{currentBScan}));
	displaySegmentation(BScanSeg{currentBScan},1:fileHeader.SizeX);
	mask{currentBScan}= maskSave;
end


% segment BScan
function segmentBScan_Callback(hObject,eventdata)
	if strcmp(get(hObject,'Tag'),'segmentAllBScans')
		toSegment = 1:numBScans;
	elseif strcmp(get(hObject,'Tag'),'segmentBScan')
		toSegment = currentBScan;
	end
	pred = startSegmentation(struct(),struct(),toSegment);
	if exist('predictionGlobal')
		clearvars -global predictionGlobal
	end
   	columnsToPlot = collector.options.columnsPred*collector.options.clipFactor + collector.options.clipRange(1) - 1 - (collector.options.clipFactor-1);

	for i = 1:length(pred)
		numBounds = size(pred{i}.prediction{1},1); 
	    predictions{toSegment(i)}.prediction = pred{i}.prediction{1};
        predictions{toSegment(i)}.predictionInit = pred{i}.prediction_init{1};
		predictions{toSegment(i)}.appearance = pred{i}.appearanceTerms.prediction;
		predictions{toSegment(i)}.q_b = pred{i}.q_b{1};
        predictions{toSegment(i)}.funcVal = pred{i}.funcVal;
		predictions{toSegment(i)}.funcValInit = pred{i}.funcValInit;

		shapeModel(toSegment(i)) = pred{i}.modelSelect{1};
		if size(BScanSeg{toSegment(i)},1) == size(predictions{toSegment(i)}.prediction,1)
			predictions{toSegment(i)}.unsignedError = mean(mean(abs(predictions{toSegment(i)}.prediction-BScanSeg{toSegment(i)}(:,columnsToPlot))));
		end
		% if the currently displayed BSCan was segmented, display the prediction result
		if toSegment(i) == currentBScan
			% create mask and activate buttons
			if length(mask{currentBScan}) == 0
				mask{currentBScan} = zeros(size(predictions{currentBScan}.prediction)); 
			end
			% active buttons that show the prediction and other results
	        setButtonProps(handlesVec,{'Enable','on'});
			currentAppLayer = 1;
			displaySegmentation(predictions{currentBScan}.prediction,columnsToPlot);
			% display unsigned error
			if isfield(predictions{currentBScan},'unsignedError')
				set(unsignedErrorText,'String',sprintf('Error: %.2f px | %.2f mum',predictions{currentBScan}.unsignedError,predictions{currentBScan}.unsignedError*3.87));
			end
			projMode_Callback(hObject,eventdata);
		end
	end
end

function pred = startSegmentation(collectorAdd,optionsAdd,toSegment)
	files.name = fileNameScan;
	sizeDiff = (int32(fileHeader.SizeX)-model.params.X);
	if sizeDiff < 0
		updateStatus('Scan width is smaller than for the trained model, aborting segmentation');
		return
	end

	collector.options = struct('numBScansPerVolume',double(fileHeader.NumBScans),'folder_labels','','loadLabels',0,'clipRange',double([1 int32(model.params.X)]+sizeDiff/2),'folder_data',pathNameScan,'numRegionsPerVolume',1,'saveAppearanceTerms',1,'printTimings',1,'clipFactor',1,'makePredGlobal',0);
	if isfield(fileHeader,'distanceBScans')
		collector.options.distBScans = fileHeader.distanceBScans;
	end
	if numBScans > 1
		collector.options.mirrorBScan = 'OS';
	else
		collector.options.mirrorBScan = '';
	end

	collector.options.loadRoutineData = ['spectralis' upper(fileExt(2)) fileExt(3:4)];
	testFunc.name = @predVariational; testFunc.options = struct('calcFuncVal',1,'alpha',str2num(get(setAlpha,'String')),'variance',str2num(get(setVariance,'String')));
	testFunc.options.appearanceModel = model.options.appearanceModel;

	% scan format, which has double resolution in y-dimension
	if fileHeader.SizeX == 1536
		collector.options.clipRange = double([1 int32(model.params.X)*2] + (int32(fileHeader.SizeX) - model.params.X*2)/2);
		collector.options.clipFactor = 2;
	end


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

	if ~isempty(get(setLeftBound,'String')) && ~isempty(get(setRightBound,'String'))
		% set mask defined by bounds
		%collector.options.columnsPred = str2num(get(setLeftBound,'String')):2:str2num(get(setRightBound,'String'));
    	mask{currentBScan} = zeros(numBounds,length(collector.options.columnsPred));
		mask{currentBScan}(:,collector.options.columnsPred < str2num(get(setLeftBound,'String')) | collector.options.columnsPred > str2num(get(setRightBound,'String'))) = 1;
	end

	if sum(mask{currentBScan}(:)) > 0
		testFunc.options.doNotPredict = mask{currentBScan};
	end
	
	for i = 1:length(toSegment)
		collector.options.labelIDs = toSegment(i);
		collectorTest.name = @collectTestData; collectorTest.options = collector.options;

		updateStatus(sprintf('Segmenting BScan %d.',toSegment(i))); 
		set(f,'Pointer','watch');
		drawnow

		pred{i} = testFunc.name(files,collectorTest,struct(),model,testFunc.options);
		updateStatus(sprintf('Finished segmenting BScan %d: Likelihood: %.3f, Likelihood/Point: %.3f.',toSegment(i),-sum(pred{i}.funcVal.q_c_data(:)),-mean(pred{i}.funcVal.q_c_data(pred{i}.funcVal.q_c_data~=0)))); drawnow
	end
	set(f,'Pointer','arrow');
end

function showAppearance_Callback(hObject,eventdata)
	currentAppLayer = get(switchAppearance,'Value');
	displayAppearance();
end

function switchAppearance_Callback(hObject,eventdata)
	currentAppLayer = get(switchAppearance,'Value'); 
	displayAppearance();
end

function displayAppearance()
	axes(axisBScan); reset(axisBScan); cla(axisBScan);
	h1 = imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap(axisBScan,'gray'); hold on;
	app = zeros(size(BScanData{currentBScan}));
	app(:,columnsToPlot) = reshape(predictions{currentBScan}.appearance{1}(currentAppLayer,:),fileHeader.SizeZ,[]);
	% normalize the maximum of each column to one
	app(:,columnsToPlot) = app(:,columnsToPlot)./repmat(max(app(:,columnsToPlot)),fileHeader.SizeZ,1);
	% interpolate columns
	columnsToInterpolate = setdiff(columnsToPlot(1):columnsToPlot(end),columnsToPlot);

	[X Y] = meshgrid(columnsToPlot,1:size(BScanData{currentBScan},1));
	[Xq Yq] = meshgrid(columnsToInterpolate,1:size(BScanData{currentBScan},1));
	app(:,columnsToInterpolate) = interp2(X,Y,app(:,columnsToPlot),Xq,Yq);

	h2 = imagesc(app);
	set(h2,'AlphaData',app);
	colormap(axisBScan,[gray(64); jet(64)])
	set(h1,'CData',get(h1,'CData')/2)
	set(h2,'CData',(get(h2,'CData')+1)/2)
end

function h = displaySegmentation(prediction,columns)
	if size(prediction,1) == length(columns)
		prediction = prediction';
	end
	axes(axisBScan); reset(axisBScan);
	imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap(axisBScan,'gray'); hold on;
	h = [];
	cmap = lines(numBounds);
	if sum(prediction(:)==0) > 0
		for i = 1:size(prediction,1)
			starts = find(([~mask{currentBScan}(i,1:end-1) 0] -[0 ~mask{currentBScan}(i,1:end-1)])==1);
			stops = find(([~mask{currentBScan}(i,1:end-1) 0] -[0 ~mask{currentBScan}(i,1:end-1)])==-1);
			for j = 1:length(starts)
				h(end+1) = plot(columns(starts(j):stops(j)-1),prediction(i,starts(j):stops(j)-1),'Color',cmap(i,:));
			end
		end
	else
		h = plot(columns,prediction);
	end
end

function showDataQuality_Callback(hObject,eventdata)
	plotFuncVal(predictions{currentBScan}.prediction,squeeze(predictions{currentBScan}.funcVal.q_c_data),margins.mu(4,:),margins.sigma(4,:));
end

function showDataInit_Callback(hObject,eventdata)
	plotFuncVal(predictions{currentBScan}.predictionInit,squeeze(predictions{currentBScan}.funcValInit.q_c_data),margins.mu(4,:),margins.sigma(4,:));
end

function showEntropy_Callback(hObject,eventdata)
	plotFuncVal(predictions{currentBScan}.prediction,squeeze(predictions{currentBScan}.funcVal.q_c_singleton),margins.mu(1,:),margins.sigma(1,:));
end

function plotFuncVal(prediction,data,mu,sigma)
    h = displaySegmentation(prediction,columnsToPlot);
    set(h(:),'Color',[0.8 0.8 1]);
    makeQualityOverlay(prediction,data,mu,sigma,columnsToPlot);
end


function showPredInit_Callback(hObject,eventdata)
	displaySegmentation(predictions{currentBScan}.predictionInit,columnsToPlot);
	updateStatus(sprintf('Finished segmenting BScan %d: Likelihood: %.3f, Likelihood/Point: %.3f.',currentBScan,-sum(predictions{currentBScan}.funcValInit.q_c_data(:)),-mean(predictions{currentBScan}.funcValInit.q_c_data(predictions{currentBScan}.funcValInit.q_c_data~=0)))); drawnow 
end
function showPred_Callback(hObject,eventdata)
	displaySegmentation(predictions{currentBScan}.prediction,columnsToPlot);
	updateStatus(sprintf('Finished segmenting BScan %d: Likelihood: %.3f, Likelihood/Point: %.3f.',currentBScan,-sum(predictions{currentBScan}.funcVal.q_c_data(:)),-mean(predictions{currentBScan}.funcVal.q_c_data(predictions{currentBScan}.funcVal.q_c_data~=0)))); drawnow 
end

function noPred_Callback(hObject,eventdata)
	axes(axisBScan); reset(axisBScan);
  	imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap(axisBScan,'gray');
end

% load model
function loadModel_Callback(hObject,eventdata)
	[FileName,PathName,FilterIndex] = uigetfile('*.mat','Load Model',[getenv('OCT_CODE_DIR') '/datafiles']);

    if FileName ~= 0
		tmp = load([PathName FileName],'model');
		try 
			model = tmp.model;
		catch
			updateStatus('Invalid model file.');
			return
		end
		model.appearanceModel{1}
		tBScan = printDataInTable(tab3,model.params);
		tBScan.Position(3) = tBScan.Extent(3); tBScan.Position(4) = tBScan.Extent(4);           
		if length(fieldnames(model)) > 0
			updateStatus('Model loaded.');
			setButtonProps({'segmentBScan','setAlpha','setVariance','resetMode'},{'Enable','on'});

			% show first mode in the respective tab
			numModes = size(model.shapeModel(1).WML,2); z = zeros(numModes,1);
			set(normZText,'String',sprintf('||z|| = %.2f',norm(z)));
			currentMode = 1; 
			stringVals = cellstr(num2str((1:ceil(numModes/modesPerView))'));
			set(selectMode,'String',stringVals,'Value',1,'Visible','on');
			set([nextMode,setLeftBound,setRightBound],'Enable','on'); 
			set(modelNameText,'String',FileName);
			switchMode(currentMode);
			updateComposition_Callback();
			numBounds = length(model.params.EdgesTrain);
		else
			updateStatus('Loading of model failed.');
		end
	end
end

% ################# SEGMENTATION MASK CALLBACKS ##############

function showMask_Callback(hObject,eventdata)
	renderMask();
end

function addToMask_Callback(hObject,eventdata)
	rect = getrect(axisBScan);
	% select columns
	columnsSelect = find(columnsToPlot>rect(1) & columnsToPlot < rect(1)+rect(3));
	% select boundaries
	tmp = find(predictions{currentBScan}.prediction(:,columnsSelect) < rect(2) + rect(4) & predictions{currentBScan}.prediction(:,columnsSelect) > rect(2));
	mask{currentBScan}(tmp + numBounds*(columnsSelect(1)-1)) = 1;
	renderMask();
end

function removeFromMask_Callback(hObject,eventdata)
	rect = getrect(axisBScan);
	% select columns
	columnsSelect = find(columnsToPlot>rect(1) & columnsToPlot < rect(1)+rect(3));
	% select boundaries
	tmp = find(predictions{currentBScan}.prediction(:,columnsSelect) < rect(2) + rect(4) & predictions{currentBScan}.prediction(:,columnsSelect) > rect(2));
	mask{currentBScan}(tmp + numBounds*(columnsSelect(1)-1)) = 0;
	renderMask();
end

function resetMask_Callback(hObject,eventdata)
	mask{currentBScan} = zeros(size(mask{currentBScan}));
	renderMask();
end

function renderMask()
	h = displaySegmentation(predictions{currentBScan}.prediction,columnsToPlot);
	set(h(:),'Color',[0 0.8 0.1]);
	makeQualityOverlay(predictions{currentBScan}.prediction,mask{currentBScan}',zeros(1,numBounds),ones(1,numBounds)*0.01,columnsToPlot);
end

function autoMask_Callback(hObject,eventdata)
	[mask{currentBScan} E Eafter] = calcMask(predictions{currentBScan}.prediction,predictions{currentBScan}.predictionInit,squeeze(predictions{currentBScan}.funcVal.q_c_data),squeeze(predictions{currentBScan}.funcValInit.q_c_data),margins);
	mask{currentBScan}(predictions{currentBScan}.prediction==0) = 1;
	renderMask();
end

% ################# SHAPE PRIOR VIEWER CALLBACKS  #############

function addMode_Callback(hObject,eventdata)
	columnsMode = mask{currentBScan}(1,:)==0;
    numColumnsMode = sum(columnsMode);
	modeShape = cos(linspace(-pi/2,pi/2,numColumnsMode))*(numColumnsMode/2);
	addMode(modeShape);
	
	stringVals = cellstr(num2str((1:ceil(numModes/modesPerView))'));
    set(selectMode,'String',stringVals,'Value',1,'Visible','on');
	updateStatus(sprintf('Added new mode, now %d modes in total\n',numModes));
end

function addMode(modeShape)
  	newMode = zeros(size(model.shapeModel(shapeModel(currentBScan)).WML(:,1)));
    columnsMode = mask{currentBScan}(1,:)==0;
    numColumnsMode = sum(columnsMode);
    % add drusen shape to boundaries 6 to 9
    columnsToAdd = repmat(find(columnsMode),4,1) + repmat((5:8)'*length(newMode)/numBounds,1,numColumnsMode);
    newMode(columnsToAdd) = repmat(modeShape,4,1);

    model.shapeModel(shapeModel(currentBScan)).WML(:,end+1) = newMode;
	numModes = size(model.shapeModel(shapeModel(currentBScan)).WML,2);
end

function rmMode_Callback(hObject,eventdata)
	model.shapeModel(shapeModel(currentBScan)).WML(:,end) = [];

	numModes = size(model.shapeModel(shapeModel(currentBScan)).WML,2);
	stringVals = cellstr(num2str((1:ceil(numModes/modesPerView))'));
    set(selectMode,'String',stringVals,'Value',1,'Visible','on');
	updateStatus(sprintf('Removed mode, now %d modes in total\n',numModes));
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
	for j = 1:modesPerView
		axes(axisMode(j)); reset(axisMode(j)); set(selectModeComp(j),'String',0)
		if j+startMode <= numModes
			if isstruct(collector)
				plotWithMask(columnsToPlot,reshape(model.shapeModel(shapeModel(currentBScan)).WML(:,j+startMode),length(model.shapeModel(shapeModel(currentBScan)).columnsShape),[])',mask{currentBScan});
				set(gca,'XLim',[min(columnsToPlot) max(columnsToPlot)]);
			else
				plot(reshape(model.shapeModel(shapeModel(currentBScan)).WML(:,j+startMode),length(model.shapeModel(shapeModel(currentBScan)).columnsShape),[]));
			end
			set(selectModeComp(j),'String',sprintf('%.2f',z(j+startMode)));
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
	set(normZText,'String',sprintf('||z|| = %.2f',norm(z)));
	set(selectModeComp(id),'String',z(id+startMode));
	updateComposition_Callback();
end

function updateComposition_Callback(hObject,eventdata)
	startMode = (currentMode-1)*modesPerView;
	for i = 1:modesPerView
		if (i+startMode) <= numModes 
			z(i+startMode) = str2num(get(selectModeComp(i),'String'));
		end
	end	
	% calculate the composition based on the current z
	if length(mask{currentBScan}) == 0
		idx_b = 1:size(model.shapeModel(shapeModel(currentBScan)).WML,1);
	else
		idx_b = find(reshape(~mask{currentBScan}',1,[]));
	end
	comp = zeros(size(model.shapeModel(shapeModel(currentBScan)).WML,1),1);
	comp(idx_b) = model.shapeModel(shapeModel(currentBScan)).mu(idx_b) + model.shapeModel(shapeModel(currentBScan)).WML(idx_b,:)*z;
	axes(axisModeComp); reset(axisModeComp);
	if isstruct(collector)
		imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap(axisBScan,'gray'); hold on;
		plotWithMask(columnsToPlot,reshape(comp,length(model.shapeModel(shapeModel(currentBScan)).columnsShape),[])',mask{currentBScan}); set(gca,'YDir','reverse');
	else
		plot(reshape(comp,length(model.shapeModel(shapeModel(currentBScan)).columnsShape),[])); set(gca,'YDir','reverse');
	end
end

function resetMode_Callback(hObject,eventdata)
	z = zeros(numModes,1);
	set(normZText,'String',sprintf('||z|| = %.2f',norm(z))); set(selectModeComp,'String',0);
	updateComposition_Callback;
end


function projMode_Callback(hObject,eventdata)
	tag = get(hObject,'Tag');
	if isfield(predictions{currentBScan},tag)
		eval(['projModeFunc(reshape(predictions{currentBScan}.' tag ''',[],1));']);
	else
		projModeFunc(reshape(predictions{currentBScan}.prediction',[],1));
	end
end


function projModeFunc(prediction)
	if isfield(predictions{currentBScan},'prediction')
		idx_b = find(reshape(~mask{currentBScan}',1,[]));
		M = model.shapeModel(shapeModel(currentBScan)).WML(idx_b,:)'*model.shapeModel(shapeModel(currentBScan)).WML(idx_b,:) + eye(numModes)*model.shapeModel(shapeModel(currentBScan)).sigmaML;
		z = inv(M)*(model.shapeModel(shapeModel(currentBScan)).WML(idx_b,:)'*(prediction(idx_b) - model.shapeModel(shapeModel(currentBScan)).mu(idx_b)));
		% update z-vec length field
		set(normZText,'String',sprintf('||z|| = %.2f',norm(z)));
		for i = 1:modesPerView
			if (i + (currentMode-1)*modesPerView <= numModes)
				set(selectModeComp(i),'String',sprintf('%.2f',z(i + (currentMode-1)*modesPerView)));
			end
		end
		comp = zeros(prod(size(predictions{currentBScan}.prediction)),1);
		comp(idx_b) = model.shapeModel(shapeModel(currentBScan)).mu(idx_b) + model.shapeModel(shapeModel(currentBScan)).WML(idx_b,:)*z;
		axes(axisModeComp); reset(axisModeComp);
		imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap(axisBScan,'gray'); hold on;

		plotWithMask(columnsToPlot,reshape(comp,length(model.shapeModel(shapeModel(currentBScan)).columnsShape),[])',mask{currentBScan}); set(gca,'YDir','reverse');
	end
end

function zGrad_Callback(hObject,eventdata)
	idx_b = find(reshape(~mask{currentBScan}',1,[]));
    M = inv(model.shapeModel(shapeModel(currentBScan)).WML(idx_b,:)'*model.shapeModel(shapeModel(currentBScan)).WML(idx_b,:) + eye(numModes)*model.shapeModel(shapeModel(currentBScan)).sigmaML);
	tmp = reshape(predictions{currentBScan}.prediction',[],1);
	vec = tmp(idx_b) - model.shapeModel(shapeModel(currentBScan)).mu(idx_b);
	tmp = M*model.shapeModel(shapeModel(currentBScan)).WML(idx_b,:)';
	% calc gradient of norm(z)
	grad = zeros(size(model.shapeModel(shapeModel(currentBScan)).WML,1),1);
	for i = 1:size(model.shapeModel(shapeModel(currentBScan)).WML,2)
	    P = tmp(i,:)'*tmp(i,:);
	    grad(idx_b) = grad(idx_b) + P*vec;
	end
	grad = grad/norm(z);
%	axes(axisModeComp); reset(axisModeComp);
	figure;
	subplot(3,1,[1,2]); 
	imagesc(sqrt(sqrt(BScanData{currentBScan}))); colormap(axisBScan,'gray'); colormap gray; hold on;
	plotWithMask(columnsToPlot,predictions{currentBScan}.prediction,mask{currentBScan}); set(gca,'YDir','reverse');
	subplot(3,1,3);
	plot(columnsToPlot,(reshape(grad,length(model.shapeModel(shapeModel(currentBScan)).columnsShape),[])));
	set(gca,'XLim',[1 size(BScanData{currentBScan},2)]);
	h = line([1 size(BScanData{currentBScan},2)],[0 0]);
	set(h,'Color',[0.7 0.7 0.7]);
end

% ############### END SHAPE PRIOR VIEWER CALLBACKS #####################

function plotWithMask(columns,toPlot,mask)
    cmap = lines(numBounds);
	if length(mask) == 0
		mask = zeros(size(toPlot));
	end
	for i = 1:size(toPlot,1)
		starts = find(([~mask(i,1:end-1) 0] -[0 ~mask(i,1:end-1)])==1);
		stops = find(([~mask(i,1:end-1) 0] -[0 ~mask(i,1:end-1)])==-1);
		for j = 1:length(starts)
			plot(columns(starts(j):stops(j)-1),toPlot(i,starts(j):stops(j)-1),'Color',cmap(i,:)); hold on;
		end
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

%% *********** OLD CODE *****************
%function addMode_Callback(hObject,eventdata)
%	% add the difference between the new segmentation and the
%	diff = predictions{currentBScan}.prediction-predictions{currentBScan}.prediction{2};
%	for i = 1:9
%    	if sum(diff(i,:)~=0) > 5
%        	PP = splinefit(1:250,diff(i,:),round(sum(diff(i,:)~=0)/2));
%	        diffSmooth(i,:) = ppval(PP,1:250);
%    	else
%	        diffSmooth(i,:) = zeros(1,250);
%	    end
%	end
%	dist = (currentBScan - (floor(collector.options.numBScansPerVolume/2)+1))*collector.options.distBScans;
%    % find closest B-scan in the shape model
%    [~,scanID] = min(abs(model.params.BScanPositions-dist));
%	
%	model.shapeModel(scanID).WML = [model.shapeModel(scanID).WML reshape(diffSmooth',[],1)*0.5];
%	updateStatus(sprintf('Added mode to %d''th shape model.',scanID)); drawnow
%	predictions{currentBScan}.prediction(2) = []; predictions{currentBScan}.funcVal(2) = []; set(addMode,'Enable','off'); set(showPredNew,'Value',0); set(showPredNew,'Enable','off');
%	set(hObject,'Tag','segmentBScan');
%	segmentBScan_Callback(hObject,eventdata);
%end
%function resegmentBScan_Callback(hObject,eventdata)
%	toSegment = currentBScan;
%
%	q_c_data = squeeze(predictions{toSegment}.funcVal{1}.q_c_data);
%	idxHealthy =  (q_c_data < repmat(margins.mu(4,:),250,1) + repmat(margins.sigma(4,:)*5,250,1))';
%	% recalculate all boundaries 6-9
%	idxHealthy(6:9,sum(~idxHealthy(6:9,:))>0) = 0;
%	optionsAdd.segmentation = predictions{toSegment}.prediction{1};
%	optionsAdd.idxRecalc = ~idxHealthy;
%	optionsAdd.onlyInitialize = 1;
%	optionsAdd.appearance = predictions{toSegment}.appearance;
%
%	pred = startSegmentation(struct(),optionsAdd,toSegment);
%
%	predictions{toSegment}.prediction{2} = predictions{toSegment}.prediction{1};
%	predictions{toSegment}.prediction{2}(~idxHealthy) = pred{1}.predictionInit(~idxHealthy);
%	predictions{toSegment}.funcVal{2} = pred{1}.funcVal;
%	set(showPredNew,'Enable','on','Value',1); set(addMode,'Enable','on');
%	displaySegmentation(predictions{toSegment}.prediction{2},columnsToPlot);
%end
%

