function textureModeViewer(model)
% textureModeViewer - If the texture models are based on PCA, shows the modes and linear combinations of them 
% 
% Syntax:
%   textureModeViewer(model)
%
% Inputs:
%   model     - [struct] trained segmentation model
%
% See also: octGUI

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 24-Nov-2016

close all

% patch-size
pX = model.options.appearanceModel{1}.width; 
pY = model.options.appearanceModel{1}.height;

numBScanRegions = model.options.appearanceModel{1}.numRegionsPerBScan;
modesPerView = 5;
currentMode = 1;
numModes = length(model.appearanceModel{1}(1,1,1).class_mean);
z = zeros(numModes,1);

f = figure('Position', [100 100 1100 750],'Tag','mainWindow');
axisViewer = axes('Units','Pixels','Position',[300 400 300 300]); set(axisViewer,'YTickLabel',[],'XTickLabel',[]);

ModeBottom = 300; ModeLeft = 100; ModeWidth = 800; ModeHeight = 400; subModeHeight = 150; subModeWidth = 150;

for i = 1:modesPerView
    posLeft = 50 + (i-1)*(subModeWidth+50);
    posBottom = 70;
    axisMode(i) = axes('Units','Pixels','Position',[posLeft,posBottom,subModeWidth,subModeHeight]); set(axisMode(i),'YTickLabel',[],'XTickLabel',[]);
    decreaseMode(i) = uicontrol('Style','pushbutton','String','<','Position',[posLeft+subModeWidth/2-60,posBottom-55,25,25],'Enable','on','Callback',{@changeComposition_Callback},'Tag',sprintf('DecMode%d',i));
    selectModeComp(i) = uicontrol('Style','edit','Units','Pixels','Position',[posLeft+subModeWidth/2-30,posBottom-55,60,25],'String',0,'Callback',{@updateComposition_Callback});
    increaseMode(i) = uicontrol('Style','pushbutton','String','>','Position',[posLeft+subModeWidth/2+35,posBottom-55,25,25],'Enable','on','Callback',{@changeComposition_Callback},'Tag',sprintf('IncMode%d',i));
% text field showing current Mode number
end
nextMode = uicontrol('Style','pushbutton','String','>','Position',[ModeLeft+ModeWidth/2+40,ModeBottom-55,25,25],'Enable','off','Tag','nextMode','Callback',{@switchMode_Callback});
prevMode = uicontrol('Style','pushbutton','String','<','Position',[ModeLeft+ModeWidth/2-60,ModeBottom-55,25,25],'Enable','off','Tag','prevMode','Callback',{@switchMode_Callback});
% text field showing current Mode number
selectMode = uicontrol('Style','popupmenu','String','','Position',[ModeLeft+ModeWidth/2-20,ModeBottom-55,50,25],'Callback',{@switchMode_Callback},'Visible','off');

axes(axisViewer); reset(axisViewer); imagesc(reshape(model.appearanceModel{1}(1,1,1).mu,pX,pY)); colormap gray;

stringVals = cellstr(num2str((1:ceil(numModes/modesPerView))'));
set(selectMode,'String',stringVals,'Value',1,'Visible','on');
switchMode(1)

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
            imagesc(reshape(model.appearanceModel{1}(1,1,1).W(:,j+startMode),pX,pY));
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
    comp = model.appearanceModel{1}(1,1,1).mu' + model.appearanceModel{1}(1,1,1).W*z;
    axes(axisViewer); reset(axisViewer);
	imagesc(reshape(comp,pX,pY));    
end


end
