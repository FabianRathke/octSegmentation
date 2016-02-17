function labelTool(filename,collector,restore)
% labelTool - Allows the user to add ground truth to a B-Scan
%
% If restore is true, the function checks for previously created label files in collector.options.saveDir.
% If no such file is detected, the fallback position is to load labels that are provided with the dataset. The function loadLabels is called
% which then uses a user-predefined loading scheme to obtain the labels.
%
% Syntax:
%   labelTool(filename,collector,restore)
%
% Input: 
%	filename  - [string] name of the matfile that holds the data
%	collector - [struct] containts options required (see the manual for the required options)
%	restore   - [boolean] whether existing labels should be used as starting point
%
% See also: loadData, loadLabels
%
% -------------
% | Shortcuts |
% -------------
%
% Different modes for left-click:
% ----------------------------------
%   a - (a)dds marker
%   d - (d)eletes (nearest) marker
%   e - (e)dits (nearest) marker
%   f - (f)uses current line with selected marker
%
% General Controls
% ----------------------------------
%   s           - saves interpolation data and markers; generates preview
%   n           - generates new line
%   r           - flattens the scan along the current boundary
%   backspace   - deletes current marker
%   control + d - deletes current line (marked green)
%   q,w	        - cycle through boundaries
%   z,x         - cycle through markers within the current boundary
%   b           - hide/show lines
%   control + b - hide/show markers
%   c           - skip cropping
%   arrow keys  - use to move the selected marker

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 17-Dec-2013

close all;
% load image data
saveFileName = sprintf('%s%s_%d_coordinates.mat',collector.options.saveDir,strrep(filename,'*.mat',''),collector.options.labelID);

% load the B-Scan
B0 = loadData(filename,collector.options);

clear Data;
% plot normalized image
f = figure('KeyPressFcn',@keyPress);
h = imagesc(B0); colormap gray;
t = title('');
set(h,'buttondownfcn',@mouseClick);
hold on;

% holding point/line information
points = struct();
initializeStruct();
interpolations = zeros(1,size(B0,2));

% counter/history
current_line = 1;
current_point = 1;
% handles for plots
line_handles = [];
point_handles = [];
active_point_handle = [];
reference_line = [];
reference_line_idx = [];

% load old save file for edit
if restore
	if exist(saveFileName)
		Data_save = load(saveFileName);
		if isfield(Data_save,'points')
			points = Data_save.points;
			interpolations = Data_save.interpolation;
		else
			points = Data_save.points_save;
			interpolations = Data_save.interpolation;
		end
		fprintf('Restored previously created ground truth\n');

		clear Data_save;
	else
		% use segmentation lines given in the image file
		B0seg = loadLabels(filename,collector.options);
		positions = linspace(1,size(B0seg,2),length(1:5:size(B0seg,2)));
		interpolations = B0seg;
		for i = 1:size(B0seg,1)
			points(i).coordinates = [round(positions); B0seg(i,round(positions))];
			points(i).time(:,1:length(positions)) = repmat(clock',1,length(positions));
			points(i).order = 1:length(positions);
		end
		fprintf('Loaded ground truth provided with the data set\n');
	end

	% initialize boundaries
	for i = 1:length(points)
		 line_handles(i) = plot(1:size(B0,2),interpolations(i,:),'-y');
		 point_handles(i) = plot(points(i).coordinates(1,:),points(i).coordinates(2,:),'sy','MarkerFaceColor','yellow','MarkerEdgeColor','blue');
	end
	
	activateLine(current_line);
	active_point_handle = plot(points(current_line).coordinates(1,current_point),points(current_line).coordinates(2,current_point),'sr','MarkerFaceColor','red','MarkerEdgeColor','white');
	activatePoint();
end

t = title('Crop upper part of the image; use ''c'' to skip');
set(t,'FontSize',12,'FontWeight','bold');
upper_crop = 1;
lower_crop = size(B0,1);
modus = 'crop_image_upper';

% mouseClick Action 
function mouseClick(src,eventdata)
	point = get(gca,'CurrentPoint'); % button down detected
	if strcmp(modus,'add_points')
		%fprintf('point added to line %d\n',current_line);
		% first point of new line, make two points on the left and right borders of the scan	
		if isempty(points(current_line).coordinates)
			points(current_line).coordinates = [[1 point(1,2) ]' [size(B0,2) point(1,2)]'];
			points(current_line).order = [1 2];
			points(current_line).time = [repmat(clock',1,2)];
			% do interpolation
			interpolations(current_line,:) = ones(1,size(B0,2))*point(1,2);
			% draw line and marker
			if (current_line == 1) line_handles(current_line) = plot(1:size(B0,2),point(1,2)*ones(1,size(B0,2)),'-g');
			else line_handles(current_line) = plot(1:size(B0,2),point(1,2)*ones(1,size(B0,2)),'-g','Visible',get(line_handles(1),'Visible')); end
			
			point_handles(current_line) = plot(points(current_line).coordinates(1,:),points(current_line).coordinates(2,:),'sg','MarkerFaceColor','green','MarkerEdgeColor','blue');
			current_point = 1;
			if isempty(active_point_handle)
				active_point_handle = plot(points(current_line).coordinates(1,1),points(current_line).coordinates(2,1),'sr','MarkerFaceColor','red','MarkerEdgeColor','white');
			else
				activatePoint();
			end
		else
			addPoint(point(1,1:2))
		end
	elseif strcmp(modus,'crop_image_upper')
		upper_crop = round(point(1,2));
		if ~isempty(points(current_line).coordinates)
			for i = 1:length(points)
				points(i).coordinates(2,:) = points(i).coordinates(2,:) - upper_crop + 1;
			end
		end
		updateInterpolationAll();

		set(h,'CData',B0(upper_crop:end,:));
		set(gca,'YLim',[0.5 size(B0,1)-upper_crop-1.5]);
		modus = 'crop_image_lower';
		set(t,'string','''Crop lower part of the image; use ''c'' to skip''')
	elseif strcmp(modus,'crop_image_lower')
		lower_crop = round(point(1,2)) + upper_crop;
		set(h,'CData',B0(upper_crop:lower_crop,:));
		set(gca,'YLim',[0.5 lower_crop - upper_crop + 0.5])
		screen_size = get(0, 'ScreenSize');
		set(f, 'Position', [40 40 screen_size(3)-150 screen_size(4)-150]);
		set(gca,'Position',[0.02 0 0.95 0.95])
		modus = 'add_points';
		set(t,'String','Add Modus');
	else
		dist = inf;
		idx_select = [];
		% search nearst point
		for i = 1:length(points)
			if length(points(i).order) > 0
				[a idx] = min(sum(abs(points(i).coordinates-repmat(point(1,1:2)',1,length(points(i).order)))));
				if a < dist
					dist = a;
					idx_select = [i idx];
				end
			end
		end
		% edit points
		if strcmp(modus,'edit_points')
			points(idx_select(1)).coordinates(:,idx_select(2)) = point(1,1:2);
			points(idx_select(1)).time(:,idx_select(2)) = clock;
			current_line = idx_select(1);
			current_point = idx_select(2);
			updateInterpolation();
		% delete points
		elseif strcmp(modus,'select_points')
			current_point = idx_select(2);
			activatePoint();
		elseif strcmp(modus,'fuse_points')
			addPoint(points(idx_select(1)).coordinates(:,idx_select(2)));
		else
			% last two points, delete line
			if length(points(idx_select(1)).order)==2
				deleteLine(idx_select(1));
				current_line = idx_select(1)-1;
				current_point = 1;
				activatePoint();
			else				
				points(idx_select(1)).coordinates(:,idx_select(2)) = [];
				points(idx_select(1)).time(:,idx_select(2)) = [];
				points(idx_select(1)).order(idx_select(2)) = [];
				current_line = idx_select(1);
				[tmp idx] = max(points(current_line).order);
				current_point = idx;

				updateInterpolation();
			end
		end

		checkLines();

		if length(points) == 0
			initializeStruct()
		end
	end
end

% Keypress action
function keyPress(src,eventdata)
	% add new line
	if length(eventdata.Modifier) == 1
		if strcmp(eventdata.Modifier{1},'control')
			if strcmp(eventdata.Key,'d') && strcmp(modus,'add_points')
				deleteLine(current_line);	
			end

			if strcmp(eventdata.Key,'b')
    			if length(point_handles) > 0 && strcmp(get(point_handles(1),'Visible'),'off')
	                set(point_handles,'Visible','on');
	            else
	                set(point_handles,'Visible','off');
	            end
			end
		end
	else
		if strcmp(eventdata.Key,'c')
			if strcmp(modus,'crop_image_upper')
				set(t,'String','Crop Image (Below Retina); use ''c'' to skip');
				modus = 'crop_image_lower';
			else
				modus = 'add_points';
				set(t,'String','Add Modus');
			end
		end

		if strcmp(eventdata.Key,'n')
			if ~isempty(points(current_line).coordinates)
				set(line_handles(current_line),'color','yellow');
				set(point_handles(current_line),'color','yellow','MarkerFaceColor','yellow');
				checkLines();
				current_line = length(points) + 1;
				points(current_line).coordinates = [];
				points(current_line).time = [];
				points(current_line).order = [];
				deactivatePoint();

				modus = 'add_points';
				set(t,'String','Add Modus');
				fprintf('inserted line %d\n',current_line);
			end
		end

		if strcmp(eventdata.Key,'r')
			numRows = lower_crop - upper_crop + 1;
			if numRows >= size(B0,1)
				 fprintf('Adjustment only possible for cropped scans.\n');
			elseif isempty(reference_line_idx)
				% get leftest x-coordinate of the current line
				x_coord = points(current_line).coordinates(2,1);
			
				% adjust scan	
				B0Tmp = zeros(numRows,size(B0,2)*2-1);
				difference = numRows - max(round(interpolations(current_line,:))) + 30;
				for i = 1:size(B0,2)
					position = round(interpolations(current_line,i)) + upper_crop - 1;
					firstPos = (difference + position) - numRows + 1;
					if firstPos < 1
						vec = [zeros(1,abs(firstPos-1)) B0(1:position+difference,i)'];
					else
						vec = B0(firstPos:position+difference,i);
					end
					B0Tmp(1:2:length(vec)*2,i) = vec;
					% stretch the plot
					B0Tmp(2:2:(length(vec)-1)*2,i) = interp1(1:2:length(vec)*2,vec,2:2:(length(vec)-1)*2,'linear');
				end
				set(h,'CData',B0Tmp);
				set(gca,'YLim',[0.5 size(B0Tmp,1) + 0.5])

				reference_line = interpolations(current_line,:);
				reference_line_idx = current_line;

				% adjust segmentations
				difference = ((numRows - difference)*2-1) - (interpolations(current_line,:)*2 - 1);
				for i = 1:length(points)
					points(i).coordinates(2,:) = points(i).coordinates(2,:)*2 - 1;
					points(i).coordinates(2,:) = points(i).coordinates(2,:) + difference(round(points(i).coordinates(1,:)));
				end
				updateInterpolationAll();
			else
				fprintf('Boundary already adjusted, only one adjustment possible yet\n');
			end
			modus = 'add_points';
		end

		if regexp(eventdata.Key,'arrow') > 0
			directions = {'right','left','up','down'};
			directionIdx = [1, 1, 2, 2];
			change = [1, -1, -0.5, 0.5];
			string = eventdata.Key(1:regexp(eventdata.Key,'arrow')-1);
			stringidx = strcmp(string,directions);

			points(current_line).coordinates(directionIdx(stringidx),current_point) = points(current_line).coordinates(directionIdx(stringidx),current_point) + change(stringidx);
			% ommit points being outside the scan
			if points(current_line).coordinates(1,current_point)<1
				points(current_line).coordinates(1,current_point) = 1;
			end
			if points(current_line).coordinates(1,current_point) > size(B0,2)
				points(current_line).coordinates(1,current_point) = size(B0,2);
			end
			updateInterpolation();
		end	

		if strcmp(eventdata.Key,'pagedown')
			points(current_line).coordinates(2,current_point) = points(current_line).coordinates(2,current_point) + 10;
			updateInterpolation();
		end	
		if strcmp(eventdata.Key,'pageup')
			points(current_line).coordinates(2,current_point) = points(current_line).coordinates(2,current_point) - 10;
			updateInterpolation();
		end	

		if strcmp(eventdata.Key,'x')
			current_point = current_point + 1;
			if current_point > length(points(current_line).order)
				current_point = 1;
			end
			activatePoint();
		end

		if strcmp(eventdata.Key,'z')
			current_point = current_point - 1;
			if current_point < 1 
				current_point = length(points(current_line).order)
			end
			activatePoint();
		end

		if strcmp(eventdata.Key,'o')
			B0Norm = B0Norm.^0.8;
			set(h,'CData',B0Norm);
		end

		if strcmp(eventdata.Key,'i')
			B0Norm = B0Norm.^1.2;
			set(h,'CData',B0Norm);
		end

		if strcmp(eventdata.Key,'b')
			if length(line_handles) > 0 && strcmp(get(line_handles(1),'Visible'),'off')
				set(line_handles,'Visible','on');
			else
				set(line_handles,'Visible','off');
			end
		end

		% remove current point
		if strcmp(eventdata.Key,'backspace')
			% line has only two points --> delete it
			if length(points(current_line).order) == 2
				deleteLine(current_line);
			else
				% delete point from line
				

				points(current_line).coordinates(:,current_point) = [];
				points(current_line).time(:,current_point) = [];
				points(current_line).order(current_point) = [];
				[tmp idx] = max(points(current_line).order);
				current_point = idx;

				% update line spline
				updateInterpolation();
			end
			if length(points) > 0 fprintf('%d points remaining in line %d\n',length(points(current_line).order),current_line); end
			%checkLines();
		end

		% save lines
		if strcmp(eventdata.Key,'s')
	%		% transform lines from normalized into original
	%		for i = 1:length(points)
	%			if length(points(i).order) > 0
	%				a = diag(Raster{1}(floor(interpolations(i,:)),:));
	%				b = diag(Raster{1}(ceil(interpolations(i,:)),:));
	%
	%				interpolations_norm(i,:) = a + (b-a).*(ceil(interpolations(i,:))-interpolations(i,:))';	
	%
	%				a = round(points(i).coordinates);
	%				a = diag(Raster{1}(a(2,:),a(1,:)));
	%				points_norm(i).coordinates = [points(i).coordinates(1,:); a'];
	%%				interpolations_norm(i,:) = interp1(points(i).coordinates(1,:),a,1:size(B0,2),'linear');
	%			end
	%		end

			% sort lines
			[a b] = sort(interpolations(:,1));
			interpolation = interpolations(b,:);
		
			reference_line_idx_save = find(b==reference_line_idx);

			points_save = points(b);
			% remove zoom factor
			if ~isempty(reference_line_idx_save)
				if length(points) > size(interpolation,1)
					points(end) = [];
				end
				for i = 1:length(points)
					if length(points(i).order) > 0
						interpolation(i,:) = (interpolation(i,:) + ones(1,size(B0,2)))/2;
					end
				end

				diff = reference_line - interpolation(reference_line_idx_save,:);

				for i = 1:length(points)
					% project onto reference line
					interpolation(i,:) = interpolation(i,:) + diff;
				end
			end

			for i = 1:length(points)
				% add crop boundary
				interpolation(i,:) = interpolation(i,:) + upper_crop - 1;
				points_save(i).coordinates(2,:) = interpolation(i,round(points_save(i).coordinates(1,:)));
			end

			fprintf('Saved labels in %s\n',saveFileName);
			save(saveFileName,'interpolation','points_save');
			figure
			imagesc(B0); colormap gray; hold on;
			plot(1:size(B0,2),interpolation,'y-','LineWidth',1);
	%		hold on;
	%		for i = 1:length(points)
	%			plot(points_norm(i).coordinates(1,:),points_norm(i).coordinates(2,:),'sy','MarkerFaceColor','yellow','MarkerSize',3);
	%		end
		end

		% switch trough lines
		if strcmp(eventdata.Key,'w') && strcmp(modus,'add_points') && length(points)>1
			if current_line == length(points)
				deactivateLine(current_line);
				current_line = 1;
				activateLine(current_line);
			else
				deactivateLine(current_line);
				current_line = current_line + 1;
				activateLine(current_line);
			end
			[tmp idx] = max(points(current_line).order);
			current_point = idx;
			activatePoint();
		end
		
		if strcmp(eventdata.Key,'q') && strcmp(modus,'add_points') && length(points)>1
			if current_line == 1
				deactivateLine(current_line);
				current_line = length(points);
				activateLine(current_line);
			else
				deactivateLine(current_line);
				current_line = current_line - 1;
				activateLine(current_line);
			end
			[tmp idx] = max(points(current_line).order);
			current_point = idx;
			activatePoint();
		end

		% switch between modi
		if strcmp(eventdata.Key,'e')
			modus = 'edit_points';
			set(t,'String','Edit Modus');
			deactivateLine(current_line);
		end

		if strcmp(eventdata.Key,'g')
			modus = 'select_points';
			set(t,'String','Select Point Modus');
		end

		if strcmp(eventdata.Key,'a')	
			activateLine(current_line);
			modus = 'add_points';
			set(t,'String','Add Modus');
		end

		if strcmp(eventdata.Key,'f')
			activateLine(current_line);
			modus = 'fuse_points';
			set(t,'String','Fuse Modus');
		end


		if strcmp(eventdata.Key,'d')
			modus = 'delete_points';
			set(t,'String','Delete Modus');
			deactivateLine(current_line);
		end
	end
end

function addPoint(point)
	points(current_line).coordinates(:,end+1) = point;
	points(current_line).order(end+1) = max(points(current_line).order)+1;
	points(current_line).time(:,end+1) = clock;
	% resort points according to x value
	[tmp idx] = sort(points(current_line).coordinates(1,:));
	points(current_line).coordinates = points(current_line).coordinates(:,idx);
	points(current_line).order = points(current_line).order(idx);
	points(current_line).time = points(current_line).time(:,idx);

	[tmp idx] = max(points(current_line).order);
	current_point = idx;

	% do interpolation
	updateInterpolation();
end

function updateInterpolation()
	if ~isempty(points(current_line).coordinates)
		interpolations(current_line,:) = interp1(points(current_line).coordinates(1,:),points(current_line).coordinates(2,:),1:size(B0,2),'pchip');
		set(line_handles(current_line),'YData',interpolations(current_line,:));
		set(point_handles(current_line),'XData',points(current_line).coordinates(1,:),'YData',points(current_line).coordinates(2,:));
		activatePoint();
	end
end

function updateInterpolationAll()
	if ~isempty(points(current_line).coordinates)
		for i = 1:length(points)
			interpolations(i,:) = interp1(points(i).coordinates(1,:),points(i).coordinates(2,:),1:size(B0,2),'pchip');
			set(line_handles(i),'YData',interpolations(i,:));
			set(point_handles(i),'XData',points(i).coordinates(1,:),'YData',points(i).coordinates(2,:));
		end
		activatePoint();
	end
end

function activateLine(a)
	set(line_handles(a),'color','green');
	set(point_handles(a),'color','green','MarkerFaceColor','green');
end

function activatePoint()
	set(active_point_handle,'XData',points(current_line).coordinates(1,current_point));
	set(active_point_handle,'YData',points(current_line).coordinates(2,current_point),'ZData',1);

	set(t,'String',sprintf('x: %.1f y:%.1f',points(current_line).coordinates(1,current_point),points(current_line).coordinates(2,current_point)));
end

function deactivatePoint()
	set(active_point_handle,'XData',-100);
	set(active_point_handle,'YData',-100);
end

function deactivateLine(a)
	set(line_handles(a),'color','yellow');
	set(point_handles(a),'color','yellow','MarkerFaceColor','yellow');
end

function checkLines()
	for i = 1:length(points)
		if length(points(i).order)==0
			deleteLine(i);
		end
	end
end

function deleteLine(a)
	if reference_line_idx > a
		reference_line_idx = reference_line_idx - 1
	end
	if current_line > a
		current_line = current_line - 1;
	end
	points(a) = [];
	interpolations(a,:) = [];
	
	set(point_handles(a),'XData',[],'YData',[]);
	set(line_handles(a),'XData',[],'YData',[]);

	point_handles(a) = [];
	line_handles(a) = [];

	if length(points) > 0
		current_line = length(points)
		current_point = 1;
		activatePoint();
		activateLine(current_line);
	else
		initializeStruct()
	end
end

function initializeStruct()
	points(1).coordinates = [];
	points(1).time = [];
	points(1).order = [];
	current_line = 1;
end

end
