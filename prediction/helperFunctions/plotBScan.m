function plotBScan(data,prediction,columns,filename,idxSet)
% plotBScan - plots the B-Scan given in data and the predicted segmentation thereof; saves the result in filename; if no filename given, prints directly to the screen;
%  
% Syntax:
%   plotBScan(data,prediction,columns,filename)
%
% Inputs:
%   data       - [matrix] gray-valued image of the B-Scan
%   prediction - [matrix] segmentation lines
%   columns    - [array] columns of the B-Scan for which there exists predictions 
%   filename   - [string] points to the file where the plot is saved as eps; if empty directly plots the scan to the screen
%	idxSet	   - [matrix] bindary matrix deconding which positions to plot
%
% See also: predVariational

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 02-May-2014

if ~isempty(filename)
	[pathstr, name, ext] = fileparts(filename);
	if length(pathstr)>0 && ~exist(pathstr,'dir')
		mkdir(pathstr);
	end
	figure('visible','off');
else
	figure;
end

imagesc(data); colormap gray; hold on;
if nargin > 4
	cmap = lines(size(prediction,1));
	for i = 1:size(prediction,1)
		starts = find(([idxSet(i,1:end-1) 0] -[0 idxSet(i,1:end-1)])==1);
		stops = find(([idxSet(i,1:end-1) 0] -[0 idxSet(i,1:end-1)])==-1);
		for j = 1:length(starts)
			plot(columns(starts(j):stops(j)-1),prediction(i,starts(j):stops(j)-1),'Color',cmap(i,:));
		end
	end
else
	plot(columns,prediction);
end

if ~isempty(filename)
	print(filename,'-depsc2');
	close gcf
end
