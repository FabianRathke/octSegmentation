function plotBScan(data,prediction,columns,filename)
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
	imagesc(data); colormap gray; hold on;
	plot(columns,prediction);
	print(filename,'-depsc2');
	close all
else
	figure;
	imagesc(data); colormap gray; hold on;
	plot(columns,prediction);
end
