function plotBScan(data,prediction,columns,filename)
% plotBScan - plots the B-Scan given in data and the predicted segmentation thereof; saves the result in filename
%  
% Syntax:
%   plotBScan(data,prediction,columns,filename)
%
% Inputs:
%   data       - [matrix] gray-valued image of the B-Scan
%   prediction - [matrix] segmentation lines
%   columns    - [array] columns of the B-Scan for which there exists predictions 
%   filename   - [string] points to the file where the plot is saved as eps
%
% See also: predVariational

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 12-Dec-2013

[pathstr, name, ext] = fileparts(filename);
if length(pathstr)>0 && ~exist(pathstr,'dir')
	mkdir(pathstr);
end
figure('visible','off');
imagesc(data); colormap gray; hold on;
plot(columns,prediction);
print(filename,'-depsc2');
close all

