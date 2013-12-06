function interpolation = loadLabels(filename,options)
% loadLabels - loads labels (i.e. ground truth of OCT-scans)
% 
% Inputs:
%   filename - [string] filename => the name of the matfile holding the data with an additional _option.labelID_coordinates
%   options  - [struct]
%       .folder_labels - [string] points to the folder of filename_options.labelID_coordinates
%       .labelID     - [int] Spectralis B-Scans are labeled B0,B1,..., for each ID there exists a mat-file with labels
%       .clip        - [boolean] indicates whether to clip a B-Scan at the left and right border
%       .clipRange   - [array](2) defines the left and right border of the B-Scan after clipping 
%
% Outputs:
%   interpolation - [matrix] the segmented boundaries for the respective scan
%
% See also: loadData, fetchPatches, collectTrnData, collectTestData

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 05-Dec-2013

if isfield(options,'labelID')
	load([options.folder_labels filename '_' num2str(options.labelID) '_coordinates.mat']);
	if options.printInformation > 1
		fprintf('Loaded Labels for %s and region %d\n',filename,options.labelID);
	end
else
	load([options.folder_labels filename '_coordinates.mat']);
end

if options.clip
	interpolation = interpolation(:,options.clipRange(1):options.clipRange(2));
end

end

