function newModes = addModes(columnsPred,sizeWML,columnsShape,collector)
% addMode - creates a set of shape modes that can then be used to augment the global shape prior; currently used to specify drusen modes
%
% Syntax: 
%	newModes = addModes(columnsPred,sizeWML,columnsShape,collector)
%
% Inputs:
%	columnsPred	 - [array] absolute positions of columns that restrict the area of the new mode
%	sizeWML 	 - [int] size of one mode vector
%	columnsShape - [array] absolute positions of the shape prior --> used to select the columns that are affected by the boundaries defined by columnsPred
%	collector	 - [struct] used to collector some information (e.g. number of boundaries) of the shape prior
%
% Output:
%	newModes	 - [cell array] each entry is a new mode
%

% the boundaries of the new modes
minMaxColumns = [min(columnsPred) max(columnsPred)];
columnsModes = find(columnsShape>=minMaxColumns(1) & columnsShape<=minMaxColumns(2));
numColumnsModes = length(columnsModes);
numBounds = length(collector.options.EdgesPred);

modeCounter = 1;

for numModes = 1:2
%numModes = 2;
	% coordinates
	interval = floor((numColumnsModes-1)/numModes);
	height = (diff(minMaxColumns)+1)/numModes/2;
	modeShape = cos(linspace(-pi/2,pi/2,interval+1))*height;
	modeShapeShift = [modeShape(ceil(end/2):end) modeShape(2:ceil(end/2))];

	for i = 1:numModes
		newMode = zeros(1,sizeWML);
		columnsMode = columnsModes((1:(interval+1)) + (i-1)*interval);
		% add drusen shape to boundaries 6 to 9
		columnsToAdd = repmat(columnsMode,4,1) + repmat((5:8)'*sizeWML/numBounds,1,interval+1);

		newMode(columnsToAdd) = repmat(modeShape,4,1);
		% add to cell array of modes
		newModes{modeCounter} = newMode; modeCounter = modeCounter+1;
		
		newMode(columnsToAdd) = repmat(modeShapeShift,4,1);
		% add to cell array of modes
		newModes{modeCounter} = newMode; modeCounter = modeCounter+1;
	end
end

% DEBUG: plot all modes
%figure; hold on;
%for i = 1:length(newModes)
%	plot(reshape(newModes{i},sizeWML/numBounds,[]));
%end

end

