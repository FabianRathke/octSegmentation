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
boundariesAffected = 1:5; N = length(boundariesAffected);

%for numModes = 1:2
% if numModes > 1, copy the shape numModes times
for numModes = 1
	% coordinates
	interval = floor((numColumnsModes-1)/numModes);
	height = (diff(minMaxColumns)+1)/numModes/2;
%	modeShape = cos(linspace(-pi/2,pi/2,interval+1))*height;
	modeShape = [linspace(1,1.5,interval+1)*height; linspace(1.5,1.5,interval+1)*height; linspace(1.5,1,interval+1)*height];
	
	for i = 1:numModes
		for j = 1:size(modeShape,1)
			newMode = zeros(1,sizeWML); % initial zero mode
			columnsMode = columnsModes((1:(interval+1)) + (i-1)*interval); % determine column indices that correspond to mode 

			% add basic shape to all boundaries, determined by boundariesAffected
			columnsToAdd = repmat(columnsMode,N,1) + repmat((boundariesAffected-1)'*sizeWML/numBounds,1,interval+1);

			newMode(columnsToAdd) = repmat(modeShape(j,:),N,1);
			% add to cell array of modes
			newModes{modeCounter} = newMode; 
			modeCounter = modeCounter+1;
		end
	end
end

% DEBUG: plot all modes
%figure; hold on;
%for i = 1:length(newModes)
%	plot(reshape(newModes{i},sizeWML/numBounds,[]));
%end

end

