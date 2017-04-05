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

%for numModes = 1:2
% if numModes > 1, copy the shape numModes times
for numModes = 1
	% coordinates
	interval = floor((numColumnsModes-1)/numModes);
	if length(strfind(collector.options.folder_data,'2015_BOE_Chiu')) || strcmp(collector.options.modeType,'DME') % DME Dataset
		boundariesAffected{1} = 1:5;
		boundariesAffected{2} = 4:5;
		if max(columnsPred > 300)
			boundariesAffected{3} = 1;
			boundariesAffected{4} = 2;
		end
		height = 20;
		modeShape{1} = [linspace(1,1.5,interval+1)*height; linspace(1.5,1.5,interval+1)*height; linspace(1.5,1,interval+1)*height]; % cysts between layers 5 and 6
		modeShape{2} = [linspace(1,1.5,interval+1)*height/5; linspace(1.5,1.5,interval+1)*height/5; linspace(1.5,1,interval+1)*height/5]; % cysts between layers 3 and 4
		if max(columnsPred > 300) % swelling only occures at nasal region
			modeShape{3} = [linspace(1,1.5,interval+1)*height/10; linspace(1.5,1.5,interval+1)*height/10;]; % increased RNFL for DME
			modeShape{4} = -[linspace(1,1.5,interval+1)*height/10; linspace(1.5,1.5,interval+1)*height/10;]; % increased RNFL for DME
		end
		modeType = [1 1 1 1];
	elseif length(strfind(collector.options.folder_data,'drusendaten')) || length(strfind(collector.options.folder_data,'Chiu_IOVS_2011')) || strcmp(collector.options.modeType,'AMD')  % AMD Dataset
		boundariesAffected{1} = 6:9;
%		boundariesAffected{2} = 6; % GA 
%		boundariesAffected{3} = 7; % GA 
		height = (diff(minMaxColumns)+1)/numModes/2;
		modeShape{1} = cos(linspace(-pi/2,pi/2,interval+1))*height;
%		modeShape{2} = [linspace(1,1.5,interval+1)*4; linspace(1.5,1.5,interval+1)*4; linspace(1.5,1,interval+1)*4]; % cysts between layers 5 and 6
%		modeShape{3} = [linspace(1,1.5,interval+1)*4; linspace(1.5,1.5,interval+1)*4; linspace(1.5,1,interval+1)*4]; % cysts between layers 5 and 6
%		modeShape{2} = [linspace(1.5,1.5,interval+1)*height/5];
		I1 = round((interval+1)/3); I2 = interval+1-I1;
%		modeShape{2} = [[linspace(1,0,I1)*height/10 linspace(0,0,I2)*height/5]; [linspace(0,0,I2)*height/10 linspace(0,1,I1)*height/10]]; % cysts between layers 3 and 4
%		modeShape{3} = cos(linspace(-pi/2,pi/2,interval+1))*height;
%		modeShape{1} = [cos(linspace(-pi/2,pi/2,interval+1))*height; cos(linspace(-pi/2,pi/2,interval+1))*height; cos(linspace(-pi/2,pi/2,interval+1))*height-5; cos(linspace(-pi/2,pi/2,interval+1))*height-5];
		modeType = [1 1 1];
	elseif length(strfind(collector.options.folder_data,'glaukom')) % Glaukom dataset
		boundariesAffected{1} = 1;
%		boundariesAffected{2} = 2:4;
%		boundariesAffected{3} = 2:4;
%		boundariesAffected{2} = 2;
		boundariesAffected{2} = 3;
		boundariesAffected{3} = 4;
%		boundariesAffected{5} = 5;
		height = 4;
		modeShape{1} = [linspace(1.5,1.5,interval+1)*height]; % cysts between layers 5 and 6
		modeShape{2} = [linspace(1.5,1.5,interval+1)*height/4]; % cysts between layers 5 and 6
		modeShape{3} = [linspace(1.5,1.5,interval+1)*height/4]; % cysts between layers 5 and 6
%		modeShape{4} = [linspace(1.5,1.5,interval+1)*height/2]; % cysts between layers 5 and 6
%		modeShape{5} = [linspace(1.5,1.5,interval+1)*height/2]; % cysts between layers 5 and 6
%		modeShape{1} = [linspace(1,1.5,interval+1)*height; linspace(1.5,1.5,interval+1)*height; linspace(1.5,1,interval+1)*height]; % cysts between layers 5 and 6
%		modeShape{2} = [linspace(3,3,interval+1)*height/4; linspace(2,2,interval+1)*height/4; linspace(1,1,interval+1)*height/4]; % cysts between layers 5 and 6
%		modeShape{2} = [linspace(2,2,interval+1)*height/4; linspace(2,2,interval+1)*height/4; linspace(2,2,interval+1)*height/4]; % cysts between layers 5 and 6
		modeType = [1 1 1 1 1];
	else
		error('Dataset not detected --> no modes added');
	end

	for k = 1:length(modeShape)
		for i = 1:numModes
			N = length(boundariesAffected{k});
			if modeType(k)
				for j = 1:size(modeShape{k},1)
					newMode = zeros(sizeWML,1); % initial zero mode
					columnsMode = columnsModes((1:(interval+1)) + (i-1)*interval); % determine column indices that correspond to mode 

					% add the shape in modeShape to *each* boundaries in boundariesAffected: 1/2 a/b ==> 1: a/b, 2: a/b
					columnsToAdd = repmat(columnsMode,N,1) + repmat((boundariesAffected{k}-1)'*sizeWML/numBounds,1,interval+1);
					newMode(columnsToAdd) = repmat(modeShape{k}(j,:),N,1);
					% add to cell array of modes
					newModes(:,modeCounter) = newMode;
					modeCounter = modeCounter+1;
				end
			else
			    newMode = zeros(sizeWML,1); % initial zero mode
				columnsMode = columnsModes((1:(interval+1)) + (i-1)*interval); % determine column indices that correspond to mode 

				% add the shape in modeShape to *all* boundaries in boundariesAffected: 1/2 a/b ==> 1: a, 2: b 
				columnsToAdd = repmat(columnsMode,N,1) + repmat((boundariesAffected{k}-1)'*sizeWML/numBounds,1,interval+1);
				newMode(columnsToAdd) = modeShape{k};
				% add to cell array of modes
				newModes(:,modeCounter) = newMode;
				modeCounter = modeCounter+1;
			end
		end
	end
end

% DEBUG: plot all modes
%figure; hold on;
%for i = 1:length(newModes)
%	plot(reshape(newModes{i},sizeWML/numBounds,[]));
%end

end
