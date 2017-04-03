function checkPredictionGlobal(filename,labelID)

global predictionGlobal
global isStoredInGlobal;
isStoredInGlobal = 0;
% check filename and scanID (for 3D scans) of saved appearance terms
if ~isempty(predictionGlobal)
	if isfield(predictionGlobal,'filename')
		if ~strcmp(predictionGlobal.filename,hash(filename,'MD5'))
			% delete outdated global variable
			clearvars -global predictionGlobal
			% reinitialize
			global predictionGlobal;
		else
			if nargin > 1
				if isfield(predictionGlobal,'data') 
					if length(predictionGlobal.data) < labelID
						isStoredInGlobal = 0;
					elseif ~isempty(predictionGlobal.data{labelID})
						isStoredInGlobal = 1;
					end
				else
					isStoredInGlobal = 0;
				end
			else
				isStoredInGlobal = 1;
			end
		end
	else
		isStoredInGlobal = 0;
	end
end

