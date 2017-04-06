function checkPredictionGlobal(filename,labelID)

if nargin < 2
	error('Not enough input arguments provided');
end

global predictionGlobal
global isStoredInGlobal;
isStoredInGlobal.prediction = 0;
isStoredInGlobal.data = 0;
% check filename and scanID (for 3D scans) of saved appearance terms
if ~isempty(predictionGlobal)
	if isfield(predictionGlobal,'filename')
		if ~strcmp(predictionGlobal.filename,hash(filename,'MD5'))
			% delete outdated global variable
			clearvars -global predictionGlobal
			% reinitialize
			global predictionGlobal;
		else
			if isfield(predictionGlobal,'BScans') & ~isempty(predictionGlobal.BScans{labelID})
				isStoredInGlobal.data = 1
			end
			if isfield(predictionGlobal,'data') 
				if length(predictionGlobal.data) < labelID
					isStoredInGlobal.prediction = 0;
				elseif ~isempty(predictionGlobal.data{labelID})
					isStoredInGlobal.prediction = 1;
				end
			else
				isStoredInGlobal.prediction = 0;
			end
		end
	end
end

