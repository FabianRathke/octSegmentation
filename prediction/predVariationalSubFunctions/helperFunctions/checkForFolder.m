function folderName = checkForFolder(folderName,numVolRegions)
	if ~exist(folderName,'dir')
		mkdir(folderName);
		mkdir(sprintf('%s/init',folderName));
	end
end
