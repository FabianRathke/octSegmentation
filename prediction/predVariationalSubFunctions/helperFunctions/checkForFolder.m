function folderName = checkForFolder(folderName,numVolRegions)
	if ~exist(folderName,'dir')
		mkdir(folderName);
		mkdir(sprintf('%s/init',folderName));
	else
		fprintf('Delete all content in %s\n',folderName);
		rmdir(folderName,'s');
		mkdir(folderName);
		mkdir(sprintf('%s/init',folderName));
	end
end
