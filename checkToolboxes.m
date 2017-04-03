% script that checks for required toolboxes in all files

% test function from the statistics toolbox
mvnrnd(0,1);

cd(getenv('OCT_CODE_DIR'));

saveToolboxes = {};
D = subdir('*.m');
counter = 1;
for i=1:length(D)
	[folder filename] = fileparts(D(i).name);
	cd(folder);
	P = fdep([filename '.m']);
	for  j=1:length(P.toolbox)
		if sum(ismember(saveToolboxes,P.toolbox(j))) == 0
			saveToolboxes(counter) = P.toolbox(j);
			counter = counter + 1;
		end
	end
end
