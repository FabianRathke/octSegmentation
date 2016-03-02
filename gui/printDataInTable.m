function t = printDataInTable(tab,data)

absoluteMax = 4;
fieldNames = fieldnames(data);
maxLength = 1;

for i = 1:length(fieldNames)
	tmp = getfield(data,fieldNames{i});
	if isnumeric(tmp)
		if length(tmp) > maxLength
			maxLength = length(tmp);
		end
	end
end

maxLengthString = 1;
cellArray = cell(length(fieldNames),min(maxLength,absoluteMax+1));
for i = 1:length(fieldNames)
	field = getfield(data,fieldNames{i});
	if isnumeric(field)
		for j = 1:length(field)
			if (j > absoluteMax)
				cellArray{i,j} = '...';
				break
			else
				cellArray{i,j} = field(j);
			end
		end
	elseif ischar(field)
		cellArray{i,1} = field;
	elseif iscell(field)
		cellArray{i,1} = 'cell array';
	end
end
	
t = uitable('Parent',tab,'Data', cellArray, 'RowName', fieldnames(data),'ColumnWidth',{120});


