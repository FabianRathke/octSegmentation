function makeQualityOverlay(prediction,qualityData,mu,sigma,valX)

valInter = [valX-1; valX+1];
for i = 1:9
    interPData(1,:,i) = interp1(valX,prediction(i,:),valInter(1,:),'pchip');
    interPData(2,:,i) = interp1(valX,prediction(i,:),valInter(2,:),'pchip');
end

ranges = [4 8 12 16 20];
numR = length(ranges);
colors = [ones(numR-1,1) linspace(1,0,numR-1)' zeros(numR-1,1)];

for i = 2:numR
	for k = 1:9
		if i == numR
		    idx = find(qualityData(:,k) > mu(k)+sigma(k)*ranges(i));
		else
	    	idx = find(qualityData(:,k) > mu(k)+sigma(k)*ranges(i-1) & qualityData(:,k) < mu(k)+sigma(k)*ranges(i));
		end
		for j = 1:length(idx)
			plot([valInter(1,idx(j)) valX(idx(j)) valInter(2,idx(j))],[interPData(1,idx(j),k) prediction(k,idx(j)) interPData(2,idx(j),k)],'LineWidth',1.5,'Color',colors(i-1,:));
	   	end
	end
end

%for k = 1:9
%    idx = find(qualityData(:,k) > mu(k)+sigma(k)*8);
%    for j = 1:length(idx)
%        plot([valInter(1,idx(j)) valX(idx(j)) valInter(2,idx(j))],[interPData(1,idx(j),k) prediction(k,idx(j)) interPData(2,idx(j),k)],'LineWidth',1.5,'Color','red');
%    end
%    idx = find(qualityData(:,k) > mu(k)+sigma(k)*4 & qualityData(:,k) < mu(k)+sigma(k)*8);
%    for j = 1:length(idx)
%        plot([valInter(1,idx(j)) valX(idx(j)) valInter(2,idx(j))],[interPData(1,idx(j),k) prediction(k,idx(j)) interPData(2,idx(j),k)],'LineWidth',1.5,'Color','yellow');
%    end
%end

