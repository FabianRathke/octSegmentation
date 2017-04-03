function makeQualityOverlay(prediction,qualityData,mu,sigma,valX)

if sum(size(qualityData) == size(prediction)) > 0
	qualityData = qualityData';
end

diff = valX(2)-valX(1);
valInter = [valX-diff/2; valX+diff/2];
for i = 1:9
    interPData(1,:,i) = interp1(valX,prediction(i,:),valInter(1,:),'pchip');
    interPData(2,:,i) = interp1(valX,prediction(i,:),valInter(2,:),'pchip');
end

%ranges = [4 8 12 16 20];
ranges = [3:3:15];
numR = length(ranges);
colors = [ones(numR-1,1) linspace(1,0,numR-1)' zeros(numR-1,1)];

mask = prediction==0;

for i = 2:numR
	for k = 1:9
		starts = find(([~mask(k,1:end-1) 0] -[0 ~mask(k,1:end-1)])==1);
		stops = find(([~mask(k,1:end-1) 0] -[0 ~mask(k,1:end-1)])==-1);

		for j = 1:length(starts)
			if i == numR
				idx = find(qualityData(starts(j):stops(j)-1,k) > mu(k)+sigma(k)*ranges(i));
			else
				idx = find(qualityData(starts(j):stops(j)-1,k) > mu(k)+sigma(k)*ranges(i-1) & qualityData(starts(j):stops(j)-1,k) < mu(k)+sigma(k)*ranges(i));
			end
			for j = 1:length(idx)
				if j == length(idx) && idx(j) ~= length(valX)
					plot([valInter(1,idx(j)) valX(idx(j))],[interPData(1,idx(j),k) prediction(k,idx(j))],'LineWidth',1.5,'Color',colors(i-1,:));
				else
					plot([valInter(1,idx(j)) valX(idx(j)) valInter(2,idx(j))],[interPData(1,idx(j),k) prediction(k,idx(j)) interPData(2,idx(j),k)],'LineWidth',1.5,'Color',colors(i-1,:));
				end
			end
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

