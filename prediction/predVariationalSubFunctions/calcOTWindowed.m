if collector.options.printTimings
	ticCalcOT = tic;
end

% x was calculated in optQB
mu_a_b = zeros(numColumnsShape*numBounds,1);
for i = 1:numWindows
	productB = sigmaML^-1*x(idx_b_Windowed{i}) - sigmaML^-1*WMLWindowed{i}(idx_b_Windowed{i},:)*(prodWMTWindowed{i}(:,idx_b_Windowed{i})*x(idx_b_Windowed{i})) - K_jj_block(idx_b_Windowed{i},idx_b_Windowed{i})*x(idx_b_Windowed{i});
	mu_a_b(idx_b_Windowed{i}) = models.shapeModel.mu(idx_b_Windowed{i}) - K_jj_inverse_block(idx_b_Windowed{i},idx_b_Windowed{i})*productB;
end

[condQB cmin cmax] = calcCondQB(mu_a_b,sigma_tilde_squared,numRows,numBounds,numColumnsShape,numColumnsPred,int32(columnsPredShapeVec(1,:)-1),columnsPredShapeFactorVec(1,:),int32(columnsPredShapeVec(2,:)-1),columnsPredShapeFactorVec(2,:),options.thresholdAccuracy);

%tmp = X - mu_a_b(:,ones(1,numRows));
%condQB = repmat(1./sqrt(2*pi*sigma_tilde_squared),numRows,1).*exp((-0.5)*tmp.*tmp./sigma_tilde_squared(ones(1,numRows),:)')';
% implicit cast to single precision
%condQB(condQB < options.thresholdAccuracy) = 0;

if collector.options.printTimings
   	if collector.options.calcOnGPU
   		GPUsync;
	end
	fprintf('[calcOT]: %.3f s\n',toc(ticCalcOT));
end
