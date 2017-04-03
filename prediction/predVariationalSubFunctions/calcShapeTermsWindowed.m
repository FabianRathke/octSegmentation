% calc $A^j_k = (K_{jj}^-1*K_{j,\setminus{j}})_{k,*}$
if collector.options.printTimings
	ticCalcShape = tic;
end

% calc M for each window
numWindows = length(options.windowSize);
MWindowed = cell(1,numWindows);
idx_b_Windowed = cell(1,numWindows);
for i = 1:numWindows
	idx_b_Windowed{i} = reshape([repmat(find(columnsShapePred{1} >= options.windowSize{i}(1) & columnsShapePred{1} <= options.windowSize{i}(2)),numBounds,1) + repmat([0:numColumnsShape:(numBounds-1)*numColumnsShape]',1,sum(columnsShapePred{1} >= options.windowSize{i}(1) & columnsShapePred{1} <= options.windowSize{i}(2)))]',1,[]);
	if isempty(options.windowModes{i})
%		WMLWindowed{i} = WML(idx_b_Windowed{i},:);
		WMLWindowed{i} = models.shapeModel.WML;
	else
%		WMLWindowed{i} = [WML(idx_b_Windowed{i},:) options.windowModes{i}(idx_b_Windowed{i},:)];
		WMLWindowed{i} = [models.shapeModel.WML options.windowModes{i}];
%		WMLWindowed{i} = WML;
	end

	MWindowed{i} = inv(WMLWindowed{i}(idx_b_Windowed{i},:)'*WMLWindowed{i}(idx_b_Windowed{i},:) + sigmaML*eye(size(WMLWindowed{i},2)));
	prodWMWindowed{i} = WMLWindowed{i}*MWindowed{i};
	prodWMTWindowed{i} = prodWMWindowed{i}';
end
windowLength = cellfun(@length,idx_b_Windowed)/numBounds;
A_k = eval(sprintf('zeros(numBounds*numColumnsShapeTotal,%s);',collector.options.dataType));
sigma_tilde_squared = zeros(1,numBounds*numColumnsShapeTotal);
idx_diag = 1:numBounds+1:numBounds^2;
s1 = zeros(1,numColumnsShapeTotal*numBounds^2); s2 = zeros(1,numColumnsShapeTotal*numBounds^2);

% naive implementation for each image column
for volRegion = 1:numVolRegions
	for i = 1:numWindows
		% 
		idxTmp = find(columnsShapePred{1} >= options.windowSize{i}(1) & columnsShapePred{1} <= options.windowSize{i}(2));
		for j = idxTmp
			idx_j = (j:numColumnsShape(volRegion):numColumnsShape(volRegion)*numBounds) + sum(numColumnsShape(1:volRegion-1))*numBounds;
			idx_not_j = idx_b_Windowed{i}; idx_not_j(ismember(idx_not_j,idx_j)) = [];
		
			% use the CPU version of WML and M explicitly
			tmp = eye(numBounds)*sigmaML^-1 - sigmaML^-1*WMLWindowed{i}(idx_j,:)*MWindowed{i}*WMLWindowed{i}(idx_j,:)';
			K_jj_inverse = inv(tmp);
			sigma_tilde_squared(idx_j) = K_jj_inverse(idx_diag)/options.alpha;
			
			countJ = sum(numColumnsShape(1:volRegion-1)) + j;
			idxRange = (countJ-1)*numBounds^2 + 1:countJ*numBounds^2;

	        [idxI(idxRange) idxJ(idxRange)] = meshgrid(idx_j);
			s1(idxRange) = K_jj_inverse(:);
			s3(idxRange) = tmp(:);

			% calculating A_k^j: note that we do not need the sigma^-2I part, since we pull a part of K that does not include the diagonal
			A_k(idx_j,idx_not_j) = K_jj_inverse*-sigmaML^-1*WMLWindowed{i}(idx_j,:)*(MWindowed{i}*WMLWindowed{i}(idx_not_j,:)');
%			s2(idxRange) = reshape(K_jj_inverse*-sigmaML^-1*WMLWindowed{i}(idx_j,:)*MWindowed{i}*WMLWindowed{i}(idx_j,:)',1,[]);
		end
	end

end

K_jj_inverse_block = sparse(idxI,idxJ,s1);
K_jj_block = sparse(idxI,idxJ,s3);
%A_k_nonzero = sparse(idxI,idxJ,s2);
%for i = 1:numWindows
%	A_k_partial{i} = -sigmaML^-1*K_jj_inverse_block(idx_b_Windowed{i},idx_b_Windowed{i})*WMLWindowed{i}(idx_b_Windowed{i},:)*MWindowed{i};
%end

% cast to different data type if required
factor = (1./sigma_tilde_squared');

%K = zeros(size(P_mu));
%for i = 1:numWindows
%   end
%A = (P_mu + K);
%for t = logspace(0,2,3)
%	AAt{log10(t*10)} = AA*t;
%end
for i = 1:numWindows
	K = sigmaML^-1*eye(length(idx_b_Windowed{i})) - sigmaML^-1*WMLWindowed{i}(idx_b_Windowed{i},:)*MWindowed{i}*WMLWindowed{i}(idx_b_Windowed{i},:)';
	P_mu = (A_k(idx_b_Windowed{i},idx_b_Windowed{i}).*factor(idx_b_Windowed{i},ones(1,length(idx_b_Windowed{i}))))'*A_k(idx_b_Windowed{i},idx_b_Windowed{i});
	AWindowed{i} = P_mu+K; 
	AAWindowed{i} = 2*AWindowed{i}*AWindowed{i};
end

% constraint matrix B
B = zeros(numBounds*(numWindows-1)*2,numColumnsShape*numBounds);
% for all transitions between windows, impose smooth boundaries
for i = 1:numWindows-1
    l1 = sum(windowLength(1:i));
    for j = 1:numBounds
        B((i-1)*numBounds*2 + (j-1)*2+1,[l1 l1+1]+(j-1)*numColumnsShape) = [1 -1];
        B((i-1)*numBounds*2 + j*2,[l1 l1+1]+(j-1)*numColumnsShape) = [-1 1];
    end
end
B = sparse(B);


if collector.options.printTimings
	if collector.options.calcOnGPU
		GPUsync;
	end
	fprintf('[Initialized Shape Terms]: %.3fs \n',toc(ticCalcShape));
end
