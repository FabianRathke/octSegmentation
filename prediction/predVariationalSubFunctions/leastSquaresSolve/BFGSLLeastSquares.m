alpha = 10^-4; beta = 0.1;
m = 10; % number of previous gradients that are used to calculate the inverse Hessian
s_k = zeros(n,m);
y_k = zeros(n,m);
sy = zeros(1,m); syInv = zeros(1,m);
tmp = zeros(size(b));
for t = logspace(0,2,3)
	%fprintf('t: %d\n',t);
	Abt = Ab*t;
	for j = 1:numWindows
		tmp(idx_b_Windowed{j}) = t*(AAWindowed{j}*x(idx_b_Windowed{j}));
	end
	grad = tmp + Abt + ((1./(delta-B*x))'*B)';
	newtonStep = -grad;
	for j = 1:numWindows
		tmp(idx_b_Windowed{j}) = AWindowed{j}*x(idx_b_Windowed{j});
	end
	funcVal = t*norm((tmp-b))^2 - sum(log(-B*x+delta));
    for iter = 1:1000
        lambdaSq = grad'*-newtonStep;

        if lambdaSq < 10^-0
            break
        end

        xNew = x + newtonStep;
		for j = 1:numWindows
			tmp(idx_b_Windowed{j}) = AWindowed{j}*xNew(idx_b_Windowed{j});
		end
        %funcValStep = t*norm((A*xNew-b))^2 - sum(log(-B*xNew+delta));
        funcValStep = t*norm((tmp-b))^2 - sum(log(-B*xNew+delta));

        step = 1;
        while ~isreal(funcValStep) || funcValStep > funcVal - step*alpha*lambdaSq
            step = beta*step;
            if step < 10^-10
                break
            end
            xNew = x + step*newtonStep;
			for j = 1:numWindows
				tmp(idx_b_Windowed{j}) = AWindowed{j}*xNew(idx_b_Windowed{j});
			end
	        funcValStep = t*norm((tmp-b))^2 - sum(log(-B*xNew+delta));
        end

        funcVal = funcValStep;
		x = xNew;
        grad_old = grad;
%        grad = t*2*A'*(A*x-b) + sum(B.*repmat((delta-B*x).^-1,1,n))';
    	for j = 1:numWindows
        	tmp(idx_b_Windowed{j}) = t*(AAWindowed{j}*x(idx_b_Windowed{j}));
        end
		grad = tmp + Abt + ((1./(delta-B*x))'*B)';
    	
		% move and delete entries
		s_k(:,2:end) = s_k(:,1:end-1); y_k(:,2:end) = y_k(:,1:end-1); sy(2:end) = sy(1:end-1); syInv(2:end) = syInv(1:end-1);
		s_k(:,1) = step*newtonStep; 
	 	y_k(:,1) = grad-grad_old; 
		sy(1) = s_k(:,1)'*y_k(:,1); syInv(1) = 1/sy(1);
	
		% choose H0
		gamma_k = sy(1)/(y_k(:,1)'*y_k(:,1));
		H0 =  gamma_k;

		% first for-loop
		q = grad;
		for l = 1:min([m,iter,length(b)]);
			alphaBFGS(l) = s_k(:,l)'*q*syInv(l);
			q = q - alphaBFGS(l)*y_k(:,l);
		end

		r = H0.*q;
		% second for-loop
		for l = min([m,iter,length(b)]):-1:1
			betaBFGS = y_k(:,l)'*r*syInv(l);
			r = r + s_k(:,l)*(alphaBFGS(l)-betaBFGS);
		end
		newtonStep = -r;
    end
	%fprintf('%d: %.3e, %.3e, %.3e \n',iter,funcValStep,norm((A*xNew-b))^2,lambdaSq);
end

