alpha = 10^-4; beta = 0.1;
for t = logspace(0,2,3)
	AA = t*2*A'*A; 
    fprintf('t: %d\n',t);
    for iter = 1:1000
        funcVal = t*norm((A*x-b))^2 - sum(log(-B*x+delta));
        grad = t*2*A'*(A*x-b) + sum(B.*repmat((delta-B*x).^-1,1,n))';
        Hessian = sparse(AA + B'*sparse(diag((delta-B*x).^-2))*B);

        newtonStep = Hessian\-grad;
        lambdaSq = grad'*-newtonStep;

        if lambdaSq < 10^-1
            break
        end

        xNew = x + newtonStep;
        funcValStep = t*norm((A*xNew-b))^2 - sum(log(-B*xNew+delta));

        step = 1;
        while ~isreal(funcValStep) || funcValStep > funcVal - step*alpha*lambdaSq
            step = beta*step;
            if step < 10^-50
                break
            end
            xNew = x + step*newtonStep;
            funcValStep = t*norm((A*xNew-b))^2 - sum(log(-B*xNew+delta));
        end

        x = xNew;
        funcVal = funcValStep;
    end
end


