alpha = 10^-4; beta = 0.1;
for t = logspace(0,5,6)
    H = eye(n);
    fprintf('t: %d\n',t);
    grad = t*2*A'*(A*x-b) + sum(B.*repmat((delta-B*x).^-1,1,n))';
    funcVal = t*norm((A*x-b))^2 - sum(log(-B*x+delta));
    for iter = 1:1000
           %Hessian = t*2*A'*A + B'*sparse(diag((delta-B*x).^-2))*B;
%        newtonStep = Hessian\-grad;

        newtonStep = -H*grad;
        lambdaSq = grad'*-newtonStep;

        if lambdaSq < 10^-5
            break
        end

        xNew = x + newtonStep;
        funcValStep = t*norm((A*xNew-b))^2 - sum(log(-B*xNew+delta));

        step = 1;
        while ~isreal(funcValStep) || funcValStep > funcVal - step*alpha*lambdaSq
            step = beta*step;
            if step < 10^-10
                break
            end
            xNew = x + step*newtonStep;
            funcValStep = t*norm((A*xNew-b))^2 - sum(log(-B*xNew+delta));
        end

        x = xNew;
        grad_old = grad;
        grad = t*2*A'*(A*x-b) + sum(B.*repmat((delta-B*x).^-1,1,n))';
        funcVal = funcValStep;

        s = step*newtonStep;
        y = grad - grad_old;
        p = 1/(s'*y);
        if iter == 1
            H = (y'*s/(y'*y))*eye(n);
        end
%       Bs = B*s; sy = s'*y;
%       B = B + y*y'/(sy) - Bs*Bs'/(s'*Bs);
        H = (eye(n)-p*s*y')*H*(eye(n)-p*y*s') + p*s*s';
        fprintf('%d: %.3e, %.3e, %.3e \n',iter,funcValStep,norm((A*xNew-b))^2,lambdaSq);
    end
end

