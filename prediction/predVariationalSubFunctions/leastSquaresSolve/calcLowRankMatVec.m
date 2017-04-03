for i = 1:numWindows
    idxBB = idx_b_Windowed{i};
    p = xx(idxBB);
	y = (A_k_partial{i}*(WMLWindowed{i}(idxBB,:)'*p) - A_k_nonzero(idxBB,idxBB)*p)./sigma_tilde_squared(idxBB)';
	Ap(idxBB) = WMLWindowed{i}(idxBB,:)*(A_k_partial{i}'*y) - A_k_nonzero(idxBB,idxBB)'*y + sigmaML^-1*p - sigmaML^-1*WMLWindowed{i}(idxBB,:)*MWindowed{i}*(WMLWindowed{i}(idxBB,:)'*p);
end


