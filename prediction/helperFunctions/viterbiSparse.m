function maxState = viterbiSparse(pStart,pTrans,pObs)

warning off

[K N] = size(pObs);
w = zeros(K,N); 

% initialize w
w(:,1) = pStart.*pObs(:,1)'; 
save_max = zeros(K,N);

% do the forward loop
for i = 2:N
    [val idx] = max(w(:,(i-1)*ones(1,K)).*pTrans{i},[],1);
    w(:,i) = val'.*pObs(:,i);
	save_max(:,i) = idx;
end

% backtrack, find the sequence corresponding to the hightest probability
maxState = zeros(1,N);
[tmp,ind] = max(w(:,N)); maxState(N) = ind(1);

for j = N-1:-1:1
    %[tmp,ind] = max(log_w(:,j)+log_pTrans(:,maxState(j+1),i)); 
    maxState(j) = save_max(maxState(j+1),j+1); 
end

warning on
