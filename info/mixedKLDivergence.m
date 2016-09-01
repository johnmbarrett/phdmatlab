function K = mixedKLDivergence(P,Q,degenerate,k,method,informativeSingletons)
    if nargin < 6
        informativeSingletons = true;
    end
    
    if ~iscell(P)
        P = {P};
    end
    
    if ~iscell(Q)
        Q = {Q};
    end
    
    nVars = numel(P);
    
    assert(numel(Q) == nVars);
    assert(numel(degenerate) == nVars);
    
    P = reshape(P,1,nVars);
    Q = reshape(Q,1,nVars);
    
    Kcontinuous = zeros(1,nVars);
    
    for ii = setdiff(1:nVars,find(logical(degenerate)))
        Kcontinuous(ii) = continuousKLDivergence(P{ii},Q{ii},k,method,informativeSingletons);
    end
    
    rowSizeFun = @(A) size(A,1);
    
    n = cellfun(rowSizeFun,P);
    p = n/sum(n);
    
    m = cellfun(rowSizeFun,Q);
    q = m/sum(m);
    
    pNonZero = p > 0;
    
    if any(pNonZero & q == 0)
        K = Inf;
        return;
    end
    
    Kdiscrete = p(pNonZero).*(log(p(pNonZero)) - log(q(pNonZero)));
    
    K = sum(Kdiscrete) + sum(p.*Kcontinuous);
end