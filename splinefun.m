function fun = splinefun(P,k)
    if isvector(P) && ismatrix(k)
        temp = P;
        P = k;
        k = temp;
    end

    k = reshape(k,numel(k),1);
    k = sort(k);
    
    assert(ismatrix(P),'Constitutent polynomials must be a specified as a matrix');
    
    nKnots = numel(k);
    assert(ismember(nKnots-1,size(P)),'There must be one fewer polynomial than there are knots');
    
    knotIndex = find(ismember(size(P),nKnots-1));
    degreeIndex = 3-knotIndex;
    P = permute(P,[degreeIndex,knotIndex]);
    
    function y = fun_(x)
        y = nan(size(x));
        
        for ii = 1:nKnots-1
            p = P(:,ii);
            
            if ii == 1
                inP = x < k(2);
            elseif ii == nKnots-1
                inP = x >= k(nKnots-1);
            else
                inP = k(ii) <= x & x < k(ii+1);
            end
            
            y(inP) = polyval(p,x(inP) - k(ii));
        end
    end

    fun = @fun_;
end