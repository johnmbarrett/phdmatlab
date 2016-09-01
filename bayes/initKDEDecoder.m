function params = initKDEDecoder(s,R,lb,ub)
    assert(all(R(:) >= lb),'Empirical response probability must be zero below the support of the response distribution');
%     assert(all(R(:) <= ub),'Empirical response probability must be zero above the support of the response distribution');
        
    n = 2^12;
    k = size(R,2);
    
    if n == k % this would break makeSubDecoder but is also extremely unlikely to happen
        n = 2^13;
    end
    
    m = max(s);
        
    bw = zeros(k,m);
    
    for ii = 1:k
        r = R(:,ii);
        
        for jj = 1:m
            rs = r(s == jj);
        
            [bw(ii,jj),~,xmesh] = kde(rs(isfinite(rs)),n,lb,ub);
            
            if ii == 1 && jj == 1 || (~all(isfinite(x)) && all(isfinite(xmesh)))
                x = xmesh;
            end
        end
    end
    
    params = struct('bw',bw,'lb',lb,'ub',ub,'x',x);
end