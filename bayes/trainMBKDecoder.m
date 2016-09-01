function model = trainMBKDecoder(trainIndices,s,R,params)
    n = 2^12;
    k = size(R,2);
    
    if n == k % this would break makeSubDecoder but is also extremely unlikely to happen
        n = 2^13;
    end
    
    m = max(s);
        
    prks = zeros(n,k,m);
    qrks = zeros(k,m); % log q, mu, sigma
    bw = params.bw;
    lb = params.lb;
    ub = params.ub;
    x = params.x;
    
    for ii = 1:k
        r = R(trainIndices,ii);
        
        for jj = 1:m
            rs = r(s(trainIndices) == jj);
        
            isResponse = isfinite(rs);
            qrks(ii,jj) = sum(~isResponse)/numel(rs);
            
            prks(:,ii,jj) = kernelPDFEstimate(rs(isResponse),x,bw(ii,jj),lb,ub);
        end
    end
    
    model = struct('bw',bw,'prks',prks,'qrks',qrks,'x',x); 
end