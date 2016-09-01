function model = trainMBGDecoder(trainIndices,s,R,varargin)
    m = max(s);
        
    k = size(R,2);
    params = zeros(k,m,3); % log q, mu, sigma
    
    for ii = 1:k
        r = R(trainIndices,ii);
        
        for jj = 1:m
            rs = r(s(trainIndices) == jj);
        
            isResponse = isfinite(rs);
            params(ii,jj,1) = sum(~isResponse)/numel(rs);

            rs = rs(isResponse);

            params(ii,jj,2) = mean(rs);
            params(ii,jj,3) = std(rs);
        end
    end
    
    model = struct('params',params); 
end