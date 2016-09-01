function model = trainWFSDecoder(trainIndices,s,R,varargin)
    k = size(R,2);
    m = max(s);
    
    latencyMatrix = zeros(k,k,m);
    
    if islogical(trainIndices)
        trainIndices = find(trainIndices);
    end
    
    for ii = 1:m
        rs = R(s(trainIndices) == ii,:);
        ns = size(rs,1);
        
        Li = repmat(rs,[1 1 k]);
        Lj = repmat(reshape(rs,[ns 1 k]),[1 k 1]);
        
        latencyMatrix(:,:,ii) = sum(Li > Lj,1)/ns; % TODO : check normalisation
    end
    
    model = struct('latencyMatrix',latencyMatrix);
end