function logLikelihood = getWFSLikelihood(response,model,varargin)
    latencyMatrix = model.latencyMatrix;
    
    k = size(response,2);
    
    ri = permute(repmat(response,[1 1 k]),[2 3 1]);
    rj = permute(repmat(reshape(response,[1 1 k]),[1 k 1]),[2 3 1]);
    
    % p(r_i,j|s)    =   p(r_i > r_j|s)  if r_i >  r_j
    %               = 1-p(r_i > r_j|s)  if r_i <= r_j
    mask = repmat(2*(ri > rj)-1,[1 1 8]);
    shift = repmat(ri <= rj,[1 1 8]);
    
    likelihood = mask.*latencyMatrix+shift;
    logLikelihood = permute(mean(log(likelihood),2),[1 3 2]);
end