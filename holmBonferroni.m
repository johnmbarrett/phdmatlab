function [hs,k,as] = holmBonferroni(ps,alpha)
    if nargin < 2
        alpha = 0.05;
    end
    
    m = numel(ps);
    as = alpha./(m+1-(1:m))';
    
    [qs,sortIndices] = sort(ps(:));
    
    k = find(qs > as,1);
    
    hs = [ones(k-1,1); zeros(m-k+1,1)];
    
    assert(numel(hs) == numel(qs));
    
    hs(sortIndices) = hs;
    as(sortIndices) = as;
end