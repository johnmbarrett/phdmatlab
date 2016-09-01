function [ hs, pk, k ] = fdrcorrect( ps, alpha )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    if nargin < 2
        alpha = 0.05;
    end
    
    m = numel(ps);
    qs = sort(reshape(ps,1,numel(ps)));
    as = alpha*(1:m)/m;
    k = find(qs <= as,1,'last');
    
    if ~isempty(k)
        pk = qs(k);
        hs = ps <= pk;
        return;
    end
    
    k = 0;
    pk = 0;
    hs = false(size(ps));
end

