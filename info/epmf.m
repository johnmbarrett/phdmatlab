function p = epmf(x)
    n = max(x);
    f = zeros(n,1);
    
    for ii = 1:numel(x)
        f(x(ii)) = f(x(ii)) + 1;
    end
    
    p = f./sum(f);
end