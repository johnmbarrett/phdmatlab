function [p,nx,ny] = ejpmf(x,y)
    nx = max(x);
    ny = max(y);
    
    f = zeros(nx,ny);
    n = numel(x);
    
    for ii = 1:n
        f(x(ii),y(ii)) = f(x(ii),y(ii)) + 1;
    end
    
    p = f./sum(sum(f));
end