function d = infDist(x,y,dim)
    assert(isequal(size(x),size(y)),'X and Y must be the same size');
    
    if nargin < 3
        if isrow(x)
            dim = 2;
        else
            dim = find(ndims(x) > 1,1);
        end
    end
    
    d = x-y;
    d(isnan(d)) = 0;
    d = sqrt(sum(d.^2,dim));
end
    