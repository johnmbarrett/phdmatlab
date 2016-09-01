function p = randpartition(n,m,c)
    if nargin < 2
        m = 1;
    end
    
    if nargin < 3
        c = 0;
    end
    
    if numel(n) == 1
        p = rand(n,1);
    else
        p = rand(n);
    end
    
    p = m*p/sum(p(:))+c;
end