function r = generalChR2RGCLNModel(I,x,y,s,p,n)
    assert(nargin >= 1,'You must provide at the very least an input stimulus');
    
    if nargin < 6
        n = @(x,p) p(1)./(1+exp(-(x-p(3))/p(4)))+p(2);
    end
    
    if nargin < 5
        p = [1 0 0 1];
    end
    
    if nargin < 4
        s = 60;
    end
    
    szI = size(I);
    
    if nargin < 3
        y = szI(1)/2;
    end
    
    if nargin < 2
        x = szI(2)/2;
    end
    
    [Y,X] = ndgrid((1:szI(1))-y,(1:szI(2))-x);
    
    k = gauss2d(X,Y,s,s,0,true);
    
    L = k.*I;
    
    R = n(L,p);
    
    r = mean(R(:));
    
    r = max(0,r);
end