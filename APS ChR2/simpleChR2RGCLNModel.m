function r = simpleChR2RGCLNModel(I,x,y,s,a,b,t,g)
    assert(nargin >= 1,'You must provide at the very least an input stimulus');
    
    if nargin < 8
        g = 1;
    end
    
    if nargin < 7
        t = 0;
    end
    
    if nargin < 6
        b = 0;
    end
    
    if nargin < 5
        a = 1;
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
    
    n = @(l) a*(1./(1+exp(-(l-t)/g)))+b;
    
    R = n(L);
    
    r = mean(R(:));
    
    r = max(0,r);
end