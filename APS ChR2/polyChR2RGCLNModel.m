function r = polyChR2RGCLNModel(I,x,y,s,p)
    assert(nargin >= 1,'You must provide at the very least an input stimulus');
    
    if nargin < 5
        p = [1 0];
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
    
    R = max(0,polyval(p,L));
    
    r = mean(R(:));
end