function r = rectifyingChR2RGCLNModel(I,x,y,s,a,b,t,g)
    assert(nargin >= 1,'You must provide at the very least an input stimulus');
    
    if nargin < 5
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
    
    R = a*log(1+exp(L-t)).^g+b;
    
    r = mean(R(:));
end