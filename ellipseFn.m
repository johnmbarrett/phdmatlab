function [x,y] = ellipseFn(x0,y0,a,b,phi)
    if nargin < 5
        phi = 0;
    end
    
    if nargin < 4
        b = 1;
    end
    
    if nargin < 3
        a = 1;
    end
    
    if nargin < 2
        y0 = 0;
    end
    
    if nargin < 1
        x0 = 0;
    end
    
    x = @(t) x0+a*cos(t)*cos(phi)-b*sin(t)*sin(phi);
    y = @(t) y0+a*cos(t)*sin(phi)+b*sin(t)*cos(phi);
end