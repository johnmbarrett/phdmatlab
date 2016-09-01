function [F, S, C] = gaussNoise(x, y, t, sigma, maxZ, filename, fps, format)
    if nargin < 1
        x = 256;
    end
    
    if nargin < 2
        y = 256;
    end
    
    if nargin < 3
        t = 300;
    end
    
    if nargin < 4
        sigma = 1;
    end
    
    if nargin < 5
        maxZ = 3;
    end

    S = normrnd(0,sigma,x,y,t);
    S = S/(2*maxZ*sigma) + 0.5;
    Z = zeros(x,y,t);
    I = ones(x,y,t);
    S = min(max(S,Z),I);
    
    clear I;
    clear Z;
    
    C = 255*S;
    
    f = figure;
    set(f,'Position',[100 100 x y])
    a = axes;
    axis off;
    set(a,'Position',[0 0 1 1])
    
    for ii = 1:t
        image(C(:,:,ii));
        colormap('gray')
        F(ii) = getframe;
    end
    
    close(f);
end