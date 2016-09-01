function [F,S,G] = smoothBinaryNoise(M,N,T,sigmaX,sigmaY)
    if nargin < 1
        M = 16;
    end
    
    if nargin < 2
        N = M;
    end
    
    if nargin < 3
        T = 256;
    end

    if nargin < 4
        sigmaX = 1;
    end
    
    if nargin < 5
        sigmaY = sigmaX;
    end
    
    [X,Y] = meshgrid((1:M)-M/2,(1:N)-N/2);
    kernel = exp(-(X.^2/(2*sigmaX^2) + Y.^2/(2*sigmaY^2)));
    
    f = figure;
    set(f,'Position',[100 100 M N]);
    a = axes;
    axis off;
    set(a,'Position',[0 0 1 1]);
    
    if nargout > 1
        S = zeros(M,N,T);
    end
    
    for tt = 1:T
        data = randi(2,2*M-1,2*N-1)-1;
        smoothed = conv2(data,kernel,'valid');
        maxS = max(max(smoothed));
        minS = min(min(smoothed));
        smoothed = 255*(smoothed-minS)/(maxS-minS);
        
        if nargout > 1
            S(:,:,tt) = smoothed;
        end
        
        image(smoothed);
        colormap(gray(256));
        F(tt) = getframe;
    end
    
    close(f)
    
    if nargout > 2
        G = gifify(S);
    end
end