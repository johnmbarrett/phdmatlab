function [F, S, C, video] = gaussNoise(x, y, t, sigma, filename, fps, format)
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

    minC = zeros(x,y);
    maxC = 255*ones(x,y);
    
    f = figure;
    set(f,'Position',[100 100 x y])
    a = axes;
    axis off;
    set(a,'Position',[0 0 1 1])
    
    if nargout > 1
        S = zeros(x,y,t);
    end
    
    if nargout > 2
        C = zeros(x,y,t);
    end
    
    if nargin > 4
        if nargin < 7
            format = 'Archival';
        end
        
        video = VideoWriter(filename,format);
        
        if nargin > 5
            video.FrameRate = fps;
        end
        
        open(video);    
    end
    
    for ii = 1:t
        s = normrnd(0,sigma,x,y);
        %s = s/(2*maxZ*sigma) + 0.5;
        
        if nargout > 1
            S(:,:,t) = s;
        end
    
        c = 255*s+128;
        c = min(max(c,minC),maxC);
        
        if nargout > 2
            C(:,:,t) = c;
        end
        
        image(c);
        colormap(gray(256))
        frame = getframe;
        
        if nargout > 0
            F(ii) = frame;
        end
        
        if nargin > 4
            writeVideo(video,frame);
        end
    end
    
    close(f);
    
    if nargin < 5
        return
    end
    
    close(video)
end