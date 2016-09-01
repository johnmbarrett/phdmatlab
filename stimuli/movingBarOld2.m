function [pixels,extraParams] = movingBar(T,X,Y,window,extraParams)
    if ~isfield(extraParams,'skipArgCheck')
        if ~isstruct(extraParams)
            extraParams = struct();
        end
        
        if ~isfield(extraParams,'angle') || ~isnumeric(extraParams.angle) || ~isequal(size(extraParams.angle),[1 1]) || extraParams.angle < 0 || extraParams.angle > 2*pi
            extraParams.angle = 0;
        end
        
        if ~isfield(extraParams,'colour') || ~isnumeric(extraParams.colour) || ~(isequal(size(extraParams.colour),[1 1]) || isequal(size(extraParams.colour),[1 3])) || min(extraParams.colour < 0) || max(extraParams.colour > 255)
            extraParams.colour = 255;
        end
        
        if ~isfield(extraParams,'length') || ~isnumeric(extraParams.length) || ~isequal(size(extraParams.length),[1 1]) || extraParams.length < 1
            extraParams.length = abs(X*sin(extraParams.angle)) + abs(Y*cos(extraParams.angle));
        end
        
        if ~isfield(extraParams,'maxT') || ~isnumeric(extraParams.maxT) || ~isequal(size(extraParams.maxT),[1 1]) || extraParams.maxT < 1
            extraParams.maxT = Inf;
        end
        
        if ~isfield(extraParams,'originX') || ~isnumeric(extraParams.originX) || ~isequal(size(extraParams.originX),[1 1]) || extraParams.originX < 1
            extraParams.originX = X/2;
        end
        
        if ~isfield(extraParams,'originY') || ~isnumeric(extraParams.originY) || ~isequal(size(extraParams.originY),[1 1]) || extraParams.originY < 1
            extraParams.originY = Y/2;
        end
        
        if ~isfield(extraParams,'speed') || ~isnumeric(extraParams.length) || ~isequal(size(extraParams.length),[1 1]) || extraParams.length < 1
            extraParams.speed = 1;
        end
        
        if ~isfield(extraParams,'width') || ~isnumeric(extraParams.width) || ~isequal(size(extraParams.width),[1 1]) || extraParams.width < 1
            extraParams.width = 1;
        end
        
        extraParams.skipArgCheck = true;
    end
    
    if T > extraParams.maxT
        pixels = NaN;
        return
    end
    
    angle = extraParams.angle;
    speed = extraParams.speed;
    radius = (T-1)*speed;
    
    centreX = round(mod(radius*sin(angle)+extraParams.originX-1,X)+1);
    centreY = round(mod(radius*cos(angle)+extraParams.originY-1,Y)+1);
    
    len = extraParams.length;
    
    radii = [ceil(-len/2)-mod(len,2)+1 floor(len/2)];
    lineXs = round(radii*sin(angle+pi/2))+centreX;
    lineYs = round(radii*cos(angle+pi/2))+centreY;
    
    Screen('DrawLines', window, [lineXs; lineYs], extraParams.width, extraParams.colour);
    
    pixels = [];
end