function [pixels,extraParams] = movingBar(T,X,Y,window,extraParams)
    if ~isfield(extraParams,'skipArgCheck')
        if ~isstruct(extraParams)
            extraParams = struct();
        end
        
        if ~isfield(extraParams,'angle') || ~isnumeric(extraParams.angle) || ~isequal(size(extraParams.angle),[1 1]) || extraParams.angle < 0 || extraParams.angle > 2*pi
            extraParams.angle = pi/2;
        end
        
        if ~isfield(extraParams,'colour') || ~isnumeric(extraParams.colour) || ~(isequal(size(extraParams.colour),[1 1]) || isequal(size(extraParams.colour),[1 3])) || min(extraParams.colour < 0) || max(extraParams.colour > 255)
            extraParams.colour = 255;
        end
        
        if ~isfield(extraParams,'finishX') || ~isnumeric(extraParams.finishX) || ~isequal(size(extraParams.finishX),[1 1]) || extraParams.finishX < 1
            extraParams.finishX = X;
        end
        
        if ~isfield(extraParams,'finishY') || ~isnumeric(extraParams.finishY) || ~isequal(size(extraParams.finishY),[1 1]) || extraParams.finishY < 1
            extraParams.finishY = Y/2;
        end
        
        if ~isfield(extraParams,'length') || ~isnumeric(extraParams.length) || ~isequal(size(extraParams.length),[1 1]) || extraParams.length < 1
            extraParams.length = 256; %abs(X*sin(extraParams.angle)) + abs(Y*cos(extraParams.angle));
        end
        
        if ~isfield(extraParams,'maxT') || ~isnumeric(extraParams.maxT) || ~isequal(size(extraParams.maxT),[1 1]) || extraParams.maxT < 1
            extraParams.maxT = Inf;
        end
        
        if ~isfield(extraParams,'startX') || ~isnumeric(extraParams.startX) || ~isequal(size(extraParams.startX),[1 1]) || extraParams.startX < 1
            extraParams.startX = 0;
        end
        
        if ~isfield(extraParams,'startY') || ~isnumeric(extraParams.startY) || ~isequal(size(extraParams.startY),[1 1]) || extraParams.startY < 1
            extraParams.startY = Y/2;
        end
        
        if ~isfield(extraParams,'speed') || ~isnumeric(extraParams.length) || ~isequal(size(extraParams.length),[1 1]) || extraParams.length < 1
            extraParams.speed = 1;
        end
        
        if ~isfield(extraParams,'width') || ~isnumeric(extraParams.width) || ~isequal(size(extraParams.width),[1 1]) || extraParams.width < 1
            extraParams.width = 1;
        end
        
        extraParams.trajectoryLength = sqrt((extraParams.finishX-extraParams.startX)^2 + (extraParams.finishY-extraParams.startY)^2);
        
%         extraParams.slope = (extraParams.finishY-extraParams.startY)/(extraParams.finishX-extraParams.startX);
        extraParams.theta = atan2(extraParams.finishY-extraParams.startY,extraParams.finishX-extraParams.startX);
        
        extraParams.skipArgCheck = true;
    end
    
    if T > extraParams.maxT
        pixels = NaN;
        return
    end
    
    pixels = [];
    
    theta = extraParams.theta;
    Xt = round(mod((T-1)*extraParams.speed,extraParams.trajectoryLength)*cos(theta))+extraParams.startX;
    Yt = round(mod((T-1)*extraParams.speed,extraParams.trajectoryLength)*sin(theta))+extraParams.startY;
    
    phi = extraParams.angle;
    
    if extraParams.width < 7.5 % TODO : detect max width
        lineXs = [1 -1]*round(extraParams.length*cos(theta+phi)/2) + Xt;
        lineYs = [1 -1]*round(extraParams.length*sin(theta+phi)/2) + Yt;
    
        Screen('DrawLines', window, [lineXs; lineYs], extraParams.width, extraParams.colour);
        
        return
    end
    
    psi = atan2(extraParams.width/2,extraParams.length/2);
    diag = sqrt(extraParams.width^2 + extraParams.length^2)/2;
    
    angles = [theta+phi-psi; theta+phi+psi; theta+pi+phi-psi; theta+pi+phi+psi];
    vertexXs = ones(4,1).*round(diag*cos(angles)) + Xt;
    vertexYs = ones(4,1).*round(diag*sin(angles)) + Yt;
    
    Screen('FillPoly', window, extraParams.colour, [vertexXs vertexYs], 1);
end