function [pixels,extraParams] = movingBar(T,X,Y,extraParams)
    if ~isfield(extraParams,'skipArgCheck')
        if ~isstruct(extraParams)
            extraParams = struct();
        end
        
        if ~isfield(extraParams,'angle') || ~isnumeric(extraParams.angle) || ~isequal(size(extraParams.angle),[1 1]) || extraParams.angle < 0 || extraParams.angle > 2*pi
            extraParams.angle = 0;
        end
        
        if ~isfield(extraParams,'bgColour') || ~isnumeric(extraParams.bgColour) || ~(isequal(size(extraParams.bgColour),[1 1]) || isequal(size(extraParams.bgColour),[1 3])) || min(extraParams.bgColour < 0) || max(extraParams.bgColour > 255)
            extraParams.bgColour = 0;
        end
        
        if ~isfield(extraParams,'colour') || ~isnumeric(extraParams.colour) || ~(isequal(size(extraParams.colour),[1 1]) || isequal(size(extraParams.colour),[1 3])) || min(extraParams.colour < 0) || max(extraParams.colour > 255)
            extraParams.colour = 255;
        end
        
        if ~isequal(size(extraParams.bgColour),size(extraParams.colour))
            extraParams.bgColour = 255-extraParams.colour;
        end
        
        if ~isfield(extraParams,'length') || ~isnumeric(extraParams.length) || ~isequal(size(extraParams.length),[1 1]) || extraParams.length < 1
            extraParams.length = abs(X*sin(extraParams.angle)) + abs(Y*cos(extraParams.angle));
        end
        
        if ~isfield(extraParams,'maxT') || ~isnumeric(extraParams.maxT) || ~isequal(size(extraParams.maxT),[1 1]) || extraParams.maxT < 1
            extraParams.maxT = Inf;
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
    
    if length(extraParams.bgColour) == 3
        pixels = repmat(reshape(extraParams.bgColour,1,1,3),[X Y 1]);
    else
        pixels = extraParams.bgColour*ones(X,Y);
    end
    
    angle = extraParams.angle;
    speed = extraParams.speed;
    radius = (T-1)*speed;
    
    originX = round(mod(radius*sin(angle)+X/2-1,X)+1);
    originY = round(mod(radius*cos(angle)+Y/2-1,Y)+1);
    
    len = extraParams.length;
    
    radii = ceil(-len/2)-mod(len,2)+1:floor(len/2);
    lineXs = round(radii*sin(angle+pi/2))+originX;
    lineYs = round(radii*cos(angle+pi/2))+originY;
    
    validLineXs = lineXs(lineXs > 0 & lineXs <= X & lineYs > 0 & lineYs <= Y);
    validLineYs = lineYs(lineXs > 0 & lineXs <= X & lineYs > 0 & lineYs <= Y);
    
    for ii = 1:length(validLineXs)
        pixels(validLineXs(ii),validLineYs(ii),:) = extraParams.colour;
    end
end