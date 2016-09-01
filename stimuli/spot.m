function f = spot(centre,radius,colour,maxT)
    if nargin < 4 || ~isnumeric(maxT) || numel(maxT) ~= 1 || isnan(maxT)
        maxT = Inf;
    end
    
    if nargin < 3 || ~isnumeric(colour) || ~all(isfinite(colour)) || ~ismember(numel(colour),[1 3 4])
        colour = 255;
    end
    
    colour = min(255,max(0,colour));
    
    if nargin < 2 || ~isnumeric(radius) || ~all(isfinite(radius)) || numel(radius) ~= 1
        radius = 0.5;
    end
    
    if nargin < 1
        error 'You must specify a centre location';
    end
    
    if ~isnumeric(centre) || ~all(isfinite(centre)) || numel(centre) ~= 2
        error 'You must specify the centre location as an (x,y) co-ordinate pair';
    end
    
    rect = [centre-floor(radius) centre+ceil(radius)];
    
    extraParams = struct('maxT',maxT);
    extraParams.colour = colour;
    extraParams.rect = rect;
    
    f = @(T,X,Y,W,P) getSpot(T,X,Y,W,extraParams);
end

function [pixels,extraParams] = getSpot(T,X,Y,window,extraParams)
    if T > extraParams.maxT
        pixels = NaN;
        return;
    end
    
    Screen('FillOval',window,extraParams.colour,extraParams.rect);
    pixels = [];
end
    