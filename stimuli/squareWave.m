function f = squareWave(colour,spatialPeriod,temporalPeriod,angle,maxT,bgColour,phase)
    if nargin < 1 || ~isnumeric(colour) || ~ismember(numel(colour), [1 3 4]) || any(colour < 0 | colour > 255 | ~isfinite(colour))
        colour = 255;
    end
    
    if nargin < 2 || ~isnumeric(spatialPeriod) || numel(spatialPeriod) ~= 1 || isnan(spatialPeriod) || spatialPeriod <= 0
        spatialPeriod = 1;
    end
    
    if nargin < 3 || ~isnumeric(temporalPeriod) || numel(temporalPeriod) ~= 1 || isnan(temporalPeriod) || temporalPeriod <= 0
        temporalPeriod = Inf;
    end
    
    if nargin < 4 || ~isnumeric(angle) || numel(angle) ~= 1 || ~isfinite(angle)
        angle = 0;
    end
    
    if nargin < 5 || ~isnumeric(maxT) || numel(maxT) ~= 1 || isnan(maxT) || maxT < 1
        maxT = Inf;
    end
    
    if nargin < 6 || ~isnumeric(bgColour) || ~ismember(numel(bgColour), [1 3 4]) || any(bgColour < 0 | bgColour > 255 | ~isfinite(bgColour))
        bgColour = zeros(size(colour));
    end
    
    if nargin < 7 || ~isnumeric(phase) || numel(phase) ~= 1 || ~isfinite(phase) || abs(phase) > 1
        phase = 0;
    end
    
    extraParams = struct('angle',angle,'spatialPeriod',spatialPeriod,'temporalPeriod',temporalPeriod,'maxT',maxT,'phase',phase);
    extraParams.colour = colour;
    extraParams.bgColour = bgColour;
    
    f = @(T,X,Y,W,P) getSquareWave(T,X,Y,W,extraParams);
end

function [pixels,extraParams,repeats,markerIndex] = getSquareWave(T,X,Y,window,extraParams)
    markerIndex = NaN;
    
    if ismember(extraParams.temporalPeriod,[0 Inf])
        repeats = extraParams.maxT;
    else
        repeats = 1;
    end
    
    if T > extraParams.maxT
        pixels = NaN;
        return;
    end
    
    [Xgrid,Ygrid] = meshgrid(0:X-1,0:Y-1);
    thetaGrid = Xgrid*cos(extraParams.angle)+Ygrid*sin(extraParams.angle);
    grating = floor(mod(2*(thetaGrid/extraParams.spatialPeriod+extraParams.phase)+2*(T-1)/extraParams.temporalPeriod,2));
    
    colourDepth = numel(extraParams.colour);
    pixels = zeros(Y,X,colourDepth);
    
    for ii = 1:colourDepth
        pixels(:,:,ii) = extraParams.colour(ii)*grating + extraParams.bgColour(ii)*(1-grating);
    end
end