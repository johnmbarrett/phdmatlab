function [pixels,extraParams] = grating(T,X,Y,window,extraParams)
    if ~isfield(extraParams,'skipArgCheck')
        if ~isfield(extraParams,'angle') || ~isnumeric(extraParams.angle) || ~isequal(size(extraParams.angle),[1 1]) || extraParams.angle < 0 || extraParams.angle > 2*pi
            extraParams.angle = 0;
        end
        
        if ~isfield(extraParams,'phaseOffset') || ~isnumeric(extraParams.phaseOffset) || ~isequal(size(extraParams.phaseOffset),[1 1]) || extraParams.phaseOffset < 0 || extraParams.phaseOffset > 2*pi
            extraParams.phaseOffset = 0;
        end
        
        if ~isfield(extraParams,'contrast') || ~isnumeric(extraParams.contrast) || ~isequal(size(extraParams.contrast),[1 1]) || extraParams.contrast < 0 || extraParams.contrast > 100
            extraParams.contrast = 100;
        end
        
        if ~isfield(extraParams,'luminance') || ~isnumeric(extraParams.luminance) || ~isequal(size(extraParams.luminance),[1 1]) || extraParams.luminance < 0 || extraParams.luminance > 255
            extraParams.luminance = 127.5;
        end
        
        if extraParams.contrast/100 > 255/extraParams.luminance - 1;
            extraParams.contrast = 100*(255/extraParams.luminance - 1);
        end
        
        if ~isfield(extraParams,'maxT') || ~isnumeric(extraParams.maxT) || numel(extraParams.maxT) ~= 1 || any(isnan(extraParams.maxT))
            extraParams.maxT = Inf;
        end
        
        if ~isfield(extraParams,'spatialPeriod') || ~isnumeric(extraParams.spatialPeriod) || ~isequal(size(extraParams.spatialPeriod),[1 1]) || extraParams.spatialPeriod < 1
            extraParams.spatialPeriod = 32;
        end
        
        if ~isfield(extraParams,'temporalPeriod') || ~isnumeric(extraParams.temporalPeriod) || ~isequal(size(extraParams.temporalPeriod),[1 1]) || extraParams.temporalPeriod < 1
            extraParams.temporalPeriod = Inf;
        end
        
        extraParams.skipArgCheck = true;
    end
    
    if T > extraParams.maxT
        pixels = NaN;
        return;
    end
    
    [Xgrid,Ygrid] = meshgrid(1:X,1:Y);
    
    thetaGrid = Xgrid*cos(extraParams.angle)+Ygrid*sin(extraParams.angle);
    
    pixels = extraParams.luminance*extraParams.contrast*sin(2*pi*(thetaGrid/extraParams.spatialPeriod+(T-1)/extraParams.temporalPeriod)+extraParams.phaseOffset)/100+extraParams.luminance;
end