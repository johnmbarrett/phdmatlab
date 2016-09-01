function [pixels,extraParams] = defaultGetPixels(T,X,Y,window,extraParams) %#ok<INUSL>
    if T > extraParams.maxT
        pixels = NaN;
        return
    end
    
    if ~isfield(extraParams,'pixelSize') || ~isnumeric(extraParams.pixelSize) || ~isfinite(extraParams.pixelSize)
        pixelSize = 1;
    else
        pixelSize = extraParams.pixelSize;
    end
    
    X = ceil(X/pixelSize);
    Y = ceil(Y/pixelSize);
    
    if isfield(extraParams,'isRandiDefined') && extraParams.isRandiDefined
        pix = (randi(2,X,Y)-1);
    else
        pix = (ceil(2*rand(X,Y))-1);
    end
    
    pix = kron(pix,ones(pixelSize,pixelSize));
    
    if ~isfield(extraParams,'colour') || ~isnumeric(extraParams.colour)
        pixels = 255*pix;
    elseif numel(extraParams.colour) == 1
        pixels = extraParams.colour*pix;
    else
        n = numel(extraParams.colour);
        pixels = zeros([size(pix) n]);
        for ii = 1:n
            pixels(:,:,ii) = extraParams.colour(ii)*pix;
        end
    end
end