function f = whiteNoiseWithFrameMarker(textureRect,markerRect,maxT,colour,pixelSize,extraParams)
    if nargin < 6 || ~isstruct(extraParams)
        extraParams = struct('maxT',maxT,'colour',colour,'pixelSize',pixelSize);
    else
        extraParams.maxT = maxT;
        extraParams.colour = colour;
        extraParams.pixelSize = pixelSize;
    end
    
    markers(:,:,:,1) = repmat(reshape(colour,1,1,numel(colour)),[diff(markerRect([2 4])) diff(markerRect([1 3])) 1]);
    markers(:,:,:,2) = zeros(size(markers));
    
    f = @(T,X,Y,W,P) whiteNoiseWithFrameMarkerHelper(T,X,Y,textureRect,markerRect,markers,extraParams);
end

function [pixels,extraParams,repeats,markerIndex] = whiteNoiseWithFrameMarkerHelper(T,X,Y,textureRect,markerRect,markers,extraParams)
    markerIndex = NaN;    
    repeats = 1;

    if T > extraParams.maxT
        pixels = NaN;
        return
    end
    
    colour = extraParams.colour;
    
    noise = defaultGetPixels(T,diff(textureRect([2 4])),diff(textureRect([1 3])),NaN,extraParams);
    
    pixels = zeros(Y,X,numel(colour));
    
    if isfield(extraParams,'bgColour')
        for ii = 1:numel(colour)
            pixels(:,:,ii) = extraParams.bgColour(ii);
        end
    end
    
    pixels(textureRect(2)+1:textureRect(4),textureRect(1)+1:textureRect(3),:) = noise;
    pixels(markerRect(2)+1:markerRect(4),markerRect(1)+1:markerRect(3),:) = squeeze(markers(:,:,:,mod(T-1,2)+1));
    
    extraParams.textureRect = textureRect;
    extraParams.markerRect = markerRect;
end