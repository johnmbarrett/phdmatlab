function [pixels,extraParams,stimulus,frame] = getNextStimulus(frame,stimulus,getPixels,getExtraParams,extraParams,textureWidth,textureHeight)
    if iscell(getPixels)
        getPixel = getPixels{stimulus,1};
        [pixels,extraParams] = getPixel(frame,textureWidth,textureHeight,NaN,getExtraParams{stimulus});
    else
        [pixels,extraParams] = getPixels(frame,textureWidth,textureHeight,NaN,extraParams);
    end

    if any(isnan(pixels))
        if iscell(getPixels) && stimulus < size(getPixels,1)
            [pixels,extraParams,stimulus,frame] = getNextStimulus(frame - 1,stimulus + 1,getPixels,getExtraParams,extraParams,textureWidth,textureHeight);
        end
        
        return;
    end

    if isfield(extraParams,'textureRect') && isnumeric(extraParams.textureRect) && all(isfinite(extraParams.textureRect)) && numel(extraParams.textureRect) == 4
        textureRect = extraParams.textureRect;
        pixels = pixels(textureRect(2)+1:textureRect(4),textureRect(1)+1:textureRect(3),:);
    else
%         pixels = reshape(pixels,textureWidth*textureHeight,1);
%         pixels = (pixels/127.5)-1;
    end
end