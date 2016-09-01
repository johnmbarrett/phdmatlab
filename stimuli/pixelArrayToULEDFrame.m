function commands = pixelArrayToULEDFrame(pixels)
    if ~isnumeric(pixels) || size(pixels,1) ~= 16 || size(pixels,2) ~= 16 || ndims(pixels) > 4
        error 'Pixels must be a 16x16xN matrix';
    end
    
    maxValue = max(max(max(pixels)));
    minValue = min(min(min(pixels)));
    
    if maxValue == minValue
        if maxValue ~= 0
            pixels = ones(size(pixels));
        end
    else
        pixels = round((pixels-minValue)./(maxValue-minValue));
    end
    
    commands = repmat(char(0),size(pixels,3),64);
    
    for t = 1:size(pixels,3)
        frame = [];

        for y = (1:16)+mod(1:16,2)-mod(0:15,2)
            for x = 0:4:12
                frame(end+1) = sum(pixels(y,x+(1:4),t).*[8 4 2 1]); %#ok<AGROW>
            end
        end
        
        commands(t,:) = sprintf('%x',frame);
    end
end