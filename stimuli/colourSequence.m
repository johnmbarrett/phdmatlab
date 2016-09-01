function f = colourSequence(colours,rect,period,maxT,useMarker,markerOffset)
    if nargin < 5 || ~all(logical(useMarker))
        markerOffset = NaN;
    elseif nargin < 6 || ~isscalar(markerOffset) || ~isnumeric(markerOffset) || ~isfinite(markerOffset)
        markerOffset = 0;
    end

    if nargin < 4
        maxT = Inf;
    end
    
    if nargin < 3
        period = Inf;
    end
    
    if nargin < 2
        rect = [0 0 640 480];
    end
    
    if nargin < 1
        colours = {0};
    end
    
    if ~iscell(colours)
        colours = {colours};
    end
    
    for ii = 1:length(colours)
        colour = colours{ii};
        
        if ~isnumeric(colour) || ~ismember(numel(colour), [1 3 4]) || any(colour < 0 | colour > 255 | ~isfinite(colour))
            error(['Colour ' num2str(ii) ' is not a valid colour']);
        end
    end
    
    extraParams = struct('period',period,'maxT',maxT);
    extraParams.colours = colours;
    extraParams.markerOffset = markerOffset;
    extraParams.rect = rect;
    
    f = @(T,X,Y,W,P) getColourSequence(T,X,Y,W,extraParams);
end
    
function [pixels,extraParams,repeats,markerIndex] = getColourSequence(T,X,Y,window,extraParams)
    repeats = min(min(extraParams.period),extraParams.maxT);

    if T > extraParams.maxT
        markerIndex = NaN;
        pixels = NaN;
        return
    end

    if numel(extraParams.period) == 1
        index = mod(floor((T-1)/extraParams.period),numel(extraParams.colours))+1;
    else
        index = find(mod(T-1,sum(extraParams.period))+1 <= cumsum(extraParams.period),1);
        repeats = extraParams.period(index);
    end
    
    if isfield(extraParams,'markerOffset')
        markerOffset = extraParams.markerOffset;
    else
        markerOffset = 0;
    end
    
    markerIndex = mod(index-1+markerOffset,2)+1;
    
    colour = extraParams.colours{index};
    
    if size(extraParams.rect,1) == 1
        rect = extraParams.rect;
    else
        rect = extraParams.rect(index,:);
    end
    
    Screen('FillRect', window, colour, rect);
    
    pixels = [];
end