function f = getColour(colour)
    f = @(T,X,Y,W,P) getColourHelper(T,X,Y,W,struct('colour',colour'));
end
    
function [pixels,extraParams] = getColourHelper(T,X,Y,W,extraParams)
    pixels = extraParams.colour*ones(X,Y);
end