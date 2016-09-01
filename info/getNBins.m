function [n,values,isInteger] = getNBins(x)
    if any(mod(x,1) ~= 0)
        values = unique(x);
        n = numel(values);
        isInteger = false;
    else
        n = max(x)+1;
        values = 0:max(x);
        isInteger = true;
    end
end