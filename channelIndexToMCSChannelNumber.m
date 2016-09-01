function channels = channelIndexToMCSChannelNumber(indices)
    channels = nan(size(indices));
    
    channels(indices < 7) = indices(indices < 7) + 11;
    channels(indices > 54) = indices(indices > 54) + 27;
    
    middle = indices > 6 & indices < 55;
    middleIndices = indices(middle);
    
    channels(middle) = 10*(floor((middleIndices-7)/8)+2)+mod(middleIndices-7,8)+1;
end