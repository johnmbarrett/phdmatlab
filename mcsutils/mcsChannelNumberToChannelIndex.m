function index = mcsChannelNumberToChannelIndex(number)
    index = zeros(size(number));
    
    index(number < 18) = number(number < 18) - 11;
    index(number > 80) = number(number > 80) - 27;
    
    middle = number > 20 & number < 79;
    middleChannels = number(middle);
    index(middle) = 8*(floor((middleChannels-21)/10))+mod(middleChannels,10)+6;
    
    index(number < 12 ...
        | number > 87 ... 
        | ismember(mod(number,10),[0 9]) ...
        | ismember(number,[18 81]) ...
        ) = nan;
end