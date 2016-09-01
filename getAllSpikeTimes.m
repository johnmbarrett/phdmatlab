function [allSpikeTimes,channels,clusters] = getAllSpikeTimes(recording,ignoreNoise,flat)
    if flat
        allSpikeTimes = {};
        channels = '';
        clusters = [];
        n = 1;
    else
        allSpikeTimes = cell(60,1);
    end
    
    function fn(spikeTimes,channelIndex,channelLabel,cluster,varargin)
        if flat
            allSpikeTimes{n} = spikeTimes;
            channels(n,:) = channelLabel;
            clusters(n) = cluster;
            n = n + 1;
        else
            allSpikeTimes{channelIndex}{cluster+~ignoreNoise} = spikeTimes;
        end
    end

    forEachChannel(recording,[],ignoreNoise,@fn);
    
    if ~flat
        return;
    end
end