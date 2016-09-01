function mungeWaveClusSpikesForBurstAnalysis(recordings)
    % see waves codes/readSpikes
    x = kron(0:9,ones(1,10));
    y = repmat(-1:8,1,10);
    
    e = channelIndexToMCSChannelNumber(1:60);
    function fn(spikeTimes,channelIndex,~,~,~,spikes)
        t{channelIndex} = kron(spikeTimes(:)*1000,ones(75,1))+repmat((0:74)'/25,numel(spikeTimes),1);
        peaks{channelIndex} = reshape(spikes',numel(spikes),1)*1e6;
    end

    for ii = 1:numel(recordings)
        t = cell(1,60);
        peaks = cell(1,60);
        e = zeros(1,60);
        recording = recordings{ii};
    
        forEachChannel(recording,[],true,@fn,false,false)
        
        save(sprintf('%s_spikes_for_burst_analysis.mat',recording),'t','e','x','y','peaks');
    end
end
        