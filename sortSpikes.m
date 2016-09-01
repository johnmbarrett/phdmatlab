function sortSpikes(recording,varargin)
    filedir = getAnalysisOutputDir(recording);
    
    options = getopt('channel chlabel forceclustered=''no'' ignorenoise=''no''',varargin{:});
    channelInfo = getMetaInfo(recording.index,'spont',true);
    channelIndices = parseChannelOptions(recording.index,'channel',options.channel,'chlabel',options.chlabel,'spont',true);
    
    nSpikes = 0;
    spikesSoFar = zeros(length(channelIndices),1);
    spikeTimesByChannel = cell(length(channelIndices),1);
    for ii = 1:length(channelIndices)
        spikeFile = [filedir '\times_' recording.dataFile '_channel_' channelInfo(channelIndices(ii)).label '_MCD_spikes.mat'];
        
        if ~exist(spikeFile,'file')
            continue;
        end
        
        load(spikeFile,'cluster_class')
        spikeTimesByChannel{ii} = cluster_class;
        spikesSoFar(ii) = nSpikes;
        nSpikes = nSpikes + size(cluster_class,1);
    end
    
    allSpikeTimes = zeros(nSpikes,1);
    allSpikeClusters = zeros(nSpikes,1);
    allSpikeChannels = zeros(nSpikes,1);
    
    for ii = 1:length(channelIndices)
        spikes = spikeTimesByChannel{ii};
        
        if isempty(spikes)
            continue;
        end
        
        indices = (1:length(spikes))+spikesSoFar(ii);
        
        allSpikeChannels(indices) = str2double(channelInfo(channelIndices(ii)).label)*ones(length(spikes),1);
        allSpikeClusters(indices) = spikes(:,1);
        allSpikeTimes(indices) = spikes(:,2)/100; % wave_clus saves spikes times in 100ths of seconds for some reason
    end
    
    [sortedSpikeTimes,sortedIndices] = sort(allSpikeTimes);
    sortedSpikeClusters = allSpikeClusters(sortedIndices);
    sortedSpikeChannels = allSpikeChannels(sortedIndices);
    
    cells = unique(sortrows([allSpikeChannels allSpikeClusters]),'rows');
    nCells = size(cells,1);
    cellIndices = zeros(max(cells(:,1)),max(cells(:,2))+1);
    
    for ii = 1:nCells
        cellIndices(cells(ii,1),cells(ii,2)+1) = ii;
    end
end