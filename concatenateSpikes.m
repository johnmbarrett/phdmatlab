function [sortedSpikeTimes,sortedSpikeClusters,sortedSpikeChannels,cells,nCells,cellIndices,clustered,allSpikeTimes,allSpikeClusters,allSpikeChannels] = concatenateSpikes(recording,varargin)
    [filedir,filename] = getAnalysisOutputDir(recording);
    
    options = getopt('channel chlabel forceclustered=''no'' ignorenoise=''no'' overwrite=''no''',varargin{:});
    
    saveFile = [filedir '\' filename '_spikes_concat_forceclustered_' options.forceclustered '_ignorenoise_' options.ignorenoise '.mat'];
    
    if exist(saveFile,'file') && ~strcmpi(options.overwrite,'yes')
        load(saveFile);
        return;
    end
    
    if isstruct(recording) && isfield(recording,'aps') && all(logical(recording.aps))
        [channelInfo,channelIndices] = getAPSChannelInfo;
    else
        try
            channelInfo = getMetaInfo(recording.index,'spont',true);
            channelIndices = parseChannelOptions(recording.index,'channel',options.channel,'chlabel',options.chlabel,'spont',true);
        catch err %#ok<NASGU>
            channelInfo = struct('label',cellstr(num2str(channelIndexToMCSChannelNumber(1:60)')));
            channelIndices = 1:60;
        end
    end
    
    nSpikes = 0;
    spikesSoFar = zeros(length(channelIndices),1);
    spikeTimesByChannel = cell(length(channelIndices),1);
    clustered = nan(87,1);
    
    spikeFileFn = @(channelLabel,prefix,infix) [filedir '\' prefix filename '_channel_' channelLabel '_MCD_' infix 'spikes.mat'];
    
    for ii = 1:length(channelIndices)
        channelLabel = channelInfo(channelIndices(ii)).label;
        clusterSpikeFile = spikeFileFn(channelLabel,'times_','');
        trimmedSpikeFile = spikeFileFn(channelLabel,'times_','trimmed_');
        allSpikeFile = spikeFileFn(channelLabel,'','');
        channelIndex = str2double(channelLabel);
        
        if exist(clusterSpikeFile,'file') || exist(trimmedSpikeFile,'file')
            if exist(clusterSpikeFile,'file')
                load(clusterSpikeFile,'cluster_class');
            elseif exist(trimmedSpikeFile,'file')
                load(trimmedSpikeFile,'cluster_class');
            end
            
            if ~isnan(channelIndex)
                clustered(channelIndex) = true;
            end
            
            if strcmpi(options.ignorenoise,'yes');
                cluster_class = cluster_class(cluster_class(:,1) > 0,:);
            end
        elseif ~strcmpi(options.forceclustered,'yes') && exist(allSpikeFile,'file')
            load(allSpikeFile,'index');
            
            if ~isnan(channelIndex)
                clustered(channelIndex) = false;
            end
            
            cluster_class = [ones(numel(index),1) index']; 
        else
            warning('No extracted spike timings found for channel %s (%d/%d) in file %s\n',channelLabel,ii,length(channelIndices),filename); %#ok<WNTAG>
            continue;
        end
        
        spikeTimesByChannel{ii} = cluster_class;
        spikesSoFar(ii) = nSpikes;
        nSpikes = nSpikes + size(cluster_class,1);
    end
    
    allSpikeTimes = zeros(nSpikes,1);
    allSpikeClusters = zeros(nSpikes,1);
    allSpikeChannels = zeros(nSpikes,2);
    
    for ii = 1:length(channelIndices)
        spikes = spikeTimesByChannel{ii};
        
        if isempty(spikes)
            continue;
        end
        
        indices = (1:size(spikes,1))+spikesSoFar(ii);
        
        allSpikeChannels(indices,:) = repmat(channelInfo(channelIndices(ii)).label,size(spikes,1),1);
        allSpikeClusters(indices) = spikes(:,1);
        allSpikeTimes(indices) = spikes(:,2)/100; % wave_clus saves spikes times in 100ths of seconds for some reason
    end
    
    [sortedSpikeTimes,sortedIndices] = sort(allSpikeTimes); 
    sortedSpikeClusters = allSpikeClusters(sortedIndices); 
    sortedSpikeChannels = allSpikeChannels(sortedIndices,:); 
    
    cells = unique(sortrows([allSpikeChannels allSpikeClusters]),'rows');
    nCells = size(cells,1);
    cellIndices = zeros(max(cells(:,1)),max(cells(:,2)),max(cells(:,3))+1);
    
    for ii = 1:nCells
        cellIndices(cells(ii,1),cells(ii,2),cells(ii,3)+1) = ii;
    end
    
    save(saveFile,'allSpikeTimes','allSpikeClusters','allSpikeChannels','sortedSpikeTimes','sortedSpikeClusters','sortedSpikeChannels','cells','nCells','cellIndices','clustered');
end