function calculateSpikeTriggeredAverage(recording,varargin)
    if ~isfield(recording,'spont') || ~all(recording.spont)
        error('calculateSpikeTriggeredAverage only works for continuous MC_Rack data files');
    end

    filedir = getAnalysisOutputDir(recording);
    
    options = getopt('channel chlabel',varargin{:});
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
    
    load([filedir '\' recording.dataFile '_vsync_times.mat']);
    load(recording.stimulusFile,'getExtraParams','getPixels','seed','textureRect','timeOffset','vbls','version');
    
    if ~exist('textureRect','var')
        if isfield(recording,'textureRect')
            textureRect = recording.textureRect;
        else
            error(['Texture dimensions unknown.  Stimulus reconstruction will be impossible for ' recording.dataFile]);
        end
    end
    
    if ~isstruct(getPixels) && ~isa(getPixels,'function_handle') && ~iscell(getPixels)
        error(['Unsupported stimulus type for ' recording.dataFile]);
    end
    
    if isa(getExtraParams,'function_handle')
        extraParams = getExtraParams();
    else
        extraParams = getExtraParams;
    end
    
    vsyncTimes = vsyncTimes-timeOffset; %#ok<NODEF>
    sortedSpikeTimes = sortedSpikeTimes+recordingStartTime-timeOffset;
    
    vsyncsBefore = zeros(numel(vbls),1);
    vsyncsAfter = zeros(numel(vbls),1);
    
    % need a faster way of doing this.  Can't just find the first and then
    % chop vsyncTimes appropriately because vsyncTimes and vbls aren't in
    % one to one correspondence due to missed frames
    for ii = 1:length(vbls)
        vsyncsBefore(ii) = find(vsyncTimes < vbls(ii),1,'last');
        vsyncsAfter(ii) = find(vsyncTimes >= vbls(ii),1);
    end
    
    textureWidth = diff(textureRect([1 3]));
    textureHeight = diff(textureRect([2 4]));
    
    stass = cell(4,1);
    spikeCountss = cell(4,1);
    for ii = 1:4
        stass{ii} = zeros(textureWidth,textureHeight,nCells);
        spikeCountss{ii} = zeros(nCells,1);
    end
    
    if str2double(version) < 7.7
        rng(struct('Seed',0,'State',seed,'Type','twister'));
    else
        rng(seed);
    end
    
    stimulus = 1;
    for ii = 1:length(vbls)
        tic;
        intervals = {vsyncTimes(vsyncsBefore(ii)+[0 1]) vsyncTimes(vsyncsAfter(ii)+[0 1])};
        
        if iscell(getPixels)
            getPixel = getPixels{stimulus};
            [pixels,extraParams] = getPixel(ii,textureWidth,textureHeight,NaN,getExtraParams{stimulus});
        else
            [pixels,extraParams] = getPixels(ii,textureWidth,textureHeight,NaN,extraParams);
        end
        
        if any(isnan(pixels))
            if iscell(getPixels)
                stimulus = stimulus + 1;
                ii = ii - 1; %#ok<FXSET>
                continue;
            else
                break;
            end
        end
        
        for jj = 1:2
            stasAllSpikes = stass{2*jj-1};
            stasFirstSpikes = stass{2*jj};
            allSpikesCounts = spikeCountss{2*jj-1};
            firstSpikesCounts = spikeCountss{2*jj};
            interval = intervals{jj};
            
            spikes = find(sortedSpikeTimes > interval(1) & sortedSpikeTimes < interval(2));
            
            if isempty(spikes)
                continue;
            end
            
            seen = [];
            for kk = 1:length(spikes)
                spike = spikes(kk);

                channel = sortedSpikeChannels(spike);
                cluster = sortedSpikeClusters(spike);
                cellIndex = cellIndices(channel,cluster+1);
                
                stasAllSpikes(:,:,cellIndex) = stasAllSpikes(:,:,cellIndex) + pixels;
                allSpikesCounts(cellIndex) = allSpikesCounts(cellIndex) + 1;
                
                if ismember(cellIndex,seen)
                    continue;
                end
                
                stasFirstSpikes(:,:,cellIndex) = stasFirstSpikes(:,:,cellIndex) + pixels;
                firstSpikesCounts(cellIndex) = firstSpikesCounts(cellIndex) + 1;
                seen = [seen;cellIndex]; %#ok<AGROW>
            end
            
            stass{2*jj-1} = stasAllSpikes;
            stass{2*jj} = stasFirstSpikes;
            spikeCountss{2*jj-1} = allSpikesCounts;
            spikeCountss{2*jj} = firstSpikesCounts;
        end
        toc;
    end
    
    timingMethodFileSuffixes = {'vsync_before' 'vsync_after'};
    timingMethodTitleSuffixes = {'VSync before' 'VSync after'};
    averagingMethodFileSuffixes = {'all_spikes' 'first_spike'};
    averagingMethodTitleSuffixes = {'All Spikes' 'First Spike'};
    
    figure;
    for ii = 1:4
        spikeCounts = spikeCountss{ii};
        stass{ii} = stass{ii}./repmat(reshape(spikeCounts,1,1,nCells),[textureWidth textureHeight 1]);
        stas = stass{ii};
        
        for jj = 1:nCells
            sta = squeeze(stas(:,:,jj));
            minSTA = min(min(sta));
            maxSTA = max(max(sta));
            sta = 255*(sta-minSTA)/(maxSTA-minSTA);
            image(sta);
            colormap(gray(256));
            channel = num2str(cells(jj,1));
            cluster = num2str(cells(jj,2));
            title(['STA Channel ' channel ' cluster ' cluster ' (' timingMethodTitleSuffixes{ceil(ii/2)} ', ' averagingMethodTitleSuffixes{2-mod(ii,2)} ', spike count: ' num2str(spikeCounts(jj)) ')']);
            I = getframe(gcf);
            imwrite(I.cdata,[filedir '\sta_ ' recording.dataFile '_channel_' channel '_cluster_' cluster '_' timingMethodFileSuffixes{ceil(ii/2)} '_' averagingMethodFileSuffixes{2-mod(ii,2)} '.png']);
        end
    end
end