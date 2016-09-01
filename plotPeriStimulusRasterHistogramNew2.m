function plotPeriStimulusRasterHistogram(recording,varargin)
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
        channelLabel = channelInfo(channelIndices(ii)).label;
        spikeFileSuffix = [recording.dataFile '_channel_' channelLabel '_spikes.mat'];
        clusterSpikeFile = [filedir '\times_' spikeFileSuffix];
        allSpikeFile = [filedir '\' spikeFileSuffix];
        
        if exist(clusterSpikeFile,'file')
            load(clusterSpikeFile,'cluster_class');
        elseif exist(allSpikeFile,'file')
            load(allSpikeFile,'index');
            cluster_class = [ones(numel(index),1) index']; %#ok<NODEF>
        else
            warning('No extracted spike timings found for channel %d (%d/%d) in file %d (%s)',channelLabel,ii,length(channelIndices),recording.index,recording.dataFile); %#ok<WNTAG>
            continue;
        end
        
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
        
    if ~iscell(extraParams)
        extraParams = {extraParams};
    end
    
    if ~isfield(recording,'groupStimuli') || ~isnumeric(recording.groupStimuli) || numel(recording.groupStimuli) ~= 1 || ~isfinite(recording.groupStimuli) || recording.groupStimuli < 1
        groupStimuli = 1;
    else
        groupStimuli = recording.groupStimuli;
    end
    
    nStimuli = numel(extraParams)/groupStimuli;
    
    hash = @(values) ['value_' regexprep(num2str(horzcat(values)),'[. ]','_')];
    
    if ~isfield(recording,'factors') || ~(iscell(recording.factors) && ~ischar(recording.factors))
        factors = NaN;
        dictionary = cell(0,1);
        levels = 1;
        nFactors = 1;
        valuess = {NaN};
    else
        if ischar(recording.factors)
            factors = {recording.factors};
        else
            factors = recording.factors;
        end
        
        nFactors = numel(recording.factors);
        dictionary = cell(nFactors,1);
        valuess = cell(nFactors,1);
        
        levels = zeros(1,nFactors);

        for ii = 1:nFactors
            factor = recording.factors{ii};
            values = [];

            for jj = 1:groupStimuli:nStimuli
                params = getExtraParams{jj};
                
                if ~isfield(params,'maxT')
                    error(['Not all stimuli specify a presentation time, raster construction will be impossible for ' recording.dataFile]);
                end
                
                value = getFactorValue(params,factor);
                
                if ismember(value,values,'rows')
                    continue;
                end

                values = [values;value];
                
                dictionary{ii}.(hash(value)) = size(values,1);
                levels(ii) = levels(ii) + 1;
            end
            
            valuess{ii} = values;
        end
    end
    
    isCumulativeMaxT = false;
    if isfield(recording,'isCumulativeMaxT')
        isCumulativeMaxT = logical(recording.isCumulativeMaxT);
    end
    
    recordingStartTime = recordingStartDate*24*3600-timeOffset;
    vsyncTimes = vsyncDates*24*3600-timeOffset;
    sortedSpikeTimes = sortedSpikeTimes+recordingStartTime;
    
    vsyncsBefore = zeros(numel(vbls),1);
    vsyncsAfter = zeros(numel(vbls),1);
    
    % need a faster way of doing this.  Can't just find the first and then
    % chop vsyncTimes appropriately because vsyncTimes and vbls aren't in
    % one to one correspondence due to missed frames
    for ii = 1:length(vbls)
        vsyncsBefore(ii) = find(vsyncTimes < vbls(ii),1,'last');
        vsyncsAfter(ii) = find(vsyncTimes >= vbls(ii),1);
    end
    
    vsyncIndicess = {vsyncsBefore vsyncsAfter};
    
    rasters = cell([nCells levels 2]);
    stimulusLengths = zeros([levels 2]);
    nRasters = numel(rasters);
    
    for ii = 1:nRasters
        rasters{ii} = {};
    end
    
    tt = 0;
    for ii = 1:groupStimuli:nStimuli
        params = extraParams{ii};
        
        index = zeros(1,nFactors);
        if ~iscell(factors)
            index = 1;
        else
            for jj = 1:nFactors
                index(jj) = dictionary{jj}.(hash(getFactorValue(params,factors{jj})));
            end
        end
        
        for jj = 1:2
            vsyncIndices = vsyncIndicess{jj};
            
            stimulusOnset = tt+1;
            
            if isCumulativeMaxT
            	stimulusOffset = extraParams{ii+groupStimuli-1}.maxT;
            else
                stimulusOffset = tt;
                for kk = 0:groupStimuli-1
                    stimulusOffset = stimulusOffset + extraParams{ii+kk}.maxT;
                end
            end
            
            intervalIndices = [vsyncIndices(stimulusOnset) NaN];
            if stimulusOffset > length(vsyncIndices)
                intervalIndices(2) = vsyncIndices(end)+1;
            else
                intervalIndices(2) = vsyncIndices(stimulusOffset);
            end
            
            interval = vsyncTimes(intervalIndices)+[-0.5; 0.5];
            spikeIndices = find(sortedSpikeTimes >= interval(1) & sortedSpikeTimes < interval(2));
            spikes = sortedSpikeTimes(spikeIndices) - interval(1);
            channels = sortedSpikeChannels(spikeIndices);
            clusters = sortedSpikeClusters(spikeIndices);
            
            stimulusLengthSubscript = num2cell([index jj]);
            stimulusLengthIndex = sub2ind(size(stimulusLengths),stimulusLengthSubscript{:});
            stimulusLengths(stimulusLengthIndex) = max(stimulusLengths(stimulusLengthIndex),diff(interval));
            
            for kk = 1:nCells
                channel = cells(kk,1);
                cluster = cells(kk,2);
                
                cellSpikes = spikes(channels == channel & clusters == cluster);
                
                rasterSubscript = num2cell([kk index jj]);
                rasterIndex = sub2ind(size(rasters),rasterSubscript{:});
                
                raster = rasters{rasterIndex};
                raster{end+1} = cellSpikes; %#ok<AGROW>
                rasters{rasterIndex} = raster;
            end
        end
        
        if isCumulativeMaxT
            tt = params.maxT;
        else
            tt = tt + params.maxT;
        end
    end
    
    histograms = cell(size(rasters));
    edgess = cell(size(rasters));
    repeats = zeros(size(rasters));
    
    for ii = 1:nRasters
        raster = rasters{ii};
        nLines = numel(raster);
        
        subscript = cell(1,ndims(repeats));
        [subscript{:}] = ind2sub(size(repeats),ii);
        repeats(ii) = nLines;
        
        stimulusLengthSubscript = subscript(2:end);
        stimulusLength = stimulusLengths(sub2ind(size(stimulusLengths),stimulusLengthSubscript{:}));
        stimulusLengthMS = ceil(stimulusLength*10);
        edges = linspace(0,stimulusLengthMS/10,stimulusLengthMS+1); % 10ms bins
        edgess{ii} = edges;
        
        for jj = 1:nLines
            histogram = histc(raster{jj},edges)/nLines;
            
            if isempty(histograms{ii})
                histograms{ii} = histogram;
            else
                histograms{ii} = histograms{ii} + reshape(histogram,size(histograms{ii}));
            end
        end
    end
    
    save([filedir '\' recording.dataFile '_psrh.mat'],'rasters','histograms','edgess','repeats','stimulusLengths','factors','levels','valuess');
    
    for ii = 1:nCells
        subscript = cell(1,ndims(rasters));
        for jj = 1:prod(levels)
            levelSubscript = cell(1,ndims(levels));
            [levelSubscript{:}] = ind2sub(size(levels),jj);
            
            for kk = 1:2
                subscript = 
end

function value = getFactorValue(params,factor)
    if iscell(factor)
        value = zeros(1,numel(factor));

        for ii = 1:numel(factor)
            value(ii) = params.(factor{ii});
        end
    else
        value = params.(factor);
    end
end