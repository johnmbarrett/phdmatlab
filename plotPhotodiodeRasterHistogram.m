function plotPhotodiodeRasterHistogram(recording,varargin)
    fileDir = getAnalysisOutputDir(recording);
    
    channelInfo = getMetaInfo(recording.index);
    opt = getopt('channel chlabel forceclustered ignorenoise',varargin{:});
    channels = parseChannelOptions(recording.index,'channel',opt.channel,'chlabel',opt.chlabel);
    
    photodiodeChannel = find(strcmp({channelInfo.label},'A1')); % TODO : specify in initRecordings
    
    % TODO : error handling
    [~,file] = ns_OpenFile([recording.dataFile '.mcd']);
    [~,~,photodiodeData] = ns_GetAnalogData(file,channelInfo(photodiodeChannel).numInFile,1,channelInfo(photodiodeChannel).segments);
    
    baseline = mode(photodiodeData);
    
    t1s = find(photodiodeData(1:end-1) == baseline & photodiodeData(2:end) < baseline);
    t2s = find(photodiodeData(1:end-1) < baseline & photodiodeData(2:end) == baseline);
    
    % TODO : nothing is right with this
    stimulusOnsets = t1s(1);
    tnext = stimulusOnsets(1)+10*25000;
    while tnext < t1s(end)
        stimulusOnsets(end+1) = t1s(find(t1s >= tnext,1)); %#ok<AGROW>
        tnext = stimulusOnsets(end)+10*25000;
    end
    
    for hh = 1:length(channels)
        channel = channels(hh);
        channelLabel = channelInfo(channel).label;
        
        % TODO : unclustered data?
        load([fileDir '\times_' recording.dataFile '_channel_' channelLabel '_MCD_spikes.mat'],'cluster_class');
        
        minClusters = min(cluster_class(:,1)); %#ok<NODEF>
        maxClusters = max(cluster_class(:,1));
        
        % TODO : extract sampling rate programmatically
        for jj = minClusters:maxClusters + 1;
            figure;
            hold on;
        
            for ii = 1:length(stimulusOnsets)
                stimulusChanges = find(t2s > stimulusOnsets(ii) & t2s < stimulusOnsets(ii)+15*25000);
            
                for kk = 1:length(stimulusChanges)
                    ll = stimulusChanges(kk);
                    fill(([t1s(ll) t2s(ll) t2s(ll) t1s(ll)] - stimulusOnsets(ii))/25000,[ii-1 ii-1 ii ii],[0.9 0.9 0.9],'LineStyle','none');
                end

                spikeIndices = find(cluster_class(:,2)/100 >= stimulusOnsets(ii)/25000-5 & cluster_class(:,2)/100 < stimulusOnsets(ii)/25000 + 10);
                
                if jj <= maxClusters
                    spikeIndices = spikeIndices(cluster_class(spikeIndices,1) == jj);
                end
                    
                spikeTimes = cluster_class(spikeIndices,2)'/100 - stimulusOnsets(ii)/25000;

                if isempty(spikeTimes)
                    continue;
                end

                line(repmat(spikeTimes,2,1),[ii-1;ii]*ones(size(spikeTimes)),'Color','k');
            end
        end
    end