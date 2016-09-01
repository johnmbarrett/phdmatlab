function [spiketimestamps,cells,vsyncTimes] = concatenateSpikeTimesAndStimulusTimes(recordings,outputFilePrefix,parsedTriggerss,stimFileIndices)
    if ischar(recordings)
        recordings = {recordings};
    end
    
    allVsyncTimes = [];
    allSpikeTimes = cell(60,1);
    nRecordings = numel(recordings);
    channelLabels = channelIndexToMCSChannelNumber(1:60);
    
    recordingStartTime = 0;
    
    for ii = 1:nRecordings
        recording = recordings{ii};
        
        vsyncFile = sprintf('%s\\%s_vsync_times.mat',recording,recording);
        
        nExpectedTriggers = sum(parsedTriggerss(stimFileIndices{ii}));
        
        if ~exist(vsyncFile,'file');
            warning('Trigger times file missing for recording %d (%s), filling in with NaNs\n',ii,recording);
            vsyncTimes = nan(nExpectedTriggers,1);
        else
            load(vsyncFile);
        end
        
        nRecordedTriggers = numel(vsyncTimes);
        
        if nRecordedTriggers ~= nExpectedTriggers
            warning('Mismatch between recorded and parsed numbers of triggers for recording %d (%s), replacing with NaNs\n',ii,recording);
            vsyncTimes = nan(nExpectedTriggers,1);
        end
        
        if ii == 1
            firstRecordingStartTime = recordingStartTime;
        end
        
        allVsyncTimes = [allVsyncTimes; vsyncTimes-firstRecordingStartTime]; %#ok<AGROW>
        
        offset = recordingStartTime - firstRecordingStartTime;
        
        for jj = 1:60
            spikesFile = sprintf('%s\\times_%s_channel_%d_MCD_trimmed_spikes.mat',recording,recording,channelLabels(jj));
            
            if ~exist(spikesFile,'file')
                continue;
            end
            
            load(spikesFile);
            
            [clusters,~,clusterIndices] = unique(cluster_class(:,1)); %#ok<NODEF>
            maxCluster = max(clusters);
            
            if ii == 1
                allSpikeTimes{jj} = cell(maxCluster+1,1);
            elseif maxCluster+1 > numel(allSpikeTimes{jj})
                allSpikeTimes{jj} = [allSpikeTimes{jj}; cell(maxCluster-numel(allSpikeTimes{jj})+1,1)];
            end
            
            for kk = 1:numel(clusters)
                cluster = clusters(kk);
                spikeTimes = cluster_class(clusterIndices == kk,2)/100+offset;
                allSpikeTimes{jj}{cluster+1} = [allSpikeTimes{jj}{cluster+1}; spikeTimes];
            end
        end
    end
    
    assert(issorted(allVsyncTimes(~isnan(allVsyncTimes))),'Stimulus times are not monotonically increasing.  Are you sure you have the files in the right order?');
    
    vsyncTimes = allVsyncTimes;
    recordingStartTime = 0; %#ok<NASGU>
    
    save(sprintf('%s_vsync_times.mat',outputFilePrefix),'vsyncTimes','recordingStartTime');
    
    spiketimestamps = vertcat(allSpikeTimes{:});
    
    assert(all(cellfun(@issorted,spiketimestamps)),'Spike times are not monotonically increasing.  Are you sure you have the files in the right order?');
    
    channels = cell2mat(arrayfun(@(n,ii) channelLabels(ii)*ones(n,1),cellfun(@numel,allSpikeTimes),(1:60)','UniformOutput',false));
    clusters = cell2mat(cellfun(@(C) (0:numel(C)-1)',allSpikeTimes,'UniformOutput',false));
    cells = [channels clusters];
    
    save(sprintf('%s_spiketimestamps.mat',outputFilePrefix),'spiketimestamps','cells');
    
    getClusterLabel = @(cl) char((cl == 0)*'U' + (cl > 0)*('a'+cl-1));
    
    channelNames = cell(9,size(cells,1));
    channelNames(1,:) = arrayfun(@(ch,cl) sprintf('Ch%02d%c',ch,getClusterLabel(cl)),channels',clusters','UniformOutput',false);
    channelNames(6,:) = num2cell(1*(clusters' == 0)); %#ok<NASGU>
    
    save(sprintf('%s_channelNames.mat',outputFilePrefix),'channelNames');
end