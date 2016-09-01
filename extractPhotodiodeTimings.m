function extractPhotodiodeTimings(recording,varargin)
    channelInfo = getMetaInfo(recording.index,'spont',true);
    
    [~,file] = ns_OpenFile([recording.dataFile '.mcd']);
    
    photodiodeChannelIndex = find(strcmp({channelInfo.label},recording.photodiodeChannel));
    
    [~,entityInfo] = ns_GetEntityInfo(file,photodiodeChannelIndex);
    [~,~,rawData] = ns_GetAnalogData(file,photodiodeChannelIndex,1,entityInfo.ItemCount);
    
    [b,a]=ellip(2,0.1,40,50*2/25000);
    filteredData = filtfilt(b,a,rawData);
    
    baseline = filteredData(1:1e5); % TODO : what if stim came on before 4s?
    threshold = mean(baseline) + 5*std(baseline);
    
    [~,analogInfo] = ns_GetAnalogInfo(file,photodiodeChannelIndex);
    
    thresholdCrossings = find(diff(filteredData <= threshold));
    recordingStartTime = getMCDStartTime(file,recording.dataFile);
    
    % Have to increment thresholdCrossings due to offset introduced by diff
    % but then have to decrement if we assume sample 1 happens at
    % recordingStartTime. +1-1 = 0, hence:
    stimulusTimes = recordingStartTime + thresholdCrossings/analogInfo.SampleRate; %#ok<NASGU>
    
    save([getAnalysisOutputDir(recording) '\' recording.dataFile '_photodiode_timings'],'recordingStartTime','stimulusTimes');
end