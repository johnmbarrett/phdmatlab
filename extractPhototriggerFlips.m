function extractPhototriggerFlips(recording,varargin)
    [fileDir,filename] = getAnalysisOutputDir(recording);
    [~,file] = ns_OpenFile([filename '.mcd']);
    [~,fileInfo] = ns_GetFileInfo(file);
    
    sampleRate = round(1/fileInfo.TimeStampResolution);
    defaultMinInterval = ceil(sampleRate/1000)+1;
    
    options = getopt(['triggerprefix=''digi'' outputsuffix=''_photodiode_timings'' mininterval=' num2str(defaultMinInterval)],varargin{:});
    
    [~,entities] = ns_GetEntityInfo(file,1:fileInfo.EntityCount);
    triggerChannel = find([entities.EntityType] == 1 & [entities.ItemCount] > 0 & strncmp(options.triggerprefix,{entities.EntityLabel},4));
    
    if numel(triggerChannel) < 1
        error('Could not find phototrigger channel for file %s',filename);
    elseif numel(triggerChannel) > 1
        warning('More than one phototrigger channel found for file %s, taking only the first',filename); %#ok<WNTAG>
        triggerChannel = triggerChannel(1);
    end
    
    [err,times] = ns_GetEventData(file,triggerChannel,1:entities(triggerChannel).ItemCount);
    
    if err
        digitalChannel = find([entities.EntityType] == 2 & strncmp(options.triggerprefix,{entities.EntityLabel},4));
       
        if numel(digitalChannel) < 1
            error('Could not find phototrigger channel for file %s',filename);
        elseif numel(digitalChannel) > 1
            warning('More than one phototrigger channel found for file %s, taking only the first',filename); %#ok<WNTAG>
            digitalChannel = digitalChannel(1);
        end
        
        [err,~,data] = ns_GetAnalogData(file,digitalChannel,1,entities(digitalChannel).ItemCount);
        
        if err
            error('Could not extract phototrigger data from file %s',filename);
        end
        
        times = find(diff(data))+1;
    else
        times = round(times*sampleRate);
    end
    
    flipTimes = times([1; find(diff(times) >= options.mininterval) + 1])/sampleRate;
    
    recordingStartTime = getMCDStartTime(fileInfo,filename);
%     flipTimes = recordingStartTime+flipTimes; %#ok<NASGU>
%     
%     save([getAnalysisOutputDir(recording) '\' recording.dataFile '_flip_times'],'flipTimes');
    stimulusTimes = recordingStartTime + flipTimes; %#ok<NASGU>
    
    save([fileDir '\' filename options.outputsuffix],'recordingStartTime','stimulusTimes');
end