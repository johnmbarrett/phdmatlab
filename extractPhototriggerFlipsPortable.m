function [stimulusTimes,recordingStartTime] = extractPhototriggerFlipsPortable(dataFile)
    [~,file] = ns_OpenFile(dataFile);
    [~,fileInfo] = ns_GetFileInfo(file);
    
    sampleRate = round(1/fileInfo.TimeStampResolution);
    
    [~,entities] = ns_GetEntityInfo(file,1:fileInfo.EntityCount);
    triggerChannel = find([entities.EntityType] == 1 & [entities.ItemCount] > 0 & strncmp('digi',{entities.EntityLabel},4));
    
    if numel(triggerChannel) < 1
        error('Could not find phototrigger channel for file %s',dataFile);
    elseif numel(triggerChannel) > 1
        warning('More than one phototrigger channel found for file %s, taking only the first',dataFile); %#ok<WNTAG>
        triggerChannel = triggerChannel(1);
    end
    
    [err,times] = ns_GetEventData(file,triggerChannel,1:entities(triggerChannel).ItemCount);
    
    if err
        digitalChannel = find([entities.EntityType] == 2 & strncmp('digi',{entities.EntityLabel},4));
       
        if numel(digitalChannel) < 1
            error('Could not find phototrigger channel for file %s',dataFile);
        elseif numel(digitalChannel) > 1
            warning('More than one phototrigger channel found for file %s, taking only the first',dataFile); %#ok<WNTAG>
            digitalChannel = digitalChannel(1);
        end
        
        [err,~,data] = ns_GetAnalogData(file,digitalChannel,1,entities(digitalChannel).ItemCount);
        
        if err
            error('Could not extract phototrigger data from file %s',dataFile);
        end
        
        times = find(diff(data))+1;
    else
        times = round(times*sampleRate);
    end
    
    stimulusTimes = times([1; find(diff(times) > 25) + 1])/sampleRate;
    
    recordingStartDate = datenum([ ...
        fileInfo.Time_Year ...
        fileInfo.Time_Month ...
        fileInfo.Time_Day ...
        fileInfo.Time_Hour ...
        fileInfo.Time_Min ...
        fileInfo.Time_Sec + fileInfo.Time_MilliSec/1000]);
    recordingStartTime = recordingStartDate*24*60*60;
end