function [datas,analogLabels,times,spikes,segmentLabels] = extractMCDSpikesExample(filename)
    % obviously you'll need to change this path to reflect wherever you
    % choose to install nsMCDLibrary
    ns_SetLibrary('H:\phd\matlab\nsMCDLibrary\Matlab\Matlab-Import-Filter\Matlab_Interface\nsMCDLibrary64.dll');
    [~,file] = ns_OpenFile(filename);
    [~,fileInfo] = ns_GetFileInfo(file);
    [~,entities] = ns_GetEntityInfo(file,1:fileInfo.EntityCount);
    entities = entities([entities.EntityType] == 2); % 2 corresponds to analog entities
    nChannels = numel(entities);
    
    datas = cell(nChannels,1);
    analogLabels = cell(nChannels,1);
    
    for ii = 1:nChannels
        entity = entities(ii);
        analogLabels{ii} = entity.EntityLabel;
        [~,~,data] = ns_GetAnalogData(file,ii,1,entity.ItemCount);
        datas{ii} = data;
    end
    
    entities = entities([entities.EntityType] == 3); % 3 corresponds to segment entities
    nChannels = numel(entities);
    
    times = cell(nChannels,1);
    spikes = cell(nChannels,1);
    segmentLabels = cell(nChannels,1);
    
    for ii = 1:nChannels
        entity = entities(ii);
        segmentLabels{ii} = entity.EntityLabel;
        
        nSpikes = entity.ItemCount;
        
        if isempty(nSpikes) || isnan(nSpikes) || isinf(nSpikes) || nSpikes <= 0
            warning('No spikes found on channel %s\n',label); %#ok<WNTAG>
            continue;
        end
        
        [~,time,spikes] = ns_GetSegmentData(file,ii,1:nSpikes);
        
        times{ii} = time;
        spikes{ii} = spikes;
    end
end