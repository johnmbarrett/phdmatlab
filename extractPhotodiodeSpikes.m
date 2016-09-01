function extractPhotodiodeSpikes(recording,varargin)
    fileDir = getAnalysisOutputDir(recording);
    
    channelInfo = getMetaInfo(recording.index,'spont',true);
    
    photodiodeChannel = find(strcmp({channelInfo.label},'A1')); % TODO : specify in initRecordings
    
    % TODO : error handling
    [~,file] = ns_OpenFile([recording.dataFile '.mcd']);
    [~,~,photodiodeData] = ns_GetAnalogData(file,channelInfo(photodiodeChannel).numInFile,1,channelInfo(photodiodeChannel).segments);
    
    baseline = mode(photodiodeData);
    
    index = find((photodiodeData(1:end-1) == baseline & photodiodeData(2:end) < baseline) ...
               | (photodiodeData(1:end-1) < baseline & photodiodeData(2:end) == baseline));
           
    [x,y] = ndgrid(index,-25:49);
    spikes = photodiodeData(x+y); %#ok<NASGU>
    
    index = index'/250; %#ok<NASGU>
    
    save([fileDir '\' recording.dataFile '_channel_A1_MCD_spikes.mat'],'index','spikes');
end