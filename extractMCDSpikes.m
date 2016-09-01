function extractMCDSpikes(recording,varargin)
    [filedir,filename] = getAnalysisOutputDir(recording);

    [~,file] = ns_OpenFile([filename ' spikes.mcd']);
    [~,fileInfo] = ns_GetFileInfo(file);
    [~,entities] = ns_GetEntityInfo(file,1:fileInfo.EntityCount);
    entities = entities([entities.EntityType] == 3);
    
    saveFileFn = @(label,infix) [filedir '\' filename '_channel_' label '_MCD_' infix 'spikes.mat'];
    
    try
        channels = parseChannelOptions(recording.index,'spikes',true,varargin{:});
    catch err %#ok<NASGU>
        channels = 1:60;
    end
    
    nChannels = numel(channels);
    
    for ii = 1:nChannels
        tic;
        entity = entities(channels(ii));
        label = entity.EntityLabel(end-1:end);
        
%         if ~(strcmp(label,'73') || strcmp(label,'77'))
%             continue;
%         end
        
        nSpikes = entity.ItemCount;
        
        if isempty(nSpikes) || isnan(nSpikes) || isinf(nSpikes) || nSpikes <= 0
            warning('No spikes found on channel %s (%d/%d) in file %s\n',label,ii,nChannels,filename); %#ok<WNTAG>
            continue;
        end
        
        [~,time,spikes] = ns_GetSegmentData(file,ii,1:nSpikes);
        
        if size(spikes,2) ~= 75 % TODO : figure out number of samples per spike and test that instead
            spikes = spikes';
        end
        
        index = time'*100; %#ok<NASGU>
        
        save(saveFileFn(label,''),'index','spikes');
        
        spikes = spikes(:,10:57); %#ok<NASGU>
        
        save(saveFileFn(label,'trimmed_'),'index','spikes');
        
        fprintf('Finished extracting MCD spikes for channel %s (%d/%d) in file %s in %f seconds\n',label,ii,nChannels,filename,toc);
    end
    
    mungeWaveClusToNex(recording,channels);
end