zips = dir('*.zip');
zips = {zips.name};

for ii = 1:numel(zips)
    zipFile = zips{ii};
    
    [~,tempDir] = fileparts(zipFile);
    
    [matchStart,~,~,~,tokens] = regexp(zipFile,'(.*)_channel_([0-9]+)_stas\.zip');
    
    if isempty(matchStart)
        continue;
    end
    
    dataFile = tokens{1}{1};
    
    if exist('./initRecordings.mat','file');
        load('./initRecordings.mat','recordings');
        
        recording = recordings(strcmp(dataFile,{recordings.dataFile}));
        
        fileDir = getAnalysisOutputDir(recording);
    else
        fileDir = getAnalysisOutputDir(dataFile);
    end
    
    channelIndex = str2double(tokens{1}{2});
    
    mkdir(tempDir);
    
    try
        unzip(zipFile,tempDir);
    catch err
        rmdir(tempDir);
        continue;
    end
    
    gifs = dir(sprintf('%s\\*.gif',tempDir));
    
    success = true;
    for jj = 1:numel(gifs)
        gif = gifs(jj);
        source = sprintf('%s\\%s',tempDir,gif.name);
        sourceChannel = sprintf('channel_%d',channelIndex);
        targetChannel = sprintf('channel_%d',channelIndexToMCSChannelNumber(channelIndex+1));
        target = sprintf('%s\\%s',fileDir,regexprep(gif.name,sourceChannel,targetChannel));
        success = success && movefile(source,target);
    end
    
    delete(sprintf('%s\\*',tempDir));
    rmdir(tempDir);
    
    if success
        delete(zipFile);
    end
end