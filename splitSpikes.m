function splitSpikes(allSpikesRecording,recordings)
    [allDir,allFile] = getAnalysisOutputDir(allSpikesRecording);
    
    % TODO : remove dependence on getMetaInfo - can get 
    totalTimeSpan = NaN;
    
    safeLoadMCDLibrary;
    [err,file] = ns_OpenFile([allFile '.mcd']);
        
    if err
        warning('Combined recording %s missing, can not verify total time span',allFile); %#ok<WNTAG>
    else
        [err,fileInfo] = ns_GetFileInfo(file);

        if err
            warning('Unable to extract file info from %s, can not verify total time span',allFile); %#ok<WNTAG>
        else
            totalTimeSpan = fileInfo.TimeSpan;
        end
    end
        
    recordingStartTimes = zeros(numel(recordings),1);
    recordingTimeSpans = zeros(numel(recordings),1);
    outDirs = cell(size(recordingStartTimes));
    outFiles = cell(size(recordingStartTimes));
    for ii = 1:numel(recordings)
        if iscell(recordings)
            recording = recordings{ii};
        else
            recording = recordings(ii);
        end
        
        [outDir,outFile] = getAnalysisOutputDir(recording);
        outDirs{ii} = outDir;
        outFiles{ii} = outFile;

        [err,file] = ns_OpenFile([outFile '.mcd']);
        
        if err
            if ii == numel(recordings)
                recordingTimeSpans(ii) = Inf;
                continue;
            else
                error('Missing data file %i (%s), unable to split spikes',ii,outFile);
            end
        end

        closeFile = onCleanup(@() ns_CloseFile(file));

        [~,fileInfo] = ns_GetFileInfo(file);
        recordingTimeSpans(ii) = fileInfo.TimeSpan;
        recordingStartTimes(ii) = getMCDStartTime(fileInfo,outFile);
    end
    
    [recordingStartTimes,sortIndices] = sort(recordingStartTimes);
    recordings = recordings(sortIndices);
    recordingTimeSpans = recordingTimeSpans(sortIndices);
    outDirs = outDirs(sortIndices);
    outFiles = outFiles(sortIndices);
    
    cumulativeTimeSpans = [0; cumsum(recordingTimeSpans)];
    
    if ~isnan(totalTimeSpan) && round(cumulativeTimeSpans(end)*100) ~= round(totalTimeSpan*100)
        error('Length of concatenated file (%f) does not match combined length of individual files (%f)',totalTimeSpan,cumulativeTimeSpans(end));
    end
    
    allPhotodiodeFile = [allDir '\' allFile '_photodiode_timings.mat'];
    
    if exist(allPhotodiodeFile,'file')
        load(allPhotodiodeFile);
        allRecordingStartTime = recordingStartTime; %#ok<NODEF>
        allStimulusTimes = stimulusTimes-allRecordingStartTime; %#ok<NODEF>
        
        for jj = 1:numel(recordings)
            outDir = outDirs{jj};
            outFile = outFiles{jj};
            photodiodeFile = strrep(strrep(allPhotodiodeFile,allDir,outDir),allFile,outFile);
            
            recordingStartTime = recordingStartTimes(jj);
            stimulusTimes = allStimulusTimes(allStimulusTimes >= cumulativeTimeSpans(jj) & allStimulusTimes < cumulativeTimeSpans(jj+1));
            stimulusTimes = stimulusTimes - cumulativeTimeSpans(jj) + recordingStartTime; %#ok<NASGU>
            
            save(photodiodeFile,'stimulusTimes','recordingStartTime');
        end
    end
    
    cumulativeTimeSpans = cumulativeTimeSpans*100;
    recordingStartIndices = zeros(numel(recordings),1);
    channels = num2str(channelIndexToMCSChannelNumber(1:60)');
    
    for ii = 1:60
        channel = channels(ii,:);
        infix = [allFile '_channel_' channel '_MCD_'];
        spikeFiles = { ...
            [allDir '\' infix 'spikes.mat'] ...
            [allDir '\' infix 'trimmed_spikes.mat'] ...
            [allDir '\times_' infix 'trimmed_spikes.mat']};
        
        for jj = 1:numel(spikeFiles)
            allSpikeFile = spikeFiles{jj};
            
            if ~exist(allSpikeFile,'file')
                continue;
            end
            
            load(allSpikeFile);
            
            isClustered = exist('cluster_class','var');
            isWaveClus = exist('inspk','var');
            
            if isClustered
                allIndex = cluster_class(:,2)';
            elseif ~exist('index','var')
                warning('Unrecognised spike file format for file %s, ignoring...\n',allSpikeFile); %#ok<WNTAG>
                continue;
            else
                allIndex = index;
            end
            
            allSpikes = spikes;
            
            if isClustered
                allClusterClass = cluster_class;
                
                if isWaveClus
                    allInspk = inspk;
                    allIPermut = ipermut;
                end
            end
            
            for kk = 1:numel(recordings)
                outDir = outDirs{kk};
                outFile = outFiles{kk};
                spikeFile = strrep(strrep(allSpikeFile,allDir,outDir),allFile,outFile);
                
                validSpikes = find(allIndex > cumulativeTimeSpans(kk) & allIndex <= cumulativeTimeSpans(kk+1));
                
                if jj == 1
                    if isempty(validSpikes)
                        recordingStartIndices(kk) = nan;
                    else
                        recordingStartIndices(kk) = validSpikes(1);
                    end
                end
                
                spikes = allSpikes(validSpikes,:);
                
                if isClustered
                    cluster_class = allClusterClass(validSpikes,:);
                    cluster_class(:,2) = cluster_class(:,2) - cumulativeTimeSpans(kk);
                    
                    saveVars = {'spikes' 'cluster_class'};
                    
                    if isWaveClus
                        saveVars(end+(1:3)) = {'ipermut' 'inspk' 'par'};
                        ipermut = allIPermut(allIPermut < max(validSpikes));
                        allIPermut = setdiff(allIPermut,ipermut);
                        inspk = allInspk(validSpikes,:);
                    end
                    
                    save(spikeFile,saveVars{:});
                    clear spikes cluster_class ipermut inspk;
                else
                    index = allIndex(validSpikes);
                    index = index - cumulativeTimeSpans(kk);
                    save(spikeFile,'spikes','index');
                    clear spikes index;
                end
            end
            
            if isClustered
                clear par;
            end
        end
        
        save([allDir '\' allFile '_channel_' channel '_start_times.mat'],'recordingStartTimes','recordingStartIndices');
    end
end