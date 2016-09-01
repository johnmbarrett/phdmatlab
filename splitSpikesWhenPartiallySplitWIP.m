function splitSpikes(allSpikesRecording,recordings)
    allDir = getAnalysisOutputDir(allSpikesRecording);
    
    % TODO : remove dependence on getMetaInfo - can get 
    totalTimeSpan = NaN;
    
    safeLoadMCDLibrary;
    [err,file] = ns_OpenFile([allSpikesRecording '.mcd']);
        
    if err
        warning('Combined recording %s missing, can not verify total time span',allSpikesRecording); %#ok<WNTAG>
    else
        [err,fileInfo] = ns_GetFileInfo(file);

        if err
            warning('Unable to extract file info from %s, can not verify total time span',allSpikesRecording); %#ok<WNTAG>
        else
            totalTimeSpan = fileInfo.TimeSpan;
        end
    end
        
    recordingStartTimes = nan(numel(recordings),1);
    cumulativeTimeSpans = zeros(numel(recordings)+1,1);
    outDirs = cell(size(cumulativeTimeSpans));
    for ii = 1:numel(recordings)
        recording = recordings{ii};
        outDirs{ii} = getAnalysisOutputDir(recording);

        [err,file] = ns_OpenFile([recording '.mcd']);
        
        if err
            if ii == numel(recordings)
                cumulativeTimeSpans(ii+1) = Inf;
                continue;
            else
                warning('Missing data file %i (%s), maybe unable to split spikes',ii,recording);
                break
            end
        end

        closeFile = onCleanup(@() ns_CloseFile(file));

        [~,fileInfo] = ns_GetFileInfo(file);
        cumulativeTimeSpans(ii+1) = cumulativeTimeSpans(ii) + fileInfo.TimeSpan;
        recordingStartTimes(ii) = getMCDStartTime(fileInfo,recording);
    end
    
    startTimeFile = [allDir '\' allSpikesRecording '_channel_' channel '_start_times.mat'];
    
    if ~any(isnan(recordingStartTimes))
        if ~isnan(totalTimeSpan) && round(cumulativeTimeSpans(end)*100) ~= round(totalTimeSpan*100)
            error('Length of concatenated file (%f) does not match combined length of individual files (%f)',totalTimeSpan,cumulativeTimeSpans(end));
        end
    
        % TODO : split photodiode file when file partially split
        photodiodeFile = [allDir '\' allSpikesRecording '_photodiode_timings.mat'];

        if exist(photodiodeFile,'file')
            load(photodiodeFile);
            allRecordingStartTime = recordingStartTime; %#ok<NODEF>
            allStimulusTimes = stimulusTimes-allRecordingStartTime; %#ok<NODEF>

            for jj = 1:numel(recordings)
                recording = recordings{jj};
                outDir = outDirs{jj};
                outFile = strrep(strrep(photodiodeFile,allDir,outDir),allSpikesRecording,recording);

                recordingStartTime = recordingStartTimes(jj);
                stimulusTimes = allStimulusTimes(allStimulusTimes >= cumulativeTimeSpans(jj) & allStimulusTimes < cumulativeTimeSpans(jj+1));
                stimulusTimes = stimulusTimes - cumulativeTimeSpans(jj) + recordingStartTime; %#ok<NASGU>

                save(outFile,'stimulusTimes','recordingStartTime');
            end
        end
        
        cumulativeTimeSpans = cumulativeTimeSpans*100;
        recordingStartIndices = zeros(numel(recordings),1);
    else
        load(startTimeFile);
        cumulativeTimeSpans = cumsum([0 recordingStartTimes]);
    end
    
    channels = num2str(channelIndexToMCSChannelNumber(1:60)');
    
    for ii = 1:60
        channel = channels(ii,:);
        infix = [allSpikesRecording '_channel_' channel '_MCD_'];
        spikeFiles = { ...
            [allDir '\' infix 'spikes.mat'] ...
            [allDir '\' infix 'trimmed_spikes.mat'] ...
            [allDir '\times_' infix 'trimmed_spikes.mat']};
        
        for jj = 1:numel(spikeFiles)
            spikeFile = spikeFiles{jj};
            
            if ~exist(spikeFile,'file')
                continue;
            end
            
            load(spikeFile);
            
            isClustered = exist('cluster_class','var');
            isWaveClus = exist('inspk','var');
            
            if isClustered
                allIndex = cluster_class(:,2)';
            elseif ~exist('index','var')
                warning('Unrecognised spike file format for file %s, ignoring...\n',spikeFile); %#ok<WNTAG>
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
                recording = recordings{kk};
                outDir = outDirs{kk};
                outFile = strrep(strrep(spikeFile,allDir,outDir),allSpikesRecording,recording);
                
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
                    
                    save(outFile,saveVars{:});
                    clear spikes cluster_class ipermut inspk;
                else
                    index = allIndex(validSpikes);
                    index = index - cumulativeTimeSpans(kk);
                    save(outFile,'spikes','index');
                    clear spikes index;
                end
            end
            
            if isClustered
                clear par;
            end
        end
        
        save(startTimeFile,'recordingStartTimes','recordingStartIndices');
    end
end