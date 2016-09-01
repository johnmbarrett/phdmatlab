function extractVSyncTimings(recording,varargin)
    [fileDir,filename] = getAnalysisOutputDir(recording);

    [result,file] = ns_OpenFile([filename '.mcd']);
            
    if result
        error(['Could not open file: ' filename]);
    end

    [result,fileInfo] = ns_GetFileInfo(file);

    if result
        error(['Could not read file info from file: ' filename]);
    end

    recordingStartTime = getMCDStartTime(fileInfo,filename);

    [result,entities] = ns_GetEntityInfo(file,1:fileInfo.EntityCount);

    if result
        error(['Could not read entity info from file: ' filename]);
    end

    vSyncChannel = NaN;
    
    for ii = 1:vertcat(fileInfo.EntityCount);
        if strcmp(entities(ii).EntityLabel(end-1:end),'A2') %recording.vSyncChannel)
            vSyncChannel = ii;
            break;
        end
    end
    
    if isnan(vSyncChannel)
        warning('No or invalid vsync channel specificed for file %s, ignoring',filename); %#ok<WNTAG>
        return;
    end

    if entities(vSyncChannel).EntityType ~= 2
        error('VSync Channel is not an analog data channel');
    end

    [result,analogInfo] = ns_GetAnalogInfo(file,vSyncChannel);

    if result
        error(['Could not load analog info from file: ' filename']);
    end

    samplingRate = analogInfo.SampleRate;

    [result,~,analogData] = ns_GetAnalogData(file,vSyncChannel,1,entities(vSyncChannel).ItemCount);

    if result
        error(['Could not load analog data from file: ' filename']);
    end
    
    options = getopt('method=''baseline'' threshold=NaN mininterval=0 startsample=1 phase=''both''',varargin{:});
    
    phases = [];
    
    if strncmpi(options.phase,'bot',3)
        phases = [1 -1];
    elseif strncmpi(options.phase,'neg',3)
        phases = -1;
    elseif strncmpi(options.phase,'pos',3)
        phases = 1;
    end
    
    if strcmpi(options.method,'baseline')
        baseline = mode(analogData);
        vsyncIndices = find(analogData(1:end-1) == baseline & analogData(2:end) < baseline);
    elseif strcmpi(options.method,'threshold')
        vsyncIndices = find(ismember(diff(analogData(options.startsample:end) > options.threshold),phases));
    else
        error('Unknown extraction method: %s',options.method);
    end
    
    vsyncIndices = vsyncIndices([1; find(diff(vsyncIndices) >= options.mininterval*samplingRate/1000) + 1]);

    % Each value in vsyncIndices gives the index of a sample
    % immediately before a vSync pulse happens, so we need to add
    % one to get the right sample.  However, if we assume that the
    % first sample was taken at the timestamp given in fileInfo and
    % an event we were interested happened in the ith sample
    % (starting from 1, per Matlab conventions), then that event
    % happened (i-1)*samplesPerDay days after the first sample (let
    % the event happen in the first sample, so i = 1, then induce).
    % The plus one and the minus one cancel out, so pulseIndices *
    % samplesPerDay gives us the date of each vSync pulse.
%     vsyncDates = recordingStartDate + vsyncIndices/(samplingRate*60*60*24);
    vsyncTimes = recordingStartTime + vsyncIndices/samplingRate;
    
    savefile = [fileDir '\' filename '_vsync_times.mat'];
    
    if ~exist([filename '.mat'],'file')
        save(savefile,'vsyncIndices','vsyncTimes','recordingStartTime');
        return;
    end
    
    load([filename '.mat'],'timeOffset','vbls');
    
    %     N = numel(vbls);
%     vbls = vbls-recordingStartTime+timeOffset;
%     vsyncIndicess = {1:N 1:N};
    
    vsyncTimes = vsyncTimes-timeOffset;
    
    vsyncIndicesBefore = nan(numel(vbls),1);
    vsyncIndicesAfter = nan(numel(vbls),1);
    
    start = 1;
    
    if vbls(1) < vsyncTimes(1)
        start = find(vbls > vsyncTimes(1),1);
    end
    
%     % need a faster way of doing this.  Can't just find the first and then
%     % chop vsyncTimes appropriately because vsyncTimes and vbls aren't in
%     % one to one correspondence due to missed frames
    for ii = start:length(vbls)
        before = find(vsyncTimes < vbls(ii),1,'last');
        after = find(vsyncTimes >= vbls(ii),1);
%         
        if isempty(before) || isempty(after)
            break;
        end
        
        vsyncIndicesBefore(ii) = before;
        vsyncIndicesAfter(ii) = after;
%         vsyncsBefore = before+0:length(vbls)-1;
%         vsyncsAfter = after+0:length(vbls)-1;
    end

%     missed = find(misses >= 0);
%     frameInterval = mean(diff(vsyncTimes));
%     dvbls = diff([vbls(1)-frameInterval; vbls]);
%     dvbls2 = dvbls+[dvbls(2:end); dvbls(end)+frameInterval];
%     fudgeFactor = 2.15;
%     toffeeFactor = 0.02395;
%     
%     vsyncsBefore = nan(numel(vbls),1);
%     vsyncsBefore(2) = find(vsyncTimes <= vbls(2),1,'last');
%     vsyncsBefore(1) = vsyncsBefore(2)-1;
%     
%     for ii = 3:numel(vbls)
%         increment = 1;
%         
%         if ismember(ii,missed) && dvbls(ii)+dvbls(ii+1) > fudgeFactor*frameInterval
%             increment = ceil(dvbls(ii)/frameInterval);
%         elseif dvbls2(ii) < toffeeFactor
%             increment = 0;
%         end
%         
%         vsyncsBefore(ii) = vsyncsBefore(ii-1)+increment;
%     end
%     
%     vsyncsAfter = vsyncsBefore;

%     before = find(vsyncTimes < vbls(ii),1,'last');
%     after = find(vsyncTimes >= vbls(ii),1);
%     N = numel(vbls);
%         
%     vsyncsBefore = before+(0:N-1);
%     vsyncsAfter = after+(0:N-1);
%     
%     frameInterval = (vbls(end)-vbls(1))/N;
%     vsyncTimes = vsyncTimes(before)+((1:numel(vsyncTimes))'-before)*frameInterval;
    
%     missed = find(misses >= 0);
%     
%     if diff(vbls(1:2)) < 0
%         previousMiss = 2;
%     else
%         previousMiss = 1;
%     end
%     
%     for ii = 1:numel(missed)+1
%         if ii > numel(missed)
%             nextMiss = numel(vbls)+1;
%         else
%             nextMiss = missed(ii);
%         end
%         
%         nextSync = vbls(previousMiss);
%         vsyncBefore = find(vsyncTimes < nextSync,1,'last')-1;
%         vsyncAfter = find(vsyncTimes >= nextSync,1)-1;
%         vsyncsBefore(previousMiss:nextMiss-1) = vsyncBefore + (0:(nextMiss - previousMiss - 1));
%         vsyncsAfter(previousMiss:nextMiss-1) = vsyncAfter + (0:(nextMiss - previousMiss - 1));
%         previousMiss = nextMiss;
%     end
%     
%     if isnan(vsyncsBefore(1))
%         vsyncsBefore(1) = vsyncsBefore(2)-1;
%     end
%     
%     if isnan(vsyncsAfter(1))
%         vsyncsAfter(1) = vsyncsAfter(2)-1;
%     end
    
    vsyncIndicess = {vsyncIndicesBefore vsyncIndicesAfter}; %#ok<NASGU>

    save(savefile,'vsyncIndicess','vsyncTimes','recordingStartTime');
end