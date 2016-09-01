function calculateSpikeTriggeredAverage(recording,varargin)
    if ~isstruct(recording)
        recordings = initRecordings;
        recording = recordings(str2double(recording));
    end
    
    if ~isfield(recording,'spont') || ~all(recording.spont)
        error('calculateSpikeTriggeredAverage only works for continuous MC_Rack data files');
    end
    
    load(recording.stimulusFile,'getExtraParams','getPixels','seed','textureRect','timeOffset','vbls','version');
    
    if isfield(recording,'textureRect')
        textureRect = recording.textureRect;
    elseif ~exist('textureRect','var')
        error(['Texture dimensions unknown.  Stimulus reconstruction will be impossible for ' recording.dataFile]);
    end

    extraParams = validateStimulusParams(getPixels,getExtraParams);

    textureWidth = diff(textureRect([1 3]));
    textureHeight = diff(textureRect([2 4]));
    
    opts = getopt('skipstimuli=0',varargin{:});
    frames = 30;

    stimulusOnly = false;
    try
        filedir = getAnalysisOutputDir(recording);
        
        [sortedSpikeTimes,sortedSpikeClusters,sortedSpikeChannels,cells,nCells,cellIndices,clustered] = concatenateSpikes(recording,varargin{:});

        photoTriggerFile = [filedir '\' recording.dataFile '_photodiode_timings.mat'];
        usePhotodiodeTimings = logical(exist(photoTriggerFile,'file'));
        if usePhotodiodeTimings
            load(photoTriggerFile);

    %         if numel(stimulusTimes) < numel(vbls)
    %             vblStep = round(numel(vbls)/numel(stimulusTimes));
    %         else
    %             vblStep = 1;
    %         end
        else
            error('Could not find photodiode timings for file %i (%s)',recording.index,recording.dataFile);
        end

        stimulusTimes = stimulusTimes(opts.skipstimuli+1:end) - recordingStartTime; %#ok<NODEF>

    %     load([filedir '\' recording.dataFile '_vsync_times.mat']);
%         sortedSpikeTimes = sortedSpikeTimes + recordingStartTime;

    %     vsyncsBefore = vsyncIndicess{1}; %#ok<USENS>
    %     vsyncsAfter = vsyncIndicess{2};

    %     vsyncTimes = vsyncTimes + timeOffset; %#ok<NODEF>

        

    %     integrationTime = 0.5;
    %     timeResolution = 0.025;
    %     
    %     ts = linspace(-integrationTime+timeResolution/2,-timeResolution/2,integrationTime/timeResolution);
        

    %     nDimensions = textureWidth*textureHeight;
    %     Cp = eye(nDimensions);
        nSpikes = zeros(numel(stimulusTimes),nCells);



        for ii = 1:numel(stimulusTimes)-1
            tic;
    %         intervals = {vsyncTimes(vsyncsBefore(ii)+[0 1]) vsyncTimes(vsyncsAfter(ii)+[0 1]) NaN};
    %         
    %         if usePhotodiodeTimings
    %             if ii == numel(vbls)
    %                 intervals{3} = stimulusTimes(ii) + [0 mean(diff(stimulusTimes))];
    %             else
    %                 intervals{3} = stimulusTimes([ii ii+1]);
    %             end
    %         end
    %         
    %         for jj = 1 % 1:2+usePhotodiodeTimings
    %             interval = intervals{jj};

                interval = stimulusTimes(ii+[0 1]);

                spikes = find(sortedSpikeTimes > interval(1) & sortedSpikeTimes <= interval(2));

                if isempty(spikes)
                    continue;
                end

                for kk = 1:length(spikes)
                    spike = spikes(kk);

    %                 time = sortedSpikeTimes(spike)-interval(1);
    %                 bin = ceil(time/timeResolution);
                    channel = sortedSpikeChannels(spike,:);
                    cluster = sortedSpikeClusters(spike);
                    cellIndex = cellIndices(channel(1),channel(2),cluster+1);

                    nSpikes(ii,cellIndex) = nSpikes(ii,cellIndex) + 1;
                end
    %         end
    %         end

            toc;
        end

    %     stas = zeros(textureHeight,textureWidth,frames,nCells,3); %nDimensions,nCells);
    catch err
        stimulusTimes = 1:54000;
        stimulusOnly = true;
    end
    
    writer = VideoWriter('Stim 2.mj2','Archival');
    % TODO : frame rate
    open(writer);
    
    stas = NaN;
    stimulus = resetStimuli(version,seed,opts.skipstimuli);
    for ii = frames:numel(stimulusTimes)-1
        tic;
        
        [pixels,extraParams,stimulus,ii] = getNextStimulus(ii,stimulus,getPixels,getExtraParams,extraParams,textureWidth,textureHeight); %#ok<FXSET>
        
        if any(isnan(pixels))
            break;
        end
        
        writeVideo(writer,uint8(pixels(1:extraParams.pixelSize:end,1:extraParams.pixelSize:end)));
        
        if stimulusOnly
            continue;
        end
        
%         if ~isequal([size(stas,1) size(stas,2)],[size(pixels,1) size(pixels,2)])
        if isnan(stas)
            stas = zeros([size(pixels,1) size(pixels,2) frames nCells]);
        end
        
%         for kk = 1 %:2+usePhotodiodeTimings
%             for tt = 1:numel(ts)
%                 for jj = 1:nCells
                    stas = stas + repmat(reshape(nSpikes(ii:-1:ii-frames+1,:),1,1,frames,nCells),[size(pixels,1) size(pixels,2) 1]).*repmat(pixels,[1 1 frames nCells]);
%                 end
%             end
%         end
        
        toc;
    end
    
    close(writer);
    
    NSpikes = sum(nSpikes);
    stas = stas./repmat(reshape(NSpikes,1,1,1,nCells),[size(stas,1) size(stas,2) frames 1 ]);
    stas(isnan(stas)) = 0;
%     stas = stas./repmat(reshape(NSpikes,1,1,numel(ts),nCells,3),[size(stas,1) size(stas,2) 1 1 1]);
%     stas(isnan(stas)) = 0;
%     stasT = stas';
%     Csta = zeros(nDimensions,nDimensions,nCells);
%     
%     for ii = 1:nCells
%         Csta(:,:,ii) = stas(:,ii)*stasT(ii,:);
%     end
%     
%     Cs = zeros(nDimensions,nDimensions,nCells);
%     
%     if str2double(version) < 7.7
%         rng(struct('Seed',0,'State',seed,'Type','twister'));
%     else
%         rng(seed);
%     end
%     
%     stimulus = resetStimuli(version,seed);
%     for ii = 1:numel(vbls)
%         tic;
%         
%         [pixels,extraParams,stimulus,ii] = getNextStimulus(ii,stimulus,getPixels,getExtraParams,extraParams,textureWidth,textureHeight); %#ok<FXSET>
%         
%         if any(isnan(pixels))
%             break;
%         end
%         
%         pixelsT = pixels';
%         Cpixels = pixels*pixelsT;
%         
%         for jj = 1:nCells
%             Cs(:,:,jj) = Cs(:,:,jj) + nSpikes(ii,jj)*(Cpixels - stas(:,jj)*pixelsT - pixels*stasT(jj,:) + Csta(:,:,jj));
%         end
%         
%         toc;
%     end
%     
%     Cs = Cs./repmat(reshape(NSpikes-1,1,1,nCells),[nDimensions nDimensions 1]);
    
    timingMethodFileSuffixes = {'vsync_before' 'vsync_after' 'phototrigger'};
    timingMethodTitleSuffixes = {'VSync before' 'VSync after' 'phototrigger'};
    averagingMethodFileSuffixes = {'all_spikes' 'first_spike'};
    averagingMethodTitleSuffixes = {'All Spikes' 'First Spike'};
    
%     figure;
%     for ii = 1 %:2+usePhotodiodeTimings
        for jj = 1:nCells
            sta = squeeze(stas(:,:,:,jj));
            minSTA = min(min(min(sta)));
            maxSTA = max(max(max(sta)));
            sta = 255*(sta-minSTA)/(maxSTA-minSTA);
            sta = reshape(sta,size(sta,1),size(sta,2),1,size(sta,3));
            
            channel = char(cells(jj,1:2));
            cluster = num2str(cells(jj,3));
            infix = '';
            
            if clustered(str2double(channel))
                infix = '_clustered';
            end
            
%             title(['STA Channel ' channel ' cluster ' cluster ' (' timingMethodTitleSuffixes{3} ', ' averagingMethodTitleSuffixes{1} ', spike count: ' num2str(NSpikes(jj)) ')']);
%             title(['STA Channel ' channel ' cluster ' cluster ' (' timingMethodTitleSuffixes{ceil(ii/2)} ', ' averagingMethodTitleSuffixes{2-mod(ii,2)} ', spike count: ' num2str(spikeCounts(jj)) ')']);

%             cdata = zeros(size(stas,1),size(stas,2),1,frames);
%             
%             for kk = 1:frames
%                 image(squeeze(sta(:,:,kk)));
%                 colormap(gray(256));
%                 I = getframe(gcf);
%                 cdata(:,:,:,kk) = I.cdata;
%             end
            
%             cdata = reshape(cdata,[size(cdata,1) size(cdata,2) 1 nCells]);
            imwrite(sta,[filedir '\sta_ ' recording.dataFile '_channel_' channel '_cluster_' cluster infix '_' timingMethodFileSuffixes{3} '_' averagingMethodFileSuffixes{1} '.gif'],'LoopCount',1,'DelayTime',0.1);
        end
%     end
    
    makeSTABrowser(recording);
end