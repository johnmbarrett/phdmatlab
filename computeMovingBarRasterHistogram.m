function computeMovingBarRasterHistogram(recording,varargin)
    if ~isfield(recording,'spont') || ~all(recording.spont)
        error('calculateSpikeTriggeredAverage only works for continuous MC_Rack data files');
    end

    filedir = getAnalysisOutputDir(recording);
    
    [sortedSpikeTimes,sortedSpikeClusters,sortedSpikeChannels,cells,nCells] = concatenateSpikes(recording,varargin);
    
    filePrefix = [filedir '\' recording.dataFile];
    vsyncTimesFile = [filePrefix '_vsync_times.mat'];
    
    if exist(vsyncTimesFile,'file');
        load(vsyncTimesFile);
    else
        error('Could not find timings file for file %i (%s), raster construction will be impossible',recording.index,recording.dataFile);
    end
    
    load(recording.stimulusFile,'getExtraParams','stimuli','timeOffset','vbls');
    
    extraParams = validateStimulusParams(struct([]),getExtraParams); %#ok<USENS>
    
    if ~iscell(extraParams)
        extraParams = {extraParams};
    end
    
    if ~isfield(recording,'stepStimuli') || ~isnumeric(recording.stepStimuli) || numel(recording.stepStimuli) ~= 1 || ~isfinite(recording.stepStimuli) || recording.stepStimuli < 1
        stepStimuli = 1;
    else
        stepStimuli = recording.stepStimuli;
    end
    
    nStimuli = numel(extraParams)/stepStimuli;
    
    hash = @(values) ['value_' regexprep(num2str(horzcat(values)),'[. ]','_')];
    
    factors = {{'startX' 'startY'} 'length' 'width' 'speed'};
    nFactors = 4;
    
    dictionary = cell(nFactors,1);
    valuess = cell(nFactors,1);

    levels = zeros(1,nFactors);

    for ii = 1:nFactors
        factor = factors{ii};
        values = [];

        for jj = 1:stepStimuli:numel(extraParams)
            params = getExtraParams{jj};

            if ~isfield(params,'maxT')
                error(['Not all stimuli specify a presentation time, raster construction will be impossible for ' recording.dataFile]);
            end

            try 
                value = getFactorValue(params,factor);
            catch err %#ok<NASGU>
                if iscell(factor)
                    factorName = [factor{:}];
                else
                    factorName = factor;
                end
                
                error('Missing parameter %s for stimulus %d in file %i (%s), most likely the recording does not represent a moving bar experiment',factorName,jj,recording.index,recording.dataFile);
            end

            if ismember(value,values,'rows')
                continue;
            end

            values = [values;value]; %#ok<AGROW>

            dictionary{ii}.(hash(value)) = size(values,1);
            levels(ii) = levels(ii) + 1;
        end

        valuess{ii} = values;
    end
    
    isCumulativeMaxT = false;
    if isfield(recording,'isCumulativeMaxT') && numel(recording.isCumulativeMaxT) == 1
        isCumulativeMaxT = logical(recording.isCumulativeMaxT);
    end

    sortedSpikeTimes = sortedSpikeTimes+recordingStartTime-timeOffset;

%     N = numel(vbls);
%     firstFrame = find(vsyncTimes < vbls(1),1,'last');
%     lastFrame = find(vsyncTimes < vbls(end),1,'last'); %#ok<COLND>
%     
%     vsyncTimes = vsyncTimes(firstFrame) + (0:N-1)'*(vsyncTimes(lastFrame)-vsyncTimes(firstFrame))/(N-1);
%     vsyncIndicess = {1:N 1:N};
    
    rasters = cell([nCells levels 2]);
    
    starts = valuess{1};
    centre = mean(starts);
    nStarts = size(starts,1);
    
    widths = valuess{3};
    nWidths = numel(widths);
    
    stimulusDistances = zeros(nStarts,1);
    subStimulusLengths = zeros([nCells levels 3 2]);
   
    recordings = initRecordings;
    
    if ~isfield(recording,'rfFileIndex') || numel(recording.rfFileIndex) ~= 1 || ~isnumeric(recording.rfFileIndex) || ~isfinite(recording.rfFileIndex) || recording.rfFileIndex < 1 || recording.rfFileIndex > numel(recordings)
        isRFKnown = false;
    else
        isRFKnown = true;
        rfRecording = recordings(recording.rfFileIndex);
        rfDir = getAnalysisOutputDir(rfRecording); 
    end
    
    for ii = 1:nStarts
        start = starts(ii,:);
        finish = start+2*(centre-start);
        
        distance = round(sqrt(sum((finish-start).^2)));
        stimulusDistances(ii) = distance;
        
        if ~isRFKnown
            subStimulusLengths(:,ii,:,:,:,1,:) = distance;
            continue;
        end
        
        translation = eye(3);
        translation(1:2,3) = -start';
        
        theta = atan2(finish(2)-start(2),finish(1)-start(1));
        
        rotation = eye(3);
        rotation(1:2,1:2) = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        for jj = 1:nCells
            channel = cells(jj,1:2);
            cluster = cells(jj,3);
            rfFile = [rfDir '\rf_' rfRecording.dataFile '_channel_' channel '_cluster_' num2str(cluster) '_clustered_phototrigger_all_spikes.mat'];
    
            if ~exist(rfFile,'file')
                subStimulusLengths(jj,ii,:,:,:,1,:) = distance;
                continue;
            end
            
            load(rfFile,'centreX','centreY','sd');
            
            rfCentre = rotation*(translation*[centreX; centreY; 1]);
            
            nearEdge = rfCentre(1)-2*sd;
            
            % TODO : length if length isn't enough to cover retina
            for kk = 1:nWidths
                width = widths(kk);
                entry = max(0,nearEdge-width/2);
                crossing = width+4*sd;
                
                if entry+crossing > distance
                    exit = distance-entry;
                else
                    exit = distance-entry-crossing;
                end
                
                subStimulusLengths(jj,ii,:,kk,:,1,:) = entry;
                subStimulusLengths(jj,ii,:,kk,:,2,:) = crossing;
                subStimulusLengths(jj,ii,:,kk,:,3,:) = exit;
            end
        end
    end
    
    nRasters = numel(rasters);
    
    for ii = 1:nRasters
        rasters{ii} = {};
    end
    
    tt = 0;
    stimulusIndices = 1:stepStimuli:numel(extraParams);
    for ii = 1:nStimuli
        stimulusIndex = stimulusIndices(ii);
        params = extraParams{stimulusIndex};
        
        index = zeros(1,nFactors);
        if ~iscell(factors)
            index = 1;
        else
            for jj = 1:nFactors
                index(jj) = dictionary{jj}.(hash(getFactorValue(params,factors{jj})));
            end
        end
        
        for jj = 1:2
            vsyncIndices = vsyncIndicess{jj}; %#ok<USENS>

            if exist('stimuli','var')
                stimulusOnset = find(stimuli == stimulusIndex,1);
                stimulusOffset = find(stimuli == stimulusIndex,1,'last') + 1;
            else
                stimulusOnset = tt+1;
                stimulusOffset = (1-isCumulativeMaxT)*tt + extraParams{stimulusIndex}.maxT + 1;
            end
            
            raster = zeros(0,4);

            for kk = stimulusOnset:stimulusOffset-1
                barCentreLocation = stimulusDistances(index(1))*([kk; kk+1]-stimulusOnset)/(stimulusOffset-stimulusOnset);

                intervalIndices = [vsyncIndices(kk); NaN];                   

                if kk == numel(vsyncIndices)
                    intervalIndices(2) = vsyncIndices(kk)+1;
                else
                    intervalIndices(2) = vsyncIndices(kk+1);
                end

                interval = vsyncTimes(intervalIndices);

                spikeIndices = find(sortedSpikeTimes >= interval(1) & sortedSpikeTimes < interval(2));
                
                if isempty(spikeIndices)
                    continue;
                end
                
                distances = barCentreLocation(1)+(barCentreLocation(2)-barCentreLocation(1))*(sortedSpikeTimes(spikeIndices)-interval(1))/diff(interval);
                channels = sortedSpikeChannels(spikeIndices,:);
                clusters = sortedSpikeClusters(spikeIndices);
                
                raster = [raster; distances channels clusters]; %#ok<AGROW>
            end
            
            for kk = 1:nCells
                channel = cells(kk,1:2);
                cluster = cells(kk,3);

                if isempty(raster)
                    cellRaster = [];
                else
                    cellRaster = raster(raster(:,2) == channel(1) & raster(:,3) == channel(2) & raster(:,4) == cluster,1);
                end
                
                cellRasterIndex = num2cell([kk index jj]);
                cellRasters = rasters{cellRasterIndex{:}};
                cellRasters{end+1} = cellRaster; %#ok<AGROW>
                rasters{cellRasterIndex{:}} = cellRasters;
            end
        end
        
        if isCumulativeMaxT
            tt = extraParams{stimulusIndex+stepStimuli-1}.maxT;
        else
            for jj = 0:stepStimuli-1
                tt = tt + extraParams{stimulusIndex+jj}.maxT;
            end
        end
    end
    
    histograms = cell(size(rasters));
    edgess = cell(size(rasters));
    repeats = zeros(size(rasters));
    
    for i1 = 1:nCells
        for i2 = 1:nStarts
            for i3 = 1:size(rasters,3)
                for i4 = 1:nWidths
                    for i5 = 1:size(rasters,5)
                        for i6 = 1:2
                            index = sub2ind(size(rasters),i1,i2,i3,i4,i5,i6);
                            raster = rasters{index};
                            nLines = numel(raster);
        
                            repeats(index) = nLines;
                            edges = 0:stimulusDistances(i2)+1;
                            edgess{index} = edges;
        
                            for jj = 1:nLines
                                histogram = histc(raster{jj},edges);

                                if isempty(histogram)
                                    continue;
                                end

                                histogram = reshape(histogram,numel(histogram),1)/nLines;

                                if isempty(histograms{index})
                                    histograms{index} = histogram;
                                else
                                    histograms{index} = histograms{index} + histogram;
                                end
                            end
                            
                            if isempty(histograms{index})
                                histograms{index} = zeros(size(edges));
                            end
                        end
                    end
                end
            end
        end
    end
    
    save([filedir '\' recording.dataFile '_mbrh.mat'],'rasters','histograms','edgess','repeats','stimulusDistances','subStimulusLengths','factors','levels','valuess');
end

function value = getFactorValue(params,factor)
    if iscell(factor)
        value = zeros(1,numel(factor));

        for ii = 1:numel(factor)
            value(ii) = params.(factor{ii});
        end
    elseif isnan(factor)
        value = NaN;
    else
        rawValue = params.(factor);
        value = reshape(rawValue,1,numel(rawValue));
    end
end