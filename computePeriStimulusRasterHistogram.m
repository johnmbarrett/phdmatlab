function computePeriStimulusRasterHistogram(recording,varargin)
    if isstruct(recording) && (~isfield(recording,'spont') || ~all(recording.spont))
        error('calculateSpikeTriggeredAverage only works for continuous MC_Rack data files');
    end
    
    options = getopt('proportionalspiketimes=false rasterfilesuffix=false bw=0.1',varargin);

    filedir = getAnalysisOutputDir(recording);
    
    [sortedSpikeTimes,sortedSpikeClusters,sortedSpikeChannels,cells,nCells] = concatenateSpikes(recording,varargin);
    
    filePrefix = [filedir '\' recording.dataFile];
    vsyncTimesFile = [filePrefix '_vsync_times.mat'];
    photodiodeTimesFile = [filePrefix '_photodiode_timings.mat'];
    
    if exist(photodiodeTimesFile,'file')
        load(photodiodeTimesFile);
        usePhotodiodeTimes = true;
    elseif exist(vsyncTimesFile,'file');
        load(vsyncTimesFile);
        usePhotodiodeTimes = false;
    else
        error('Could not find timings file for file %i (%s), raster construction will be impossible',recording.index,recording.dataFile);
    end
    
    load(recording.stimulusFile,'getExtraParams','stimuli','timeOffset','vbls');
    
    extraParams = validateStimulusParams(struct([]),getExtraParams);
    
    if ~iscell(extraParams)
        extraParams = {extraParams};
    end
    
    if usePhotodiodeTimes && (numel(extraParams) + mod(numel(extraParams),2)) ~= numel(stimulusTimes)
        warning('Number of photodiode level changes does not match number of stimuli for file %i (%s), rasters may be nonsensical',recording.index,recording.dataFile); %#ok<WNTAG>
    end
    
    hash = @(values) ['value_' regexprep(num2str(horzcat(values)),'[. ]','_')];
    
    isFactorsDefined = true;
    if ~isfield(recording,'factors') || ~(iscell(recording.factors) || ischar(recording.factors))
        factors = {NaN};
        dictionary = {struct(hash(NaN),1)};
        levels = 1;
        nFactors = 1;
        valuess = {NaN}; %#ok<NASGU>
        isFactorsDefined = false;
    elseif isfield(recording,'superStimulus') && isnumeric(recording.superStimulus) && numel(recording.superStimulus) == 1 && isfinite(recording.superStimulus) && recording.superStimulus > 0 && recording.superStimulus <= numel(extraParams)
        if numel(recording.factors) ~= 1
            error 'superstimulus option only supported for single factor experiments';
        end

        factor = recording.factors;

        if iscell(recording.factors)
            factor = factor{1};
        end
            
        for hh = 1:recording.superStimulus
            superStimulus = extraParams{hh};
            superFactor = superStimulus.(factor);
            
            if hh < recording.superStimulus
                stimulusTimes = stimulusTimes((size(superFactor,1)+1):end);
                continue;
            end
            
            extraParams = cell(size(superFactor,1),1);

            for ii = 1:size(superFactor,1)
                params = struct(superStimulus);
                params.(factor) = superFactor(ii,1);
                extraParams{ii} = params;
            end
        end
    end
    
    if ~isfield(recording,'stepStimuli') || ~isnumeric(recording.stepStimuli) || numel(recording.stepStimuli) ~= 1 || ~isfinite(recording.stepStimuli) || recording.stepStimuli < 1
        stepStimuli = 1;
    else
        stepStimuli = recording.stepStimuli;
    end
    
    if ~isfield(recording,'groupStimuli') || ~isnumeric(recording.groupStimuli) || numel(recording.groupStimuli) ~= 1 || ~isfinite(recording.groupStimuli) || recording.groupStimuli < 1
        groupStimuli = 1;
    else
        groupStimuli = recording.groupStimuli;
    end
    
    if ~isfield(recording,'skipStimuli') || ~isnumeric(recording.skipStimuli) || numel(recording.skipStimuli) ~= 1 || ~isfinite(recording.skipStimuli) || recording.skipStimuli < 0
        skipStimuli = 0;
    else
        skipStimuli = recording.skipStimuli;
    end
    
    if ~isfield(recording,'skipFrames') || ~isnumeric(recording.skipFrames) || numel(recording.skipFrames) ~= 1 || ~isfinite(recording.skipFrames) || recording.skipFrames < 0
        skipFrames = skipStimuli;
    else
        skipFrames = recording.skipFrames;
    end
    
    if ~isfield(recording,'dropStimuli') || ~isnumeric(recording.dropStimuli) || numel(recording.dropStimuli) ~= 1 || ~isfinite(recording.dropStimuli) || recording.dropStimuli < 0
        dropStimuli = 0;
    else
        dropStimuli = recording.dropStimuli;
    end
    
    if ~isfield(recording,'dropFrames') || ~isnumeric(recording.dropFrames) || numel(recording.dropFrames) ~= 1 || ~isfinite(recording.dropFrames) || recording.dropFrames < 0
        dropFrames = dropStimuli;
    else
        dropFrames = recording.dropFrames;
    end
    
    stimulusIndices = (1+skipStimuli):stepStimuli:(size(extraParams,1)-dropStimuli);
    frameIndices = (1+skipFrames):stepStimuli:(size(stimulusTimes,1)-dropFrames);
    nStimuli = numel(stimulusIndices);
    
    if isFactorsDefined
        if ischar(recording.factors)
            factors = {recording.factors};
        else
            factors = recording.factors;
        end
        
        nFactors = numel(factors);
        dictionary = cell(nFactors,1);
        valuess = cell(nFactors,1);
        
        levels = zeros(1,nFactors);

        for ii = 1:nFactors
            factor = factors{ii};
            values = [];

            for jj = 1:nStimuli
                stimulusIndex = stimulusIndices(jj);
                params = extraParams{stimulusIndex};
                
                if isempty(params)
                    warning('Missing params for stimulus %d in file %d (%s), assuming it and all following stimuli were never present.',stimulusIndices(jj),recording.index,recording.dataFile); %#ok<WNTAG>
                    stimulusIndices = stimulusIndices(1:jj-1);
                    nStimuli = jj-1;
                    break;
                end
                
                if ~isfield(params,'maxT')
                    error(['Not all stimuli specify a presentation time, raster construction will be impossible for ' recording.dataFile]);
                end
                
                value = getFactorValue(params,factor);
           
                for kk = 1:size(value,1)
                    if ~isempty(values) && ismember(value(kk,:),values,'rows')
                        continue;
                    end

                    values = [values;value(kk,:)]; %#ok<AGROW>

                    dictionary{ii}.(hash(value(kk,:))) = size(values,1);
                    levels(ii) = levels(ii) + 1;
                end
            end
            
            valuess{ii} = values;
        end
    end
    
    
    isCumulativeMaxT = false;
    if isfield(recording,'isCumulativeMaxT') && numel(recording.isCumulativeMaxT) == 1
        isCumulativeMaxT = logical(recording.isCumulativeMaxT);
    end

    sortedSpikeTimes = sortedSpikeTimes+recordingStartTime;
    
    if ~usePhotodiodeTimes
        sortedSpikeTimes = sortedSpikeTimes-timeOffset;
    end

%     N = numel(vbls);
%     firstFrame = find(vsyncTimes < vbls(1),1,'last');
%     lastFrame = find(vsyncTimes < vbls(end),1,'last'); %#ok<COLND>
%     
%     vsyncTimes = vsyncTimes(firstFrame) + (0:N-1)'*(vsyncTimes(lastFrame)-vsyncTimes(firstFrame))/(N-1);
%     vsyncIndicess = {1:N 1:N};
    
    rasters = cell([nCells levels 2]);
    stimulusTimings = cell([levels 2]);
    stimulusLengths = zeros([levels 2]);
    subStimulusLengths = zeros([levels groupStimuli 2]);
    
    nRasters = numel(rasters);
    
    for ii = 1:nRasters
        rasters{ii} = {};
    end
    
    tt = 0;
    escape = false;
    for ii = 1:nStimuli
        stimulusIndex = stimulusIndices(ii);
        params = extraParams{stimulusIndex};
        
        frameIndex = frameIndices(ii);
        
        index = zeros(1,nFactors);
        if ~iscell(factors)
            index = 1;
        else
            for jj = 1:nFactors
                index(jj) = dictionary{jj}.(hash(getFactorValue(params,factors{jj})));
            end
        end
        
        if options.proportionalspiketimes
            stimulusLength = 0;
            
            for jj = 1:groupStimuli
                subStimulus = extraParams{stimulusIndex+jj-1};
                
                if isempty(subStimulus) || ~isfield(subStimulus,'maxT')
                    continue;
                end
                
                subStimulusLength = subStimulus.maxT/60;
                
                for kk = 1:2
                    subStimulusLengthsSubscript = num2cell([index jj kk]);
                    subStimulusLengthsIndex = sub2ind(size(subStimulusLengths),subStimulusLengthsSubscript{:});
                    
                    if subStimulusLengths(subStimulusLengthsIndex) == 0
                        subStimulusLengths(subStimulusLengthsIndex) = subStimulusLength;
                    elseif subStimulusLengths(subStimulusLengthsIndex) ~= subStimulusLength
                        warning('Same combination of stimulus parameters presented for different lengths of time in file %d (%s), try using maxT as a factor',recording.index,recording.dataFile); %#ok<WNTAG>
                    end
                end
                
                stimulusLength = stimulusLength + subStimulusLength;
            end
            
            for jj = 1:2
                stimulusLengthsSubscript = num2cell([index jj]);
                stimulusLengthsIndex = sub2ind(size(stimulusLengths),stimulusLengthsSubscript{:});
                stimulusLengths(stimulusLengthsIndex) = stimulusLength;
            end
        end
        
        for jj = 1:2
            if usePhotodiodeTimes
                indices = frameIndex+(0:groupStimuli);
                stimulusTiming = nan(size(indices));
                valid = find(indices <= numel(stimulusTimes));
                
                if isempty(valid)
                    warning('Stimulus %d appears to have been presented after the end of the recording for file %d (%s), halting raster calculation',indices(1),recording.index,recording.dataFile); %#ok<WNTAG>
                    escape = true;
                    break;
                end
                
                stimulusTiming(valid) = stimulusTimes(indices(valid));
                
                for kk = valid(end)+1:numel(indices)
                    startTimeIndices = kk-1:groupStimuli:numel(stimulusTimes);
                    endTimeIndices = kk:groupStimuli:numel(stimulusTimes);
                    stimulusTiming(kk) = stimulusTiming(kk-1)+median(diff(stimulusTimes([startTimeIndices(1:numel(endTimeIndices)); endTimeIndices])));
                end
                
                subIntervals = stimulusTiming([1:end-1;2:end]);
                
                if size(subIntervals,1) ~= 2
                    subIntervals = subIntervals';
                end
                
                interval = stimulusTiming([1; end]);
            else
                vsyncIndices = vsyncIndicess{jj}; %#ok<USENS>

                subStimulusOnsets = nan(groupStimuli,1);
                subStimulusOffsets = nan(groupStimuli,1);

                if exist('stimuli','var')
                    for kk = 1:groupStimuli
                        subStimulusOnset = find(stimuli == frameIndex+kk-1,1);

                        if ~isempty(subStimulusOnset)
                            subStimulusOnsets(kk) = subStimulusOnset;
                        end

                        subStimulusOffset = find(stimuli == frameIndex+kk-1,1,'last');

                        if ~isempty(subStimulusOffset)
                            subStimulusOffsets(kk) = subStimulusOffset;
                        end
                    end
                else
                    subStimulusOnsets(1) = tt+1;

                    if isCumulativeMaxT
                        for kk = 1:groupStimuli
                            subStimulusOffsets(kk) = extraParams{stimulusIndex+kk-1}.maxT;
                        end
                    else
                        subStimulusOffsets(1) = tt + extraParams{stimulusIndex}.maxT;

                        for kk = 2:groupStimuli
                            subStimulusOffsets(kk) = subStimulusOffsets(kk-1) + extraParams{stimulusIndex+kk-1}.maxT;
                        end
                    end

                    for kk = 2:groupStimuli
                        subStimulusOnsets(kk) = subStimulusOffsets(kk-1)+1;
                    end
                end

                subStimulusOnsets = subStimulusOnsets(~isnan(subStimulusOnsets));
                subStimulusOffsets = subStimulusOffsets(~isnan(subStimulusOffsets));

                stimulusOnset = subStimulusOnsets(1);
                stimulusOffset = subStimulusOffsets(end);

                if stimulusOnset > length(vsyncIndices)
                    warning('Stimulus presentation appears to go on longer than recording for file %d (%s), # vbls: %d\t#vsyncs %d\n',recording.index,recording.dataFile,numel(vbls),numel(vsyncTimes)); %#ok<WNTAG>
                    escape = true;
                    break;
                end

                intervalIndices = [vsyncIndices(stimulusOnset); NaN];

                if stimulusOffset > length(vsyncIndices)
                    intervalIndices(2) = vsyncIndices(end)+1;
                else
                    intervalIndices(2) = vsyncIndices(stimulusOffset);
                end

                if any(isnan(intervalIndices))
                    warning('Stimulus #%d appears to have been presented before recording started for file %d (%s), tt: %d\tvbls(1): %8.6f\tvsyncTimes(1): %8.6f\n',stimulusIndex,recording.index,recording.dataFile,vbls(1),vsyncTimes(1)); %#ok<WNTAG>
                    continue;
                end
                
                subIntervalIndices = nan(2,numel(subStimulusOnsets));
            
                for kk = 1:numel(subStimulusOnsets)
                    if subStimulusOnsets(kk) <= length(vsyncIndices)
                        subIntervalIndices(1,kk) = vsyncIndices(subStimulusOnsets(kk));
                    end

                    if subStimulusOffsets(kk) <= length(vsyncIndices)
                        subIntervalIndices(2,kk) = vsyncIndices(subStimulusOffsets(kk));
                    end
                end

                if ~any(isnan(subIntervalIndices))
                    subIntervals = vsyncTimes(subIntervalIndices);

                    if ~options.proportionalspiketimes
                        for kk = 1:numel(subStimulusOnsets)
                            subStimulusLengthSubscript = num2cell([index kk jj]);
                            subStimulusLengthIndex = sub2ind(size(subStimulusLengths),subStimulusLengthSubscript{:});
                            subStimulusLengths(subStimulusLengthIndex) = max(subStimulusLengths(subStimulusLengthIndex),diff(subIntervals(:,kk)));
                        end
                    end
                elseif options.proportionalspiketimes
                    error('I can''t even remember in what circumstances this is supposed to happen');
                end   
                
                interval = vsyncTimes(intervalIndices);
                % interval = vbls(intervalIndices)+[-1; 1];
            end
            
            spikeIndices = find(sortedSpikeTimes >= interval(1)-1 & sortedSpikeTimes < interval(2)+1);
            spikes = sortedSpikeTimes(spikeIndices);
            stimulusLengthSubscript = num2cell([index jj]);
            stimulusLengthIndex = sub2ind(size(stimulusLengths),stimulusLengthSubscript{:});
            
            if options.proportionalspiketimes
                stimulusLength = 0;
                
                spikes(spikes < interval(1)) = spikes(spikes < interval(1)) - interval(1);
                actualStimulusLength = diff(interval);
                
                for kk = 1:size(subIntervals,2)
                    subInterval = subIntervals(:,kk);
                    subSpikes = find(spikes >= subInterval(1) & spikes < subInterval(2));
                    subStimulusLengthSubscript = num2cell([index kk jj]);
                    subStimulusLengthIndex = sub2ind(size(subStimulusLengths),subStimulusLengthSubscript{:});
                    subStimulusLength = subStimulusLengths(subStimulusLengthIndex);
                    spikes(subSpikes) = stimulusLength + subStimulusLength*(spikes(subSpikes) - subInterval(1))/diff(subInterval);
                    subIntervals(:,kk) = [stimulusLength; stimulusLength + subStimulusLength];
                    stimulusLength = subIntervals(2,kk);
                end
                
                spikes(spikes >= interval(2)) = spikes(spikes >= interval(2))-interval(1)-actualStimulusLength+stimulusLength;
            else
                spikes = spikes - interval(1);
                stimulusLengths(stimulusLengthIndex) = max(stimulusLengths(stimulusLengthIndex),diff(interval));
            end
                        
            channels = sortedSpikeChannels(spikeIndices,:);
            clusters = sortedSpikeClusters(spikeIndices);
            
            if isempty(stimulusTimings{stimulusLengthIndex})
                stimulusTimings{stimulusLengthIndex}(:,:,1) = subIntervals;
            else
                stimulusTimings{stimulusLengthIndex}(:,:,end+1) = subIntervals;
            end
            
            for kk = 1:nCells
                channel = cells(kk,1:2);
                cluster = cells(kk,3);
                
                cellSpikes = spikes(channels(:,1) == channel(1) & channels(:,2) == channel(2) & clusters == cluster);
                
                rasterSubscript = num2cell([kk index jj]);
                rasterIndex = sub2ind(size(rasters),rasterSubscript{:});
                
                raster = rasters{rasterIndex};
                raster{end+1} = cellSpikes; %#ok<AGROW>
                rasters{rasterIndex} = raster;
            end
        end
        
        if escape
            break;
        end
        
        if isCumulativeMaxT
            tt = extraParams{stimulusIndex+stepStimuli-1}.maxT;
        else
            for jj = 0:stepStimuli-1
                tt = tt + extraParams{stimulusIndex+jj}.maxT;
            end
        end
    end
    
    [histograms,edgess,repeats] = rastersToPSRH(rasters,stimulusLengths,options.bw); %#ok<NASGU,ASGLU>
    
    saveFilePrefix = [filedir '\' recording.dataFile '_psrh'];
    
    if ischar(options.rasterfilesuffix)
        saveFileSuffix = options.rasterfilesuffix;
    elseif options.proportionalspiketimes
        saveFileSuffix = 'rel';
    else
        saveFileSuffix = 'abs';
    end
    
    save([saveFilePrefix '_' saveFileSuffix '.mat'],'rasters','histograms','edgess','repeats','stimulusLengths','stimulusTimings','subStimulusLengths','factors','levels','valuess','options','cells');
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
        
        if isvector(rawValue)
            value = reshape(rawValue,1,numel(rawValue));
        else
            value = rawValue;
        end
    end
end