function apsDataToPSRH(outputFile,spikeTimess,stimulusTimess,sequenceFile,conditions,blankStimulus,channelNames,factors,idx)
    if nargin < 9
        idx = 1:numel(stimulusTimess);
    end
    
    fin = fopen(sequenceFile,'r');
    
    for ii = 1:9
        fgetl(fin);
    end
    
    sequence = fscanf(fin,'%d\n');
    
    fclose(fin);
    
    firstStimulus = find(sequence ~= blankStimulus,1);
    sequence = sequence(firstStimulus:end);
    
    stimulusChangeIndices = [1; find(diff(sequence))+1];
    isBlank = sequence(stimulusChangeIndices) == blankStimulus;
    blankStimuli = find(isBlank);
    realStimuli = find(~isBlank);
    nStimuli = numel(realStimuli);
    
%     assert(isequal(size(conditions),[nStimuli numel(factors)]));
    
    options = struct('getopt__version',1,'proportionalspiketimes',false,'rasterfilesuffix',false,'bw',0.1); %#ok<NASGU>
    
    % 'rasters','stimulusTimings','stimulusLengths','subStimulusLengths','histograms','edgess','repeats'
    
    channelNames = vertcat(channelNames{1,idx});
    nCells = size(channelNames,1);
%     channels = unique(channelNames(:,1:end-1),'rows');
%     nChannels = size(channels,1);
    
    cells = [cellstr(channelNames(:,1:end-1)) num2cell(channelNames(:,end)-'a'+1)]; %#ok<NASGU>
    
    nFactors = numel(factors);
    levels = zeros(1,nFactors);
    valuess = cell(1,nFactors);
    
    for ii = 1:nFactors
        values = unique(conditions(:,ii),'rows');
        valuess{ii} = values;
        levels(ii) = numel(values);
    end
    
    if ~iscell(stimulusTimess)
        stimulusTimess = {stimulusTimess};
    end

    ns = zeros(levels);
    stimulusTimings = cell(levels);
    stimulusLengths = zeros(levels);
    subStimulusLengths = zeros([levels 2]);
    rasters = cell([nCells levels]);
    
    nBlocks = numel(stimulusTimess);
    nRepeats = nBlocks*nStimuli/prod(levels);

    for ii = 1:numel(stimulusTimings)
        stimulusTimings{ii} = zeros(2,2,nRepeats);
    end
    
    for ii = 1:numel(rasters)
        rasters{ii} = cell(1,nRepeats);
    end
    
    for hh = 1:nBlocks
        stimulusTimes = stimulusTimess{hh}(firstStimulus:end);
        stimulusTimes = stimulusTimes(stimulusChangeIndices);

        sumBlankTimes = 0;
        for ii = 1:nStimuli
            stimulus = sequence(stimulusChangeIndices(realStimuli(ii)));
            subscript = cell(1,nFactors);

            for jj = 1:nFactors
                subscript{jj} = find(ismember(valuess{jj},conditions(stimulus,jj)));
            end

            index = sub2ind(levels,subscript{:});
            ns(index) = ns(index)+1;

            timings = [stimulusTimes([realStimuli(ii) blankStimuli(ii)]); NaN];

            if ii == nStimuli
                timings(3) = timings(2) + sumBlankTimes/(nStimuli-1);
            else
                timings(3) = stimulusTimes(realStimuli(ii+1));
            end

            timings = timings([1 2; 2 3]);

            stimulusTimings{index}(:,:,ns(index)) = timings;

            subStimulusLength = diff(timings)';
            subStimulusLengths(subscript{:},:) = max(squeeze(subStimulusLengths(subscript{:},:)),subStimulusLength);

            stimulusLength = sum(subStimulusLength);
            stimulusLengths(index) = max(stimulusLength,stimulusLengths(index));

            sumBlankTimes = sumBlankTimes+subStimulusLength(2);

            for jj = 1:nCells
                spikeTimes = spikeTimess{idx(jj)};
                rasters{jj,subscript{:}}{ns(index)} = spikeTimes(spikeTimes >= timings(1,1)-1 & spikeTimes < timings(2,2)+1)-timings(1);
            end
        end
    end

    reps = [ones(size(levels)) 2];
    stimulusTimings = repmat(stimulusTimings,reps); %#ok<NASGU>
    stimulusLengths = repmat(stimulusLengths,reps);
    subStimulusLengths = repmat(subStimulusLengths,[ones(size(levels)) 1 2]); %#ok<NASGU>
    rasters = repmat(rasters,[1 ones(size(levels)) 2]);
    
    [histograms,edgess,repeats] = rastersToPSRH(rasters,stimulusLengths,0.1); %#ok<NASGU,ASGLU>
    
    save(outputFile,'-v7.3','options','cells','factors','levels','valuess','rasters','stimulusTimings','stimulusLengths','subStimulusLengths','histograms','edgess','repeats');
end