function apsToStimulusSpikeTrainPairs(outputFile,spikeTimess,stimulusTimess,sequenceFile,blankStimulus,conditions,channelNames)
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
    blankStimuli = stimulusChangeIndices(isBlank);
    realStimuli = stimulusChangeIndices(~isBlank);
    nStimuli = numel(realStimuli);
    
    cells = channelNames(1,:)';
    responseIndices = find([channelNames{6,:}]' == 0);
    nCells = numel(responseIndices);
    responsiveCells = cells(responseIndices)'; %#ok<NASGU>
    spikeTimess = spikeTimess(responseIndices);
    
    widths = unique(conditions(:,1));
    nWidths = numel(widths);
    
    stimuli = zeros(nStimuli/nWidths,nWidths);
    
    for ii = 1:nWidths
        stimuli(:,ii) = conditions(sequence(realStimuli(conditions(sequence(realStimuli),1) == widths(ii))),2);
    end
    
    nBlocks = numel(stimulusTimess);
    stimuli = repmat(stimuli,nBlocks,1);
    
    responses = cell(nCells,size(stimuli,1),2,4);
    
    for ii = 1:nBlocks
        stimulusTimes = stimulusTimess{ii}(firstStimulus:end);
        
        sumBlankTimes = 0;
        ns = zeros(nWidths,1);
        
        for jj = 1:nStimuli
            widthIndex = find(ismember(widths,conditions(sequence(realStimuli(jj)),1)));
            ns(widthIndex) = ns(widthIndex) + 1;
            
            timings = [stimulusTimes([realStimuli(jj) blankStimuli(jj)]); NaN];

            if jj == nStimuli
                timings(3) = timings(2) + sumBlankTimes/(nStimuli-1);
            else
                timings(3) = stimulusTimes(realStimuli(jj+1));
            end
            
            for kk = 1:2
                interval = timings([kk kk+1]);
                
                for ll = 1:nCells
                    spikeTimes = spikeTimess{ll};
                    spikesInTrial = spikeTimes(spikeTimes >= interval(1) & spikeTimes < interval(2))-interval(1);
                    responses{ll,(nStimuli/nWidths)*(ii-1)+ns(widthIndex),kk,widthIndex} = spikesInTrial;
                end
            end
            
            sumBlankTimes = sumBlankTimes + diff(timings([2 3]));
        end
    end
    
    save(outputFile,'cells','responseIndices','responses','responsiveCells','stimuli','widths','-v7.3');
end