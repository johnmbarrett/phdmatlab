function shuffleResponses(infile,outfile,shuffleStyle,blocks)
    if nargin < 4
        blocks = NaN;
        nBlocks = NaN;
    elseif isscalar(blocks)
        nBlocks = blocks;
        blocks = 1:blocks;
    else
        assert(isnumeric(blocks),'blocks should be numeric');
        blocks = blocks(:);
        nBlocks = numel(blocks);
    end
    
    if nargin < 3
        shuffleStyle = 'shift';
    end
    
    assert(isequal(infile(end-3:end),'.mat'),'Input must be a MAT file');
    
    if nargin < 2
        outfile = strrep(infile,'.mat','_shuffle.mat');
    end
    
    load(infile);
    
    assert(exist('stimuli','var') & exist('responses','var'),'Input file must be include stimulus and response variables'); % TODO : specify var names
    
    stimSize = size(stimuli); %#ok<NODEF>
    nTrials = stimSize(1); % TODO : specify dim
    nLevels = prod(stimSize(2:end));
    nStimDims = ndims(stimuli);
    
    if iscell(responses) && isvector(responses) %#ok<NODEF>
        assert(all(cellfun(@(R) testResponseSize(R,nTrials,stimSize(2:end)),responses)));
        nCells = numel(responses);
    else
        nRespDims = ndims(responses);
        colons = repmat({':'},1,nRespDims-1);
        assert(testResponseSize(responses(colons{:},1),nTrials,stimSize(2:end)));
        nCells = size(responses,nRespDims);
    end
    
    for ii = 1:nLevels
        [uniqueStimuli,~,stimulusIndices] = unique(stimuli(:,ii));
        
        nStimuli = numel(uniqueStimuli);
        
        % TODO : limit blocks, other shuffling styles
        for jj = 1:nStimuli
            thisStimulusIndices = find(stimulusIndices == jj);
            
            if isnan(nBlocks)
                nBlocks = numel(thisStimulusIndices);
                blocks = 1:nBlocks;
            end
            
            if strcmpi(shuffleStyle,'shift')
                shuffleIndices = blocks([(2:nBlocks)';1]);
            elseif strncmpi(shuffleStyle,'rand',4)
                shuffleIndices = blocks(randperm(nBlocks));
            end
            
            if iscell(responses) && isvector(responses)
                for kk = 1:nCells
                    nRespDims = ndims(responses{kk});
                    colons = repmat({':'},1,nRespDims-nStimDims);
                    responses{kk}(thisStimulusIndices(blocks),colons{:},ii) = responses{kk}(thisStimulusIndices(shuffleIndices),colons{:},ii); %#ok<AGROW>
                end
            else
                colons = repmat({':'},1,nRespDims-nStimDims-1);
                responses(thisStimulusIndices(blocks),colons{:},ii,:) = responses(thisStimulusIndices(shuffleIndices),colons{:},ii,:); %#ok<AGROW>
            end
        end
    end
    
    save(outfile,'cells','responseIndices','responses','responsiveCells','stimuli','widths');
end

function isCorrectSize = testResponseSize(R,nBlocks,stimDims)
    sizeR = size(R);
    
    isCorrectSize = size(R,1) == nBlocks && isequal(sizeR(end-numel(stimDims)+1:end),stimDims);
end
    
    
    