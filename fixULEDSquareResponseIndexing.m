function fixULEDSquareResponseIndexing(recording,responsiveCellFile)
    [fileDir,filename] = getAnalysisOutputDir(recording);
    
    saveFile = sprintf('%s\\%s_uled_square_responses_newmethod.mat',fileDir,filename);
    
    load(saveFile);
    
    allSNR = snr(:,:); %#ok<NODEF>
    someNSpikes = nSpikes; %#ok<NODEF>
    somePSpikes = pSpikes; %#ok<NODEF>
    someThresholds = thresholds; %#ok<NODEF>
    
    if exist('responseIndices','var')
        warning('Recording %s appears to have already been fixed; please verify\n',recording);
        return;
    end
    
    load(responsiveCellFile,'responsiveCells');
    nResponses = size(responsiveCells,1); %#ok<NODEF>
    
    if size(nSpikes,3) == nResponses
        responseIndices = (1:nResponses)'; %#ok<NASGU>
        cellIndices = zeros(nResponses,1);
        
        for ii = 1:nResponses
            cellIndices(ii) = find(ismember([channels clusters],responsiveCells(ii,:),'rows'));
        end
        
        snr = snr(:,cellIndices); %#ok<NASGU>
        
        save(saveFile,'responseIndices','cellIndices','snr','allSNR','-append')
        
        return;
    end
    
    cellIndices = [];
    responseIndices = [];
    
    for ii = 1:nResponses
        index = find(ismember([channels clusters],responsiveCells(ii,:),'rows'));
        
        if ~isempty(index)
            cellIndices(end+1,1) = index; %#ok<AGROW>
            responseIndices(end+1,1) = ii; %#ok<AGROW>
        end
    end
    
    if size(someNSpikes,3) < numel(responseIndices)
        warning('Recording %s appears not to have been reanalysed with all responsive cells.\n',recording);
        return;
    end
    
    snr = zeros(6,nResponses);
    snr(:,responseIndices) = allSNR(:,cellIndices); %#ok<NASGU>
    
    nSpikes = zeros(10,6,nResponses);
    nSpikes(:,:,responseIndices) = someNSpikes; %#ok<NASGU>
    
    pSpikes = zeros(6,nResponses);
    pSpikes(:,responseIndices) = somePSpikes; %#ok<NASGU>
    
    thresholds = inf(nResponses,1);
    thresholds(responseIndices) = someThresholds; %#ok<NASGU>
    
    save(saveFile,'responseIndices','cellIndices','snr','nSpikes','pSpikes','thresholds','allSNR','-append');
end