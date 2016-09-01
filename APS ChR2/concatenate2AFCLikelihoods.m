function concatenate2AFCLikelihoods(independentVariable,decoder)
    filePrefix = sprintf('gratings_perf_vs_%s_decoder_%s',independentVariable,decoder);
    
    load(sprintf('%s_init.mat',filePrefix));
    
    likelihoodsFilePrefix = sprintf('%s_likelihoods',filePrefix);
    
    nCells = cellfun(@(R) size(R,2),Rs);
    
    gratingLikelihoods = arrayfun(@(n) zeros(nStimuli,n),nCells,'UniformOutput',false);
    maskLikelihoods = arrayfun(@(n) zeros(nStimuli,n),nCells,'UniformOutput',false);
    
    allFiles = cell(nStimuli,1);
    
    for jj = 1:nStimuli
        likelihoodsFile = sprintf('%s_trials_%d_control_cells_1-%d_drug_cells_1-%d.mat',likelihoodsFilePrefix,jj,nCells(1),nCells(2));

        likelihoods = load(likelihoodsFile);
        
        allFiles{jj} = likelihoodsFile;

        for ii = 1:2
            gratingLikelihoods{ii}(jj,:) = likelihoods.gratingLikelihoods{ii};
            maskLikelihoods{ii}(jj,:) = likelihoods.maskLikelihoods{ii};
        end
    end
    
    stimulusIter = 1:nStimuli; %#ok<NASGU>
    conditionIter = 1:2; %#ok<NASGU>
    cellIter = arrayfun(@(n) 1:n,nCells,'UniformOutput',false); %#ok<NASGU>
    
    save(sprintf('%s.mat',likelihoodsFilePrefix),'stimulusIter','conditionIter','cellIter','gratingLikelihoods','maskLikelihoods');
    
    delete(allFiles{:});
end