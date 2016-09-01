currDir = pwd;

badBlocks = { ...
%     56  4   49;             ...
    57  4   [1:4 16:25];    ...
    };

letterExperiments = [56:59 61 62 63];

contGratingExperiments = 57:61; % 62 is where I started on the luminance gratings

for ii = 1:numel(contGratingExperiments)
    expt = contGratingExperiments(ii);
    
    isLetters = ismember(expt,letterExperiments);
    
    blocksIndex = ismember([badBlocks{:,1}],expt);
    
    if any(blocksIndex) && ismember(badBlocks{blocksIndex,2},[1 3+isLetters])
        blocks = badBlocks{blocksIndex,3};
    else
        blocks = 25;
    end
    
    try
        cd(sprintf('JBOG%04d',expt));

        grating2AFCAnalysis('contrast',1,3+isLetters,false,25,'jacobskde',struct('section','finish'));
        
        cd(currDir);
    catch err
        logMatlabError(err);
        cd(currDir);
    end
end

freqGratingExperiments = 57:62; % 63 was crap

for ii = 1:numel(freqGratingExperiments)
    expt = freqGratingExperiments(ii);
    
    isLetters = ismember(expt,letterExperiments);
    
    blocksIndex = ismember([badBlocks{:,1}],expt);
    
    if any(blocksIndex) && ismember(badBlocks{blocksIndex,2},[2 4+isLetters])
        blocks = badBlocks{blocksIndex,3};
    else
        blocks = 25;
    end
    
    try
        cd(sprintf('JBOG%04d',expt));

        grating2AFCAnalysis('frequency',2,4+isLetters,false,25,'jacobskde',struct('section','finish'));
        
        cd(currDir);
    catch err
        logMatlabError(err);
        cd(currDir);
    end
end
