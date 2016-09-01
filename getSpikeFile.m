function [filepath,isClustered] = getSpikeFile(recording,channel,useClustered,isTrimmed)
    if nargin < 4
        isTrimmed = true;
    end
    
    % useClustered is trinary: 1 = force clustered, 0 = clustered if
    % available, -1 = force unclustered.
    % TODO : the else statement here is backwards compatability for the
    % legacy syntax, assuming that if it's unspecified then we want
    % clustered if available, if it's specified false then we definitely
    % want unclustered and if it's specified true then we definitely want
    % clustered.  Should really go through all usages and check, but MATLAB
    % is not a good IDE.
    if nargin < 3
        useClustered = 0;
    elseif islogical(useClustered)
        useClustered = 2*useClustered-1;
    end
    
    if ~ischar(channel)
        channel = num2str(channel);
    end
    
    [fileDir,dataFile] = getAnalysisOutputDir(recording);
    
    if isTrimmed
        infix = '_trimmed';
    else
        infix = '';
    end
    
    filename = sprintf('%s_channel_%s_MCD%s_spikes.mat', dataFile, channel, infix);
    isClustered = false;
    
    if useClustered > -1
        clusteredFilename = sprintf('times_%s',filename);
        
        if useClustered == 1 || exist(sprintf('%s\\%s',fileDir,clusteredFilename),'file')
            filename = clusteredFilename;
            isClustered = true;
        end
    end
    
    filepath = sprintf('%s\\%s',fileDir,filename);
end
    