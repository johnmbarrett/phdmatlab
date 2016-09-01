function [outputdir,dataFile,parentDir] = getAnalysisOutputDir(recording)
    parentDir = nan;
    
    if ischar(recording)
        [~,outputdir] = fileparts(recording);
        dataFile = outputdir;
    elseif isstruct(recording)
        if isfield(recording,'fileDir')
            outputdir = recording.fileDir;
        else
            outputdir = sprintf('%03d_%s',recording.index,recording.dataFile);
        end
        
        if isfield(recording,'parentDir')
            parentDir = recording.parentDir;
        end
        
        dataFile = recording.dataFile;
    elseif iscell(recording)
        [outputdir,dataFile,parentDir] = getAnalysisOutputDir(recording{1});
    else
        disp(recording);
        error('Recording must be a struct or a string');
    end
    
    if ~exist(['./' outputdir],'dir')
        mkdir(pwd,outputdir);
    end
end

