function [file,fileInfo,filename] = openMCDFile(recording)
    if isstruct(recording) && isfield(recording,'file') && isfield(recording,'fileInfo') && isfield(recording,'filename')
        file = recording.file;
        fileInfo = recording.fileInfo;
        filename = recording.filename;
        return;
    elseif ~ischar(recording)
        error('Recording must be a path to an MCD file or an MCD file handle and fileInfo struct\n');
    end
    
    filename = recording;
    
    if ~strcmp(filename(end-3:end),'.mcd')
        filename = sprintf('%s.mcd',filename);
    end
    
    safeLoadMCDLibrary;
    
    [err,file] = ns_OpenFile(filename);
    
    if err
        error('Unable to open file %s\n',recording);
    end
    
    [err,fileInfo] = ns_GetFileInfo(file);
    
    if err
        error('Unable to read file info for %s\n',recording);
    end
end