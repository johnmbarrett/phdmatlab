function batchFileWrapper(analysisFunction,functionDescription,recordingIndices,varargin)
    err = safeLoadMCDLibrary;
    
    if err
        error('Failed to load valid MCD Neuroshare library');
    end
    
    recordings = initRecordings;
    
    try
        for ii = 1:numel(recordingIndices)
            recordingIndex = recordingIndices(ii);
            recording = recordings(recordingIndex);

            tic;
            analysisFunction(recording,varargin{:});
            fprintf('Finished %s for recording %s (%d/%d) in %f seconds\n',functionDescription,recording.dataFile,ii,numel(recordingIndices),toc);
        end
    catch err
        logMatlabError(err);
    end
end